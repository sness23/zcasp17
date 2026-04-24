#!/usr/bin/env python3
"""Extract binding site residues for ligands in CASP15 PDB structures."""

import os
import math
from collections import defaultdict

PDB_DIR = '/home/sness/data/vaults/casp/casp15_ligands/pdb'
OUTPUT_FILE = '/home/sness/data/vaults/casp/casp15_ligands/binding_sites.md'

# Distance cutoff for binding site (Angstroms)
DISTANCE_CUTOFF = 4.5

# Ligands to analyze (exclude solvents/buffers)
SOLVENTS = {'HOH', 'WAT', 'DOD', 'H2O', 'SOL'}
BUFFERS = {'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'MPD', 'DMS', 'ACT', 'IMD', 'CIT',
           'TRS', 'MES', 'EPE', 'BOG', 'BME', 'DTT', 'CL', 'NA', 'K'}

# Standard amino acids
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M', 'SEC': 'C'  # Modified residues
}

# Target info with ligand names
TARGETS = {
    '7UTD': {'target': 'H1114', 'protein': '[NiFe]-hydrogenase Huc',
             'key_ligands': ['F3S', 'MQ9', '3NI', 'FCO']},
    '8B43': {'target': 'T1118v1', 'protein': 'FoxA transporter',
             'key_ligands': ['OX8']},
    '7UX8': {'target': 'T1124', 'protein': 'MfnG methyltransferase',
             'key_ligands': ['SAH', 'TYR']},
    '7Z8Y': {'target': 'H1135', 'protein': 'SUN1-KASH6 complex',
             'key_ligands': []},  # Protein-protein interface
    '8SX8': {'target': 'T1158v1', 'protein': 'MRP4 + PGE1',
             'key_ligands': ['XPG']},
    '8SXB': {'target': 'T1158v2', 'protein': 'MRP4 + PGE2',
             'key_ligands': ['P2E']},
    '8SX7': {'target': 'T1158v3', 'protein': 'MRP4 + DHEAS',
             'key_ligands': ['ZWY']},
    '8SWN': {'target': 'T1158v4', 'protein': 'MRP4 + ATP',
             'key_ligands': ['ATP']},
    '7PBR': {'target': 'T1170', 'protein': 'RuvB hexamer',
             'key_ligands': ['ADP', 'AGS']},
    '7PBL': {'target': 'H1171', 'protein': 'RuvB-RuvA complex',
             'key_ligands': ['ADP', 'AGS']},
    '7PBP': {'target': 'H1172', 'protein': 'RuvB-2xRuvA complex',
             'key_ligands': ['ADP', 'AGS']},
    '8AD2': {'target': 'T1187', 'protein': 'Nictaba lectin',
             'key_ligands': ['NAG']},
    '8C6Z': {'target': 'T1188', 'protein': 'Chitinase ChiB',
             'key_ligands': ['DW0']},
    '8XBP': {'target': 'T1127', 'protein': 'AtNATA1',
             'key_ligands': ['ACO']}
}

def distance(coord1, coord2):
    """Calculate Euclidean distance between two 3D coordinates."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

def parse_pdb_atoms(filepath):
    """Parse PDB file and return protein atoms and ligand atoms."""
    protein_atoms = []  # (chain, resname, resnum, atomname, x, y, z)
    ligand_atoms = defaultdict(list)  # ligand_id -> list of (atomname, x, y, z)

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                try:
                    atomname = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21]
                    resnum = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    protein_atoms.append((chain, resname, resnum, atomname, x, y, z))
                except (ValueError, IndexError):
                    continue

            elif line.startswith('HETATM'):
                try:
                    atomname = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21]
                    resnum = line[22:26].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    if resname not in SOLVENTS and resname not in BUFFERS:
                        lig_id = f"{resname}_{chain}_{resnum}"
                        ligand_atoms[lig_id].append((atomname, x, y, z, resname))
                except (ValueError, IndexError):
                    continue

    return protein_atoms, ligand_atoms

def find_binding_site(protein_atoms, ligand_atoms, cutoff=DISTANCE_CUTOFF):
    """Find protein residues within cutoff distance of ligand atoms."""
    binding_residues = set()  # (chain, resname, resnum)
    contacts = defaultdict(list)  # residue -> list of (lig_atom, prot_atom, distance)

    for lig_atom in ligand_atoms:
        lig_name, lx, ly, lz, _ = lig_atom
        lig_coord = (lx, ly, lz)

        for prot_atom in protein_atoms:
            chain, resname, resnum, atomname, px, py, pz = prot_atom
            prot_coord = (px, py, pz)

            d = distance(lig_coord, prot_coord)
            if d <= cutoff:
                res_id = (chain, resname, resnum)
                binding_residues.add(res_id)
                contacts[res_id].append((lig_name, atomname, round(d, 2)))

    return binding_residues, contacts

def format_residue(chain, resname, resnum):
    """Format residue as string."""
    aa = AA_3TO1.get(resname, resname)
    return f"{aa}{resnum}({chain})"

def main():
    output = []
    output.append("# CASP15 Ligand Binding Site Residues\n")
    output.append(f"Residues within {DISTANCE_CUTOFF}Å of ligand atoms.\n")
    output.append("=" * 100 + "\n\n")

    pdb_files = sorted([f for f in os.listdir(PDB_DIR) if f.endswith('.pdb')])

    for pdb_file in pdb_files:
        pdb_code = pdb_file.replace('.pdb', '')
        if pdb_code not in TARGETS:
            continue

        filepath = os.path.join(PDB_DIR, pdb_file)
        target_info = TARGETS[pdb_code]

        print(f"Processing {pdb_code} ({target_info['target']})...")

        protein_atoms, all_ligands = parse_pdb_atoms(filepath)

        output.append(f"## {target_info['target']} - {pdb_code}\n")
        output.append(f"**Protein**: {target_info['protein']}\n\n")

        # Focus on key ligands
        key_ligands = target_info['key_ligands']
        if not key_ligands:
            output.append("*Protein-protein complex - no small molecule ligands*\n\n")
            output.append("-" * 100 + "\n\n")
            continue

        for lig_name in key_ligands:
            # Find all instances of this ligand
            matching_ligs = {k: v for k, v in all_ligands.items() if k.startswith(lig_name + '_')}

            if not matching_ligs:
                output.append(f"### {lig_name}\n")
                output.append("*Ligand not found in structure*\n\n")
                continue

            # Analyze first instance (or combine all for same binding site)
            for lig_id, lig_atoms in list(matching_ligs.items())[:1]:  # First instance
                output.append(f"### {lig_name} ({lig_id})\n\n")

                binding_res, contacts = find_binding_site(protein_atoms, lig_atoms)

                if not binding_res:
                    output.append("*No protein residues within cutoff distance*\n\n")
                    continue

                # Sort by chain and residue number
                sorted_res = sorted(binding_res, key=lambda x: (x[0], x[2]))

                # Group by chain
                by_chain = defaultdict(list)
                for chain, resname, resnum in sorted_res:
                    by_chain[chain].append((resname, resnum))

                output.append(f"**Binding site residues** ({len(sorted_res)} residues):\n\n")

                for chain in sorted(by_chain.keys()):
                    residues = by_chain[chain]
                    res_str = ', '.join([f"{AA_3TO1.get(rn, rn)}{rnum}" for rn, rnum in sorted(residues, key=lambda x: x[1])])
                    output.append(f"- **Chain {chain}**: {res_str}\n")

                # Detailed contacts
                output.append(f"\n**Key contacts** (closest atoms):\n\n")
                output.append("| Residue | Protein Atom | Ligand Atom | Distance (Å) |\n")
                output.append("|---------|--------------|-------------|-------------:|\n")

                # Get closest contact for each residue
                for res_id in sorted_res[:15]:  # Limit to top 15
                    chain, resname, resnum = res_id
                    res_contacts = contacts[res_id]
                    if res_contacts:
                        # Find closest contact
                        closest = min(res_contacts, key=lambda x: x[2])
                        lig_atom, prot_atom, dist = closest
                        res_fmt = f"{AA_3TO1.get(resname, resname)}{resnum}({chain})"
                        output.append(f"| {res_fmt} | {prot_atom} | {lig_atom} | {dist} |\n")

                output.append("\n")

        output.append("-" * 100 + "\n\n")

    # Write output
    with open(OUTPUT_FILE, 'w') as f:
        f.write(''.join(output))

    print(f"\nSaved to {OUTPUT_FILE}")

if __name__ == '__main__':
    main()
