#!/usr/bin/env python3
"""Analyze PDB structures and extract ligand binding information for comparison with papers."""

import os
import re
from collections import defaultdict

PDB_DIR = '/home/sness/data/vaults/casp/casp15_ligands/pdb'
OUTPUT_FILE = '/home/sness/data/vaults/casp/casp15_ligands/pdb_ligand_analysis.md'

# Common solvent/buffer molecules to exclude
SOLVENTS = {'HOH', 'WAT', 'DOD', 'H2O', 'SOL'}
BUFFERS = {'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'MPD', 'DMS', 'ACT', 'IMD', 'CIT',
           'TRS', 'MES', 'EPE', 'BOG', 'BME', 'DTT', 'CL', 'NA', 'K', 'MG', 'CA',
           'ZN', 'MN', 'FE', 'CU', 'NI', 'CO'}

# Target to PDB mapping with expected ligands from papers
TARGETS = {
    '7UTD': {
        'target': 'H1114',
        'protein': '[NiFe]-hydrogenase Huc',
        'expected_ligands': ['NiFe cluster', '[3Fe-4S] clusters', 'menaquinone'],
        'paper_binding': 'Narrow hydrophobic gas channels for H2; [3Fe-4S] clusters modulate properties; menaquinone transported 94Å'
    },
    '8B43': {
        'target': 'T1118v1',
        'protein': 'FoxA TonB-dependent transporter',
        'expected_ligands': ['siderophores', 'ferrioxamine', 'thiocillin'],
        'paper_binding': 'Outer membrane transporter; binds siderophores and antibiotics'
    },
    '7UX8': {
        'target': 'T1124',
        'protein': 'MfnG tyrosine O-methyltransferase',
        'expected_ligands': ['SAM', 'SAH', 'tyrosine'],
        'paper_binding': 'SAM cofactor binding; tyrosine substrate pocket'
    },
    '7Z8Y': {
        'target': 'H1135',
        'protein': 'SUN1-KASH6 LINC complex',
        'expected_ligands': ['KASH peptide'],
        'paper_binding': 'Asymmetric 9:6 complex; two KASH conformations'
    },
    '8SX8': {
        'target': 'T1158v1',
        'protein': 'MRP4 transporter + PGE1',
        'expected_ligands': ['PGE1', 'ATP'],
        'paper_binding': 'Prostaglandin E1 in substrate binding cavity'
    },
    '8SXB': {
        'target': 'T1158v2',
        'protein': 'MRP4 transporter + PGE2',
        'expected_ligands': ['PGE2', 'ATP'],
        'paper_binding': 'Prostaglandin E2 in substrate binding cavity'
    },
    '8SX7': {
        'target': 'T1158v3',
        'protein': 'MRP4 transporter + DHEAS',
        'expected_ligands': ['DHEAS', 'ATP'],
        'paper_binding': 'DHEA-sulfate in substrate binding cavity'
    },
    '8SWN': {
        'target': 'T1158v4',
        'protein': 'MRP4 transporter + ATP',
        'expected_ligands': ['ATP', 'Mg'],
        'paper_binding': 'ATP-Mg2+ bound; outward-occluded conformation'
    },
    '7PBR': {
        'target': 'T1170',
        'protein': 'RuvB hexamer',
        'expected_ligands': ['ATP', 'ADP', 'DNA'],
        'paper_binding': 'AAA+ ATPase; nucleotide binding; DNA substrate'
    },
    '7PBL': {
        'target': 'H1171',
        'protein': 'RuvB-RuvA complex',
        'expected_ligands': ['ATP', 'ADP', 'DNA'],
        'paper_binding': 'RuvB with RuvA domain 3; Holliday junction processing'
    },
    '7PBP': {
        'target': 'H1172',
        'protein': 'RuvB-2xRuvA complex',
        'expected_ligands': ['ATP', 'ADP', 'DNA'],
        'paper_binding': 'Full RuvAB complex with Holliday junction'
    },
    '8AD2': {
        'target': 'T1187',
        'protein': 'Nictaba lectin',
        'expected_ligands': ['chitotriose', 'GlcNAc'],
        'paper_binding': 'Jelly-roll fold; chitotriose binding centered on central GlcNAc'
    },
    '8C6Z': {
        'target': 'T1188',
        'protein': 'Clostridium chitinase ChiB',
        'expected_ligands': ['chitin', 'GlcNAc', 'Bisdionin C'],
        'paper_binding': 'GH18 chitinase; prefers linear chitin substrates'
    },
    '8XBP': {
        'target': 'T1127',
        'protein': 'AtNATA1',
        'expected_ligands': ['Acetyl-CoA'],
        'paper_binding': 'N-acetyltransferase bound to Acetyl-CoA'
    }
}

def parse_pdb(filepath):
    """Parse PDB file and extract structure information."""
    info = {
        'title': '',
        'resolution': None,
        'method': '',
        'chains': set(),
        'ligands': defaultdict(list),  # ligand_name -> list of residue info
        'metals': [],
        'modified_residues': [],
        'hetatm_count': 0,
        'atom_count': 0,
        'sequence_length': defaultdict(int),
    }

    seen_residues = set()

    with open(filepath, 'r') as f:
        for line in f:
            # Title
            if line.startswith('TITLE'):
                info['title'] += line[10:].strip() + ' '

            # Resolution
            if line.startswith('REMARK   2 RESOLUTION'):
                match = re.search(r'(\d+\.\d+)', line)
                if match:
                    info['resolution'] = float(match.group(1))

            # Method
            if line.startswith('EXPDTA'):
                info['method'] = line[10:].strip()

            # Atoms
            if line.startswith('ATOM'):
                info['atom_count'] += 1
                chain = line[21]
                info['chains'].add(chain)
                resname = line[17:20].strip()
                resnum = line[22:26].strip()
                info['sequence_length'][chain] = max(info['sequence_length'][chain], int(resnum) if resnum.isdigit() else 0)

            # Heteroatoms (ligands, metals, etc.)
            if line.startswith('HETATM'):
                info['hetatm_count'] += 1
                resname = line[17:20].strip()
                chain = line[21]
                resnum = line[22:26].strip()

                res_id = f"{resname}_{chain}_{resnum}"
                if res_id not in seen_residues:
                    seen_residues.add(res_id)

                    if resname in SOLVENTS:
                        continue
                    elif resname in BUFFERS:
                        if resname in ['ZN', 'MN', 'FE', 'CU', 'NI', 'CO', 'MG', 'CA']:
                            info['metals'].append({'name': resname, 'chain': chain, 'resnum': resnum})
                    else:
                        info['ligands'][resname].append({'chain': chain, 'resnum': resnum})

            # Modified residues
            if line.startswith('MODRES'):
                info['modified_residues'].append(line[12:27].strip())

    info['title'] = info['title'].strip()
    return info

def main():
    output = []
    output.append("# CASP15 Ligand Targets - PDB Structure Analysis\n")
    output.append("Comparison of paper descriptions with actual PDB ligand content.\n")
    output.append("=" * 100 + "\n\n")

    # Get all PDB files
    pdb_files = sorted([f for f in os.listdir(PDB_DIR) if f.endswith('.pdb')])

    for pdb_file in pdb_files:
        pdb_code = pdb_file.replace('.pdb', '')
        filepath = os.path.join(PDB_DIR, pdb_file)

        if pdb_code not in TARGETS:
            continue

        target_info = TARGETS[pdb_code]
        print(f"Analyzing {pdb_code} ({target_info['target']})...")

        pdb_info = parse_pdb(filepath)

        output.append(f"## {target_info['target']} - {pdb_code}\n")
        output.append(f"**Protein**: {target_info['protein']}\n\n")

        # Structure info
        output.append("### Structure Information\n")
        output.append(f"- **Title**: {pdb_info['title']}\n")
        output.append(f"- **Method**: {pdb_info['method']}\n")
        if pdb_info['resolution']:
            output.append(f"- **Resolution**: {pdb_info['resolution']} Å\n")
        output.append(f"- **Chains**: {', '.join(sorted(pdb_info['chains']))}\n")
        output.append(f"- **Atoms**: {pdb_info['atom_count']:,} protein, {pdb_info['hetatm_count']:,} heteroatoms\n\n")

        # Ligands found
        output.append("### Ligands Found in PDB\n")
        if pdb_info['ligands']:
            for lig_name, instances in sorted(pdb_info['ligands'].items()):
                output.append(f"- **{lig_name}**: {len(instances)} instance(s)")
                if len(instances) <= 5:
                    chains = [f"{i['chain']}{i['resnum']}" for i in instances]
                    output.append(f" at {', '.join(chains)}")
                output.append("\n")
        else:
            output.append("- No non-solvent ligands found\n")

        # Metals
        if pdb_info['metals']:
            output.append("\n### Metal Ions\n")
            metal_counts = defaultdict(int)
            for m in pdb_info['metals']:
                metal_counts[m['name']] += 1
            for metal, count in sorted(metal_counts.items()):
                output.append(f"- **{metal}**: {count} ion(s)\n")

        # Comparison with paper
        output.append("\n### Comparison with Paper Description\n")
        output.append(f"**Expected from paper**: {', '.join(target_info['expected_ligands'])}\n\n")
        output.append(f"**Paper binding description**: {target_info['paper_binding']}\n\n")

        # Analysis
        found_ligands = set(pdb_info['ligands'].keys())
        expected_keywords = [e.lower() for e in target_info['expected_ligands']]

        matches = []
        missing = []
        unexpected = []

        for exp in target_info['expected_ligands']:
            exp_lower = exp.lower()
            found = False
            for lig in found_ligands:
                lig_lower = lig.lower()
                if exp_lower in lig_lower or lig_lower in exp_lower:
                    matches.append(f"{exp} → {lig}")
                    found = True
                    break
                # Check common abbreviations
                if exp_lower == 'atp' and lig in ['ATP', 'ANP', 'ADP', 'AMP']:
                    matches.append(f"{exp} → {lig}")
                    found = True
                    break
                if 'glcnac' in exp_lower and 'NAG' in lig:
                    matches.append(f"{exp} → {lig}")
                    found = True
                    break
            if not found:
                missing.append(exp)

        if matches:
            output.append("**✓ Matched ligands**:\n")
            for m in matches:
                output.append(f"  - {m}\n")

        if missing:
            output.append("\n**? Expected but not found** (may be different name or not co-crystallized):\n")
            for m in missing:
                output.append(f"  - {m}\n")

        if found_ligands - {m.split(' → ')[1] for m in matches if ' → ' in m}:
            other = found_ligands - {m.split(' → ')[1] for m in matches if ' → ' in m}
            output.append(f"\n**Additional ligands in structure**: {', '.join(sorted(other))}\n")

        output.append("\n" + "-" * 100 + "\n\n")

    # Write output
    with open(OUTPUT_FILE, 'w') as f:
        f.write('\n'.join(output))

    print(f"\nSaved to {OUTPUT_FILE}")

if __name__ == '__main__':
    main()
