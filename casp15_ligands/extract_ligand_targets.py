#!/usr/bin/env python3
"""Extract CASP15 protein-ligand targets and create mapping with PDB codes."""

import re
import csv
import json

# Read the CSV file
input_file = '/home/sness/data/vaults/casp/casp15_original_inputs/casp15_targetlist.csv'

ligand_targets = []

with open(input_file, 'r') as f:
    reader = csv.reader(f, delimiter=';')
    header = next(reader)

    for row in reader:
        if len(row) < 10:
            continue

        target_id = row[0].strip()
        target_type = row[1].strip()
        residues = row[2].strip()
        oligo_state = row[3].strip()
        description = row[9].strip() if len(row) > 9 else ''

        # Skip RNA targets (start with R)
        if target_id.startswith('R'):
            continue

        # Check if it's a ligand target
        if 'Ligand' not in target_type and 'ligand' not in description.lower():
            continue

        # Skip canceled targets without structures
        if 'Canceled' in description and 'no structure' in description.lower():
            continue
        if 'Canceled (no ligand)' in description:
            continue

        # Extract PDB code (4 character alphanumeric)
        pdb_match = re.search(r'\b([0-9][a-zA-Z0-9]{3})\b', description)
        pdb_code = pdb_match.group(1).upper() if pdb_match else None

        # Clean description - remove HTML tags
        clean_desc = re.sub(r'<[^>]+>', '', description).strip()
        # Remove multiple spaces
        clean_desc = re.sub(r'\s+', ' ', clean_desc)

        # Determine if auxiliary
        is_auxiliary = 'auxiliary' in description.lower() or 'Not a TS target' in description

        # Determine status
        status = 'active'
        if 'Canceled' in description:
            status = 'canceled'

        ligand_targets.append({
            'target_id': target_id,
            'pdb_code': pdb_code,
            'residues': residues,
            'oligo_state': oligo_state,
            'description': clean_desc,
            'is_auxiliary': is_auxiliary,
            'status': status
        })

# Print summary
print(f"Found {len(ligand_targets)} protein-ligand targets")
print(f"With PDB codes: {len([t for t in ligand_targets if t['pdb_code']])}")
print(f"Missing PDB codes: {len([t for t in ligand_targets if not t['pdb_code']])}")
print()

# Print table
print("Target ID\tPDB\tResidues\tAuxiliary\tDescription")
print("-" * 100)
for t in ligand_targets:
    aux = "aux" if t['is_auxiliary'] else ""
    pdb = t['pdb_code'] or "---"
    desc = t['description'][:60] + "..." if len(t['description']) > 60 else t['description']
    print(f"{t['target_id']}\t{pdb}\t{t['residues']}\t{aux}\t{desc}")

# Save to JSON
output_file = '/home/sness/data/vaults/casp/casp15_ligands/ligand_targets.json'
with open(output_file, 'w') as f:
    json.dump(ligand_targets, f, indent=2)
print(f"\nSaved to {output_file}")

# Save targets with PDB codes for downloading
pdb_targets = [t for t in ligand_targets if t['pdb_code']]
pdb_file = '/home/sness/data/vaults/casp/casp15_ligands/pdb_codes.txt'
with open(pdb_file, 'w') as f:
    for t in pdb_targets:
        f.write(f"{t['pdb_code']}\n")
print(f"PDB codes saved to {pdb_file}")
