#!/usr/bin/env python3
"""Download PDB structures and fetch publication DOIs for CASP15 ligand targets."""

import csv
import json
import os
import urllib.request
import urllib.error
import time

OUTPUT_DIR = '/home/sness/data/vaults/casp/casp15_ligands'
PDB_DIR = os.path.join(OUTPUT_DIR, 'pdb')
os.makedirs(PDB_DIR, exist_ok=True)

# Read the mapping file
mapping_file = os.path.join(OUTPUT_DIR, 'ligand_targets_mapping.csv')
targets = []

with open(mapping_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        targets.append(row)

print(f"Found {len(targets)} ligand targets")
print(f"With PDB codes: {len([t for t in targets if t['pdb_code']])}")
print()

# Download PDB files and get metadata
results = []

for target in targets:
    pdb_code = target['pdb_code'].strip().lower() if target['pdb_code'] else None
    target_id = target['target_id']

    if not pdb_code:
        print(f"[{target_id}] No PDB code - skipping")
        results.append({
            **target,
            'pdb_downloaded': False,
            'doi': None,
            'title': None,
            'authors': None,
            'journal': None,
            'year': None
        })
        continue

    # Download PDB file
    pdb_url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    pdb_file = os.path.join(PDB_DIR, f"{pdb_code.upper()}.pdb")

    pdb_downloaded = False
    if not os.path.exists(pdb_file):
        try:
            print(f"[{target_id}] Downloading {pdb_code.upper()}...", end=" ")
            urllib.request.urlretrieve(pdb_url, pdb_file)
            print("OK")
            pdb_downloaded = True
        except urllib.error.HTTPError as e:
            print(f"FAILED ({e.code})")
    else:
        print(f"[{target_id}] {pdb_code.upper()} already exists")
        pdb_downloaded = True

    # Get publication info from RCSB API
    doi = None
    title = None
    authors = None
    journal = None
    year = None

    try:
        api_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_code}"
        with urllib.request.urlopen(api_url) as response:
            data = json.loads(response.read().decode())

            # Get citation info
            if 'rcsb_primary_citation' in data:
                citation = data['rcsb_primary_citation']
                doi = citation.get('pdbx_database_id_doi')
                title = citation.get('title')
                journal = citation.get('journal_abbrev')
                year = citation.get('year')

                # Get authors
                if 'rcsb_authors' in citation:
                    authors = '; '.join(citation['rcsb_authors'][:5])
                    if len(citation['rcsb_authors']) > 5:
                        authors += ' et al.'

            # If no primary citation, try to get entry info
            if not title and 'struct' in data:
                title = data['struct'].get('title')

    except Exception as e:
        print(f"  [API Error] {e}")

    results.append({
        **target,
        'pdb_downloaded': pdb_downloaded,
        'doi': doi,
        'title': title,
        'authors': authors,
        'journal': journal,
        'year': year
    })

    time.sleep(0.2)  # Be nice to RCSB

# Save results
output_file = os.path.join(OUTPUT_DIR, 'ligand_targets_with_dois.csv')
fieldnames = ['target_id', 'pdb_code', 'residues', 'oligo_state', 'protein_name',
              'ligand_info', 'is_auxiliary', 'status', 'notes',
              'pdb_downloaded', 'doi', 'title', 'authors', 'journal', 'year']

with open(output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for r in results:
        writer.writerow(r)

print(f"\nSaved to {output_file}")

# Print summary table
print("\n" + "="*120)
print("CASP15 Ligand Targets Summary")
print("="*120)
print(f"{'Target':<10} {'PDB':<6} {'DOI':<40} {'Title':<50}")
print("-"*120)

for r in results:
    target_id = r['target_id']
    pdb = r['pdb_code'] if r['pdb_code'] else '---'
    doi = r['doi'] if r['doi'] else 'No DOI'
    title = (r['title'][:47] + '...') if r['title'] and len(r['title']) > 50 else (r['title'] or 'No title')
    print(f"{target_id:<10} {pdb:<6} {doi:<40} {title}")

# Also save as JSON for easy programmatic access
json_file = os.path.join(OUTPUT_DIR, 'ligand_targets_with_dois.json')
with open(json_file, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nJSON saved to {json_file}")

# Count statistics
with_doi = len([r for r in results if r['doi']])
with_pdb = len([r for r in results if r['pdb_downloaded']])
print(f"\nStatistics:")
print(f"  Targets with PDB downloaded: {with_pdb}/{len(results)}")
print(f"  Targets with DOI: {with_doi}/{len(results)}")
