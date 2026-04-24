#!/usr/bin/env python3
"""Extract protein structure and ligand binding information from CASP15 papers."""

import subprocess
import re
import os

PAPERS_DIR = '/home/sness/data/vaults/casp/casp15_ligands/papers'
OUTPUT_FILE = '/home/sness/data/vaults/casp/casp15_ligands/binding_descriptions.md'

# Paper metadata
papers = [
    {
        'file': 'H1114_7UTD_nature_hydrogenase.pdf',
        'target': 'H1114',
        'pdb': '7UTD',
        'protein': '[NiFe]-hydrogenase Huc',
        'ligand': 'NiFe cofactor, menaquinone',
        'keywords': ['NiFe', 'active site', 'catalytic', 'cofactor', 'hydrogen', 'binding', 'menaquinone', 'H2', 'cluster', 'channel']
    },
    {
        'file': 'T1118v1_8B43_pnas_foxA.pdf',
        'target': 'T1118v1',
        'pdb': '8B43',
        'protein': 'TonB-dependent transporter FoxA',
        'ligand': 'siderophores, antibiotics',
        'keywords': ['binding', 'siderophore', 'ferrioxamine', 'antibiotic', 'transport', 'ligand', 'pocket', 'site', 'interaction']
    },
    {
        'file': 'T1124_7UX8_protsci_mfnG.pdf',
        'target': 'T1124',
        'pdb': '7UX8',
        'protein': 'MfnG tyrosine O-methyltransferase',
        'ligand': 'SAM/SAH (S-adenosylmethionine)',
        'keywords': ['SAM', 'SAH', 'methyl', 'binding', 'active site', 'catalytic', 'substrate', 'tyrosine', 'cofactor']
    },
    {
        'file': 'H1135_7Z8Y_commbiol_sun1kash6.pdf',
        'target': 'H1135',
        'pdb': '7Z8Y',
        'protein': 'SUN1-KASH6 LINC complex',
        'ligand': 'KASH peptide',
        'keywords': ['binding', 'KASH', 'SUN', 'interface', 'interaction', 'complex', 'groove', 'domain']
    },
    {
        'file': 'T1158_8SX_nsmb_mrp4.pdf',
        'target': 'T1158',
        'pdb': '8SX8/8SXB/8SX7/8SWN',
        'protein': 'MRP4 (ABCC4) transporter',
        'ligand': 'prostaglandins (PGE1, PGE2), DHEAS, ATP',
        'keywords': ['binding', 'prostaglandin', 'PGE', 'ATP', 'substrate', 'pocket', 'cavity', 'transport', 'efflux']
    },
    {
        'file': 'T1170_7PBR_nature_ruvab.pdf',
        'target': 'T1170',
        'pdb': '7PBR/7PBL/7PBP',
        'protein': 'RuvAB-Holliday junction complex',
        'ligand': 'ATP/ADP, DNA',
        'keywords': ['ATP', 'ADP', 'binding', 'nucleotide', 'active site', 'catalytic', 'DNA', 'Holliday', 'junction']
    },
    {
        'file': 'T1187_8AD2_glycobiol_nictaba.pdf',
        'target': 'T1187',
        'pdb': '8AD2',
        'protein': 'Nictaba lectin',
        'ligand': 'carbohydrate/GlcNAc',
        'keywords': ['binding', 'carbohydrate', 'sugar', 'GlcNAc', 'lectin', 'site', 'interaction', 'glycan']
    },
    {
        'file': 'T1188_8C6Z_plos_chitinases.pdf',
        'target': 'T1188',
        'pdb': '8C6Z',
        'protein': 'Clostridium perfringens chitinase',
        'ligand': 'chitin/GlcNAc oligomers',
        'keywords': ['binding', 'chitin', 'substrate', 'active site', 'catalytic', 'GlcNAc', 'cleavage', 'glycoside']
    }
]

def extract_text(pdf_path):
    """Extract text from PDF using pdftotext."""
    try:
        result = subprocess.run(
            ['pdftotext', '-layout', pdf_path, '-'],
            capture_output=True, text=True, timeout=60
        )
        return result.stdout
    except Exception as e:
        return f"Error extracting text: {e}"

def find_relevant_sections(text, keywords, context_lines=10):
    """Find paragraphs containing keywords."""
    lines = text.split('\n')
    relevant = []
    seen = set()

    for i, line in enumerate(lines):
        line_lower = line.lower()
        for kw in keywords:
            if kw.lower() in line_lower:
                # Get context
                start = max(0, i - 2)
                end = min(len(lines), i + context_lines)
                chunk = '\n'.join(lines[start:end])

                # Deduplicate by checking first 100 chars
                chunk_key = chunk[:100]
                if chunk_key not in seen and len(chunk.strip()) > 50:
                    seen.add(chunk_key)
                    relevant.append(chunk)
                break

    return relevant

def extract_abstract(text):
    """Try to extract the abstract."""
    # Common abstract patterns
    patterns = [
        r'Abstract[:\s]*\n(.*?)(?:Introduction|Keywords|Background|\n\n\n)',
        r'ABSTRACT[:\s]*\n(.*?)(?:INTRODUCTION|KEYWORDS|\n\n\n)',
    ]
    for pattern in patterns:
        match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
        if match:
            abstract = match.group(1).strip()
            # Clean up
            abstract = re.sub(r'\s+', ' ', abstract)
            if len(abstract) > 100:
                return abstract[:2000]
    return None

def main():
    output = []
    output.append("# CASP15 Ligand Targets - Protein & Binding Site Descriptions\n")
    output.append("Extracted from original publications for comparison with PDB structures.\n")
    output.append("=" * 80 + "\n\n")

    for paper in papers:
        pdf_path = os.path.join(PAPERS_DIR, paper['file'])
        if not os.path.exists(pdf_path):
            output.append(f"## {paper['target']} - {paper['protein']}\n")
            output.append(f"**FILE NOT FOUND**: {paper['file']}\n\n")
            continue

        print(f"Processing {paper['target']}...")
        text = extract_text(pdf_path)

        output.append(f"## {paper['target']} ({paper['pdb']}) - {paper['protein']}\n")
        output.append(f"**Ligand(s)**: {paper['ligand']}\n\n")

        # Extract abstract
        abstract = extract_abstract(text)
        if abstract:
            output.append("### Abstract\n")
            output.append(f"{abstract}\n\n")

        # Extract binding-related sections
        output.append("### Key Binding/Structure Information\n")
        relevant = find_relevant_sections(text, paper['keywords'])

        if relevant:
            for i, section in enumerate(relevant[:15]):  # Limit to 15 most relevant
                # Clean up the text
                section = re.sub(r'\s+', ' ', section).strip()
                if len(section) > 50:
                    output.append(f"- {section[:500]}{'...' if len(section) > 500 else ''}\n\n")
        else:
            output.append("No specific binding information extracted.\n\n")

        output.append("-" * 80 + "\n\n")

    # Write output
    with open(OUTPUT_FILE, 'w') as f:
        f.write('\n'.join(output))

    print(f"\nSaved to {OUTPUT_FILE}")

if __name__ == '__main__':
    main()
