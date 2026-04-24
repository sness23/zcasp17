# CASP15 Ligand Binding Site ChimeraX Scripts

Visualization scripts for CASP15 protein-ligand competition targets.

## Usage

Open ChimeraX and run a script:
```
open /home/sness/data/vaults/casp/casp15_ligands/chimerax/SCRIPT_NAME.cxc
```

Or from command line:
```bash
chimerax --script /home/sness/data/vaults/casp/casp15_ligands/chimerax/SCRIPT_NAME.cxc
```

## Scripts

| Script | Target | PDB | Protein | Ligand(s) |
|--------|--------|-----|---------|-----------|
| `7PBL_H1171_ruvab.cxc` | H1171 | 7PBL | RuvB-RuvA complex | ADP, AGS |
| `7PBP_H1172_ruvab2.cxc` | H1172 | 7PBP | RuvB-2xRuvA complex | ADP, AGS |
| `7PBR_T1170_ruvb.cxc` | T1170 | 7PBR | RuvB hexamer | ADP, AGS |
| `7UTD_H1114_hydrogenase.cxc` | H1114 | 7UTD | [NiFe]-hydrogenase Huc | F3S, MQ9, 3NI, FCO |
| `7UX8_T1124_mfng.cxc` | T1124 | 7UX8 | MfnG methyltransferase | SAH, TYR |
| `8AD2_T1187_nictaba.cxc` | T1187 | 8AD2 | Nictaba lectin | NAG |
| `8B43_T1118_foxa.cxc` | T1118v1 | 8B43 | FoxA transporter | OX8 (ferrioxamine) |
| `8C6Z_T1188_chitinase.cxc` | T1188 | 8C6Z | Chitinase ChiB | DW0 (Bisdionin C) |
| `8SWN_T1158v4_mrp4_atp.cxc` | T1158v4 | 8SWN | MRP4 + ATP | ATP |
| `8SX7_T1158v3_mrp4_dheas.cxc` | T1158v3 | 8SX7 | MRP4 + DHEAS | ZWY |
| `8SX8_T1158v1_mrp4_pge1.cxc` | T1158v1 | 8SX8 | MRP4 + PGE1 | XPG |
| `8SXB_T1158v2_mrp4_pge2.cxc` | T1158v2 | 8SXB | MRP4 + PGE2 | P2E |
| `8XBP_T1127_atnata1.cxc` | T1127 | 8XBP | AtNATA1 | ACO (Acetyl-CoA) |
| `load_all.cxc` | All | All | Load all structures | — |

## Color Scheme

- **Gold**: Primary binding site residues
- **Cornflower blue**: Secondary binding site / cofactor site
- **Magenta**: Ligands (ball style)
- **Red**: Key catalytic residues
- **Cyan**: Key polar contact residues
- **Green**: Metal ions (Mg, Ca)
- **Orange**: Fe ions

## Features

Each script:
1. Opens the PDB structure from local files
2. Shows protein as semi-transparent cartoon
3. Highlights binding site residues as sticks
4. Shows ligand in ball-and-stick representation
5. Adds transparent surface around binding pocket
6. Labels ligands
7. Centers view on the binding site
8. Saves a session file (.cxs)

## Binding Site Cutoff

Residues within **4.5 Å** of any ligand atom are shown as binding site residues.
