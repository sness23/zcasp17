# Conventions

Naming, structural, and metadata conventions used across the canonical
inputs tree. Enforced by `build_canonical.py` and relied on by all
adapters in `adapters/`.

## Target IDs

- **CASP15**: `T<4-digit>` for protein monomers/complexes (T1124),
  `H<4-digit>` for hard multi-chain (H1135), `R<4-digit>` for
  RNA-containing targets (R1126). We preserve the official CASP
  naming — no prefixes or transformations.
- **CASP16**: `L<set><3-digit>` for pharma-ligand targets
  (L1001–L1017 for Chymase, L2001/L2002 for Cathepsin G,
  L3001–L3231 for Autotaxin, L4001–L4028 for Mpro — with some gaps).
  Main-track CASP16 targets (D1273, R1261, T1214, etc.) are not
  included here because the ligand pipeline only handles pharma track.
- **CASP17**: TBD once targets are released.

## Chain IDs inside `target.json`

Synthetic chain IDs are used for deduped entities to avoid ambiguity:

- **Single chain per target:** `id` = the source PDB's chain ID (e.g. `"A"`).
- **Multi-chain target:** `id` = `chain<SRC_ID>_<LEN>aa` (e.g.
  `"chainA_306aa"`). The length suffix keeps near-identical chains
  distinguishable even if their single-letter IDs are the same.
- **True homodimer** (identical sequences): a single entity with
  `count: 2`, `id` = same synthetic style. The `source_chains` field
  lists both original chain IDs.

## Homodimers and near-homodimers

Per [schema.md](schema.md), chains are deduped by **exact sequence
identity**. Rationale:

1. **Exact homodimer** (both chains = 226 aa, identical sequence) → one
   entity with `count: 2`. This is what Protenix's native dimer
   notation handles best.
2. **Near-homodimer** (226 aa vs 227 aa, one residue difference) → two
   separate entities with `count: 1` each. Forcing these into a
   `count: 2` would lose information and misrepresent the structure.
3. **Concatenation is never allowed.** A previous bug in our CASP16
   builder concatenated near-homodimer chains into a single long
   sequence; this caused catastrophic ligand misplacement (LDDT_pli
   0.00, RMSD 50+ Å). See
   [`../../casp16_ligands_protenix/docs/03_discussion.md §1`](../../casp16_ligands_protenix/docs/03_discussion.md#1-the-homodimer-bug)
   for the full story.

## Ligand IDs

- **Stable internal ID:** `lig_NN` where NN is a zero-padded index
  within the target (`lig_01`, `lig_02`, ...).
- **Groupings:** ligands with the same CCD are merged into one entity
  with `count > 1`. Novel pharma SMILES (no CCD) get their own entity
  each (`count: 1` usually).
- **SDF filenames:** Original SDF filenames from the source refs are
  preserved, not renamed. A ligand entity's `sdf_paths` lists them
  in the order they came from disk.

## Receptor PDB

- Contains **protein + RNA + DNA only** — all HETATM records stripped.
  This distinguishes the reference receptor from the "combined"
  PDB in some CASP bundles that has ligand HETATMs in the same file.
- Used for **scoring only** (as the `-r` argument to OST's
  `compare-ligand-structures`). Predictors never consume this file;
  they get sequence-only inputs.

## SMILES resolution

- **CASP15:** ligands are CCD codes (3-letter PDB chemical component
  codes). SMILES is fetched from RCSB's Chemical Component Dictionary
  via `https://data.rcsb.org/rest/v1/core/chemcomp/<CCD>`. Results
  cached locally in `.smiles_cache.json` at the `inputs/` root.
- **CASP16:** pharma-track ligands have no CCD assignment; the SMILES
  comes from the CASP16 target TSV (`casp16_ligands/<SET>/<T>.tsv`,
  column `SMILES`).
- **CASP17:** TBD once targets release.

## Sequence conventions

- **Protein:** one-letter amino acid code, uppercase, no gaps.
  Non-standard residues are dropped (mirrors
  `gemmi.find_tabulated_residue(r.name).one_letter_code.isalpha()`).
- **RNA:** A/U/G/C, uppercase, no gaps.
- **DNA:** A/T/G/C, uppercase, no gaps. The `D` prefix in `DA/DT/DG/DC`
  residue names is stripped at extraction time.

## Canonical tree is authoritative

- Pipelines should never modify the canonical tree; they read it via
  adapters.
- The canonical tree itself can be regenerated from source data by
  re-running `build_canonical.py --all`. This will overwrite existing
  target directories but preserves the structure.
- If edits are needed (e.g. fixing a bad ligand SMILES), edit the
  source data AND re-run `build_canonical.py`, not the canonical tree
  directly. This keeps derivations reproducible.
