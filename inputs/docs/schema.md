# `target.json` schema

Canonical representation of a single CASP target.

## Top-level fields

```json
{
  "target_id": "string, CASP target identifier (e.g. 'T1124', 'L4001')",
  "casp_year": "int, 15 / 16 / 17",
  "set": "string or null — CASP16 pharma-ligand set id ('L1000', 'L2000', 'L3000', 'L4000') or null",
  "chains": [...],
  "ligands": [...],
  "receptor_pdb": "string, relative path to crystal receptor PDB",
  "notes": "string, freeform — anything not captured above"
}
```

## `chains`

An array of chain entities. Each entity has exactly one sequence, with
`count` ≥ 1 if that sequence is represented multiple times:

```json
{
  "id": "string, chain letter or synthetic ID ('A', 'B', 'chain1')",
  "type": "one of 'protein', 'rna', 'dna'",
  "sequence": "string, one-letter code",
  "count": "int, number of copies with this exact sequence",
  "source_chains": ["A", "B"]   // optional — original chain IDs from source PDB
}
```

**Dedup rule:** chains with identical sequences are merged into a single
entity with `count > 1`. Chains that differ even by a single residue
(truncation variants, sequence polymorphisms) are kept as separate
entities with `count: 1` each. This is critical for multi-chain
receptors — see [conventions.md](conventions.md#homodimers-and-near-homodimers)
for why.

Example: L2001 Cathepsin G has chain A (226 aa) and chain B (227 aa,
one residue longer). Emit two entities, not one concatenated 453-mer
and not a `count: 2` homodimer.

## `ligands`

An array of ligand entities:

```json
{
  "id": "string, stable internal ID like 'lig_01'",
  "ccd": "string or null, 3-letter PDB CCD code if available (e.g. 'SAH')",
  "smiles": "string, canonical SMILES from CCD lookup or source TSV",
  "count": "int, number of crystallographic instances",
  "sdf_paths": ["ligands/lig_01_SAH.sdf", "ligands/lig_02_SAH.sdf"]
}
```

**CCD vs SMILES:**
- CASP15 ligands are almost always CCD (SAH, ATP, NAG, EPE, etc.).
- CASP16 pharma ligands are novel drugs with no CCD assignment —
  only `smiles` is populated, `ccd` is null.
- When both are present, SMILES takes priority in adapters.

**Multi-copy ligands:** a crystallographic ligand present twice
(e.g. T1124 with 2 × SAH) yields a single ligand entity with
`count: 2` and two SDF paths. The SDFs are separate files.

## `receptor_pdb`

Relative path to a clean protein-only PDB (no HETATM ligands). Derived
from the CASP source PDB by stripping ligand atoms. Used for ground-
truth scoring, not fed into predictors (predictors get the sequence,
not the crystal coords).

## Examples

### CASP15 mono-ligand

```json
{
  "target_id": "T1146",
  "casp_year": 15,
  "set": null,
  "chains": [
    {"id": "A", "type": "protein", "sequence": "MKDQR...", "count": 1}
  ],
  "ligands": [
    {
      "id": "lig_01",
      "ccd": "NAG",
      "smiles": "CC(=O)N[C@@H]1[C@H]([C@@H]([C@@H](O[C@@H]1O)CO)O)O",
      "count": 1,
      "sdf_paths": ["ligands/lig_01_NAG.sdf"]
    }
  ],
  "receptor_pdb": "receptor.pdb",
  "notes": ""
}
```

### CASP15 multi-chain with multi-copy ligand

```json
{
  "target_id": "T1124",
  "casp_year": 15,
  "set": null,
  "chains": [
    {"id": "chainA_379aa", "type": "protein", "sequence": "NVSL...YFQY", "count": 1},
    {"id": "chainB_362aa", "type": "protein", "sequence": "NVSL...RATY", "count": 1}
  ],
  "ligands": [
    {
      "id": "lig_01",
      "ccd": "SAH",
      "smiles": "c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CSCC[C@@H](C(=O)O)N)O)O)N",
      "count": 2,
      "sdf_paths": ["ligands/lig_01_SAH.sdf", "ligands/lig_02_SAH.sdf"]
    }
  ],
  "receptor_pdb": "receptor.pdb",
  "notes": "Methyltransferase-like fold; two SAH cofactor copies, one per subunit"
}
```

### CASP16 pharma (homodimer, novel SMILES)

```json
{
  "target_id": "L4001",
  "casp_year": 16,
  "set": "L4000",
  "chains": [
    {"id": "chainA_306aa", "type": "protein", "sequence": "SGF...Q", "count": 1},
    {"id": "chainB_304aa", "type": "protein", "sequence": "SGF...T", "count": 1}
  ],
  "ligands": [
    {
      "id": "lig_01",
      "ccd": null,
      "smiles": "O=C(Cn1ccnn1)Nc1ccc(Br)c([N+](=O)[O-])c1",
      "count": 1,
      "sdf_paths": ["ligands/lig_01_ligand_LIG_A_401.sdf"]
    }
  ],
  "receptor_pdb": "receptor.pdb",
  "notes": "Mpro drug candidate; homodimer with 2aa-different chains"
}
```

## Validation

A valid `target.json` passes these checks (run automatically by
`build_canonical.py`):

1. `target_id` is non-empty.
2. `casp_year` in {15, 16, 17}.
3. `chains` has ≥ 1 entity, all with valid one-letter sequences.
4. `ligands` has ≥ 0 entities (0 is legal for protein-only CASP targets,
   though those aren't the focus here).
5. Every `sdf_paths[*]` file exists on disk.
6. `receptor_pdb` file exists on disk and parses as PDB.
7. For each ligand entity, `len(sdf_paths) == count`.
8. `ccd` is 3 characters uppercase, or null.

Validation failures abort the build and print a target-specific error.
