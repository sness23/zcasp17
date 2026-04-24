# Structural Templates in AlphaFold 3 / Protenix

A cheat-sheet on what "templates" are, how Protenix consumes them, and what
actually happens when you pass `--use_template true`.

## The idea in one sentence

A *template* is a known 3D structure of a homolog of your query sequence,
injected as structural prior. MSAs tell the model "these residues coevolve
→ probably contact"; templates tell it "in a related protein, these
residues were arranged *like so*".

## Why they help ligand docking specifically

The CASP16 pharma targets (Chymase, Cathepsin G, Autotaxin, Mpro) are all
well-characterized drug-discovery receptors with many deposited
protein-ligand complexes in the PDB. Turning on templates means the model
starts from a structurally-informed prior for the receptor scaffold and
can spend its capacity on the ligand pose rather than re-deriving the
fold.

Protenix-v1 and -v2 both support templates; chai-lab does not consume
them explicitly in this pipeline. This is one of the concrete knobs we
expect Protenix to win on.

## How the pipeline builds them

`protenix prep` runs a three-stage search, all against the PDB sequence
database (`pdb_seqres_2022_09_28.fasta`), using HMMER:

```
pairing.a3m + non_pairing.a3m          # already built by the MSA step
          │
          ▼  hmmbuild
      profile.hmm                      # HMM profile of the query family
          │
          ▼  hmmsearch  -o hits        # search pdb_seqres with the profile
      hits with PDB-chain identifiers
          │
          ▼  kalign                    # realign hits back onto the query
      hmmsearch.a3m                    # merged per-chain template alignment
```

`hmmsearch.a3m` is what the Protenix `TemplateEmbedder` actually consumes.
It is referenced from the per-chain `templatesPath` field of the input
JSON — `protenix prep` writes that path in for you.

Alternative format: `.hhr` files (HH-suite output) are also accepted.
See `examples/examples_with_template/example_mgyp004658859411.json` in
the Protenix repo for the `.hhr` shape.

## What the model does with them

Inside `protenix/model/modules/pairformer.py::TemplateEmbedder`:

1. Each template contributes a pair representation — distance bins between
   template Cα atoms plus amino-acid identity features.
2. A small stack of pairformer-style blocks (2 blocks in `protenix-v2`,
   controlled by `model.template_embedder.n_blocks`) refines this.
3. The refined pair features are **added** to the main pair
   representation before the full pairformer runs.

Net effect: the pair stack starts from a biased state that reflects the
template geometry, rather than a zero initialization.

## Required dependencies

Templates need two external binaries on the `PATH`:

| Binary | Package | Purpose |
|---|---|---|
| `hmmbuild` / `hmmsearch` | `hmmer` | Build & run HMM profile search |
| `kalign` | `kalign` | Realign hits to the query |

On Debian/Ubuntu pods: `apt-get install -y hmmer kalign`. These are
already in the official Protenix Docker image.

You can also pass explicit binary paths to `protenix prep` if the tools
live outside `PATH`:

```bash
protenix prep \
    --hmmbuild_binary_path  /opt/hmmer/bin/hmmbuild \
    --hmmsearch_binary_path /opt/hmmer/bin/hmmsearch \
    --kalign_binary_path    /opt/kalign/kalign \
    --input jsons/L1001.json --out_dir jsons/
```

## Database

The template database (`pdb_seqres_2022_09_28.fasta`) auto-downloads from
ByteDance TOS into `$PROTENIX_ROOT_DIR/search_database/` the first time
`protenix prep` (or `protenix mt`) runs. It's a few hundred MB, not the
giant wwPDB bundle.

**Training-data cutoff implications**: this FASTA is a 2022-09-28 snapshot.
Templates deposited after that date won't be searchable. For CASP16
pharma ligands (Chymase/Cathepsin G/Autotaxin/Mpro), this is a
non-issue — all four families have ample pre-2022 PDB coverage.

## Turning it off

Pass `--use_template false` (or omit `templatesPath` from the JSON). The
`TemplateEmbedder` is then bypassed and the pair rep is initialized from
MSA + pair embeddings only. Variants without template support (v0.5.0,
mini, tiny) require this.

## How to inspect what matched

After `protenix prep`, open `jsons/<TARGET>.prep.json` and look at the
`templatesPath` the tool wrote in. The `.a3m` file at that path has a
header per hit of the form:

```
>pdb|<PDBID>_<CHAIN>
<aligned sequence, gaps and insertions encoded a3m-style>
```

Counting the hits (`grep -c '^>' hmmsearch.a3m`) gives a quick sanity
check — a well-studied target like Mpro should return dozens to hundreds.
