# What `protenix prep` Actually Does

Deep dive into the phase we're currently running — the gap between
"I have a JSON with a sequence and a SMILES" and "the GPU can start
predicting". This is what `prep_local.sh` calls on each target.

## Big picture

Protenix-v2 (and v1, and AlphaFold3) don't predict structure from
sequence alone. They consume **evolutionary** and **structural** priors:

- **Evolutionary**: Multiple Sequence Alignments (MSAs) — many related
  protein sequences stacked up. Columns that co-vary encode contacts.
- **Structural**: Templates — known 3D structures of homologs. Provide
  a geometric starting point for the pair representation.

`protenix prep` produces both. Without it, inference still runs but
quality drops meaningfully, especially on hard targets.

```
  jsons/<T>.json                  (protein sequence + ligand SMILES)
         │
         ▼  protenix prep
         │
         │  Stage 1:  MSA search    (remote  — minutes, queues)
         │  Stage 2:  Template search (local — seconds, hmmer + kalign)
         │  Stage 3:  JSON rewrite  (local — instant)
         │
         ▼
  jsons_prepped/<T>.json           (protein + ligand + MSA paths + template path)
  jsons_prepped/<T>/msa/0/
      pairing.a3m                  (taxonomically-paired MSA)
      non_pairing.a3m              (unpaired MSA)
      hmmsearch.a3m                (template alignment)
```

---

## Stage 1 — MSA search (the slow part)

This is what you watch sit on `PENDING` for minutes.

### What's being computed

For the query protein, the server runs **MMseqs2** against two large
databases:
- **UniRef30** — a clustered version of UniRef (sequence universe)
- **ColabFoldDB** / **env_nr** — metagenomic environmental sequences

The result: up to thousands of hits aligned to your query. Two flavors
are returned:

| File | What it is | Why |
|---|---|---|
| `pairing.a3m` | Hits with **taxonomy tags** in the header | Used to pair chains of a complex by species for inter-chain co-evolution signal |
| `non_pairing.a3m` | All hits, no taxonomy filtering | Used for intra-chain evolution within a single chain |

Protenix-v2 consumes both: the MSAModule attends over `non_pairing.a3m`
for per-chain features, and the pairing logic uses `pairing.a3m` when
multiple chains need to be stacked.

### The network call

```
POST https://protenix-server.com/api/msa/ticket/msa     ← submit (non-pairing)
POST https://protenix-server.com/api/msa/ticket/pair    ← submit (pairing)
GET  https://protenix-server.com/api/msa/ticket/{ID}    ← poll status
GET  https://protenix-server.com/api/msa/result/download/{ID}  ← fetch tar.gz
```

The tar.gz unpacks to:
- `0.a3m` → the raw MMseqs2 alignment
- `uniref_tax.m8` → tab-separated UniRef hits with taxonomy
- `pdb70_220313_db.m8` → hits against a PDB-sequence database (for templates later)
- `msa.sh` → provenance script the server ran
- `pairing.a3m` (from the pair ticket)

The parser in `colab_request_parser.py` reshapes these into the final
`pairing.a3m` + `non_pairing.a3m` layout.

### The `SUBMIT → PENDING → RUNNING → COMPLETE` states

What you see in the progress bar:

| State | What the server is doing |
|---|---|
| `SUBMIT` | Receiving your request |
| `PENDING` | Queued — waiting for a worker |
| `RUNNING: N%` | Worker is searching — N% through the MMseqs2 pipeline |
| `COMPLETE` | Result ready for download |

### Why timings vary so wildly

What we're seeing (9 s vs 600 s vs 2000 s) is **all queue + cache
variance on the server side**. The actual search usually takes ~30 s
when it runs. The rest is:

- **Cache hits**: if your exact protein sequence was searched recently,
  the server returns the cached tar immediately (~9 s total round-trip).
  That's what L1001/L1002/L1004 etc. hit.
- **Cache miss + light queue**: a real search runs — a few minutes.
- **Cache miss + heavy queue**: you wait behind other users. L4008 at
  1953 s was mostly queue, not compute.

**Implication**: the 17 unique Mpro sequences in L4000 each trigger a
fresh search. The clone-by-sequence trick (`clone_prepped.py`) avoids
duplicated searches within the batch, but can't help unique sequences.

### Alternatives we discussed

- `--msa_server_mode colabfold` + `MMSEQS_SERVICE_HOST_URL=https://api.colabfold.com`
  → hits the Söding lab's public server instead. Sometimes less queued,
  different schema though — both flags must be set together.
- `--use_msa false` → skip MSA entirely. Fast, worse predictions.
- Self-host the MMseqs2 server → ~$900/mo infra, only pays off at
  hundreds of queries/month.

---

## Stage 2 — Template search (the fast part)

This runs locally. Tools: `hmmbuild`, `hmmsearch`, `kalign` (already
installed in the conda env).

### What it does

1. Build an HMM profile from the MSAs we just fetched:
   ```
   hmmbuild --informat stockholm --amino --hand profile.hmm (MSA input)
   ```
2. Scan the PDB sequence database for structural homologs:
   ```
   hmmsearch --noali --cpu 8 -E 100 --incE 100 \
       -A hits.sto profile.hmm pdb_seqres_2022_09_28.fasta
   ```
3. Realign hits to the query and write `hmmsearch.a3m`.

### The database

`pdb_seqres_2022_09_28.fasta` is a ~300 MB snapshot of all PDB chain
sequences as of Sept 2022. It auto-downloads to
`$PROTENIX_ROOT_DIR/search_database/` on first prep run.

Because the cutoff is 2022-09-28, any PDB deposited after that date is
invisible to template search. For CASP16 pharma ligands (Chymase /
Cathepsin G / Autotaxin / Mpro), there's dense pre-2022 coverage, so
this isn't a practical limitation.

### Why it's fast

All the work is on your local CPU. hmmbuild takes ~0.2 s, hmmsearch
takes ~7 s against a 300 MB DB. No queues, no network.

### What it produces

`hmmsearch.a3m` — an a3m-format alignment of template hits, e.g.:
```
>query
IIGGRESRPHSRPYMAYLQIQSP...
>pdb|1T32_A
--GGKESKQHSRPYMAFLLISQS...
>pdb|1T31_B
...
```

Protenix-v2's `TemplateEmbedder` reads this, extracts Cα-Cα distance
bins and residue identities from each hit's PDB structure, and fuses
them into the pair representation. Details are in
`docs/templates.md`.

---

## Stage 3 — JSON rewrite

The trivial step. Protenix writes two output files:

| File | Role |
|---|---|
| `jsons/<T>-update-msa.json` | Intermediate — has MSA paths but no template path. `prep_local.sh` deletes this. |
| `jsons/<T>-final-updated.json` | Final — has MSA paths **and** template path. `prep_local.sh` moves it to `jsons_prepped/<T>.json`. |

The final JSON's `proteinChain` block gains three keys:

```json
{
  "proteinChain": {
    "sequence": "IIGG...RSFKL",
    "count": 1,
    "pairedMsaPath":   "/.../jsons_prepped/<T>/msa/0/pairing.a3m",
    "unpairedMsaPath": "/.../jsons_prepped/<T>/msa/0/non_pairing.a3m",
    "templatesPath":   "/.../jsons_prepped/<T>/msa/0/hmmsearch.a3m"
  }
}
```

These paths are **absolute**, which is why the planned migration to
`/workspace/casp16_ligands_protenix/` needs either a filesystem at the
same prefix on both sides (user's `/workspace` mount handles this) or a
sed rewrite.

---

## What happens next (on the pod)

When we scp everything to the pod and run `protenix pred` in
`run_batch.sh`:

1. `run_batch.sh` sees `jsons_prepped/<T>.json` exists → skips the prep
   call entirely.
2. `protenix pred` reads the JSON, resolves the three `*Path` entries,
   loads the MSA / template files.
3. Data pipeline builds feature tensors:
   - MSA feature: [N_seq × N_residues × 49] encoded alignment
   - Template feature: [N_template × N_residues × N_residues × N_pair_feats]
     distance bins + identity
4. Those feed the Pairformer + Diffusion + Confidence heads on GPU.

Prep is 100% CPU/network; pred is GPU-bound. That's why it's worth
running prep locally on a machine that's not costing $0.86/hr —
the pod is idle during queue waits otherwise.

---

## Troubleshooting cheat-sheet

| Symptom | Cause | Fix |
|---|---|---|
| `PENDING` for >30 min | Server overloaded | Wait, or switch to colabfold endpoint (both `--msa_server_mode` and `MMSEQS_SERVICE_HOST_URL`) |
| `FileNotFoundError: bfd.mgnify30.metaeuk30.smag30.a3m` | Mode/host mismatch — colabfold parser against protenix server | Set `MMSEQS_SERVICE_HOST_URL=https://api.colabfold.com` too |
| `ImportError: cannot import name 'FastaBatchedDataset'` | `fair-esm` vs `esm 3.x` collision | Fresh conda env; don't mix the two packages |
| `hmmsearch: command not found` | Missing deps | `conda install -c bioconda hmmer` + `apt install kalign` |
| Prep reports success but `jsons_prepped/<T>.json` doesn't exist | Output name mismatch | Our wrapper handles this: it looks for `<T>-final-updated.json` in the source dir and moves it |
| Different targets, same sequence, both doing fresh searches | Cache not shared | `clone_prepped.py` — copies one prepped JSON to all siblings with matching protein sequence |
