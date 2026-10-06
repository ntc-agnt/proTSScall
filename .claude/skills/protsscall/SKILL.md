---
name: protsscall
description: >
  Use when picking the single dominant, active TSS per gene from PRO-seq 3' (1-nt)
  bedGraphs plus a TSS/TES gene list — calling active vs inactive TSSs by TSS-proximal
  read counts (TSS to +150), choosing one dominant TSS per gene, and emitting the
  dominant/non-dominant/inactive/gene-body formakeheatmap files — or when a task in a
  related NTC repo touches these. Upstream inputs come from pipeline2 PRO-seq
  (3' bedGraphs) and a formakeheatmap TSS/TES list (GeneAnnotationScripts /
  Genome_Onboarding); downstream `_dominant_formakeheatmap.txt` feeds
  metagene_heatmap_tool (`-proTSS` arm) and pauseIndex. NOT a fork of TSScall and does
  not run it — naming similarity only. Aliases: proTSScall, proTSScall.R.
---

# proTSScall

## Purpose
proTSScall answers a single question: given PRO-seq transcription activity, **which one
TSS is the dominant active start site for each gene?** From per-sample 1-nt 3' PRO-seq
bedGraphs and a candidate TSS/TES list, it measures TSS-proximal read activity (TSS to
+150 bp), splits TSSs into active/inactive, then selects exactly one dominant TSS per gene
and deduplicates co-located calls. Repo: `proTSScall` (single script `proTSScall.R`);
deployed under `/n/data1/cores/ntc/scripts` or the Adelman tree. Maintainer: NTC /
sgoldman101.

## Architecture
**Standalone R script — NOT a smartSlurm/rAP pipeline** (no `runAsPipeline`, no `#@`
annotations, no `ssbatch`, no flag folder, no resume). `proTSScall.R` runs top-to-bottom
in one process, orchestrating two external NIH_scripts binaries via `system()`. It is
run inside an sbatch wrapper (README: `-p short -t 0-3 --mem=50G -c 6`,
`module load gcc/6.2.0 R/4.0.1`). The `-c 6` matters: `make_heatmap` is invoked with
`-t 6`.

Pipeline inside the script:
1. In `bedGraphs/`, make `forward/`+`reverse/`, symlink (`ln -rs`) every
   `*_3pr_{forward,reverse}.bedGraph`, merge each strand via `bedgraphs2stdBedGraph` →
   `<outPrefix>_{F,R}.bedGraph`.
2. `make_heatmap` over merged tracks, ±2 kb around each TSS in 25-bp bins (160 bins) →
   `<P>_<tssPrefix>_25mer_+-2kb.txt` matrix.
3. Sum bins `X0.24`..`X125.149` (= TSS to +150) per TSS → activity score; `<=threshold`
   inactive, `>threshold` active.
4. Pick one dominant TSS per gene (by `V7` gene ID), dedup co-located calls, write outputs
   to `pro_tss_tmp/` then `mv` to `pro_tss/`, drop `complete.txt`.

Full algorithm (tie-breaks, dedup, bin math, column layout, output table):
`reference/algorithm-and-io.md`.

## I/O contract
**Invocation** (CWD **must** be the pipeline outdir containing `bedGraphs/`):
```
Rscript /path/to/proTSScall.R /full/path/to/<genome>.<annot>_formakeheatmap.txt 4
```
`arg1` = full path to the TSS/TES list (formakeheatmap format); `arg2` = numeric
read_threshold (README uses 4).

**Inputs:**
- `<outdir>/bedGraphs/*_3pr_{forward,reverse}.bedGraph` — per-sample 1-nt 3' PRO-seq
  tracks (from pipeline2); **all are merged together** per strand.
- TSS/TES list: tab-separated, **no header**, read with `row.names=1`. Load-bearing
  columns: col1 = transcript ID (row name), `V3` chr, `V4` start, `V5` end, `V6` strand,
  `V7` gene ID, `V10` transcript length.

**Outputs** (in `<outdir>/pro_tss/`, all keep the input tab layout, no header):
`*_25mer_+-2kb.txt` (matrix), `*_inactive_formakeheatmap.txt`,
`*_non-dominant_formakeheatmap.txt`, **`*_dominant_formakeheatmap.txt`** (the key file —
one TSS per gene), `*_+250.TES.min400_formakeheatmap.txt` (gene bodies), `complete.txt`.
Full schema: `reference/algorithm-and-io.md`.

## Dependencies & environment
- **R 4.0.1** under `gcc/6.2.0` (README), base R only — no CRAN/Bioc packages used. The
  `ntc_basic` conda env (`ntc_env` skill) is the ecosystem standard but is **not stated in
  this repo** (TODO verify which env is actually used).
- External binaries from **NIH_scripts** (AdelmanLab), via **hardcoded** `/n/data1` paths
  resolved from the invoked script's own path:
  `.../AdelmanLab/NIH_scripts/bedgraphs2stdBedGraph/bedgraphs2stdBedGraph` and
  `.../make_heatmap/make_heatmap`.
- `scriptsPath` is chosen by where `proTSScall.R` itself lives: `/n/data1/cores/ntc` →
  `ntc/scripts`; `/n/data1/hms/bcmp/adelman` → `adelman/Scripts/ntc`; else warns and
  defaults to the ntc tree.
- Reference data: the TSS/TES formakeheatmap list (genome/annotation specific).

## Modifying it safely (gotchas & assumptions)
- **Hardcoded O2 paths** — must be run from the ntc or adelman script tree, or NIH_scripts
  won't resolve. Off-cluster it cannot run unmodified.
- **Not resume-safe.** `mkdir pro_tss` has no `-p` (fails if it exists), and `system()`
  exit codes for the heatmap/merge steps are ignored — a partial run leaves broken state.
  The `pro_tss_tmp/` → `mv` → `pro_tss/` staging + `complete.txt` marker (added 2024) is
  the only completeness guard; re-run from a clean outdir.
- **Symlinks are recreated** in `forward/`/`reverse/` each run (`ln -rs`); stale links from
  a prior run can linger.
- **Rigid column names.** Hard dependence on `V3`–`V7` and `V10` positions — any TSS list
  not in exact formakeheatmap column order silently miscomputes.
- **ENSG tie-break is human-specific.** Final dedup tie is broken by the lowest numeric
  part of the gene ID (ENSG number); assumes Ensembl-style IDs.
- **`_+250.TES.min400` naming is misleading.** Despite "TES", it only shifts the *start*
  inward by 250 bp (on `-` strand, the end by −250) for transcripts ≥400 bp — it does not
  recompute a TES.

## Integration points
- **Upstream:** pipeline2 PRO-seq (`bedGraphs/*_3pr_{forward,reverse}.bedGraph`; a commit
  explicitly made it "compatible with pipeline2"); NIH_scripts (`bedgraphs2stdBedGraph`,
  `make_heatmap`); the TSS/TES formakeheatmap list from GeneAnnotationScripts /
  Genome_Onboarding (TODO confirm exact producer).
- **Downstream:** `*_dominant_formakeheatmap.txt` feeds metagene_heatmap_tool (`-proTSS`
  arm) and pauseIndex.
- **Shared formats:** the `formakeheatmap` gene-list format in and out.
- **Not related to TSScall** beyond the name — proTSScall never invokes TSScall and is not
  a fork of it (see `gene-annotations` / `tsscall` skills for the actual TSScall path).

## Reference files
- `reference/algorithm-and-io.md` — full activity/selection algorithm, tie-breaks, dedup,
  bin math, input column layout, output file table, evolution & gotchas.
- `.claude/reference/runAsPipeline-grammar.md` — smartSlurm grammar (proTSScall is *not*
  rAP, but the surrounding ecosystem is; context only).
