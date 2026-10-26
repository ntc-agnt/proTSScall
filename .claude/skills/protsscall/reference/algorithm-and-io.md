# proTSScall — algorithm, I/O schema, gotchas

`proTSScall.R` is a single-process R script. It takes two positional args parsed from
`commandArgs()`:

- `arg1` (`tssPath`) — full path to the TSS/TES formakeheatmap list.
- `arg2` (`readThreshold`) — numeric TSS-proximal read cutoff (README uses `4`).

It also inspects its own invoked path (`--file=`) to pick `scriptsPath`
(`/n/data1/cores/ntc/scripts` for the ntc tree, `/n/data1/hms/bcmp/adelman/Scripts/ntc`
for the adelman tree, else warn + default to ntc). CWD must be the pipeline outdir that
contains `bedGraphs/`.

## Full algorithm

### 1. Merge bedGraphs
- `setwd("bedGraphs")`; `mkdir -p forward reverse`.
- Symlink all `*_3pr_forward.bedGraph` into `forward/` and all `*_3pr_reverse.bedGraph`
  into `reverse/` (`ln -rs`, relative symlinks).
- `outPrefix` = the name of the outdir (second-to-last path element of the pre-`setwd`
  working dir).
- In each strand dir run `bedgraphs2stdBedGraph <outPrefix>_{F|R}` → a single merged
  `<outPrefix>_F.bedGraph` / `<outPrefix>_R.bedGraph`. **All per-sample tracks are merged
  into one** — proTSScall operates on the pooled signal, not per sample.

### 2. make_heatmap matrix (±2 kb, 25-bp bins)
- `tssPrefix` = the TSS list basename with `_formakeheatmap.txt` stripped.
- `mkdir pro_tss` (**no `-p`**) and `mkdir pro_tss_tmp`.
- Invocation (cores = 6):
  ```
  make_heatmap -t 6 -l s -s s --nohead \
    -p bedGraphs/forward/<P>_F.bedGraph -m bedGraphs/reverse/<P>_R.bedGraph \
    -- <tssPath> pro_tss/<P>_<tssPrefix>_25mer_+-2kb.txt -2000 25 160
  ```
  `-2000 25 160` = window start −2000, 25-bp bins, 160 bins → spans −2000..+2000 (±2 kb).
  `-l s -s s` select sense length/strand; `-p` = plus/forward track, `-m` = minus/reverse.
- The matrix is read back with `header=T, row.names=1`; row names are the transcript IDs,
  columns are bins named like `X<start>.<end>` (e.g. `X0.24`, `X25.49`, …).

### 3. Activity score and active/inactive split
- Per-TSS activity = `rowSums` over bin columns **`X0.24` through `X125.149`** — i.e. the
  6 bins covering **TSS to +150 bp** (positions 0–149, 6×25 bp).
- `tssReads <= readThreshold` → **inactive**; `tssReads > readThreshold` → **active**.
- Inactive rows written to `<...>_inactive_formakeheatmap.txt` (into `pro_tss_tmp/`).

### 4. One dominant TSS per gene
Iterate active TSSs; group by gene ID (`V7`). For each gene keep the current best TSS and
replace it when the candidate is better:
1. **More TSS-proximal reads** wins (`tssReads[i]` greater).
2. **Read tie → more upstream (strand-aware):** on `+` strand the smaller `V4` (start)
   wins; on `-` strand the larger `V5` (end) wins. Any other strand symbol → `stop()`.

Active TSSs not selected as any gene's dominant → `<...>_non-dominant_formakeheatmap.txt`.

### 5. Deduplicate co-located dominant TSSs
Two genes' dominant calls can share a start site. Key each dominant TSS by its start
coordinate + strand (`V3.V4.+` on `+`, `V3.V5.-` on `-`). When two dominants collide:
1. **Longer transcript wins** (larger `V10`).
2. **Length tie → lower ENSG number:** extract the first run of digits from the gene ID
   (`V7`) of each and keep the smaller numeric value. **Human/Ensembl-specific.**
The surviving unique dominants → `<...>_dominant_formakeheatmap.txt`.

### 6. Gene-body file + completion
- From the dominant file, keep transcripts with `V10 >= 400` and shift the start inward by
  250 bp (`+` strand: `$4+250`; `-` strand: `$5-250`), via an `awk` one-liner →
  `<...>_+250.TES.min400_formakeheatmap.txt`. Despite "TES" in the name it only moves the
  start; it does not compute a true TES.
- `mv pro_tss_tmp/* pro_tss/`, then `rmdir pro_tss_tmp` (fails loudly and `quit(status=1)`
  if the tmp dir is non-empty).
- Write `pro_tss/complete.txt` with a timestamp as the success marker.

## Input TSS/TES list — column layout
Tab-separated, **no header**, read with `row.names=1` (so col1 is the row name, not `V1`):

| Position | Name in script | Meaning |
|--|--|--|
| col 1 | (row name) | transcript ID (e.g. ENST) |
| col 3 | `V3` | chromosome |
| col 4 | `V4` | start coordinate |
| col 5 | `V5` | end coordinate |
| col 6 | `V6` | strand (`+` / `-`) |
| col 7 | `V7` | gene ID (e.g. ENSG) — grouping key |
| col 10 | `V10` | transcript length |

Columns V2, V8, V9 are carried through but not used by the logic. File naming is typically
`<genome>.<annot>_formakeheatmap.txt` (e.g. `hg38.basic_formakeheatmap.txt`).

## Output files (`pro_tss/`)
All keep the input's tab layout with no header (written via `write.table(..., quote=F,
sep="\t", col.names=F)`). `<P>` = outPrefix, `<tssPrefix>` = TSS-list prefix.

| File | Contents | Consumer |
|--|--|--|
| `<P>_<tssPrefix>_25mer_+-2kb.txt` | make_heatmap matrix, ±2 kb / 25-bp bins | internal (activity scoring) |
| `<P>_<tssPrefix>_inactive_formakeheatmap.txt` | TSSs with activity ≤ threshold | QC / audit |
| `<P>_<tssPrefix>_non-dominant_formakeheatmap.txt` | active but not dominant TSSs | QC / audit |
| **`<P>_<tssPrefix>_dominant_formakeheatmap.txt`** | **one dominant TSS per gene (key output)** | metagene_heatmap_tool `-proTSS`, pauseIndex |
| `<P>_<tssPrefix>_+250.TES.min400_formakeheatmap.txt` | dominant gene bodies, start shifted +250, `V10≥400` | gene-body / metagene analyses |
| `complete.txt` | timestamp success marker | run-completeness check |

## Evolution, rationale & gotchas
- Core script written 2021–22 (gmnelson). In 2024 Brian Lawney added the `pro_tss_tmp/`
  staging dir and the `complete.txt` marker — the only safeguards against half-written
  outputs populating `pro_tss/`. Merged in 2025-02 (sgoldman101); otherwise stable.
- A commit explicitly made it "compatible with pipeline2" — the `*_3pr_{forward,reverse}`
  naming convention is that interface; renaming PRO-seq tracks breaks the symlink glob.
- **Not resume-safe:** `mkdir pro_tss` lacks `-p` (re-running in a used outdir errors), and
  the heatmap/merge `system()` calls ignore exit codes, so an external-tool failure can
  pass silently until downstream rows are empty. Clean the outdir and re-run.
- **Hardcoded `/n/data1` paths**, resolved from the script's own location — it only works
  from the ntc or adelman tree; elsewhere it warns and defaults to the ntc path, which will
  then fail to find the binaries.
- **Rigid positional columns** (`V3`–`V7`, `V10`) assume exact formakeheatmap order; a
  reordered list miscomputes without error.
- The **ENSG numeric tie-break** and the **`+250.TES` misnomer** are the two most
  surprising behaviors to anyone reading output — see steps 5–6 above.
- Despite the name, proTSScall has **no relationship to TSScall**: it does not call,
  vendor, or fork it. The shared "TSS call" wording is coincidental.
