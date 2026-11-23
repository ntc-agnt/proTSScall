#!/usr/bin/env bash
# downsample_chr21.sh -- chr21-only copies of 3' PRO-seq bedGraphs for in-sandbox
# proTSScall testing. Preserves the *_3pr_{forward,reverse}.bedGraph names that
# proTSScall globs on, and keeps any leading `track` header line.
#
# Usage: downsample_chr21.sh <src_outdir> <dest_dir> [chrom]
#   <src_outdir> : a pipeline2 outdir that contains bedGraphs/
#   <dest_dir>   : destination; a bedGraphs/ subdir is created under it
#   [chrom]      : chromosome to keep (default chr21)
set -euo pipefail
SRC=${1:?usage: downsample_chr21.sh <src_outdir> <dest_dir> [chrom]}
DEST=${2:?usage: downsample_chr21.sh <src_outdir> <dest_dir> [chrom]}
CHR=${3:-chr21}

[[ -d "$SRC/bedGraphs" ]] || { echo "ERROR: $SRC/bedGraphs not found" >&2; exit 1; }
mkdir -p "$DEST/bedGraphs"

shopt -s nullglob
n=0
for f in "$SRC"/bedGraphs/*_3pr_forward.bedGraph "$SRC"/bedGraphs/*_3pr_reverse.bedGraph; do
  awk -v c="$CHR" '/^track/ || $1==c' "$f" > "$DEST/bedGraphs/$(basename "$f")"
  n=$((n+1))
done
[[ $n -gt 0 ]] || { echo "ERROR: no *_3pr_{forward,reverse}.bedGraph in $SRC/bedGraphs" >&2; exit 1; }
echo "wrote $n $CHR bedGraph(s) to $DEST/bedGraphs/"
ls -l "$DEST/bedGraphs/"
