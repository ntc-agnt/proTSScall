# proTSScall_optimized.R
#
# Optimized, hardened variant of proTSScall.R. The original is kept untouched;
# this is the parallel file (TSScall-style). Intended to eventually replace the
# original once validated superior.
#
# Changes vs proTSScall.R, by category:
#
# BYTE-IDENTICAL (same outputs on valid input; proven by
#   test/equivalence_dominant_dedup.R, 20000 randomized trials, 0 mismatches):
#   * Lever A -- the two per-row R loops (dominant-TSS-per-gene and co-located
#     dedup) are replaced by vectorized group operations. The originals index the
#     data.frame by character row name every iteration (O(active x transcripts));
#     the vectorized form is a single grouped pass. Tie-breaks are reproduced
#     exactly: dominant selection keeps the FIRST transcript on a full tie (reads
#     then strand-aware upstream), dedup keeps the LAST on a full tie (V10 then
#     ENSG-digit; the original's which.min ties to the candidate). The
#     write.table subsetting expressions are kept verbatim so row ORDER -- hence
#     the bytes -- is identical.
#
# BEHAVIOR CHANGES (only on BAD input / re-runs; identical on a valid fresh run):
#   * NTC-001 robustness: validate args + input list; check exit codes on every
#     external step (merges, make_heatmap, awk, mv) so a failed tool stops the run
#     instead of silently producing empty/wrong output; `mkdir -p` and `ln -rsf`
#     for resume-safety; assert the make_heatmap matrix has the TSS+150 activity
#     columns before scoring.
#   * NTC-015 chr-prefix guard: before make_heatmap, assert the TSS-list chroms
#     and the merged bedGraph chroms share naming; a mismatch (chr1 vs 1) would
#     otherwise yield all-zero coverage and call everything inactive, exit 0.
#
# NOT changed here (deferred): NTC-003 config-file portability -- the invoked-path
#   -> scriptsPath detection block below is kept as-is for now.
#
# Byte-identity is defined over the five data outputs (the matrix + 4
# formakeheatmap files). `complete.txt` carries a Sys.time() stamp and is
# intentionally not reproducible.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Run an external command and stop (non-zero exit) with the failing command, so
# a broken merge / make_heatmap can't pass silently into downstream steps.
run <- function(cmd) {
  rc <- system(cmd)
  if (rc != 0L) stop(sprintf("external command failed (exit %d): %s", rc, cmd))
}

# First run of digits in an ID, as numeric (NA if none) -- vectorized form of the
# original's per-pair regexpr/substr ENSG-number extraction.
ensg_digits <- function(ids) {
  m <- regexpr("[[:digit:]]+", ids)
  as.numeric(substr(ids, as.integer(m), as.integer(m) + attr(m, "match.length") - 1))
}

# One dominant TSS per gene. Byte-identical to proTSScall.R's dominant loop:
# per gene (V7), the winner has the most TSS-proximal reads; ties go to the more
# upstream TSS (smaller V4 on +, larger V5 on -); a full tie keeps the first TSS
# in active (matrix-row) order. Returns the dominant transcript IDs in
# gene-first-appearance order (== unlist(domTss)).
pick_dominant <- function(tssActive, tssReads, tssList) {
  if (length(tssActive) == 0) return(character(0))
  reads  <- tssReads[tssActive]
  gene   <- tssList[tssActive, "V7"]
  strand <- tssList[tssActive, "V6"]
  if (any(strand != "+" & strand != "-")) stop("Unexpected strand symbol")
  v4 <- tssList[tssActive, "V4"]
  v5 <- tssList[tssActive, "V5"]
  # strand-aware upstream key: ascending == more upstream on both strands.
  upstream <- ifelse(strand == "+", v4, -v5)
  idx <- seq_along(tssActive)
  ord <- order(gene, -reads, upstream, idx)          # max reads, most upstream, first-seen
  firstPick <- !duplicated(gene[ord])
  win_tx   <- tssActive[ord][firstPick]
  win_gene <- gene[ord][firstPick]
  gene_first <- tapply(idx, gene, min)               # gene-insertion order
  win_tx[order(gene_first[win_gene])]
}

# Deduplicate co-located dominant TSSs. Byte-identical to proTSScall.R's dedup
# loop: among dominants sharing a start key (V3.start.strand), the winner has the
# longest transcript (V10); ties go to the lower ENSG number; a full tie keeps the
# LAST in dominant-list order (the original replaces on a tie). Returns the
# surviving transcript IDs (a set; the dominant file is written in tssList order
# via %in%).
dedup_dominant <- function(domList, tssList) {
  if (length(domList) == 0) return(character(0))
  strand <- tssList[domList, "V6"]
  if (any(strand != "+" & strand != "-")) stop("Unexpected strand symbol")
  coord <- ifelse(strand == "+", tssList[domList, "V4"], tssList[domList, "V5"])
  key <- paste(tssList[domList, "V3"], coord, strand, sep = ".")
  len <- tssList[domList, "V10"]
  dig <- ensg_digits(tssList[domList, "V7"])
  idx <- seq_along(domList)
  ord <- order(key, -len, dig, -idx)                 # max V10, min ENSG-digit, last-seen
  domList[ord][!duplicated(key[ord])]
}

# ---------------------------------------------------------------------------
# Arguments + input list
# ---------------------------------------------------------------------------
cArgs <- commandArgs()
fileArg <- which(startsWith(cArgs, "--file="))
argsArg <- which(startsWith(cArgs, "--args"))
runningScript <- normalizePath(sapply(strsplit(cArgs[fileArg], "--file="), "[", 2))

if (length(argsArg) == 0L || length(cArgs) < argsArg + 2L) {
  stop("usage: Rscript proTSScall_optimized.R <tss.tes_list> <read_threshold>")
}
tssPath <- cArgs[argsArg + 1]
readThreshold <- as.numeric(cArgs[argsArg + 2])
if (!is.finite(readThreshold)) {
  stop(sprintf("read_threshold must be a finite number; got '%s'", cArgs[argsArg + 2]))
}
if (!file.exists(tssPath)) {
  stop(sprintf("TSS/TES list not found: %s", tssPath))
}

tssList <- read.table(tssPath, sep = "\t", row.names = 1, stringsAsFactors = F)
# logic needs columns V3-V7 and V10 (file columns 3-7, 10 -> >= 9 data columns
# after the row-name column is consumed).
if (ncol(tssList) < 9L) {
  stop(sprintf("TSS/TES list needs >= 10 tab columns (V3-V10 are used); %s has %d data columns",
               tssPath, ncol(tssList)))
}

# ---------------------------------------------------------------------------
# scriptsPath (NTC-003 will replace this detection with a config file)
# ---------------------------------------------------------------------------
if (startsWith(runningScript, "/n/data1/cores/ntc")) {
  scriptsPath <- "/n/data1/cores/ntc/scripts"
} else if (startsWith(runningScript, "/n/data1/hms/bcmp/adelman")) {
  scriptsPath <- "/n/data1/hms/bcmp/adelman/Scripts/ntc"
} else {
  message("Warning: invoked script not on 'ntc' or 'adelman' paths")
  scriptsPath <- "/n/data1/cores/ntc/scripts"
}

# ---------------------------------------------------------------------------
# Merge per-sample 3' bedGraphs into one forward / one reverse track
# ---------------------------------------------------------------------------
if (!dir.exists("bedGraphs")) {
  stop("expected a 'bedGraphs/' directory in the working dir -- run from the pipeline outdir")
}
setwd("bedGraphs")
run("mkdir -p forward reverse")
if (length(Sys.glob("*_3pr_forward.bedGraph")) == 0L ||
    length(Sys.glob("*_3pr_reverse.bedGraph")) == 0L) {
  stop("no *_3pr_{forward,reverse}.bedGraph in bedGraphs/ -- wrong outdir, or tracks not named *_3pr_*")
}
run("ln -rsf *_3pr_forward.bedGraph forward")
run("ln -rsf *_3pr_reverse.bedGraph reverse")
final_pro_tss_dir <- 'pro_tss'

outPrefix <- sapply(strsplit(getwd(), "/"), tail, 2)[1]
setwd("forward")
run(paste(paste(scriptsPath, "AdelmanLab/NIH_scripts/bedgraphs2stdBedGraph/bedgraphs2stdBedGraph", sep = "/"), paste(outPrefix, "F", sep = "_")))
setwd("../reverse")
run(paste(paste(scriptsPath, "AdelmanLab/NIH_scripts/bedgraphs2stdBedGraph/bedgraphs2stdBedGraph", sep = "/"), paste(outPrefix, "R", sep = "_")))
setwd("../../")
tssPrefix <- sapply(strsplit(sapply(strsplit(tssPath, "/"), "[", length(strsplit(tssPath, "/")[[1]])), "_formakeheatmap.txt"), "[", 1)
run(sprintf("mkdir -p %s", final_pro_tss_dir))

# a temporary directory which will prevent incomplete files from populating the final `pro_tss` directory:
tmp_pro_tss_dir <- 'pro_tss_tmp'
run(sprintf('mkdir -p %s', tmp_pro_tss_dir))

# ---------------------------------------------------------------------------
# NTC-015: chromosome-naming guard (before make_heatmap)
# A chr-prefix mismatch between the TSS list and the bedGraphs makes make_heatmap
# score all-zero coverage -> everything "inactive", exit 0. Fail loudly instead.
# ---------------------------------------------------------------------------
fwdMerged <- file.path("bedGraphs", "forward", paste0(outPrefix, "_F.bedGraph"))
bgChroms <- system(sprintf("awk '$1 !~ /^track/ {print $1}' %s | sort -u", shQuote(fwdMerged)),
                   intern = TRUE)
listChroms <- unique(as.character(tssList$V3))
if (length(bgChroms) > 0L && length(intersect(listChroms, bgChroms)) == 0L) {
  stop(sprintf(paste0("chromosome-naming mismatch (NTC-015): TSS-list chroms e.g. {%s} share no ",
                      "name with bedGraph chroms e.g. {%s}; make_heatmap would score all-zero ",
                      "coverage and call every TSS inactive."),
               paste(head(listChroms, 3), collapse = ", "),
               paste(head(bgChroms, 3), collapse = ", ")))
}

# ---------------------------------------------------------------------------
# make_heatmap matrix (+/- 2 kb, 25-bp bins)
# ---------------------------------------------------------------------------
run(paste(scriptsPath, "/AdelmanLab/NIH_scripts/make_heatmap/make_heatmap -t 6 -l s -s s --nohead -p bedGraphs/forward/", outPrefix, "_F.bedGraph -m bedGraphs/reverse/", outPrefix, "_R.bedGraph -- ", tssPath, " ", final_pro_tss_dir, "/", outPrefix, "_", tssPrefix, "_25mer_+-2kb.txt -2000 25 160", sep = ""))

tssMat <- read.table(paste(final_pro_tss_dir, "/", outPrefix, "_", tssPrefix, "_25mer_+-2kb.txt", sep = ""), header = T, sep = "\t", row.names = 1)
# the activity score sums TSS..+150 (bins X0.24 .. X125.149); fail if absent.
if (!all(c("X0.24", "X125.149") %in% colnames(tssMat))) {
  stop(sprintf("make_heatmap matrix is missing the TSS+150 activity columns X0.24..X125.149; got columns starting: %s",
               paste(head(colnames(tssMat)), collapse = ", ")))
}

# get TSS.+150 counts
tssReads <- rowSums(tssMat[, which(colnames(tssMat) == "X0.24"):which(colnames(tssMat) == "X125.149")])
tssInactive <- names(tssReads)[tssReads <= readThreshold]

# print inactive TSSs (into tmp dir)
write.table(tssList[tssInactive, ], paste(tmp_pro_tss_dir, "/", outPrefix, "_", tssPrefix, "_inactive_formakeheatmap.txt", sep = ""), quote = F, sep = "\t", col.names = F)
tssActive <- names(tssReads)[tssReads > readThreshold]

# identify 1 dominant TSS per active gene (vectorized; byte-identical to the loop)
domList <- pick_dominant(tssActive, tssReads, tssList)

tssNonDom <- tssActive[!(tssActive %in% domList)]
# print active, non-dominant TSSs (into tmp dir)
write.table(tssList[tssNonDom, ], paste(tmp_pro_tss_dir, "/", outPrefix, "_", tssPrefix, "_non-dominant_formakeheatmap.txt", sep = ""), quote = F, sep = "\t", col.names = F)

# deduplicate co-located dominant TSSs (vectorized; byte-identical to the loop)
uniqueTss <- dedup_dominant(domList, tssList)

# print dominant TSSs to tmp
write.table(tssList[rownames(tssList) %in% uniqueTss, ], paste0(tmp_pro_tss_dir, "/", outPrefix, "_", tssPrefix, "_dominant_formakeheatmap.txt"), quote = F, sep = "\t", col.names = F)

# print dominant gene bodies (reading the file in the tmp dir and directing stdout to that same tmp dir)
run(paste("awk -F'\t' -v OFS='\t' '{if($10>=400) {if($6==\"+\") print $1, $2, $3, $4+250, $5, $6, $7, $8, $9, $10; else print $1, $2, $3, $4, $5-250, $6, $7, $8, $9, $10}}'",
          paste0(tmp_pro_tss_dir, "/", outPrefix, "_", tssPrefix, "_dominant_formakeheatmap.txt"), ">",
          paste0(tmp_pro_tss_dir, "/", outPrefix, "_", tssPrefix, "_+250.TES.min400_formakeheatmap.txt")))

# If we have reached here, we can move the files to their final destination:
run(sprintf('mv %s/* %s/', tmp_pro_tss_dir, final_pro_tss_dir))

# attempt to remove the tmp directory. Note that this only works if the directory is empty.
# If it's not empty, then we catch that and error
exit_code <- system(sprintf('rmdir %s', tmp_pro_tss_dir), intern = FALSE)
if (exit_code != 0) {
  message('Could not clean the tmp directory. Files were still present. Exit.')
  quit(status = 1)
}

# Create a job completion file with a timestamp.
cat(format(Sys.time(), "%h %d, %Y:%H:%M:%S"), file = paste0(final_pro_tss_dir, "/complete.txt"))
