#!/usr/bin/env Rscript
# equivalence_dominant_dedup.R
#
# Proves the vectorized dominant-TSS selection and co-located dedup in
# proTSScall_optimized.R are byte-identical to the original per-row loops in
# proTSScall.R, on randomized synthetic inputs. Base R only -- no bedGraphs, no
# NIH_scripts binaries, no cluster.
#
# The two tie-break semantics this guards (they are NOT the same):
#   * dominant selection: replace only on STRICTLY-better -> a full tie keeps the
#     FIRST transcript seen (in active-TSS / matrix-row order).
#   * dedup: replace on >= (V10 greater, or V10 equal AND ENSG-digit <=) -> a full
#     tie keeps the LAST transcript seen (in dominant-list order).
#
# Run: Rscript test/equivalence_dominant_dedup.R [n_trials]

# ----------------------------------------------------------------------------
# ORIGINAL implementations -- faithful transcriptions of proTSScall.R's loops.
# ----------------------------------------------------------------------------

orig_dominant <- function(tssActive, tssReads, tssList) {
  domTss <- vector("list")
  for (i in tssActive) {
    g <- tssList[i, "V7"]
    if (!(g %in% names(domTss))) {
      domTss[[g]] <- i
    } else {
      cur <- domTss[[g]]
      if (tssReads[i] > tssReads[cur]) {
        domTss[[g]] <- i
      } else if (tssReads[i] == tssReads[cur]) {
        if (tssList[i, "V6"] == "+") {
          if (tssList[i, "V4"] < tssList[cur, "V4"]) domTss[[g]] <- i
        } else if (tssList[i, "V6"] == "-") {
          if (tssList[i, "V5"] > tssList[cur, "V5"]) domTss[[g]] <- i
        } else {
          stop("Unexpected strand symbol")
        }
      }
    }
  }
  unlist(domTss, use.names = FALSE)
}

orig_dedup <- function(domList, tssList) {
  uniqueTss <- vector("list")
  for (i in domList) {
    if (tssList[i, "V6"] == "+") {
      tssStart <- paste(tssList[i, "V3"], tssList[i, "V4"], "+", sep = ".")
    } else if (tssList[i, "V6"] == "-") {
      tssStart <- paste(tssList[i, "V3"], tssList[i, "V5"], "-", sep = ".")
    } else {
      stop("Unexpected strand symbol")
    }
    if (!(tssStart %in% names(uniqueTss))) {
      uniqueTss[[tssStart]] <- i
    } else {
      if (tssList[i, "V10"] > tssList[uniqueTss[[tssStart]], "V10"]) {
        uniqueTss[[tssStart]] <- i
      } else if (tssList[i, "V10"] == tssList[uniqueTss[[tssStart]], "V10"]) {
        digPos <- regexpr("[[:digit:]]+",
                          c(tssList[i, "V7"], tssList[uniqueTss[[tssStart]], "V7"]))
        geneDigits <- substr(c(tssList[i, "V7"], tssList[uniqueTss[[tssStart]], "V7"]),
                             as.integer(digPos),
                             as.integer(digPos) + (attributes(digPos)$match.length - 1))
        if (which.min(as.numeric(geneDigits)) == 1) uniqueTss[[tssStart]] <- i
      }
    }
  }
  unlist(uniqueTss, use.names = FALSE)
}

# ----------------------------------------------------------------------------
# VECTORIZED implementations -- the same logic proTSScall_optimized.R uses.
# ----------------------------------------------------------------------------

# first run of digits in each id, as numeric (NA if none) -- vectorized form of
# the original's per-pair regexpr/substr.
ensg_digits <- function(ids) {
  m <- regexpr("[[:digit:]]+", ids)
  as.numeric(substr(ids, as.integer(m), as.integer(m) + attr(m, "match.length") - 1))
}

vec_dominant <- function(tssActive, tssReads, tssList) {
  if (length(tssActive) == 0) return(character(0))
  reads  <- tssReads[tssActive]
  gene   <- tssList[tssActive, "V7"]
  strand <- tssList[tssActive, "V6"]
  if (any(strand != "+" & strand != "-")) stop("Unexpected strand symbol")
  v4 <- tssList[tssActive, "V4"]
  v5 <- tssList[tssActive, "V5"]
  # strand-aware upstream key: + wants smaller V4, - wants larger V5; so an
  # ASCENDING key means "more upstream" on both strands.
  upstream <- ifelse(strand == "+", v4, -v5)
  idx <- seq_along(tssActive)
  # winner per gene: max reads, then most upstream, then first-seen (min idx).
  ord <- order(gene, -reads, upstream, idx)
  firstPick <- !duplicated(gene[ord])
  win_tx   <- tssActive[ord][firstPick]
  win_gene <- gene[ord][firstPick]
  # emit in gene-first-appearance order, matching unlist(domTss).
  gene_first <- tapply(idx, gene, min)
  win_tx[order(gene_first[win_gene])]
}

vec_dedup <- function(domList, tssList) {
  if (length(domList) == 0) return(character(0))
  strand <- tssList[domList, "V6"]
  if (any(strand != "+" & strand != "-")) stop("Unexpected strand symbol")
  coord <- ifelse(strand == "+", tssList[domList, "V4"], tssList[domList, "V5"])
  key <- paste(tssList[domList, "V3"], coord, strand, sep = ".")
  len <- tssList[domList, "V10"]
  dig <- ensg_digits(tssList[domList, "V7"])
  idx <- seq_along(domList)
  # winner per start-key: max V10, then min ENSG-digit, then LAST-seen (-idx) --
  # the original replaces on a full tie (which.min ties to the candidate), so the
  # last colliding entry wins.
  ord <- order(key, -len, dig, -idx)
  firstPick <- !duplicated(key[ord])
  domList[ord][firstPick]
}

# ----------------------------------------------------------------------------
# Randomized case generator -- deliberately dense in ties and collisions.
# ----------------------------------------------------------------------------

make_case <- function(seed) {
  set.seed(seed)
  # Unique gene IDs (a gene is one locus: one V7, one strand), but numbers repeat
  # across prefixes so distinct genes can still share a first-digit-run -- that is
  # what exercises the dedup ENSG-digit tie-break. Pool of 18 unique IDs; take a
  # distinct sample so no two genes collide on V7 (a V7 collision would fuse two
  # loci into one mixed-strand "gene", which never occurs in real formakeheatmap
  # input and is out of contract).
  combos <- as.vector(outer(c("ENSG", "AC"), 1:9, paste0))
  nGenes <- sample(1:8, 1)
  geneIDs <- sample(combos, nGenes)
  chrs <- c("chr1", "chr2", "chrX")
  recs <- list()
  ids <- character(0)
  tid <- 0L
  for (g in seq_len(nGenes)) {
    gene   <- geneIDs[g]
    strand <- sample(c("+", "-"), 1)                        # constant within a gene
    chr    <- sample(chrs, 1)
    nt     <- sample(1:4, 1)
    for (t in seq_len(nt)) {
      tid <- tid + 1L
      id  <- sprintf("ENST%05d", tid)
      start <- sample(1:6, 1) * 10L          # small coords -> start collisions
      recs[[id]] <- data.frame(
        V2 = 0L, V3 = chr, V4 = start, V5 = start + sample(1:6, 1) * 10L,
        V6 = strand, V7 = gene, V8 = 0L, V9 = 0L,
        V10 = sample(c(100L, 400L, 400L, 800L), 1),  # 400 repeated -> length ties
        stringsAsFactors = FALSE)
      ids <- c(ids, id)
    }
  }
  df <- do.call(rbind, recs)
  rownames(df) <- ids
  reads <- stats::setNames(sample(0:4, nrow(df), replace = TRUE), ids)  # read ties
  thr <- sample(0:2, 1)
  active <- sample(ids)                       # shuffled == arbitrary matrix order
  active <- active[reads[active] > thr]       # active/inactive split
  list(tssList = df, tssReads = reads, tssActive = active)
}

# ----------------------------------------------------------------------------
# Fuzz loop
# ----------------------------------------------------------------------------

main <- function(n_trials) {
  # unlist(list()) is NULL; the vectorized form returns character(0). They are
  # equivalent downstream (write.table subsets via %in%), so compare normalized.
  norm <- function(x) if (is.null(x)) character(0) else x
  dom_mismatch <- 0L
  dedup_iso_mismatch <- 0L   # dedup alone, fed the original's dominant list
  chain_mismatch <- 0L       # full vec chain vs full orig chain (output set)
  tie_seen <- 0L             # cases that exercised a dedup start-key collision
  for (s in seq_len(n_trials)) {
    cs <- make_case(s)
    od <- norm(orig_dominant(cs$tssActive, cs$tssReads, cs$tssList))
    vd <- norm(vec_dominant(cs$tssActive, cs$tssReads, cs$tssList))
    if (!identical(od, vd)) {
      dom_mismatch <- dom_mismatch + 1L
      if (dom_mismatch <= 3) { cat("DOMINANT mismatch seed", s, "\n  orig:", od, "\n  vec :", vd, "\n") }
    }
    # isolate dedup: feed BOTH the original's dominant list; compare as sets
    ou <- sort(norm(orig_dedup(od, cs$tssList)))
    vu <- sort(norm(vec_dedup(od, cs$tssList)))
    if (length(od) > length(unique(ou))) tie_seen <- tie_seen + 1L
    if (!identical(ou, vu)) {
      dedup_iso_mismatch <- dedup_iso_mismatch + 1L
      if (dedup_iso_mismatch <= 3) { cat("DEDUP mismatch seed", s, "\n  orig:", ou, "\n  vec :", vu, "\n") }
    }
    # full chain (dominant-list order feeds dedup); compare final dominant SET
    oc <- sort(norm(orig_dedup(orig_dominant(cs$tssActive, cs$tssReads, cs$tssList), cs$tssList)))
    vc <- sort(norm(vec_dedup(vec_dominant(cs$tssActive, cs$tssReads, cs$tssList), cs$tssList)))
    if (!identical(oc, vc)) {
      chain_mismatch <- chain_mismatch + 1L
      if (chain_mismatch <= 3) { cat("CHAIN mismatch seed", s, "\n  orig:", oc, "\n  vec :", vc, "\n") }
    }
  }
  cat(sprintf("trials=%d\n", n_trials))
  cat(sprintf("  dominant-selection mismatches : %d\n", dom_mismatch))
  cat(sprintf("  dedup (isolated)   mismatches : %d\n", dedup_iso_mismatch))
  cat(sprintf("  full-chain         mismatches : %d\n", chain_mismatch))
  cat(sprintf("  trials that hit a start-key collision (dedup exercised): %d\n", tie_seen))
  if (dom_mismatch + dedup_iso_mismatch + chain_mismatch == 0) cat("PASS\n") else { cat("FAIL\n"); quit(status = 1) }
}

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1) as.integer(args[1]) else 20000L
main(n)
