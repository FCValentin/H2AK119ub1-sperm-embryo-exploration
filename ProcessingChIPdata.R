# =============================================================================
# ProcessingChIPdata.R
# -----------------------------------------------------------------------------
# Downstream ChIP-seq statistical analyses for H2AK119ub1 characterisation.
# Standalone script — NOT called by the shell pipeline, run independently.
#
# Sections:
#   1. Permutation tests — USP21-sensitive gene overlap with H2AK119ub1 peaks
#   2. Randomisation tests — CpG islands, repeats, genomic features
#   3. TF motif dot heatmap (HOMER findMotifs.pl output)
#   4. Peak annotation with ChIPseeker
#   5. TSS profile p-value (Wilcoxon) from deepTools plotProfile output
#   6. Proportion tests — H2A vs non-H2A / Nucl vs sub-nucl per genomic bin
#   7. Hidden Markov Model (HMM) state assignment on fragment probabilities
#
# Language : R
#
# Author   : Valentin FRANCOIS--CAMPION, PhD
# Contact  : valentin.francoiscampion@gmail.com
# GitHub   : https://github.com/FCValentin/H2AK119ub1-sperm-embryo-exploration
# Paper    : François-Campion V. et al., Nature Communications, 2025
#            DOI: 10.1038/s41467-025-58615-7
# =============================================================================


# =============================================================================
# DEPENDENCIES
# =============================================================================

.load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

for (pkg in c("GenomicFeatures", "GenomicRanges", "ChIPseeker",
              "genomation", "regioneR", "clusterProfiler",
              "patchwork", "cowplot", "viridis", "ggrepel",
              "ggplot2", "dplyr", "RHmm")) {
  .load_pkg(pkg)
}


# =============================================================================
# PARAMETERS — edit here
# =============================================================================

DATA_DIR   <- "."            # working directory (BED / TSV files)
GFF_DIR    <- "GFF"          # GTF + chromosome BED directory
FIGURES_DIR <- "Figures/ChIPseq"

N_PERM_USP21 <- 10000        # permutations for USP21 overlap tests
N_PERM_RND   <- 1000         # randomisations for genomic feature tests

# TSS profile parameters (deepTools plotProfile output)
TSS_FLANK_BINS <- c(81, 120) # BorneDown:BorneUp — 1 kb around TSS at binsize 50
TSS_ROW_COMPARE <- 37        # row index: sample to compare (MZ genes)
TSS_ROW_BASELINE <- 42       # row index: baseline (All genes)

# HMM parameters
HMM_STATES   <- 2
HMM_CHUNK    <- 4000         # chunk size for Viterbi decoding

create_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}
create_dir(FIGURES_DIR)


# =============================================================================
# I. PERMUTATION TESTS — USP21-sensitive gene overlap with H2AK119ub1 peaks
# =============================================================================

message("[1/7] Permutation tests — USP21 gene overlap...")

# Genomic universes
universe_tss       <- toGRanges(file.path(DATA_DIR, "Genes9_2TSS5Kb_chrosomose.bed"))
universe_genes     <- toGRanges(file.path(DATA_DIR, "Genes9_2_chrosomose.bed"))
universe_enhancers <- toGRanges(file.path(DATA_DIR, "Genes9_2Enhancers_chrosomose.bed"))

# USP21-sensitive gene sets
bg_tss       <- toGRanges(file.path(DATA_DIR, "USP21TSS_chr.bed"))
bg_genes     <- toGRanges(file.path(DATA_DIR, "USP21Genes.bed"))
bg_enhancers <- toGRanges(file.path(DATA_DIR, "USP21Enhancers.bed"))

# Peak set
peak_set <- toGRanges(file.path(DATA_DIR, "H2AK119ub1SpermPeaks.bed"))

pdf(file.path(FIGURES_DIR, "EnrichmentUSP21ReplicatedU21.pdf"),
    width = 10, height = 10)

pt_tss <- permTest(
  A = bg_tss, B = peak_set,
  ntimes = N_PERM_USP21, count.once = TRUE,
  randomize.function = resampleRegions,
  evaluate.function  = numOverlaps,
  universe = universe_tss,
  allow.overlaps = TRUE, verbose = TRUE
)
summary(pt_tss); plot(pt_tss)

pt_genes <- permTest(
  A = bg_genes, B = peak_set,
  ntimes = N_PERM_USP21, count.once = TRUE,
  randomize.function = resampleRegions,
  evaluate.function  = numOverlaps,
  universe = universe_genes,
  allow.overlaps = TRUE, verbose = TRUE
)
summary(pt_genes); plot(pt_genes)

pt_enhancers <- permTest(
  A = bg_enhancers, B = peak_set,
  ntimes = N_PERM_USP21, count.once = TRUE,
  randomize.function = resampleRegions,
  evaluate.function  = numOverlaps,
  universe = universe_enhancers,
  allow.overlaps = TRUE, verbose = TRUE
)
summary(pt_enhancers); plot(pt_enhancers)

dev.off()
message("  Saved: EnrichmentUSP21ReplicatedU21.pdf")


# =============================================================================
# II. RANDOMISATION TESTS — CpG islands, repeats, genomic features
# =============================================================================

message("[2/7] Randomisation tests — CpG / repeats / intergenic...")

genome_gr <- toGRanges(file.path(GFF_DIR, "XL9_Chromosomes.bed"))
peak_u21  <- toGRanges(file.path(DATA_DIR, "Egg-U21_only.bed"))
cpg       <- toGRanges(file.path(DATA_DIR, "cpgIslandExt.bed"))
repeats   <- toGRanges(file.path(DATA_DIR, "repeats.bed"))
intergenic <- toGRanges(file.path(DATA_DIR, "intergenic.bed"))

pdf(file.path(FIGURES_DIR, "EnrichmentTestFig_CpG_repeats_CO.pdf"),
    width = 10, height = 10)

for (test_cfg in list(
  list(A = peak_u21,  B = cpg,        label = "CpG islands"),
  list(A = peak_u21,  B = repeats,    label = "Repeats"),
  list(A = intergenic, B = peak_u21,  label = "Intergenic regions")
)) {
  res <- permTest(
    A = test_cfg$A, B = test_cfg$B,
    ntimes = N_PERM_RND,
    randomize.function = randomizeRegions,
    evaluate.function  = numOverlaps,
    count.once = TRUE, genome = genome_gr,
    allow.overlaps = TRUE, verbose = TRUE
  )
  summary(res); plot(res, main = test_cfg$label)
}

dev.off()
message("  Saved: EnrichmentTestFig_CpG_repeats_CO.pdf")


# =============================================================================
# III. TF MOTIF DOT HEATMAP (HOMER findMotifs.pl output)
# =============================================================================

message("[3/7] TF motif dot heatmap...")

motif_file <- file.path(DATA_DIR, "HomerLogPValueInf5.tsv")
if (file.exists(motif_file)) {
  gene_cluster <- read.table(motif_file, header = TRUE,
                              sep = "\t", stringsAsFactors = FALSE)
  p <- ggplot(gene_cluster,
              aes(x = Cluster, y = MOTIF_NAME,
                  size  = log10(1 - MinlogP.Value),
                  color = as.numeric(Percentage_motif))) +
    geom_point() +
    ylab("Gene Motif") +
    scale_colour_gradientn(
      colours = c("red", "blue", "green"),
      values  = c(0, 0.1, 1),
      limits  = c(0, 100),
      name    = "% Motif"
    ) +
    labs(size = "log10(1 - p-value)") +
    theme_bw(base_size = 8)

  pdf(file.path(FIGURES_DIR, "DotHeatmapHOMERFilter.pdf"),
      width = 100, height = 100)
  print(p)
  dev.off()
  message("  Saved: DotHeatmapHOMERFilter.pdf")
} else {
  warning("HOMER motif file not found: ", motif_file, " — skipping.")
}


# =============================================================================
# IV. PEAK ANNOTATION WITH CHIPSEEKER
# =============================================================================

message("[4/7] Peak annotation with ChIPseeker...")

txdb <- makeTxDbFromGFF(file.path(GFF_DIR, "XENLA_9.2_Xenbase.gtf"))
promoter <- getPromoters(TxDb = txdb, upstream = 5000, downstream = 5000)

peak_files <- list.files(
  path       = DATA_DIR,
  pattern    = "_0.001.bed$",
  full.names = TRUE
)

if (length(peak_files) == 0) {
  warning("No *_0.001.bed files found in: ", DATA_DIR)
} else {
  pdf(file.path(FIGURES_DIR, "ChipSeekerHMMGroups.pdf"),
      width = 10, height = 10)
  for (f in peak_files) {
    message("  Annotating: ", basename(f))
    peak_anno <- annotatePeak(f,
                               tssRegion = c(-5000, 5000),
                               TxDb = txdb)
    plotAnnoPie(peak_anno)
  }
  dev.off()
  message("  Saved: ChipSeekerHMMGroups.pdf")
}


# =============================================================================
# V. TSS PROFILE P-VALUE (Wilcoxon) from deepTools plotProfile output
# =============================================================================

message("[5/7] TSS profile Wilcoxon test...")

profile_file <- file.path(DATA_DIR, "Droso_USP21Sensitive_LadderH2Aub.tsv")
if (file.exists(profile_file)) {
  tbl <- read.table(profile_file, header = FALSE,
                    sep = "\t", stringsAsFactors = FALSE)
  # BorneDown:BorneUp selects ~1 kb around TSS (binsize 50 bp)
  bins         <- TSS_FLANK_BINS[1]:TSS_FLANK_BINS[2]
  baseline     <- as.numeric(tbl[TSS_ROW_BASELINE, bins])
  compared     <- as.numeric(tbl[TSS_ROW_COMPARE,  bins])

  my_data <- data.frame(
    group  = rep(c("AllGenes", "MZ"), each = length(baseline)),
    values = c(baseline, compared)
  )

  wt <- wilcox.test(values ~ group, data = my_data,
                    paired = TRUE, alternative = "less")
  message("  Wilcoxon p-value (MZ vs All): ", signif(wt$p.value, 4))

  stats <- group_by(my_data, group) %>%
    summarise(n      = n(),
              median = median(values, na.rm = TRUE),
              IQR    = IQR(values,    na.rm = TRUE),
              .groups = "drop")
  print(stats)
} else {
  warning("Profile TSV not found: ", profile_file, " — skipping.")
}


# =============================================================================
# VI. PROPORTION TESTS — H2A vs non-H2A / Nucl vs sub-nucl per genomic bin
# Output: .probabilities file read by H2A_enrichment.sh downstream
# =============================================================================

message("[6/7] Proportion tests (H2A / Nucl models)...")

#' Run proportion tests on a fragment count bedgraph
#'
#' @param frag_file Character. Path to 7-column bedgraph:
#'   chr, start, end, count_A, count_B, count_C, total
#' @param model     Character. "H2A" or "Nucl".
#' @return Invisible NULL. Writes .probabilities file.
run_proportion_tests <- function(frag_file, model = c("H2A", "Nucl")) {
  model <- match.arg(model)
  if (!file.exists(frag_file)) {
    warning("Fragment file not found: ", frag_file, " — skipping.")
    return(invisible(NULL))
  }
  message("  Model: ", model, " | File: ", basename(frag_file))

  frag <- read.table(frag_file, header = FALSE, sep = "\t",
                     stringsAsFactors = FALSE)
  colnames(frag) <- c("Chr", "Start", "End", "C150", "C110", "C70", "Total")
  n <- nrow(frag)

  if (model == "H2A") {
    # H2A-positive (C150 + C110) vs H2A-negative (C70)
    col_sum <- colSums(frag[, 4:7])
    f_h2a   <- 0.5 * (col_sum[1] / col_sum[4] + col_sum[2] / col_sum[4])
    f_non   <- col_sum[3] / col_sum[4]
    p_h2a   <- numeric(n)
    p_non   <- numeric(n)
    for (i in seq_len(n)) {
      if (i %% 50000 == 0) message("    Row: ", i, " / ", n)
      if (frag$Total[i] > 0) {
        p_h2a[i] <- prop.test(frag$C150[i] + frag$C110[i], frag$Total[i],
                               p = f_h2a, alternative = "greater")$p.value
        p_non[i] <- prop.test(frag$C70[i], frag$Total[i],
                               p = f_non,  alternative = "greater")$p.value
      } else {
        p_h2a[i] <- 1
        p_non[i] <- 1
      }
    }
    out <- cbind(frag,
                 (frag$C150 + frag$C110) / frag$Total,
                 frag$C70 / frag$Total,
                 p_h2a, p_non, 1)
    suffix <- "_H2A_vs_NonH2A.probabilities"

  } else {
    # Nucleosome (C150) vs sub-nucleosome (C110 + C70)
    col_sum  <- colSums(frag[, 4:7])
    f_nucl   <- col_sum[1] / col_sum[4]
    f_semi   <- 0.5 * (col_sum[2] / col_sum[4] + col_sum[3] / col_sum[4])
    p_nucl   <- numeric(n)
    p_semi   <- numeric(n)
    for (i in seq_len(n)) {
      if (i %% 50000 == 0) message("    Row: ", i, " / ", n)
      if (frag$Total[i] > 0) {
        p_nucl[i] <- prop.test(frag$C150[i], frag$Total[i],
                                p = f_nucl, alternative = "greater")$p.value
        p_semi[i] <- prop.test(frag$C110[i] + frag$C70[i], frag$Total[i],
                                p = f_semi, alternative = "greater")$p.value
      } else {
        p_nucl[i] <- 1
        p_semi[i] <- 1
      }
    }
    out <- cbind(frag,
                 frag$C150 / frag$Total,
                 (frag$C110 + frag$C70) / frag$Total,
                 p_nucl, p_semi, 1)
    suffix <- "_Nucl_vs_SemiNucl.probabilities"
  }

  out_file <- paste0(sub("\\.bedgraph$", "", frag_file), suffix)
  write.table(out, file = out_file,
              sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  message("  Saved: ", basename(out_file))
  invisible(NULL)
}

run_proportion_tests(
  file.path(DATA_DIR,
    "Fragment/Input_Sp_WithoutDupl.Fragm150-110_vs70_sum_nodup_250bpbin_slid50_above5_H2A_vs_NonH2A.full.bedgraph"),
  model = "H2A"
)
run_proportion_tests(
  file.path(DATA_DIR,
    "Fragment/Input_Sp_WithoutDupl.Fragm150_vs110-70_sum_nodup_250bpbin_slid50_above5_Nucl_vs_SemiNucl.full.bedgraph"),
  model = "Nucl"
)


# =============================================================================
# VII. HIDDEN MARKOV MODEL — chromatin state assignment on fragment probs
# =============================================================================

message("[7/7] HMM state assignment (RHmm)...")

hmm_file <- file.path(DATA_DIR,
  "Fragment/Input.Fragm150_110_70_sum_nodup_above5_slid50.full.bedgraph.probabilities.discrete")

if (file.exists(hmm_file)) {
  hmm_data <- read.table(hmm_file, stringsAsFactors = FALSE,
                          sep = "\t", header = FALSE)
  colnames(hmm_data) <- c(
    "Chr", "Begin", "End",
    "Count150", "Count110", "Count70", "TotalCount",
    "Pct150", "Pct110", "Pct70",
    "p150", "p110", "p70",
    "obs150", "obs110", "obs70"
  )
  hmm_data[is.na(hmm_data)] <- 0

  message("  Fitting HMM (", HMM_STATES, " states)...")
  hmm_fit <- HMMFit(hmm_data, nStates = HMM_STATES, dis = "DISCRETE")

  # Chunked Viterbi decoding (avoids memory issues on large datasets)
  n      <- nrow(hmm_data)
  breaks <- c(seq(1, n, by = HMM_CHUNK), n + 1)
  states <- c()
  for (i in seq_len(length(breaks) - 1)) {
    chunk  <- hmm_data[breaks[i]:(breaks[i + 1] - 1), ]
    vit    <- viterbi(hmm_fit, chunk)
    states <- c(states, vit$states)
  }

  out_hmm <- file.path(DATA_DIR,
    "Fragment/Input.Fragm150_above5_hmm_state_50bpw.full")
  write.table(states, file = out_hmm,
              sep = "\t", quote = FALSE,
              col.names = FALSE, row.names = FALSE)
  message("  HMM states saved: ", basename(out_hmm))
} else {
  warning("HMM input file not found: ", hmm_file, " — skipping.")
}

message("ProcessingChIPdata complete.")
