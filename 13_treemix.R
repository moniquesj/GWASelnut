##### TREEMIX ANALYSIS

## convert vcf to treemix input files
vcf2treemix.sh population.snps.hz.tombul.filtered.pruned.leaf.vcf.gz clust.clust

## run treemix for 5 m and with 5 iterations
for m in 0 1 2 3 4 5; do
for rep in {1..5}; do     # ≥3 reps per m so OptM is happy
treemix -i population.snps.hz.tombul.filtered.pruned.leaf.treemix.frq.gz \
-k 25 -m $m -se -seed $RANDOM -o ./treemix_output_new/treemix_m${m}_r${rep}
done
done

## Visualize treemix results

outdir <- "./treemix_output_withiterations"

ll_m0 <- list.files(outdir, pattern = "^treemix_m0_r\\d+\\.llik$", full.names = TRUE)
setNames(lapply(ll_m0, readLines), basename(ll_m0))

# Extract first numeric value from a .llik file (handles labels, blanks, sci notation)
read_llik_numeric <- function(f) {
  txt <- readLines(f, warn = FALSE)
  if (!length(txt)) return(NA_real_)
  # pick first number anywhere in the file
  m <- regexpr("-?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?", txt, perl = TRUE)
  take <- which(m > 0)[1]
  if (is.na(take)) return(NA_real_)
  val <- substr(txt[take], m[take], m[take] + attr(m, "match.length")[take] - 1L)
  suppressWarnings(as.numeric(val))
}

prefixes_for_m <- function(m) {
  verts <- list.files(outdir,
                      pattern = sprintf("^treemix_m%d_r\\d+\\.vertices\\.gz$", m),
                      full.names = TRUE)
  sub("\\.vertices\\.gz$", "", verts)
}

best_prefix_for_m <- function(m) {
  # try lliks first
  lliks <- list.files(outdir,
                      pattern = sprintf("^treemix_m%d_r\\d+\\.llik$", m),
                      full.names = TRUE)
  prefs <- prefixes_for_m(m)
  
  # if no .llik at all, fall back immediately
  if (length(lliks) == 0L) {
    if (length(prefs) == 0L) stop(sprintf("No replicates found for m=%d in %s", m, outdir))
    warning(sprintf("No .llik files for m=%d; using first replicate: %s", m, basename(prefs[1])))
    return(prefs[1])
  }
  
  vals <- vapply(lliks, read_llik_numeric, numeric(1))
  if (all(is.na(vals))) {
    # all NA — fall back to first available replicate
    if (length(prefs) == 0L) stop(sprintf("No replicates found for m=%d in %s", m, outdir))
    warning(sprintf("All .llik values NA for m=%d; using first replicate: %s", m, basename(prefs[1])))
    return(prefs[1])
  }
  
  # pick the llik file with the max numeric value
  best_llik <- lliks[ which.max(vals) ]
  sub("\\.llik$", "", best_llik)
}

ms <- 0:5

# ensure to not overwrite an open file (Windows sharing violation)
if (file.exists("treemix_trees_m0-5_best.pdf")) unlink("treemix_trees_m0-5_best.pdf", force = TRUE)

pdf("treemix_trees_m0-5_best.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(3, 3, 2, 1))
for (m in ms) {
  pref <- best_prefix_for_m(m)
  plot_tree(pref)
  title(main = sprintf("TreeMix tree (m = %d) — best rep", m), line = 0.5, cex.main = 1)
}
dev.off()

pdf("treemix_residuals_m0-5_best.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(11, 11, 3, 2))
par(cex.axis = 0.8, cex.lab = 0.9)
for (m in ms) {
  pref <- best_prefix_for_m(m)
  plot_resid(pref, "pop_order.txt")
  title(main = sprintf("TreeMix residuals (m = %d) — best rep", m), line = 0.5, cex.main = 1)
}
dev.off()

##### estimate optimal m value #####
res <- optM(folder = outdir, method = "Evanno")

# Evanno method stores the pick as 'best_m'
res$best_m   # suggested m

# Save OptM elbow plot to a PDF (as the function expects) and convert to PNG
# (adds one lightweight dependency)
if (!requireNamespace("pdftools", quietly = TRUE)) install.packages("pdftools")

plot_optM(res, pdf = "optM_elbow.pdf")  # OptM handles the plotting/saving itself

# Convert the single-page PDF to PNG
pdftools::pdf_convert("optM_elbow.pdf", filenames = "optM_elbow.png", dpi = 150)
