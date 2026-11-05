# Project 4: NKX2-1 expression in TCGA-LUAD (Primary Tumor vs Solid Tissue Normal)
# Run in the folder with gdc_sample_sheet.tsv and tcga_data/

ss <- read.delim("gdc_sample_sheet.tsv")
ss <- ss[ss$Sample.Type %in% c("Solid Tissue Normal", "Primary Tumor"), ]
ss$Sample.Type <- factor(ss$Sample.Type, levels = c("Solid Tissue Normal", "Primary Tumor"))

# Read each file as text and take only the NKX2-1 line; column 7 = tpm_unstranded
get_tpm <- function(id, f) {
  p <- file.path("tcga_data", id, f)
  if (!file.exists(p)) return(NA_real_)
  x <- readLines(p)
  as.numeric(strsplit(x[grepl("\tNKX2-1\t", x, fixed = TRUE)][1], "\t")[[1]][7])
}
cores <- if (.Platform$OS.type == "windows") 1L else 4L
ss$tpm <- unlist(parallel::mcmapply(get_tpm, ss$File.ID, ss$File.Name, mc.cores = cores))
ss <- ss[!is.na(ss$tpm), ]
ss$log2tpm <- log2(ss$tpm + 1)

# Check against the value given in the assignment (53.1718)
stopifnot(abs(ss$tpm[ss$File.ID == "e6e83fd4-2f80-44e1-8dce-588a37308dbb"] - 53.1718) < 1e-3)

# Sample counts and Wilcoxon rank-sum test
n <- table(ss$Sample.Type)
print(n)
pval <- wilcox.test(log2tpm ~ Sample.Type, data = ss, exact = FALSE)$p.value
print(pval)
p_lab <- if (pval < 0.001) sprintf("p = %.1e", pval) else sprintf("p = %.2f", pval)

# Output file name (if the old image is locked/open in a viewer, use a new name)
out_file <- "NKX2-1_boxplot.png"
if (file.exists(out_file) && !suppressWarnings(file.remove(out_file))) {
  out_file <- format(Sys.time(), "NKX2-1_boxplot_%H%M%S.png")
}

# Boxplot
png(out_file, width = 1400, height = 1000, res = 150)
par(mar = c(7, 5, 4, 2), mgp = c(3, 1.5, 0), las = 1)
ytop <- max(ss$log2tpm) * 1.08
boxplot(log2tpm ~ Sample.Type, data = ss, notch = TRUE, pch = 16, cex = 0.5,
        col = c("skyblue", "tomato"), ylim = c(0, ytop * 1.05),
        names = paste0(names(n), "\n(n = ", n, ")"),
        main = "NKX2-1 Expression in Lung Adenocarcinoma (LUAD)",
        ylab = "Log2 (TPM + 1)", xlab = "",
        cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1)
mtext("Tissue Type", side = 1, line = 5.5, cex = 1.2)
segments(1, ytop, 2, ytop)
text(1.5, ytop, paste("Wilcoxon", p_lab), pos = 3, cex = 1.1)
dev.off()
cat("Plot saved:", out_file, "\n")
