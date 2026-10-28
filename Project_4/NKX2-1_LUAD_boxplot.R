# Project 4: NKX2-1 expression in LUAD (TCGA), Primary Tumor vs Solid Tissue Normal

ss <- read.delim("gdc_sample_sheet.tsv", stringsAsFactors = FALSE)
ss <- subset(ss, Sample.Type %in% c("Primary Tumor", "Solid Tissue Normal"))

get_tpm <- function(file_id, file_name) {
  f <- file.path("tcga_data", file_id, file_name)
  if (!file.exists(f)) return(NA_real_)
  d <- read.delim(f, comment.char = "#")
  d$tpm_unstranded[d$gene_name == "NKX2-1"][1]
}

ss$tpm <- mapply(get_tpm, ss$File.ID, ss$File.Name)
ss <- ss[!is.na(ss$tpm), ]
ss$log2tpm <- log2(ss$tpm + 1)

# Sanity check: should print 53.1718
print(ss$tpm[ss$File.ID == "e6e83fd4-2f80-44e1-8dce-588a37308dbb"])

print(table(ss$Sample.Type))
pval <- wilcox.test(log2tpm ~ Sample.Type, data = ss)$p.value
print(pval)

png("NKX2-1_boxplot.png", width = 800, height = 700)
boxplot(log2tpm ~ Sample.Type, data = ss, notch = TRUE,
        main = "NKX2-1 expression in LUAD",
        xlab = "Tissue type", ylab = "log2(TPM + 1)",
        col = c("skyblue", "tomato"))
legend("topleft", legend = paste("p =", signif(pval, 2)), bty = "n")
dev.off()
