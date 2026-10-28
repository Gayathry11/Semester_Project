# Project 4: NKX2-1 Expression in TCGA LUAD

## Aim
To compare the expression of the biomarker gene NKX2-1 between Primary Tumor
and Solid Tissue Normal samples in TCGA Lung Adenocarcinoma (LUAD), using
RNA-seq data.

## Data
Downloaded from the figshare link provided with the assignment:
- `gdc_sample_sheet.tsv`: sample information (sample type for each file)
- `gdc_manifest.txt`: list of files
- `tcga_data.tar.gz`: gene expression files (one folder per sample)

The large data files are not included in this repository.

## Method
1. Read the sample sheet and keep only "Primary Tumor" and "Solid Tissue Normal" samples.
2. For each sample, open its `rna_seq.augmented_star_gene_counts.tsv` file and
   extract the `tpm_unstranded` value for NKX2-1.
3. Transform the values to log2(TPM + 1).
4. Compare the two groups with a Wilcoxon rank-sum test.
5. Draw a boxplot of the two groups.

Sanity check: the sample e6e83fd4-2f80-44e1-8dce-588a37308dbb gives a TPM of
53.1718, as given in the assignment.

## How to run
Requires R. Place `gdc_sample_sheet.tsv` and the unzipped `tcga_data` folder
in the same directory as the script, then run:

```
Rscript NKX2-1_LUAD_boxplot.R
```

## Output
- `NKX2-1_boxplot.png`: boxplot of log2(TPM + 1) for Primary Tumor and Solid Tissue Normal

## Results
- Primary Tumor: n = 539
- Solid Tissue Normal: n = 59
- Wilcoxon rank-sum test: p = 0.78

The median NKX2-1 expression is almost the same in the two groups (about 7 on
the log2(TPM + 1) scale), so no significant difference was found in this dataset.

## Files
- `NKX2-1_LUAD_boxplot.R`: analysis script
- `NKX2-1_boxplot.png`: output plot
- `README4.md`: this file
