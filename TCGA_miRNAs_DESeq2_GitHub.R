## R version 3.6.2
## installation of packages see TCGA_DESeq2_GitHub.R
library(DESeq2)
library(org.Hs.eg.db)

input_dir <- "C:/Users/gcanha/OneDrive - TUNI.fi/Documents/Secondary GBM/microRNAs/"

load(paste0(input_dir, "lgg_miRNA_hg38_raw.RData"))
load(paste0(input_dir, "gbm_miRNA_hg38_raw.RData"))

rownames(lgg_miRNA_hg38) <- lgg_miRNA_hg38[, 1]
rownames(gbm_miRNA_hg38) <- gbm_miRNA_hg38[, 1]
lgg_miRNA_matrix <- lgg_miRNA_hg38[, grep("read_count", colnames(lgg_miRNA_hg38))]
gbm_miRNA_matrix <- gbm_miRNA_hg38[, grep("read_count", colnames(gbm_miRNA_hg38))]

# read in the WHO_2021 classification from Serafiina Jaatinen
TCGA_2021_classification <- read.table("C:/Users/gcanha/OneDrive - TUNI.fi/Documents/Secondary GBM/TCGA_WHO_2021_classification.txt",
                                       sep = "\t", header = TRUE)

## filter for all cases with WHO2021 classification
# for all of the cases remove the "read_count", and the last characters
# GBMs
tcga_gbm_hg38_miRNA_short <- substr(colnames(gbm_miRNA_matrix), 12, 23)
WHO_2021_classification_gbm_miRNA <- c()
for (i in 1:length(colnames(gbm_miRNA_matrix))) {
  index <- grep(tcga_gbm_hg38_miRNA_short[i], TCGA_2021_classification[, 1])
  if (length(index) > 0) {
    WHO_2021_classification_gbm_miRNA[i] <- as.character(TCGA_2021_classification[index, 12])
  } else {
    WHO_2021_classification_gbm_miRNA[i] <- NA
  }
}
table(WHO_2021_classification_gbm_miRNA)
## no GBs with a WHO2021 classification

# LGGs
tcga_lgg_hg38_miRNA_short <- substr(colnames(lgg_miRNA_matrix), 12, 23)
WHO_2021_classification_lgg_miRNA <- c()
for (i in 1:length(colnames(lgg_miRNA_matrix))) {
  index <- grep(tcga_lgg_hg38_miRNA_short[i], TCGA_2021_classification[, 1])
  if (length(index) > 0) {
    WHO_2021_classification_lgg_miRNA[i] <- as.character(TCGA_2021_classification[index, 12])
  } else {
    WHO_2021_classification_lgg_miRNA[i] <- NA
  }
}
table(WHO_2021_classification_lgg_miRNA)
# n-counts: oligos: 82 gr2, 72 gr3
# astros: 110 gr2, gr3: 97, gr4:15
## GBs: 67

## filter for the samples with correct classification
lgg_miRNA_matrix_filtered <- lgg_miRNA_matrix[, !(is.na(WHO_2021_classification_lgg_miRNA))]
## remove the "read_count" from the colnames for the following steps to work
colnames(lgg_miRNA_matrix_filtered) <- substr(colnames(lgg_miRNA_matrix_filtered), 12, nchar(colnames(lgg_miRNA_matrix_filtered)) - 8)

## check whether the samples are primary or recurrent
primary_info <- read.delim(paste0(input_dir, "tumor_classes_recurrent_genes_expr_2209.txt"))
primary_IDHmut_tumors <- primary_info[which(primary_info$idh_codel_primary_recurrent == "IDHmut-non-codel_primary"), ]
primary_status <- as.character(primary_IDHmut_tumors[, 1])
primary_IDHmut_tumors[, 1] <- substr(primary_status, 1, nchar(primary_status) - 8)
## it seems like the last two letters do not match, try it without them
primary_IDHmut_matrix <- lgg_miRNA_matrix_filtered[, colnames(lgg_miRNA_matrix_filtered) %in% primary_IDHmut_tumors[, 1]]

## filter the primary_IDHmut_tumor to make it match the primary_IDHmut_matrix
primary_IDHmut_tumors_filtered <- primary_IDHmut_tumors[primary_IDHmut_tumors[, 1] %in% colnames(primary_IDHmut_matrix), ]

colData = primary_IDHmut_tumors_filtered[c("sample", "idh_codel_subtype", "grade", "primary_recurrent", "IDH_codel_grade", "idh_codel_primary_recurrent")]
## export the colData for making the violin plots
write.table(colData, file = paste0(input_dir, "TCGA_miRNA_classification_for_violins.txt"),
            sep = "\t", quote = FALSE)

colData$grade_group = as.character(colData$grade)
colData[colData$grade=="G3" | colData$grade=="G2", ]$grade_group = "G2_G3"
colData$grade_group <- as.factor(colData$grade_group)

## make sure the colData and the primary_IDHmut_matrix are in the same order
new_order <- colData[, 1]
primary_IDHmut_matrix <- primary_IDHmut_matrix[, order(match(colnames(primary_IDHmut_matrix), new_order))]

## using the primary information
dds = DESeqDataSetFromMatrix(countData = primary_IDHmut_matrix,
                             colData = colData,
                             design = ~grade_group)

dds = DESeq(dds)

res = results(dds)

resTable = data.frame(ensembl=res@rownames, baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, 
                      lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj)

resTable_all = data.frame(ensembl=res@rownames, baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, 
                          lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj)

resTable = resTable[!is.na(resTable$padj),]

resTable = resTable[resTable$padj<0.1 & abs(resTable$log2FoldChange)>1, ]
## 30 DE genes with 0.05, 35 with 0.1

## write to file for the Venns
save(dds, file = paste0(input_dir, "TCGA_miRNA_dds_object.RData"))
write.table(resTable_all, file = paste0(input_dir, "TCGA_miRNA_DESeq2_all_genes.txt"),
            sep = "\t", quote = FALSE)
write.table(resTable, file = paste0(input_dir, "TCGA_miRNA_DESeq2_0.1_FC_1.txt"),
            sep = "\t", quote = FALSE)

## normalise together with oligos for the violin plots
primary_info <- read.delim(paste0(input_dir, "tumor_classes_recurrent_genes_expr_2209.txt"))
### only primary tumors
## colData and filtering
colData_all <- primary_info[c("sample", "idh_codel_subtype", "grade", "primary_recurrent", "IDH_codel_grade", "idh_codel_primary_recurrent")]
colData_all[, 1] <- as.character(colData_all[, 1])
colData_all[, 1] <- substr(colData_all[, 1], 1, nchar(colData_all[, 1]) - 8)
primary_colData_all <- colData_all[which(colData_all$primary_recurrent == "primary"), ]
primary_diffuse_astro_matrix <- lgg_miRNA_matrix_filtered[, (colnames(lgg_miRNA_matrix_filtered) %in% primary_colData_all$sample)]
primary_colData_all <- primary_colData_all[primary_colData_all$sample %in% colnames(primary_diffuse_astro_matrix), ]
## sample numbers
## IDHmut codel: gr2: 81, 70 gr3
## IDHmut astro: gr2: 107, gr3: 96, gr4: 11
## IDHwt: 64

## make sure the colData and the TCGA_GBM_LGG_matrix are in the same order
setdiff(primary_colData_all$sample, colnames(primary_diffuse_astro_matrix))
new_order <- primary_colData_all[, 1]
primary_diffuse_astro_matrix <- primary_diffuse_astro_matrix[, order(match(colnames(primary_diffuse_astro_matrix), new_order))]

## export the colData for making the violin plots
write.table(colData_all, file = paste0(input_dir, "TCGA_classification_for_violins_all_primary_tumors_miRNAs.txt"),
            sep = "\t", quote = FALSE)

## the design variable doesn't matter
dds <- DESeqDataSetFromMatrix(countData = primary_diffuse_astro_matrix,
                              colData = primary_colData_all, design = ~ IDH_codel_grade)
dds <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))

write.table(normalized_counts, file = paste0(input_dir, "TCGA_miRNA_primary_norm_counts_GBMs_oligos_included.txt"),
            sep = "\t", quote = FALSE)
