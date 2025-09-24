# R version 3.6.2, Bioconductor v.3.10, org.Hs... v.3.10, DESeq2 v.1.26

BiocManager::install(version="3.10")
BiocManager::install("DESeq2", version = "3.10", force = TRUE)
BiocManager::install("org.Hs.eg.db", version = "3.10")
BiocManager::install("ggplot2", version = "3.10")
BiocManager::install("ggplot2", version = "3.10")
BiocManager::install("scales", version = "3.10", force = TRUE, type = "source")
BiocManager::install("vctrs", version = "3.10")
BiocManager::install("pillar", version = "3.10")
BiocManager::install("latticeExtra", version = "3.10")
BiocManager::install("htmlwidgets", version = "3.10")
BiocManager::install("htmlTable", version = "3.10")
BiocManager::install("knitr", version = "3.10")
BiocManager::install("glue", version = "3.10", type = "source")
BiocManager::install("stringr", version = "3.10")
BiocManager::install("checkmate", version = "3.10")
BiocManager::install("htmltools", version = "3.10", type = "source")
BiocManager::install("locfit", version = "3.10")
BiocManager::install("xml2", version = "3.10")
BiocManager::install("XML", version = "3.10", type = "binary")

library(DESeq2)
library(org.Hs.eg.db)

## the data was downloaded from TCGA via TCGA biolinks
input_dir <- "C:/Users/gcanha/OneDrive - TUNI.fi/Documents/Secondary GBM/microRNAs/"

writeLines(capture.output(sessionInfo()), paste0(input_dir, "TCGA_analysis_sessionInfo.txt"))

load(paste0(input_dir, "lgg_RNA_hg38_raw.RData"))
load(paste0(input_dir, "gbm_RNA_hg38_raw.RData"))

gbm_matrix = SummarizedExperiment::assay(gbm_gene_hg38, "unstranded")
lgg_matrix = SummarizedExperiment::assay(lgg_gene_hg38, "unstranded")

# read in the WHO_2021 classification from Serafiina Jaatinen
TCGA_2021_classification <- read.table("C:/Users/gcanha/OneDrive - TUNI.fi/Documents/Secondary GBM/TCGA_WHO_2021_classification.txt",
                                       sep = "\t", header = TRUE)

## filter for all cases with WHO2021 classification
# for all of the cases remove the last 4 characters for matching
# GBMs
tcga_gbm_hg38_RNA_short <- substr(colnames(gbm_gene_hg38), 1, 12)
WHO_2021_classification_gbm_RNA <- c()
for (i in 1:length(colnames(gbm_gene_hg38))) {
  index <- grep(tcga_gbm_hg38_RNA_short[i], TCGA_2021_classification[, 1])
  if (length(index) > 0) {
    WHO_2021_classification_gbm_RNA[i] <- as.character(TCGA_2021_classification[index, 12])
  } else {
    WHO_2021_classification_gbm_RNA[i] <- NA
  }
}
table(WHO_2021_classification_gbm_RNA)

# LGGs
tcga_lgg_hg38_RNA <- substr(colnames(lgg_gene_hg38), 1, 12)
WHO_2021_classification_lgg_RNA <- c()
for (i in 1:length(colnames(lgg_gene_hg38))) {
  index <- grep(tcga_lgg_hg38_RNA[i], TCGA_2021_classification[, 1])
  if (length(index) > 0) {
    WHO_2021_classification_lgg_RNA[i] <- as.character(TCGA_2021_classification[index, 12])
  } else {
    WHO_2021_classification_lgg_RNA[i] <- NA
  }
}
table(WHO_2021_classification_lgg_RNA)

# filter for the samples for which we have sample information
# 175 samples but 174 unique submitter IDs
# TCGA-06-0156-01A has two samples, 02R and 03R
# they are from the same sample, just two different portions
# "TCGA-06-0156-01A-02R-1849-01" is excluded
gbm_matrix_filtered <- gbm_matrix[, !(is.na(WHO_2021_classification_gbm_RNA))]
lgg_matrix_filtered <- lgg_matrix[, !(is.na(WHO_2021_classification_lgg_RNA))]
# combine the LGG and GBM matrices for normalisation
# double-check that the rownames are the same between both
sum(rownames(gbm_matrix_filtered) == rownames(lgg_matrix_filtered))
TCGA_GBM_LGG_matrix <- cbind(gbm_matrix_filtered, lgg_matrix_filtered)

#expr_log2 = log2(expr_matrix + 1)

## filter for only the primary ones
primary_info <- read.delim(paste0(input_dir, "tumor_classes_recurrent_genes_expr_2209.txt"))
primary_IDHmut_tumors <- primary_info[which(primary_info$idh_codel_primary_recurrent == "IDHmut-non-codel_primary"), ]
primary_IDHmut_matrix <- TCGA_GBM_LGG_matrix[, colnames(TCGA_GBM_LGG_matrix) %in% primary_IDHmut_tumors[, 1]]

colData = primary_IDHmut_tumors[c("sample", "idh_codel_subtype", "grade", "primary_recurrent", "IDH_codel_grade", "idh_codel_primary_recurrent")]
colData$grade_group = as.character(colData$grade)
## export the colData for making the violin plots
write.table(colData, file = paste0(input_dir, "TCGA_classification_for_violins.txt"),
            sep = "\t", quote = FALSE)
colData[colData$grade=="G3" | colData$grade=="G2", ]$grade_group = "G2_G3"
colData$grade_group <- as.factor(colData$grade_group)

## make sure the colData and the primary_IDHmut_matrix are in the same order
new_order <- colData[, 1]
primary_IDHmut_matrix <- primary_IDHmut_matrix[, order(match(colnames(primary_IDHmut_matrix), new_order))]

# TCGA astrocytoma primary grade 2-3 vs. primary grade 4

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
resTable = resTable[resTable$padj<0.05 & abs(resTable$log2FoldChange)>1, ]

## save the DESeq2 object for further usage in the lncRNAs
save(dds, file = paste0(input_dir, "TCGA_mRNA_dds_object.RData"))
write.table(resTable_all, file = paste0(input_dir, "TCGA_mRNA_DESeq2_all_genes.txt"),
            sep = "\t", quote = FALSE)

## attach the gene symbols
resTable_all$ensembl <- gsub("\\..*", "", resTable_all$ensembl)

ENSEMBL_IDs <- as.character(resTable_all$ensembl)
annots <- select(org.Hs.eg.db, keys = ENSEMBL_IDs,
                 columns = "SYMBOL", keytype = "ENSEMBL")
annots_save = annots
annots = unique(annots)

colnames(annots_save)[1] <- "ensembl"
all_genes_with_symbol <- merge(resTable_all, annots_save,
                               by.x = "ensembl")

write.table(all_genes_with_symbol, file = paste0(input_dir, "TCGA_DESeq2_all_genes_symbols.txt"),
            sep = "\t", quote = FALSE)

## check how many lncRNAs are included in the DE genes
## read in the list of all lncRNAs
all_lncRNAs <- rtracklayer::import("C:/Users/gcanha/OneDrive - TUNI.fi/Documents/Secondary GBM/gencode.v30.long_noncoding_RNAs.gtf")
all_lncRNAs_df <- as.data.frame(all_lncRNAs)

### Only primary tumors
## normalise the data together with oligos, etc. for the lncRNA violin plots
colData_all = primary_info[c("sample", "idh_codel_subtype", "grade", "primary_recurrent", "IDH_codel_grade", "idh_codel_primary_recurrent")]
## make sure the TCGA_GBM_LGG_matrix is filtered
primary_colData_all <- colData_all[which(colData_all$primary_recurrent == "primary"), ]
primary_diffuse_astro_matrix <- TCGA_GBM_LGG_matrix[, (colnames(TCGA_GBM_LGG_matrix) %in% primary_colData_all$sample)]

## make sure the colData and the TCGA_GBM_LGG_matrix are in the same order
setdiff(primary_colData_all$sample, colnames(primary_diffuse_astro_matrix))
## remove the "TCGA-DU-6392-01A-11R-1708-07", that is not in the diffuse_astro_matrix
primary_colData_all <- primary_colData_all[-(which(primary_colData_all$sample == "TCGA-DU-6392-01A-11R-1708-07")), ]
new_order <- primary_colData_all[, 1]
primary_diffuse_astro_matrix <- primary_diffuse_astro_matrix[, order(match(colnames(primary_diffuse_astro_matrix), new_order))]
## numbers:
## IDHmutcodel G2: 81, G3: 70
## IDHmut noncodel G2: 110, G3: 97, G4: 19
## IDHwt: 200

## export the colData for making the violin plots
write.table(primary_colData_all, file = paste0(input_dir, "TCGA_classification_for_violins_all_primary_tumors.txt"),
            sep = "\t", quote = FALSE)

## the design variable doesn't matter
dds <- DESeqDataSetFromMatrix(countData = primary_diffuse_astro_matrix,
                              colData = primary_colData_all, design = ~ IDH_codel_grade)
dds <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
normalized_counts$ensembl <- gsub("\\..*", "", rownames(normalized_counts))

## turn ENSEMBL IDs to gene symbols
ENSEMBL_IDs <- as.character(normalized_counts$ensembl)
annots <- select(org.Hs.eg.db, keys = ENSEMBL_IDs,
                 columns = "SYMBOL", keytype = "ENSEMBL")
annots_save = annots
annots = unique(annots)

colnames(annots_save)[1] <- "ensembl"
norm_counts_with_symbol <- merge(normalized_counts, annots_save,
                                 by.x = "ensembl")

write.table(norm_counts_with_symbol, file = paste0(input_dir, "TCGA_primary_norm_counts_GBMs_oligos_included.txt"),
            sep = "\t", quote = FALSE)
