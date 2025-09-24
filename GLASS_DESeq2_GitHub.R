## R version 3.6.2
## installation of packages see TCGA_DESeq2_GitHub.R
library(DESeq2)
library(org.Hs.eg.db)

input_dir <- "C:/Users/gcanha/OneDrive - TUNI.fi/Documents/Secondary GBM/GLASS data/"

## read in the summed up transcripts
GLASS_transcripts <- read.delim(paste0(input_dir, "GLASS_transcripts_summed_up.txt"),
                                sep = " ")
## read in the information about the cases
GLASS_cases <- read.delim(paste0(input_dir, "GLASS_RNA_seq_info_R.txt"),
                          sep = "\t")

## remove TCGA-DU-7304 as it is only progressing to grade 3
GLASS_transcripts <- GLASS_transcripts[, -(grep("TCGA.DU", colnames(GLASS_transcripts)))]
## remove the GLSS-LX-0267 (no treatment before progression)
GLASS_transcripts <- GLASS_transcripts[, -(grep("GLSS.LX.0267", colnames(GLASS_transcripts)))]
## remove the GLSS-HK-0002-TP and GLSS.SF.0003 (only TMZ before progression)
GLASS_transcripts <- GLASS_transcripts[, -(grep("GLSS.HK.0002", colnames(GLASS_transcripts)))]
GLASS_transcripts <- GLASS_transcripts[, -(grep("GLSS.SF.0003", colnames(GLASS_transcripts)))]
## remove the GLSS-LX-0357 (no primary data)
GLASS_transcripts <- GLASS_transcripts[, -(grep("GLSS.LX.0357", colnames(GLASS_transcripts)))]

## round the samples to integers for usage with DESeq2
GLASS_transcripts <- round(GLASS_transcripts)

## make the colData
colData <- GLASS_cases[, c(2, 3, 7)]
colData$grade_group <- colData$grade_new
colData[colData$grade_new=="3" | colData$grade_new=="2", ]$grade_group = "2_3"
colData_primary <- colData[-which(colData$grade_group == "2_3" & colData$surgery_number > 1), ]
## replace the "-" in the sample name with a "." for matching
colData_primary$sample_barcode <- gsub("-", ".", colData_primary$sample_barcode)
## remove the SF-0003 from here
colData_primary <- colData_primary[-c(13, 14), ]

colnames(GLASS_transcripts) <- substr(colnames(GLASS_transcripts), 1, nchar(colnames(GLASS_transcripts)[1]) - 15)
GLASS_transcripts_primary <- GLASS_transcripts[, colnames(GLASS_transcripts) %in% colData_primary$sample_barcode]

## sort according to the order of the colData
new_order <- order(colData_primary$sample_barcode)
colData_primary <- colData_primary[new_order, ]

new_order <- order(colnames(GLASS_transcripts_primary))
GLASS_transcripts_primary <- GLASS_transcripts_primary[, new_order]

## paired DESeq2
## for GLASS it's b and a because they are denoted as "R" and "TP"
condition <- factor(c(rep(c("b", "a"), 8)))
patient <- factor(substr(colData_primary$sample_barcode, 1, nchar(colData_primary$sample_barcode) - 3))
colData_paired <- data.frame(condition, patient)

dds <- DESeqDataSetFromMatrix(countData = GLASS_transcripts_primary,
                              colData = colData_paired,
                              design = ~patient + condition)
dds <- DESeq(dds)
res <- results(dds)
resTable_all <- data.frame(ensembl=res@rownames, baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, 
                          lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj)
resTable <- resTable_all[!is.na(resTable_all$padj),]
resTable_sig <- resTable[resTable$padj < 0.05 & abs(resTable$log2FoldChange) > 1, ]
## 759 significant genes

## write to file for later usage
save(dds, file = paste0(input_dir, "GLASS_RNA_dds_object.RData"))
write.table(resTable_all, file = paste0(input_dir, "GLASS_RNA_DESeq2_all_genes.txt"),
            sep = "\t", quote = FALSE)
write.table(resTable_sig, file = paste0(input_dir, "GLASS_RNA_DESeq2_sig_genes.txt"),
            sep = "\t", quote = FALSE)

ENSEMBL_IDs <- as.character(resTable_sig$ensembl)
annots <- select(org.Hs.eg.db, keys = ENSEMBL_IDs,
                 columns = "SYMBOL", keytype = "ENSEMBL")
annots_save <- annots
annots <- unique(annots)

## exclude the ones with NA and then export for enrichr
annots_filtered <- annots[- which(is.na(annots$SYMBOL) == TRUE), ]
write.table(annots_filtered[, 2], file = paste0(input_dir, "GLASS_DE_genes_symbols.txt"),
            quote = FALSE, col.names = FALSE, row.names = FALSE)

## repeat for all genes
resTable_all <- read.delim(paste0(input_dir, "GLASS_RNA_DESeq2_all_genes.txt"),
                           sep = "\t")

ENSEMBL_IDs <- as.character(resTable_all$ensembl)
annots <- select(org.Hs.eg.db, keys = ENSEMBL_IDs,
                 columns = "SYMBOL", keytype = "ENSEMBL")
annots_save <- annots
annots <- unique(annots)

colnames(annots_save)[1] <- "ensembl"
all_genes_with_symbol <- merge(resTable_all, annots_save,
                               by.x = "ensembl")

write.table(all_genes_with_symbol, file = paste0(input_dir, "GLASS_RNA_DESeq2_all_genes_symbols.txt"),
            sep = "\t", quote = FALSE)

## check how many lncRNAs are measured
## check how many lncRNAs are measured
DE_lncRNAs_all <- c("LINC00836", "ARHGEF26-AS1", "LINC00844", "GTSCR1",
                    "LINC00842", "CRNDE", "NRG3-AS1", "FAM95B1", "RMST",
                    "MIR4435-2HG", "LNX1-AS1", "LINC01241", "LINC01736",
                    "LINC01088", "LINC00648", "CYTOR", "WDR11-AS1", "SOX5-AS1",
                    "FBXL19-AS1", "COLCA1", "LINC01094", "KCNIP4-IT1",
                    "SMC2-AS1", "EDNRB-AS1", "NEAT1", "LINC00624", "FAM225B",
                    "HAGLR", "PVT1", "LINC01224", "FAM225A", "DNMBP-AS1",
                    "DISC1-IT1")

which(DE_lncRNAs_all %in% all_genes_with_symbol$SYMBOL)
## measured: COLCA1
which(all_genes_with_symbol$SYMBOL == "COLCA1") # significant

## export the gene symbols with the LFC for the Venn Diagrams
## merge the annots with the resTable
colnames(annots_save)[1] <- "ensembl"
sig_genes_with_symbol <- merge(resTable_sig, annots_save,
                               by.x = "ensembl")

write.table(sig_genes_with_symbol, file = paste0(input_dir, "GLASS_RNA_DESeq2_sig_genes_symbols.txt"),

            sep = "\t", quote = FALSE)
