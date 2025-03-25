## R version 3.6.2
## installation of packages see TCGA_DESeq2_GitHub.R
library(DESeq2)
library(org.Hs.eg.db)
library(openxlsx)

input_dir <- "C:/Users/gcanha/OneDrive - TUNI.fi/Documents/Secondary GBM/microRNA Paper/Sloan data/"

## read in the summed up transcripts
Sloan_raw_matrix <- read.xlsx(paste0(input_dir, "GSE276176_hCS_feature_counts.xlsx"),
                               sep = " ")
## remove the length
Sloan_raw_matrix <- Sloan_raw_matrix[, -(2)]

## identify the duplicate rownames first
dups <- which(duplicated(Sloan_raw_matrix$Geneid) == TRUE)
## remove the duplicates
Sloan_raw_matrix <- Sloan_raw_matrix[-(dups), ]
## make the GeneID the rownames
rownames(Sloan_raw_matrix) <- Sloan_raw_matrix$Geneid
Sloan_raw_matrix <- Sloan_raw_matrix[, -(1)]

## create a short table to include some basic sample info
sample_names <- colnames(Sloan_raw_matrix)
## C3 is male, C4 is female
colData <- data.frame(
  Sample = sample_names,
  Sex = ifelse(grepl("^C3", sample_names), "Male", "Female")
)

## design needs to be filled somehow
dds <- DESeqDataSetFromMatrix(countData = Sloan_raw_matrix,
                              colData = colData,
                              design = ~Sex)

## extract the normalised counts
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)

## export them
write.xlsx(normalized_counts, file = paste0(input_dir, "Sloan_norm_counts.xlsx"),
           rowNames = TRUE)