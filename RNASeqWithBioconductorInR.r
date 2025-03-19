install.packages(c("RColorBrewer", "pheatmap", "tidyverse"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")  # up-tp-date 3.20, AIMs 3.18, laptop 3.14

BiocManager::install("DESeq2")

library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)

head(df_rawCountMatrix, 6)
str(df_rawCountMatrix)

ggplot(df_rawCountMatrix) + 
    geom_histogram(aes(x = wt_normal1), 
        stat = "bin", 
        bins = 200) + 
    xlab("Raw Expression Counts") + 
    ylab("Number of Genes")

genotype <- c("wt", "wt", "wt", "wt", "wt", "wt", "wt")

condition <- c("normal", "fibrosis", "normal", 
    "fibrosis", "normal", "fibrosis", "fibrosis")

metadata <- data.frame(genotype, condition)

rownames(metadata) <- c("wt_normal3", "wt_fibrosis3", "wt_normal1", 
    "wt_fibrosis2", "wt_normal2", "wt_fibrosis4", "wt_fibrosis1")

object_DESeq2 <- DESeqDataSetFromMatrix(
    countData = df_rawCountMatrix, 
    colData = metadata, 
    design = ~ condition
    )


