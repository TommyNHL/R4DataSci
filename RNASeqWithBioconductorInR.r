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

