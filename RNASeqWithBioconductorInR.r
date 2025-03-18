install.packages(c("RColorBrewer", "pheatmap", "tidyverse"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")  # up-tp-date 3.20, AIMs 3.18, laptop 3.14

BiocManager::install("DESeq2")

library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)


