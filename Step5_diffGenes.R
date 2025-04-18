# the R script below is adapted from Daniel Beiting and modified by me

# Load packages -----
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(gt)
library(DT)
library(plotly)
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
targets <- read_tsv("studydesign.txt")
sampleLabels <- targets$sample
load(file = "../myDGEList1011FilterNorm")
log2.cpm.filtered.norm <- edgeR::cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

# ==============================================================================

# Set up your design matrix ----
group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# ==============================================================================

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)

# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# ==============================================================================

# Contrast matrix ----
contrast.matrix <- makeContrasts(infection = disease - healthy,
                                 levels=design)

# ==============================================================================

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)

#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
write.fit(ebFit, file="../lmfit1011_results.txt")

# ==============================================================================

# TopTable to view DEGs -----
#myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=10, sort.by="logFC")
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=10990, sort.by="logFC")

# convert to a tibble
myTopHits.df <- myTopHits %>% 
    as_tibble(rownames = "geneID")

temp <- left_join(log2.cpm.filtered.norm.df, myTopHits.df)
write.csv(temp, "../myDGEList1011log2_filtered_normalized_stat.csv")

#gt(myTopHits.df)
# TopTable (from Limma) outputs a few different stats:
# logFC, AveExpr, and P.Value should be self-explanatory
# adj.P.Val is your adjusted P value, also known as an FDR 
#     (if BH method was used for multiple testing correction)
# B statistic is the log-odds that that gene is differentially expressed. 
#     If B = 1.5, then log odds is e^1.5, where e is euler's constant (approx. 2.718).  
#     So, the odds of differential expression os about 4.8 to 1
# t statistic is ratio of the logFC to the standard error 
#     (where the error has been moderated across all genes...because of Bayesian approach)

# ==============================================================================

# Volcano Plots ----
vplot <- ggplot(myTopHits.df) + 
    aes(y=-log10(P.Value), x=logFC, text = paste("Symbol:", geneID)) + 
    geom_point(size=2) + 
    geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=0.8) + 
    geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=0.8) + 
    geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=0.8) + 
    annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.05), ymax = 6.5, alpha=.2, fill="#BE684D") + 
    annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.05), ymax = 6.5, alpha=.2, fill="#2C467A") + 
    labs(title="RNA-Seq: Liver Fibrosis & Healthy Groups") + 
    theme_bw()

EnhancedVolcano(myTopHits.df, 
                x="logFC", 
                y="P.Value", 
                xlim = c(min(-12), max(12)), 
                ylim = c(0, max(-log10(myTopHits.df[['P.Value']]), na.rm = TRUE) + 0.1), 
                xlab = bquote(~Log[2] ~ "Fold Change"), 
                ylab = bquote(~-Log[10] ~italic(P)~"-values"), axisLabSize = 18, 
                #title = "", 
                subtitle = "RNA-Seq: Liver Fibrosis & Healthy Groups", 
                lab = myTopHits.df$geneID, 
                pCutoff = 5e-2, 
                FCcutoff = 1)

# Now make the volcano plot above interactive with plotly
#ggplotly(vplot)

# ==============================================================================

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="none", p.value=0.05, lfc=1)

# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include=c("up", "down"))

# ==============================================================================

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes)

#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

#write your DEGs to a file
write_csv(diffGenes.df, "../DiffGenes1011.csv")

# ==============================================================================

# create interactive tables to display your DEGs ----
datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in Liver Fibrosis', 
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100"))) %>% 
    formatRound(columns=c(2:11), digits=2)

# ==============================================================================

# OPTIONAL: differential transcript usage (DTU) analysis ----
library(IsoformSwitchAnalyzeR)

# The IsoformSwitchAnalyzeR package looks for certain column headers in our study design
# So, the first step is to make sure our study design contains the following:
# unique sample IDs must be contained in column called 'sampleID'
# covariate(s) of interest must be in column labeled 'condition'
# remove extraneous columns
targets.mod <- targets %>%
  dplyr::rename(sampleID = sample, condition = group) %>%
  dplyr::select(sampleID, condition)

# import transcript Kallisto quant data
# using the same path variable we set way back in the step 1 script
Txi_trans <- importIsoformExpression(sampleVector = path)

# fix column headers of abundance and counts data to match sampleID in target.mod
colnames(Txi_trans$abundance) <- c("isoform_id", sampleLabels)
colnames(Txi_trans$counts) <- c("isoform_id", sampleLabels)

# import data
mySwitchList <- importRdata(
  isoformCountMatrix   = Txi_trans$counts,
  isoformRepExpression = Txi_trans$abundance,
  designMatrix         = targets.mod,
  removeNonConvensionalChr = TRUE,
  addAnnotatedORFs=TRUE,
  ignoreAfterPeriod=TRUE,
  # the files below must be from the same ensembl release (in this case release 108), and must match the reference release version that we originally mapped our reads to at the beginning of the course
  # you can find version 108 of the gtf file below here: https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/
  isoformExonAnnoation = "Homo_sapiens.GRCh38.111.chr_patch_hapl_scaff.gtf.gz",
  isoformNtFasta       = "Homo_sapiens.GRCh38.cdna.all.fa",
  showProgress = TRUE)


# We'll do the isoform analysis in one step, but there's a lot to unpack here, so you should really read the package documentation at:
# https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html
# Note that without additional manual work here (beyond the scope of this class), we'll only capture isoform annotations for 1) intron retention; 2) ORF sequence similarity; and 3) nonsense mediate decay (NMD)

#NOTE: THIS NEXT BIT COULD TAKE A WHILE!
mySwitchList <- isoformSwitchAnalysisCombined(
  switchAnalyzeRlist   = mySwitchList,
  pathToOutput = 'isoform_output') # directory must already exist

# now look at the directory that you just created above
# in case you missed the summary output from the function above
extractSwitchSummary(mySwitchList)

# extract the top n isoform switching events
extractTopSwitches(
  mySwitchList,
  filterForConsequences = TRUE, # these 'consequences' related to the annotations I reference above.
  n = 50,
  sortByQvals = FALSE) #change to TRUE if you want this list sorted by FDR-adusted Pval (a.k.a., q value)

# visualize by making a 'switch plot'
switchPlot(
  mySwitchList,
  gene='FCGR3B',
  condition1 = 'disease',
  condition2 = 'healthy',
  localTheme = theme_bw())
