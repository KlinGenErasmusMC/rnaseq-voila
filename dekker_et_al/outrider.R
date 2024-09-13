#! /usr/bin/env Rscript

# The typical workflow for analyzing RNA-seq data described below
# is mainly based on the paper by Love et al, F1000 Research 2016

# An automatic workflow to look for outliers in a dataset with Bioconductor package "OUTRIDER"

# use arguments between quotes: 'experiment' 'species' 'treatment' 'fragment' (fe '6-24' 'fib' '0' 'genes')

################################################################################

# packages/libraries
################################################################################

library("genefilter")
library("pheatmap")
library("RColorBrewer")
library("AnnotationDbi")
library("OUTRIDER")
library("labeling")
library("farver")
library("stats")
library("BiocParallel")
library("EnsDb.Hsapiens.v79")
library("tidyr")
library("Glimma")
library("ggplot2")
library("ggrepel")
library("readr")

# An OUTRIDER analysis in detail
################################################################################

# Read arguments
# 1: experiment; 2: species; 3: treatment; 4: fragment
args <- commandArgs(trailingOnly = TRUE)

if (length(args)< 4) {
  stop("Fill in the following arguments between quotes: 'experiment' 'species' 'treatment' 'fragment' (fe '6-24' 'fib' '0' 'genes')", call. = FALSE)
}

experiment <- args[1]

species <- args[2]
if (species == 'FIB'){
  species <- 'fib'
} else if (species == 'pax'){
  species <- 'PAX'
}

treatment <- args[3]

fragment <- args[4]
if (fragment == 'gene'){
  fragment <- 'genes'
} else if (fragment == 'exon'){
  fragment <- 'exons'
} else if (fragment == 'intron'){
  fragment <- 'introns'
}

# Define files
countData = paste0("/../countfiles_all_exp/exp",experiment,"_ht2_",fragment,"_counts.tsv")
metaData = paste0("/..metadata/exp",experiment,"_metadata.csv")
countTPM = paste0("/../countfiles_all_exp/exp",experiment,"_ht2_genes_TPM.tsv")

resultfile = paste0("/../res_outr_exp",experiment,"_",species,"_",treatment,"_ht2_",fragment,"_counts.csv")
plotheat = paste0("/../plots/exp",experiment,"_",species,"_",treatment,"_ht2_",fragment,"_counts_heatplots.pdf")


### OutriderDataSet

# Reading in count data (made with current pipeline)
countdata <- read.table(countData,
                        header = TRUE, 
                        stringsAsFactors = FALSE)

# make data integer
countdata[,] <- apply(countdata[,],2,function(x) as.integer(as.character(x)))


# reading sample metadata
metadata <- read.table(metaData,
                       sep = ";", 
                       header = TRUE, 
                       stringsAsFactors = TRUE)
colnames(metadata)[colnames(metadata)=="sample_id"] <- "sampleID"

# make a subset
metadata_subset <- metadata[ which (metadata$species == species & metadata$treatment == treatment),]

# filter only wanted samples in countdata
count_subset <- countdata[,(names(countdata) %in% metadata_subset$sampleID)]

# split rownames in new columns to create column gene
df <- count_subset
df$row = row.names(df)

if(fragment == 'exons') {
df_split <- separate(df, row, into=c("chr","start","end","gene-nm"), sep="_", extra="merge")
df_split <- separate(df_split, "gene-nm", into=c("gene","nm"), sep="-", extra="merge")
df_split1 <- subset(df_split, select = -c(chr, start, end, nm))
df_split <- subset(df_split, select = -c(nm))
df_split$geneID <- rownames(df_split)
} else {
df_split <- separate(df, row, into=c("chr","start","end","gene"), sep="_", extra="merge")
df_split1 <- subset(df_split, select = -c(chr, start, end))
df_split$geneID <- rownames(df_split)  
}

#filter counts based on TPM_genes with 95th percentile > 0.1875
counttpm <- read.table(countTPM,
                       header = TRUE, 
                       stringsAsFactors = FALSE)
# order countdata on column (alphabetically)
counttpm <- counttpm[,order(colnames(counttpm))]
# filter out samples
counttpm_subset <- counttpm[,(names(counttpm) %in% metadata_subset$sampleID)]
# filter only genes with percentile 95% <0.1875
counttpm_subset$perc95 = apply(counttpm_subset, 1, quantile, .95)
counttpm_subset_perc95 <- counttpm_subset[(counttpm_subset$perc95<.1875),]

# split rownames to new columns to create column gene
df_tpm <- counttpm_subset_perc95
df_tpm$row = row.names(df_tpm)
df_tpm_split <- separate(df_tpm, row, into=c("chr","start","end","gene"), sep="_", extra="merge")

# filter countdata for only those which are left in the counttpm table
count_filter <- df_split1[!(df_split1$gene %in% df_tpm_split$gene),]
# remove not necessary columns
count_filter <- subset(count_filter, select = -c(gene))

# extra filter to remove alle samples with 0 counts in all samples
if(fragment != 'genes') {
count_filter <- count_filter + 1
no.samples <- ncol(count_filter)
count_filter <- count_filter[!(rowSums(count_filter==1)==no.samples),]
}

# create an OutriderDataSet object [1]
ods <- OutriderDataSet(countData = count_filter, colData = metadata_subset)

### Preprocessing
###################
#Genes which do not have any expression in all samples are filtered out of the dataset with the function "filterExpression" [1].*
ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)

### Controlling for Confounders
#The samples vary from each other due to technical (like sequencing batches) or biological (like sex, age) differences.

# Heatmap of the sample correlation
# it can also annotate the clusters resulting from the dendrogram

# save all result plots in a pdf file
pdf(plotheat,onefile = TRUE)

# Heatmap of the sample correlation
# it can also annotate the clusters resulting from the dendrogram
ods <- plotCountCorHeatmap(ods, colGroups=c("experiment_GS", "gender_genesis"),
                           normalized=FALSE, nRowCluster=4)

# Heatmap of the gene/sample expression
ods <- plotCountGeneSampleHeatmap(ods, colGroups=c("experiment_GS", "gender_genesis"),
                                  normalized=FALSE, nRowCluster=4)

#These differences can affect the detection of aberrant features, so control them is necessary. To control for confounders present in the data  we first calculate the sizeFactors as done in DESeq2. Additionally, the "controlForConfounders" function is used. This function calls an autoencoder that automatically controls for confounders present in the data. Therefore an encoding dimension q needs to be set.The optimal value of q can be determined using the "findEncodingDim" function.Use this q to control for confounders. The default iterations is 15.[1]*
ods <- estimateSizeFactors(ods)

# calculate optimal q
register(SerialParam())
q <- findEncodingDim(ods, params = seq(2, min(100, ncol(ods) -1, nrow(ods) -1), 2), zScore = 3, BPPARAM = MulticoreParam(3))
opt_q <- q@metadata[["optimalEncDim"]]

# control for confounders (default zScore = 3)
ods <- controlForConfounders(ods, q = opt_q, iterations=15, BPPARAM = MulticoreParam(3))
ods <- fit(ods,BPPARAM = MulticoreParam(3))

#After controlling for confounders the heatmap should be plotted again. If it worked, no batches should be present and the correlations between samples should be reduced and close to zero. [1]*
ods <- plotCountCorHeatmap(ods, colGroups=c("experiment_GS", "gender_genesis"),
                           normalized=TRUE)

dev.off()

### Excluding samples from the autoencoder fit
#Outrider expects that each sample within the population is independent of all others. 
#It is recommended to exclude the possible replicates. In our dataset no replicates are used, 
#so excluding of samples is not neccessary.

### P-value calculation
#Two-sided p-values are computed. Multiple testing correction is performed using Benjamini-Yekutielis False Discovery Rate method [1].*
ods <- computePvalues(ods, alternative="two.sided", method="BY", BPPARAM=MulticoreParam(3))

### Z-score calculation
#Z-scores are computed on the log transformed counts [1].
ods <- computeZscores(ods)

## Results
#if there are no significant events, use all=TRUE to show all events (sort by pvalue).
res <- results(ods, all=TRUE)

#add sampleID, chr, start, end and gene columns to the results
res_merge <- merge(res, metadata_subset[,c(1,2,13,14)], by = "sampleID", all.res = TRUE)
split <- subset(df_split[(ncol(df_split)-5):ncol(df_split)])
res_merge <- merge(res_merge, split, by = "geneID", all.res = TRUE)

### Exporting results
# save results table in a CSV file
write_csv(res_merge, resultfile, append = FALSE, col_names = TRUE)
