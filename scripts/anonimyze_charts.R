##Anonymize patient ids in OUTRIDER AND FRASER output files to generate top 20 lists per analysis
##Code to generate plots of output data, regarding Final Master Project 

library(readr)
library(tibble)
library(plyr)
library(ggplot2)
library(pals)

##OUTRIDER output GENE LEVEL untreated samples
urlg <- "/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/downloads/protein_coding/finalrefset_cntrl/umcu_rnaseq_fib_untreated_res_outrider_gene_counts_pheno.tsv"
umc_res_g <- read_delim(urlg)

umc_res_g$sampleID <- as.factor(umc_res_g$sampleID)
levels(umc_res_g$sampleID)

old_levels <- levels(umc_res_g$sampleID)
replace_levels <- c("Individual_02", "Individual_08","Individual_14","Individual_10","Individual_12",
                    "Individual_15","Individual_17","Individual_04","Individual_13","Individual_18",
                    "Individual_09","Individual_11","Individual_22","Individual_23","Individual_05",
                    "Individual_16","Individual_19","Individual_07","Individual_20","Individual_03",
                    "Individual_01","Individual_24","Individual_21","Individual_06","Individual_25")
replace_levels
class(umc_res_g)
#umc_res_g$indiv <-  mapvalues(umc_res_g$sampleID, from=levels(umc_res_g$sampleID), to=replace_levels)
umc_res_g$indiv <-  mapvalues(umc_res_g$sampleID, from=old_levels, to=replace_levels)

#write anonymous tsv gene
genes_anonym <- umc_res_g
genes_anonym$sampleID <- paste0(genes_anonym$indiv,"_untreated")
write_csv(genes_anonym[c(1:23)], paste0("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/","umcu_rnaseq_fib_untreated_res_outrider_genes_counts_anonymous.tsv"), append=FALSE, col_names = TRUE)
#write anonymous csv gene top 100
top100 <- genes_anonym[order(genes_anonym$pValue),]
write_csv(top100[c(1:100),c(1:23)], paste0("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/top100/","umcu_rnaseq_fib_untreated_res_outrider_genes_counts_anonymous_top100.csv"), append=FALSE, col_names = TRUE)


# add a column for outliers or not
umc_res_g$Outlier <- umc_res_g$indiv
umc_res_g$Outlier[abs(umc_res_g$zScore) < 2.5 | umc_res_g$pValue > 0.01] <- "Outside boundaries"
p <- ggplot(data=umc_res_g, aes(x=zScore, y=-log10(pValue), col=Outlier)) + geom_point() + theme_minimal()
# Add lines 
p2 <- p + geom_vline(xintercept=c(-2.5, 2.5), col="#bf1932") +
  geom_hline(yintercept=-log10(0.01), col="#bf1932")
#Add colors
p3 <- p2 + scale_color_manual(values=c(unname(polychrome(25)), "#CCC")) + 
  ggtitle("Vulcanoplot gene expression untreated samples") + 
  guides(colour = guide_legend(ncol = 1, )) +
  theme(plot.title = element_text(hjust = 0.5, size = 26, margin = margin(0,0,30,0)),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        )
p3




####OUTRIDER EXON RATIO
urle <- "/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/downloads/protein_coding/finalrefset_ratiofpkm260624/umcu_rnaseq_fib_untreated_res_outrider_exons_counts_countsadded_pheno_added.tsv"
umc_res_e <- read_delim(urle)

umc_res_e$sampleID <- as.factor(umc_res_e$sampleID)
levels(umc_res_e$sampleID)

old_levels <- levels(umc_res_e$sampleID)
old_levels
replace_levels <- c("Individual_02", "Individual_08","Individual_14","Individual_10","Individual_12",
                    "Individual_15","Individual_17","Individual_04","Individual_13","Individual_18",
                    "Individual_09","Individual_11","Individual_22","Individual_23","Individual_05",
                    "Individual_16","Individual_19","Individual_07","Individual_20","Individual_03",
                    "Individual_01","Individual_24","Individual_21","Individual_06","Individual_25")
replace_levels
umc_res_e$indiv <-  mapvalues(umc_res_e$sampleID, from=old_levels, to=replace_levels)
dim(umc_res_e)


exons_anonym <- umc_res_e
exons_anonym$sampleID <- paste0(exons_anonym$indiv,"_untreated")
write_csv(exons_anonym[c(1:23)], paste0("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/","umcu_rnaseq_fib_untreated_res_outrider_exons_counts_ratio_anonymous.tsv"), append=FALSE, col_names = TRUE)
#write anonymous csv gene top 100
top100 <- exons_anonym[order(exons_anonym$pValue),]
head(top100[,c(1:5)])
write_csv(top100[c(1:100),c(1:23)], paste0("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/top100/","umcu_rnaseq_fib_untreated_res_outrider_exons_counts_ratio_anonymous_top100.csv"), append=FALSE, col_names = TRUE)



# add a column of outlier or not
umc_res_e$Outlier <- umc_res_e$indiv
umc_res_e$Outlier[umc_res_e$zScore > -2.5 | umc_res_e$pValue > 0.01] <- "Outside boundaries"
#umc_res_g$diffexpressed[abs(umc_res_g$zScore) > 2.5 & umc_res_g$pValue < 0.01] <- "DIFF"
p <- ggplot(data=umc_res_e, aes(x=zScore, y=-log10(pValue), col=Outlier)) + geom_point() + theme_minimal()
#p
# Add lines 
p2 <- p + geom_vline(xintercept=c(2.5), col="#bf1932") +
  geom_hline(yintercept=-log10(0.01), col="#bf1932")
#Add colors
# 
p3 <- p2 + scale_color_manual(values=c(unname(polychrome(25)), "#CCC")) + 
  ggtitle("Vulcanoplot exon ratio") + 
  guides(colour = guide_legend(ncol = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, margin = margin(0,0,30,0)))
p3



##filtering STEPS 1 - 4 for EXON RATIO results

dim(umc_res_g)
umc_res_g$genesample <- paste(umc_res_g$gene,umc_res_g$sampleID,sep="-")
genes_keep <- umc_res_g[umc_res_g$normcounts>200,]
dim(genes_keep)
##filtering exon ratio
#keep only those genes with multiple exons per gene
dim(umc_res_e)
umc_res_e$genesample <- paste(umc_res_e$gene,umc_res_e$sampleID,sep="-")

#duplicate_exonsexist <- unique(umc_res_e$gene[duplicated(umc_res_e$genesample)])
#keep only those genes with 1 exon per gene with a zScore < -2
exonsZ2 <- umc_res_e[umc_res_e$zScore < -2,]
duplicates <- unique(exonsZ2$genesample[duplicated(exonsZ2$genesample)])

#STEP 1
exonskeep <- umc_res_e[umc_res_e$genesample %in% genes_keep$genesample,]
dim(exonskeep)
#STEP 2
#exonskeep <- umc_res_e[!umc_res_e$genesample %in% duplicates & umc_res_e$gene %in% duplicate_exonsexist & umc_res_e$genesample %in% genes_keep2$genesample,]
exonskeep <- umc_res_e[!umc_res_e$genesample %in% duplicates & umc_res_e$genesample %in% genes_keep$genesample,]
dim(exonskeep)

## Volcano Plot step 3
p <- ggplot(data=exonskeep, aes(x=zScore, y=-log10(pValue), col=Outlier)) + geom_point() + theme_minimal()
#p
# Add lines 
p2 <- p + geom_vline(xintercept=c(-2.5), col="#bf1932") +
  geom_hline(yintercept=-log10(0.01), col="#bf1932")
#Add colors
# 
p3 <- p2 + scale_color_manual(values=c(unname(polychrome(25)), "#CCC")) + 
  ggtitle("Vulcanoplot exon ratio - step 3") + 
  guides(colour = guide_legend(ncol = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 26, margin = margin(0,0,30,0)),
      axis.title = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 12)
)
p3

#STEP 3
exonskeep2 <- exonskeep[exonskeep$zScore < -2.5 & exonskeep$pValue<0.01, ]
dim(exonskeep2)
head(exonskeep2$Outlier)

##Plot step 4
exonskeep2$Outlier[exonskeep2$meanCorrected < 100] <- "Outside boundaries"
p <- ggplot(data=exonskeep2, aes(x=log10(meanCorrected), y=-log10(pValue), col=Outlier)) + geom_point() + theme_minimal()
#p
# Add lines 
p2 <- p + geom_vline(xintercept=c(2), col="#bf1932") 
#Add colors
p3 <- p2 + scale_color_manual(values=c(unname(polychrome(25)), "#CCC")) + 
  ggtitle("Mean exon expression ratio - step 4") + 
  guides(colour = guide_legend(ncol = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 26, margin = margin(0,0,30,0)),
      axis.title = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 12)
  )
p3

#STEP 4
exons_ratiofinal <- exonskeep2[exonskeep2$meanCorrected>100,]
dim(exons_ratiofinal)


exons_anonym <- exons_ratiofinal
exons_anonym$sampleID <- paste0(exons_anonym$indiv,"_untreated")
#write_csv(exons_anonym[c(1:23)], paste0("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/","umcu_rnaseq_fib_untreated_res_outrider_exons_counts_ratio_anonymous.tsv"), append=FALSE, col_names = TRUE)
#write anonymous csv gene top 100
top100 <- exons_anonym[order(exons_anonym$pValue),]
head(top100[,c(1:5)])
write_csv(top100[c(1:100),c(1:23)], paste0("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/top100/","umcu_rnaseq_fib_untreated_res_outrider_exons_counts_ratiostep4_anonymous_top100.csv"), append=FALSE, col_names = TRUE)



## OUTRIDER output gene level - CHX treated samples

urlgchx <- "/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/downloads/protein_coding/finalrefset_chx/umcu_rnaseq_fib_chx_res_outrider_gene_counts.tsv"
umc_res_gchx <- read_delim(urlgchx)

dim(umc_res_gchx)
umc_res_gchx$sampleID <- as.factor(umc_res_gchx$sampleID)
levels(umc_res_gchx$sampleID)

old_levels <- levels(umc_res_gchx$sampleID)
old_levels
replace_levels <- c("Individual_02", "Individual_08","Individual_14","Individual_10","Individual_12",
                    "Individual_15","Individual_17","Individual_04","Individual_13","Individual_18",
                    "Individual_09","Individual_11","Individual_22","Individual_23","Individual_05",
                    "Individual_16","Individual_19","Individual_07","Individual_20","Individual_03",
                    "Individual_01","Individual_24","Individual_21","Individual_06","Individual_25")

replace_levels
#umc_res_gchx$indiv <- replace(umc_res_gchx$sampleID, umc_res_gchx$sampleID %in% old_levels, replace_levels)
umc_res_gchx$indiv <-  mapvalues(umc_res_gchx$sampleID, from=levels(umc_res_gchx$sampleID), to=replace_levels)

head(umc_res_gchx$indiv)
genes_anonym <- umc_res_gchx
genes_anonym$sampleID <- paste0(genes_anonym$indiv,"_chx")
genes_anonym$phenotype <- "outside boundaries"
head(genes_anonym[c(1:3,24)])
write_csv(genes_anonym[c(1:22,24)], paste0("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/","umcu_rnaseq_fib_chx_res_outrider_genes_counts_anonymous.tsv"), append=FALSE, col_names = TRUE)



## OUTRIDER output exon level - untreated samples

urlex <- "/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/umcu_rnaseq_fib_untreated_res_outrider_exons_counts_finalexon.tsv"
umc_res_ex <- read_delim(urlex)

dim(umc_res_ex)
umc_res_ex$sampleID <- as.factor(umc_res_ex$sampleID)
levels(umc_res_ex$sampleID)

old_levels <- levels(umc_res_ex$sampleID)
old_levels
replace_levels <- c("Individual_02", "Individual_08","Individual_14","Individual_10","Individual_12",
                    "Individual_15","Individual_17","Individual_04","Individual_13","Individual_18",
                    "Individual_09","Individual_11","Individual_22","Individual_23","Individual_05",
                    "Individual_16","Individual_19","Individual_07","Individual_20","Individual_03",
                    "Individual_01","Individual_24","Individual_21","Individual_06","Individual_25")

umc_res_ex$indiv <-  mapvalues(umc_res_ex$sampleID, from=levels(umc_res_ex$sampleID), to=replace_levels)

head(umc_res_ex$indiv)
genes_anonym <- umc_res_ex
genes_anonym$sampleID <- paste0(genes_anonym$indiv,"_untreated")
head(genes_anonym[c(1:3,23)])
write_csv(genes_anonym[c(1:23)], paste0("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/","umcu_rnaseq_fib_untreated_res_outrider_exons_counts_anonymous.tsv"), append=FALSE, col_names = TRUE)

top100 <- genes_anonym[order(genes_anonym$pValue),]
write_csv(top100[c(1:100),c(1:23)], paste0("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/top100/","umcu_rnaseq_fib_untreated_res_outrider_exons_counts_anonymous_top100.csv"), append=FALSE, col_names = TRUE)



#OUTPUT data FRASER untreated samples - anonimize and generate top 100 list and Vulcano plot
urlf <- "/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/fraser/300624/fraser_result_sample2.csv"
umc_res_f <- read_delim(urlf)

umc_res_f <- read_delim(urlf,show_col_types = TRUE,delim="\t")
problems(umc_res_f)
urlf <- "/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/top100/fraser_result_top100.csv"
umc_res_f <- read_delim(urlf)
dim(umc_res_f)

umc_res_f$sampleID <- as.factor(umc_res_f$sampleID)
levels(umc_res_f$sampleID)

old_levels <- levels(umc_res_f$sampleID)
old_levels
replace_levels <- c("Individual_02", "Individual_08","Individual_14","Individual_10","Individual_12",
                    "Individual_15","Individual_17","Individual_04","Individual_13","Individual_18",
                    "Individual_09","Individual_11","Individual_22","Individual_23","Individual_05",
                    "Individual_16","Individual_19","Individual_07","Individual_20","Individual_03",
                    "Individual_01","Individual_24","Individual_21","Individual_06","Individual_25")
replace_levels
umc_res_f$indiv <-  mapvalues(umc_res_f$sampleID, from=old_levels, to=replace_levels)

genes_anonym <- umc_res_f
genes_anonym$sampleID <- paste0(genes_anonym$indiv,"")

#write_csv(genes_anonym[,c(1:21)], paste0("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/top100/","fraser_results_anonymous.csv"), append=FALSE, col_names = TRUE)


# add a column for outliers or not
dim(umc_res_f)
umc_res_f$Outlier <- umc_res_f$indiv
umc_res_f$Outlier[abs(umc_res_f$deltaPsi) < 0.05 | umc_res_f$pValue > 0.001] <- "Outside boundaries"
dim(umc_res_f[umc_res_f$Outlier!='Outside boundaries',])
p <- ggplot(data=umc_res_f, aes(x=deltaPsi, y=-log10(pValue), col=Outlier)) + geom_point() + theme_minimal()
# Add lines 
p2 <- p + geom_vline(xintercept=c(-0.05, 0.05), col="#bf1932") +
  geom_hline(yintercept=-log10(0.001), col="#bf1932")
#Add colors
p3 <- p2 + scale_color_manual(values=c(unname(polychrome(25)), "#CCC")) + 
  ggtitle("Vulcanoplot splicing FRASER untreated samples") + 
  guides(colour = guide_legend(ncol = 1, )) +
  theme(plot.title = element_text(hjust = 0.5, size = 26, margin = margin(0,0,30,0)),
        axis.title = element_text(size = 18, margin = margin(10,10,0,0)),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
  )
p3

