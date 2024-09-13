#!/usr/bin/env Rscript
# Code to generate exon ratio files (FPKM exon / FPKM gene)*100 for each provided sample/file

# Import statements, alphabetic order of main package.
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("OUTRIDER"))
suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))

# Argument parser
parser <- ArgumentParser(description = "Process some integers")
parser$add_argument("refexon", metavar = "ref exon folder", nargs = "+",
                    help = "Files or directories containing the exon input files.")
parser$add_argument("refgene", metavar = "ref gene folder", nargs = "+",
                    help = "Files or directories containing the gene input files.")
parser$add_argument("-o","--output_path", metavar = "output_path", nargs = "+",
                    help = "Path where output of OUTRIDER will be stored.", default="./")
args <- parser$parse_args()


read_input_files <- function(input){
  ##read input files 
  sampleIDs <- c()
  count_tables <- lapply(input, function(f) {
    input_ext <- tools::file_ext(f)
    if(input_ext == "Rds") {
      rse <- readRDS(f)  # RangedSummarizedExperiment
      sampleIDs <- append(sampleIDs, colnames(rse))
      return(as_tibble(rownames_to_column(as.data.frame(assays(rse)$counts), var="rownames")))
    } else if (input_ext %in% c("txt", "tsv")) {
      count_table <- read_delim(f, show_col_types=FALSE, skip=1)
      # col 1: samples IDs, col 9: counts
      ct <- count_table[,c(1,9)]
      names(ct)[1] <- "rownames"
      sampleIDs <<- append(sampleIDs, colnames(ct)[2])
      return(ct)
    } else {
      stop("Input file extension is not supported.")
    }
  })
  return(list("sampleIDs"=sampleIDs, "count_tables"=count_tables))
}


get_input <- function(input){
  # If directories are provided
  if(all(sapply(input, function(x) dir.exists(x)))) {
    retrieved_files <- sapply(input, function(d){
      list.files(path = d, pattern = "Rds|txt|tsv", full.names = TRUE, recursive = FALSE)
    })
  } else if(all(sapply(input, function(x) file.exists(x)))) {  # If files are provided.
    retrieved_files <- input
  } else {  # If both dir and files are provided, or different tyes (character, int etc)
    stop("Input is neither dir or file.")
  }
  
  return(read_input_files(retrieved_files))
}


merge_count_tables <- function(lst_ref){
  # merge count tables together.
  all_counts <- lst_ref %>%
    Reduce(function(dtf1,dtf2) dplyr::full_join(dtf1, dtf2,by="rownames"), .)
}


get_pc_exons <- function(all_counts, refexon){
  ##GET ALL featurecounts exon-level count data in a matrix and keep only exons within exon-level filtering boundaries
  all_counts_matrix <- as.matrix(all_counts)[,-1]
  mode(all_counts_matrix) <- "integer"
  rownames(all_counts_matrix) <- all_counts$rownames
  ods <- OutriderDataSet(countData=all_counts_matrix)
  
  ##STANDARD OUTRIDER FILTERING FOR EXON LEVEL
  ct <- read_delim(paste0(refexon,sub(".bam", "", colnames(all_counts)[2]),".cntrl.exon.txt"), show_col_types=FALSE, skip=1)
  print("identical exon ids:")
  print(identical(ct$Geneid,rownames(assays(ods)$counts)))
  readD <- apply(assays(ods)$counts, 2, function(x) x / sum(x) * 10^6)
  countsFPKM <- readD / ct$Length * 10^3
  perc95e <- apply(countsFPKM, 1, function(x) quantile(x,probs=0.95))
  keep <- perc95e>1
  #mcols(ods)$basepairs <- ct$Length
  ##KEEP only exons within filtering boundaries
  keepexons <- keep[keep==TRUE]
  #length(keepexons)
  return(keepexons)
}


print_ratiofpkm_exons <- function(keepexons, sample_ids, refexon, refgene, out_path){
  ## calculate ratio FPKM exon / FPKM gene per sample file, based on exons to keep 
  dir.create(file.path(out_path))
  lapply(sample_ids, function(sample_id){
    sample <- sub(".bam", "", sample_id)
    #    urlg <- paste("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/downloads/protein_coding/gene_pc_all/",sample,".","cntrl.gene.txt",sep="")
    #    urle <- paste("/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/downloads/protein_coding/exon_pc_all/",sample,".","cntrl.exon.txt",sep="")
    urlg <- paste(refgene,sample,".","cntrl.gene.txt",sep="")
    urle <- paste(refexon,sample,".","cntrl.exon.txt",sep="")
    ctg <- read_delim(urlg, skip=1,show_col_types = FALSE)
    cte <- read_delim(urle, skip=1,show_col_types = FALSE)
    
    colnames(ctg)[c(1,9)] <- c('gene_id','counts')
    reads <- ctg$counts/ sum(ctg$counts) * 10^6
    ctg$FPKMg <- reads / (ctg$Length) * 10^3
    ctg <- ctg[,c(1,10)]
    
    colnames(cte)[9] <- 'counts'
    reads <- cte$counts/ sum(cte$counts) * 10^6
    cte$FPKMe <- reads / (cte$Length) * 10^3
    ct_merge <- merge(cte, ctg, by = "gene_id", all.res = TRUE)
    ct_merge$ratio <- round(ct_merge$FPKMe / ct_merge$FPKMg * 100 , 0)
    
    ct_merge$ratio[is.na(ct_merge$ratio)]<-0
    ct_merge$ratio[is.infinite(ct_merge$ratio)]<-round(ct_merge$FPKMe[is.infinite(ct_merge$ratio)]/mean(ctg$FPKMg) * 100 , 0)
    
    ct_merge <- ct_merge[,c(2:7,1,8,12)]
    colnames(ct_merge)[9] <- sample_id #paste(sample,".bam",sep="")
    ct_merge_expr <- ct_merge[ct_merge$Geneid %in% rownames(as.data.frame(keepexons)),]
    cat("# ratio fpkm file ",sample,"\n", file=paste0(out_path,sample,".cntrl.ratiofpkm.exon.txt"))
    write_tsv(ct_merge_expr, paste0(out_path,sample,".cntrl.ratiofpkm.exon.txt"), append=TRUE, col_names = T)
    print(paste0(sample," is added"))
  })
}

main <- function(refexon, refgene, output_path){
  ##read exon counts input files - apply standard exon level OUTRIDER filtering - calculate EXON ratio output per sample and print exon ratio file
  ref_data_exon <- get_input(refexon)
  all_counts <- merge_count_tables(ref_data_exon$count_tables)
  keepexons <- get_pc_exons(all_counts, refexon)
  print_ratiofpkm_exons(keepexons, ref_data_exon$sampleIDs, refexon, refgene, output_path)
}

main(args$refexon, args$refgene, args$output_path)
