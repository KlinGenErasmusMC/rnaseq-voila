# Code to replace OUTRIDER exon-ratio results with corresponding exon-level results for raw counts - normalized counts and meanCorrected
# OUTRIDER exon ratio results are filtered (STEP 1 - 4) based on the exon level counts results

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))

parser <- ArgumentParser(description = "Process some integers")
parser$add_argument("ratioresult", metavar = "ratio fpkm outrider result file", nargs = "+",
                    help = "The ratio fpkm outrider result file.")
parser$add_argument("exonresult", metavar = "exon outrider result file", nargs = "+",
                    help = "The exon outrider result file.")
parser$add_argument("-o","--output_path", metavar = "output_path", nargs = "+",
                    help = "Path where output file with replaced normalized and mean correct values.", default="./")
args <- parser$parse_args()

# Get the exon ratio results file
umc_res_fpkm <- read_delim(args$ratioresult,show_col_types = FALSE)

# Get the exon level results file
umc_res_e <- read_delim(args$exonresult,show_col_types = FALSE)

# Get the correct columns and merge the files by geneID and sampleID
umc_res_e <- umc_res_e[,c(1,2,7:9)]
umc_res_merge <- right_join(umc_res_e, umc_res_fpkm, by=c("geneID","sampleID"))
# SELECT the correct columns to create exon ratio output results with correct exon-level counts
res_merge_out <- umc_res_merge[,c(1:9,13:25)]
colnames(res_merge_out)[3:5] <- c("rawcounts","normcounts","meanCorrected")

# write the ratio results file with corresponding exon-level raw counts - normalized counts and meanCorrected results
write_csv(res_merge_out, paste0(args$output_path,"umcu_rnaseq_fib_untreated_res_outrider_exons_counts_countsadded.tsv"), append=FALSE, col_names = TRUE)
