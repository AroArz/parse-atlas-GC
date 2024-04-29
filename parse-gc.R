#!/usr/bin/env Rscript

library("argparse")
library("rhdf5")
library("httr")
library("readr")
library("arrow")
library("dplyr")
library("logr")

#######################
# SETTING UP ARGPARSER #
########################

parser <- ArgumentParser(description = "minimalistic script to parse atlas GC per KEGG pathway, KEGG module, GO terms\n 
					Normalize to cpm \n
					Aggregate cpm per sample")

parser$add_argument("--input", type = "character", required = TRUE, 
                    help = "Path to atlas directory")

# aarse the command line arguments
args <- parser$parse_args()

# access the input argument
input_arg <- args$input
metadata_arg <- args$metadata
tax_arg <- args$tax


#########################################
# LOADING LIBRARIES & SETTING UP LOGGER #
#########################################

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
logfile <- paste0("logfile_", timestamp, ".log")
lf = log_open(logfile)

options("logr.compact" = TRUE)
options("logr.notes" = FALSE)
put(paste("user provided input: ", input_arg, sep = ""))
put(paste("user provided metadata: ", metadata_arg, sep = ""))
put(paste("user provided taxtable: ", tax_arg, sep = ""))
options("logr.notes" = TRUE)

#####################
# LOADING FUNCTIONS #
#####################

sep("loading functions")

source("scripts/parse-gc-functions.R")

GC_output_path = "output/GC/"
check_dirs(GC_output_path)

## Read Input Files
abundance_file_path = paste(input_arg, "/Genecatalog/counts/median_coverage.h5")
eggnogGC_path = paste(input_arg, "/Genecatalog/annotations/eggNOG.parquet")
gene_coverage_stats_path = paste(input_arg, "/Genecatalog/counts/gene_coverage_stats.parquet")
sample_coverage_stats_path = paste(input_arg, "/Genecatalog/counts/sample_coverage_stats.tsv")

eggnogGC = read_parquet(eggnogGC_path)
gene_coverage_stats = read_parquet(gene_coverage_stats_path)
sample_coverage_stats = read.csv(sample_coverage_stats_path, sep = "\t", row.names = "X")


