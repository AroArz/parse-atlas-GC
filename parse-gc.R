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

parser <- ArgumentParser(description = "attempt for a minimalistic script to rarefy otu data")
parser$add_argument("--input", type = "character", required = TRUE, 
                    help = "Path OTU Count Matrix. TAXA ARE ROWS")
parser$add_argument("--metadata", type = "character", required = TRUE, 
                    help = "Path metadata")
parser$add_argument("--tax", type = "character", required = TRUE, 
                    help = "Path Tax table")

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

