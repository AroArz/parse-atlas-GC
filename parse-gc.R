#!/usr/bin/env Rscript

library("argparse")
library("rhdf5")
library("httr")
library("readr")
library("arrow")
library("dplyr")
library("stringr")
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


#########################################
# LOADING LIBRARIES & SETTING UP LOGGER #
#########################################

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
logfile <- paste0("logfile_", timestamp, ".log")
lf = log_open(logfile)

options("logr.compact" = TRUE)
options("logr.notes" = FALSE)
put(paste("user provided input: ", input_arg, sep = ""))
options("logr.notes" = TRUE)

#####################
# LOADING FUNCTIONS #
#####################

sep("loading functions")

source("scripts/parse-gc-functions.R")

GC_output_path = "output/GC/"
KEGG_pathway_output_path = "output/GC/KEGG_pathways"
KEGG_module_output_path = "output/GC/KEGG_modules"
GO_terms_output_path = "output/GC/GO_terms"

check_dirs(GC_output_path)
check_dirs(KEGG_pathway_output_path)
check_dirs(KEGG_module_output_path)
check_dirs(GO_terms_output_path)

sep("loading input files")

## Read Input Files
abundance_file_path = paste(input_arg, "/Genecatalog/counts/median_coverage.h5", sep = "")
eggnogGC_path = paste(input_arg, "/Genecatalog/annotations/eggNOG.parquet", sep = "")
gene_coverage_stats_path = paste(input_arg, "/Genecatalog/counts/gene_coverage_stats.parquet", sep = "")
sample_coverage_stats_path = paste(input_arg, "/Genecatalog/counts/sample_coverage_stats.tsv", sep = "")

eggnogGC = read_parquet(eggnogGC_path)
gene_coverage_stats = read_parquet(gene_coverage_stats_path)
sample_coverage_stats = read.csv(sample_coverage_stats_path, sep = "\t", row.names = "X")


####################
# GETTING UNIQUE PATHWAYS; MODULES; GO'S
##################


sep("getting unique pathways, modules, GO's etc ...")
## Get unique pathways, modules, GOs etc...
eggnogGC_KEGGpwy = extractUniqueColValues(eggnogGC, "KEGG_Pathway")
eggnogGC_KEGGpwy = eggnogGC_KEGGpwy %>% subset((grepl("map", eggnogGC_KEGGpwy$unique_values)))

write.csv(eggnogGC_KEGGpwy, 
          paste(KEGG_pathway_output_path, "/eggnogGC_KEGGpwy.csv", sep = ""))

#eggnogGC_KEGGmodule = extractUniqueColValues(eggnogGC, "KEGG_Module")
#eggnogGC_KEGGmodule = eggnogGC_KEGGmodule %>% subset(!(unique_values %in% c("KEGG_Module", "-")))

#write.csv(eggnogGC_KEGGmodule, 
#          paste(KEGG_module_output_path, "/eggnogGC_KEGGmodule.csv", sep = ""))

#eggnogGC_GO = extractUniqueColValues(eggnogGC, "GO_terms")
#eggnogGC_GO = eggnogGC_GO %>% subset(!(unique_values %in% c("GOs", "-")))

#write.csv(eggnogGC_GO, 
#          paste(GO_terms_output_path, "/eggnogGC_GO.csv", sep = ""))


sep("setting up variables ...")
## setting up variables
h5overview <- rhdf5::h5ls(abundance_file_path)
dim <- h5overview[1, "dim"] %>%
  stringr::str_split(" x ", simplify = T) %>%
  as.numeric()

Ngenes <- dim[1]
Nsamples <- dim[2]

total_coverage <- sample_coverage_stats[, "Sum_coverage"]
names(total_coverage) <- rownames(sample_coverage_stats)


sep("generating KEGG pathway data")
for (ID in unique(eggnogGC_KEGGpwy$unique_values)) {
    
    
    
    cat_msg = paste("summarising KEGG PATHWAY: ", ID, sep = "")
    put(cat_msg)
    
    ID_summarised = process_anno(eggnogGC, abundance_file_path, total_coverage, "KEGG_Pathway", paste(ID))
    
    if (!is.null(ID_summarised)) {
        
        cat_msg = paste("saving: ", ID, sep = "")
        put(cat_msg)

        write.csv(ID_summarised,
                  paste(KEGG_pathway_output_path, "/", ID, ".csv", sep = "")
                 )
        
    } else {
        
        put(paste("No data to save for: ", ID, "\n"))
        
    }
    
}

#sep("generating KEGG module data")
#for (ID in unique(eggnogGC_KEGGmodule$unique_values)) {
#    
#    
#    
#    cat_msg = paste("summarising KEGG MODULE: ", ID, sep = "")
#    put(cat_msg)
#    
#    ID_summarised = process_anno(eggnogGC, abundance_file_path, total_coverage, "KEGG_Module", paste(ID))
#    
#    cat_msg = paste("saving: ", ID, sep = "")
#    put(cat_msg)
#    
#    write.csv(ID_summarised,
#              paste(KEGG_module_output_path, "/", ID, ".csv", sep = ""))
#    
#}
#
#sep("generating GO TERMS data")
#for (ID in unique(eggnogGC_GO$unique_values)) {
#    
#    
#    
#    cat_msg = paste("summarising GO TERM: ", ID, sep = "")
#    put(cat_msg)
#    
#    ID_summarised = process_anno(eggnogGC, abundance_file_path, total_coverage, "GO_terms", paste(ID))
#    
#    cat_msg = paste("saving: ", ID, sep = "")
#    put(cat_msg)
#    
#    write.csv(ID_summarised,
#              paste(GO_terms_output_path, "/", ID, ".csv", sep = ""))
#    
#}