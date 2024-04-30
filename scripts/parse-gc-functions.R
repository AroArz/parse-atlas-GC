#!/usr/bin/env Rscript

#########################
# CREATE DIRS FUNCTIONS #
#########################
check_dirs <- function(subdir_path) {
    
        if (!dir.exists(subdir_path)) {
        dir.create(subdir_path, recursive = TRUE)
        put(paste("Directory", subdir_path, "created."))
    } else {
        put(paste("Directory", subdir_path, "already exists."))
    }
}

# tutorial from Silas
# https://github.com/metagenome-atlas/Tutorial/blob/master/R/Analyze_genecatalog.Rmd
# Some helper functions to convert between gene names and gene numbers

##########################
# PROCESS ANNO FUNCTIONS #
##########################

GeneNr2GeneName <- function(GeneNumbers, MaxGeneNumber = Ngenes) {
  paste0("Gene", formatC(format = "d", GeneNumbers, flag = "0", width = ceiling(log10(MaxGeneNumber))))
}

GeneName2GeneNr <- function(GeneNames) {
  as.integer(str_sub(GeneNames, start = 5))
}


getQueries <- function(genecatalog, column, query) {
    
    GC_query_genes = genecatalog %>% subset(grepl(query, genecatalog[[column]]))
    GC_query_genelist = unique(GC_query_genes[["Query"]]) %>% GeneName2GeneNr()
    
    return(GC_query_genelist)
}


# helper function to load subset of genes from hdf5 file
load_subset_of_genes <- function(abundance_file, indexes_of_genes_to_load) {
  put("Load ", length(indexes_of_genes_to_load), " genes\n")


  data <- rhdf5::h5read(
    file = abundance_file, name = "data",
    index = list(indexes_of_genes_to_load, NULL)
  )

  # add sample names
  attributes <- rhdf5::h5readAttributes(abundance_file, "data")
  colnames(data) <- attributes$sample_names
  rownames(data) <- GeneNr2GeneName(indexes_of_genes_to_load)

  return(data)
}

normalize_subset <- function(subset_count, total_coverage) {
    
    subset_cpm_norm <- subset_count %*% diag(1 / total_coverage[colnames(subset_count)]) * 1e6
    colnames(subset_cpm_norm) <- colnames(subset_count)
    
    return(subset_cpm_norm)
}

group_anno_by_sample <- function(normalized_cpm, anno) {
    
    normalized_cpm_melted = normalized_cpm %>% 
        reshape2::melt() %>%
        dplyr::rename("Gene" = "Var1", "Sample" = "Var2") %>%
        group_by(Sample) %>%
        dplyr::summarise(total_cpm_norm = sum(value))
    
    normalized_cpm_melted$annotation = anno
    
    return(normalized_cpm_melted)
}

process_anno <- function(genecatalog, abund_file_path, total_coverage, column, query) {
    
    put("getting query genelist\n")
    query_genelist = getQueries(genecatalog, column, query)
    
    put("loading cpm data for query\n")
    query_count_subset = load_subset_of_genes(abund_file_path, query_genelist)
    
    put("normalizing cpm data for query\n")
    query_cpm_subset_norm = normalize_subset(query_count_subset, total_coverage)
    
    put("grouping normalized cpm data for query per sample\n")
    query_cpm_norm_grouped = group_anno_by_sample(query_cpm_subset_norm, query)
    
    return(query_cpm_norm_grouped)
    
}

##################
# GET ANNO LISTS #
##################

extractUniqueColValues <- function(df, column_name) {
  # Ensure the column exists in the dataframe
  if (!(column_name %in% names(df))) {
    stop("Column name does not exist in the dataframe")
  }
  
  # Split the column on commas and unlist into a single vector
  split_values <- strsplit(as.character(df[[column_name]]), ",")
  unique_values <- unique(unlist(split_values))
  
  # Return a dataframe of unique values
  return(data.frame(unique_values))
}

filterByValue <- function(df, column_name, search_string) {
  # Ensure the column exists in the dataframe
  if (!(column_name %in% names(df))) {
    stop("Column name does not exist in the dataframe")
  }
  
  # Filter the dataframe based on the presence of the search string in the column
  subset_df <- df[grepl(paste0("\\b", search_string, "\\b"), df[[column_name]]), ]
  
  return(subset_df)
    
}

get_kegg_data <- function(category) {
  # gets KEGG list data
  valid_categories <- c("pathway", "module", "enzyme", "disease", "drug", "compound")
  if (!(category %in% valid_categories)) {
    stop("Invalid category. Please choose from: ", paste(valid_categories, collapse=", "))
  }
  
  url <- paste0("https://rest.kegg.jp/list/", category)
  
  response <- GET(url)
  
  if (status_code(response) == 200) {
    # Read the content of the response
    data_text <- content(response, "text")
    
    # Convert the text data to a dataframe
    data_df <- read_delim(data_text, delim = "\t", col_names = c("KEGG_ID", "Description"), show_col_types = FALSE)
    
    return(data_df)
  } else {
    stop("Failed to retrieve data: HTTP status", status_code(response))
  }
}
