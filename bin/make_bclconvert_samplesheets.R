#! /usr/bin/env Rscript

# Sam Lee
# @samleenz

library(tidyverse)

### commandline arguments ------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("
  Usage: make_bcl_sheet.R \\
    <sampletable.csv> \\
    <single_template.txt> \\
    <dual_template.txt> \\
    <single_index.json> \\
    <dual_index.json> \\
    <feature_index_type: single|dual> \\
    <gex_index_type: single|dual>
  ")
}

s_path          <- args[1] # sample table CSV
single_template <- args[2] # path to bclconvert_single_template.txt
dual_template   <- args[3] # path to bclconvert_dual_template.txt
single_index    <- args[4] # path to truseq_i7_indexes.json
dual_index      <- args[5] # path to chromium_10X_indexes.json
feature_index_type <- tolower(args[6]) # single or dual
gex_index_type     <- tolower(args[7]) # single or dual

valid_index_types <- c("single", "dual")
if (!feature_index_type %in% valid_index_types) {
  stop("feature_index_type must be one of: single, dual")
}
if (!gex_index_type %in% valid_index_types) {
  stop("gex_index_type must be one of: single, dual")
}

if(FALSE){
  # for testing purposes
  s_path          <- "data/example_sampletable.csv"
  single_template <- "data/bcl_samplesheet_templates/bclconvert_single_template.txt" # path to bclconvert_single_template.txt
  dual_template   <- "data/bcl_samplesheet_templates/bclconvert_dual_template.txt" # path to bclconvert_dual_template.txt
  single_index    <- "data/truseq_i7_indexes.json" # path to truseq_i7_indexes.json
  dual_index      <- "data/chromium_10X_indexes.json" # path to chromium_10X_indexes.json
  feature_index_type <- "single"
  gex_index_type <- "dual"
}

### custom functions -------------------------------------------------------

reverse_complement <- function(dna) {
  dna <- toupper(dna)
  
  comp_map <- c(
    A = "T",
    T = "A",
    G = "C",
    C = "G",
    N = "N"
  )
  
  dna %>%
    map_chr(\(index){
      index |>
        strsplit("") |>
        pluck(1) |>
        map_chr(\(base){ comp_map[[base]]}) |>
        rev() |>
        paste0(collapse = "")
    } 
    )
}



make_bcl_sheet <- function(sID, index, type, revcomp = TRUE,
                           single_template, dual_template,
                           single_idx, dual_idx){
  if(type == "single"){
    # print("single")
    template <- readLines(single_template)
    single_index <- jsonlite::read_json(single_idx)
    
    # get the i7 index seq(s)
    idx <- purrr::map(index, \(x){
      single_index[[x]][[1]]
    }) |>
      purrr::list_c()
    
    # create sample IDs
    # sID <- as.character(glue::glue("S00{1:length(index)}"))  
    
    # append the sample sheet lines
    bcl_sheet <- c(
      template,
      paste(sID, idx, sep=",")
    )
  } else if(type == "dual"){
    # print("dual")
    template <- readLines(dual_template)
    dual_index <- jsonlite::read_json(dual_idx)
    
    # get the dual index seq(s)
    idx1 <- purrr::map(index, \(x){
      dual_index[[x]][[1]]
    }) |>
      purrr::list_c()
    idx2 <- purrr::map(index, \(x){
      dual_index[[x]][[2]]
    }) |>
      purrr::list_c()
    
    # reverse complement the i5 index if needed
    if(revcomp){
      idx2 <- reverse_complement(idx2)
    }

    
    # create sample IDs
    # sID <- as.character(glue::glue("S00{1:length(index)}"))  
    
    
    # append the sample sheet lines
    bcl_sheet <- c(
      template,
      paste(sID, idx1, idx2, sep=",")
    )
    
  } else{
    stop("type must be single or dual")
  }
  
  bcl_sheet
}

# read in the sample table

sampletable <- read_csv(s_path)

write_lines(
  make_bcl_sheet(
    sampletable$sample_id,
    sampletable$feature_index, 
    feature_index_type,
    revcomp = TRUE,
    single_template, dual_template,
    single_index, dual_index
    ),
  "legotek_samplesheet.csv"
)

write_lines(
  make_bcl_sheet(
    sampletable$sample_id,
    sampletable$gex_index,
    gex_index_type,
    revcomp = TRUE,
    single_template, dual_template,
    single_index, dual_index
    ),
  "gex_samplesheet.csv"
)
