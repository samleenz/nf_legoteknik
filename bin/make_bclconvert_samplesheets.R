#! /usr/bin/env Rscript

# Sam Lee
# @samleenz

library(tidyverse)

# sampletable path
s_path <- commandArgs(trailingOnly=TRUE)
# s_path <- "data/example_sampletable.csv"


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



make_bcl_sheet <- function(index, type, revcomp = TRUE){
  if(type == "single"){
    # print("single")
    template <- readLines("data/bcl_samplesheet_templates/bclconvert_single_template.txt")
    single_index <- jsonlite::read_json("data/truseq_i7_indexes.json")
    
    # get the i7 index seq(s)
    idx <- purrr::map(index, \(x){
      single_index[[x]][[1]]
    }) |>
      purrr::list_c()
    
    # create sample IDs
    sID <- as.character(glue::glue("S00{1:length(index)}"))  
    
    # append the sample sheet lines
    bcl_sheet <- c(
      template,
      paste(sID, idx, sep=",")
    )
  } else if(type == "dual"){
    # print("dual")
    template <- readLines("data/bcl_samplesheet_templates/bclconvert_dual_template.txt")
    dual_index <- jsonlite::read_json("data/chromium_10X_indexes.json")
    
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
    sID <- as.character(glue::glue("S00{1:length(index)}"))  
    
    
    # append the sample sheet lines
    bcl_sheet <- c(
      template,
      paste(sID, idx, idx2, sep=",")
    )
    
  } else{
    stop("type must be single or dual")
  }
  
  bcl_sheet
}

# read in the sample table

sampletable <- read_csv(s_path)

write_lines(
  make_bcl_sheet(sampletable$feature_index, "single", revcomp = TRUE),
  "legotek_samplesheet.csv"
)

write_lines(
  make_bcl_sheet(sampletable$gex_index, "dual", revcomp = TRUE),
  "gex_samplesheet.csv"
)
