# CHiDO is a no-code platform to integrate multi-omics data to build, train and test
# linear mixed models for identifying candidates for desired GxE interactions.
#
# Copyright (C) 2024 Francisco Gonzalez, Diego Jarquin, and Julian Garcia
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

source("code/functions/utils.R")

# Load file data and return a data.frame object
load_data <- function(path = NULL) {
  file <- NULL
  
  if(!is.null(path)) {
    # RDA files
    if(grepl("\\.rda$", path)) {
      file <- get(load(path))
    }
    # CSV files
    else if (grepl("\\.csv$", path)) {
      file <- read.csv(path)
    }
    # Text files
    else if (grepl("\\.txt$", path)) {
      file <- read.delim(path, header = TRUE)
    }
    else {
      showNotification("Invalid file format: only CSV, RDA, and TXT files permitted", type="error")
      return()
    }
  }
  
  # Clean up session
  rm(path)
  # Return data.frame object
  return(data.frame(file))
}

# Create omic metadata list, no actual omic data is attached yet
create_omic_config <- function(file, label, data_type,modtype, id_col, gid_col = NULL,
                               g1id_col = NULL, g2id_col = NULL,eid_col = NULL, uid_col=NULL, link = NULL) {
  if(!is.null(file)) {
    
    config <- list(
      label = label,
      path = file$name,
      type = data_type
    )
    
    if(label=="Y") {
      if (modtype=="Genotype level"){
      config$trait_col <- as.integer(id_col)
      config$gid_col <- as.integer(gid_col)
      config$eid_col <- as.integer(eid_col)
      config$uid_col <- as.integer(uid_col)
      }else{
      config$trait_col <- as.integer(id_col)
      config$gid_col <- as.integer(gid_col)
      config$g1id_col <- as.integer(g1id_col)
      config$g2id_col <- as.integer(g2id_col)
      config$eid_col <- as.integer(eid_col)
      config$uid_col <- as.integer(uid_col)
      }
    } else {
      
      if (modtype=="Genotype level"){
      if (tolower(data_type) %in% c("genomic markers", "pedigree data", "phenomic markers")) {
        linkage_type = "Genotype / Line ID"
      } else if (tolower(data_type) == "environmental markers") {
        linkage_type = "Environment ID"
      } else {
        linkage_type = link
      }}else{
        if (tolower(data_type) %in% c("genomic markers", "pedigree data", "phenomic markers")) {
          linkage_type = link
        } 
          else if (tolower(data_type) == "environmental markers") {
          linkage_type = link
        } else {
          linkage_type = link
        }
      }
      
      config$id_col <- as.integer(id_col)
      config$linkage_type <- linkage_type
      #config$uid_col<-as.integer(uid_col)
    }
    
    # Return omic metadata
    return(config)
    
  } else {
    # File was not specifiec
    return()
  }
}

verify_omic_object <- function(config, phen = FALSE, modtype="Genotype level") {
  if(phen) {         # Data is from phenotype response (Y) file
    if(modtype=="Genotype level"){
    req_fields <- c("gid_col", "eid_col", "trait_col")
    if(config$AMMI){req_fields<-c(req_fields,"uid_col")}
    }else{
    req_fields <- c("gid_col","g1id_col","g2id_col", "eid_col", "trait_col") 
    }
  } else {           # Data is from another (X) omic file (e.g. markers)
    req_fields <- c("label", "id_col")
  }
  
  num_fields <- get_num_fields_from_list(config)
  actual_nums <- sapply(num_fields, function(field) config[[field]])
  
  # Verify that required fields are filled
  if(any(sapply(req_fields, function(field) is.null(config[[field]]) || 
                config[[field]] == "" || is.na(config[[field]])))) {
    showNotification("Please fill out required fields (marked by *)!", type="error")
    return(FALSE)
  }
  
  # Verify no numeric values are out-of-range
  if (any(actual_nums < 1 | actual_nums > 999)) {
    showNotification("Allowed column ranges are between 1 and 999", type = "error")
    return(FALSE)
  }
  
  # Check that the column values given are unique, do not overlap
  if (length(unique(actual_nums)) != length(actual_nums)) {
    showNotification("All column values should be unique", type = "error")
    return(FALSE)
  }
  
  return(TRUE)
}

update_table <- function(curr_data, config) {
  # Check if a row with same label already exists
  existing_row_idx <- which(curr_data$label == config$label)
  
  new_data <- data.frame(
    Label = config$label,
    File = config$path,
    Type = config$type,
    ID = config$id_col,
    Linkage = config$linkage_type
  )
  
  if (length(existing_row_idx) > 0) {
    # Update the existing data
    curr_data[existing_row_idx,] <- new_data
    showNotification("Available Omics table has just been updated!", type="message")
  } else {
    # Add new entry to the data table
    curr_data <- rbind(curr_data, new_data)
    showNotification("New data added to Available Omics table!", type="message")
  }
  
  return(curr_data)
}

# Function to validate the formula
validate_formula <- function(terms, labels) {
  # Get all unique values in formula
  unique_terms <- unique(unlist(lapply(terms, function(item) strsplit(item, "*", fixed = TRUE))))
  # Ensure the unique terms are actual labels
  return(all(unique_terms %in% labels))
}