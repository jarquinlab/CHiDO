# CHiDO is a no-code platform to integrate multi-omics data to build, train and test
# linear mixed models for identifying candidates for desired GxE interactions.
#
# Copyright (C) 2025 Francisco Gonzalez, Diego Jarquin, and Julian Garcia, Vitor Sagae
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

library(dplyr)
library(ggplot2)
library(BGLR)
library(gridExtra)

reformat_model <- function(item) {
  # Remove parentheses
  clean_item <- gsub("[()]", "", item)
  
  # Split by '*'
  split_item <- strsplit(clean_item, "\\*")[[1]]
  
  # If there's only one item after splitting, return it as a character. 
  # Otherwise, return the split item as is
  if(length(split_item) == 1) {
    return(as.character(split_item))
  } else {
    return(split_item)
  }
}

get_join_id <- function(curr_term, y_data) {
  if (curr_term$linkage_type == "Environment ID") {
    return(y_data$eid_col)
  } else if (curr_term$linkage_type == "Genotype / Line ID") {
    return(y_data$gid_col)
  } else if (curr_term$linkage_type == "Parent Group 1 ID") {
    return(y_data$g1id_col)
  } else if (curr_term$linkage_type == "Parent Group 2 ID") {
    return(y_data$g2id_col)
  } else if(curr_term$linkage_type == "Compound ID / UID") {
    return(y_data$uid_col)
  }
  }


# Filter genomic data sets to remove any markers that exceed NaN limit
filter_markers <- function(data, nan_freq, out_path = NULL) {
  
  # Determine maximum number of allowed NaN values per marker
  nan_limit <- nrow(data) * nan_freq / 100
  # Set bool values based on number of NaNs in a column
  nan_matrix <- colSums(is.na(data)) <= nan_limit
  
  if (!is.null(out_path)) {
    # Create report of NaNs in markers dataset
    nan_report <- data.frame(
      marker = colnames(data),
      count = colSums(is.na(data)),
      percent = round((colSums(is.na(data)) * 100 / nrow(data)), 2)
    )
    
    # Create new directory if needed
    if (!dir.exists(dirname(out_path))) {
      dir.create(dirname(out_path), recursive = TRUE)
    }
    
    # Save report
    write.table(nan_report, file=out_path, sep=",", row.names=FALSE, col.names=TRUE)
  }
  
  # Return filtered markers data frame
  filtered_data <- data[, nan_matrix]
  return(filtered_data)
}

create_z_matrix <- function(y_data, join_id) {
  ids <- factor(y_data[, join_id])
  return(as.matrix(model.matrix(~ids - 1)))
}

pre_process_data <- function(data, type, s, wht, ctr, std, prop_maf_j) {
  
  # TODO -- Ask Diego for examples to integrate 
  if(wht) {
    weight <- scan(weight.file, skip = 1)
  }
  
  # Naive imputation, Centering, standardizing and adjusting weights
  for (i in 1:ncol(data)) {
    mean_i <- mean(data[, i], na.rm = TRUE)
    data[, i] <- ifelse(is.na(data[, i]), mean_i, data[, i])
    if (ctr) { data[, i] <- data[, i] - mean_i }
    if (std) {
      sd_i <- sd(data[, i])
      # Handle zero standard deviation
      if (sd_i == 0) { data[, i] <- 0 } 
      else { data[, i] <- data[, i] / sd_i }
    }
    if (wht) { data[, i] <- data[, i] * weight[i] }
    s <- s + var(data[, i], na.rm=TRUE)
  }

  # TODO -- Ask if this is necessary for any other omic except G
  if (type != "phenomic markers") {

  }
  
  # Filter for desired MAF proportion
  if (type == "genomic markers" && !is.null(prop_maf_j)) {
    p <- colMeans(data, na.rm = T) / 2
    p <- ifelse(p <= 0.5, p, 1 - p)
    # Look through columns in data to keep 
    maf_idx <- which(p >= prop_maf_j)
    data <- data[, maf_idx]
  }

  # Centering, standardizing and adjusting weights
  #for (i in 1:ncol(data)) {
  #  mean_i <- mean(data[, i], na.rm = TRUE)
  #}
  
  # row x row -- Initial G but not final G matrix
  if (type %in% c("genomic markers")) {
    G <- tcrossprod(data) / s
  } else {
    G <- tcrossprod(data) / ncol(data)
  }
  
  return(G)
  
}

create_g_matrix <- function(o_data, y_data, wht=FALSE, ctr=FALSE, std=FALSE, 
                            nan_freq = NULL, prop_maf_j=NULL) {
  
  data <- o_data$data
  type <- o_data$type
  id_col <- o_data$id_col
  join_id <- o_data$y_col

  ids <- factor(y_data$data[, join_id])
  Z <- as.matrix(model.matrix(~ids - 1))

  # No data, create G from Y values only
  if (is.null(data) || o_data$label %in% c("E","L","L1","L2")) {
    
    d <- colSums(Z)
    #V <- Z / sqrt(d)
    # Final G and EVD values
    #EVD <- list(vectors = V, values = d)
    #G <- tcrossprod(Z)
    
    ## Adjusted
    V<-Z
    for(i in 1:ncol(Z)){ V[,i]<-V[,i]/sqrt(d[i]) } 
    EVD<-list(vectors=V,values=d)
    G<-tcrossprod(Z)
  } else {
    
    # Filter genomic data based on NaN threshold limit
    if (type %in% c("genomic markers")) {
      data <- filter_markers(data, nan_freq)
    }
    
    O <- as.matrix(data[,-id_col])
    rownames(O) <- data[[id_col]]
    
    if (type != "pedigree") {
      # Do the processing for data to reach initial G
      s <- 0
      GO <- pre_process_data(O, type, s, wht, ctr, std, prop_maf_j)
      rownames(GO)<-data[[id_col]]
      colnames(GO)<-data[[id_col]]

      # Ensure all Y IDs 
      if(!all(ids %in% rownames(GO))) {
        # There was an error, check that data is available for all genotypes
        showNotification("There was an error, check that data is available for all genotypes",
                         type = "error")
        return ()
      }
      GO <- GO[ rownames(GO) %in% unique(ids)  ,  colnames(GO) %in% unique(ids) ]
      
      ids <- factor(ids, levels = rownames(GO))
      Z <- as.matrix(model.matrix(~ids - 1))
      GO <- tcrossprod(Z, GO)
      
    } else {
      # Initial O matrix for pedigree data
      GO <- tcrossprod(Z, 2*as.matrix(O))
    }
    
    G <- tcrossprod(GO, Z)
    EVD <- eigen(G)
    rownames(EVD$vectors) <- rownames(G)
    
  }
  
  return(list(G=G, EVD=EVD))
  
}

create_interaction_matrix <- function(termvec, tmpdir) {
  # Initial interaction matrix
  GI <- NULL
  
  # Loop through components in the interaction term
  for (i in seq_along(termvec)) {
    # Check directory for term G matrix
    g_file <- file.path(tmpdir, termvec[i], "G.rda")
    
    # Notify if file doesn't exist
    if (!file.exists(g_file)) { return() }
    
    GT <- get(load(g_file))
    
    if (i == 1) { 
      GI <- GT
    } else {
      GI <- GI * GT
    }
  }
  
  EVD <- eigen(GI)
  
  return(list(G = GI, EVD = EVD))
}

get_evd_plots <- function(EVD, out_path) {
  plt1 <-
    ggplot(data = as.data.frame(EVD$vectors), aes(x = EVD$vectors[, 1], y = EVD$vectors[, 2])) +
    geom_point() +
    labs(title = "",
         xlab = "First Component",
         ylab = "Second Component")
  plt2 <-
    ggplot(data = as.data.frame(EVD$vectors), aes(x = EVD$vectors[, 1], y = EVD$vectors[, 3])) +
    geom_point() +
    labs(title = "",
         xlab = "First Component",
         ylab = "Third Component")
  plt3 <-
    ggplot(data = as.data.frame(EVD$vectors), aes(x = EVD$vectors[, 2], y = EVD$vectors[, 3])) +
    geom_point() +
    labs(title = "",
         xlab = "Second Component",
         ylab = "Third Component")
  
  plt4 <-
    ggplot(data = as.data.frame(EVD$values), aes(x = seq_along(EVD$values), y = EVD$values)) +
    geom_line() +
    labs(title = "Eigen-values", x = "Components", y = "Eigen-value") +
    theme_minimal()
  
  plt5 <-
    ggplot(data = as.data.frame(EVD$values), aes(
      x = seq_along(EVD$values),
      y = cumsum(EVD$values) / sum(EVD$values)
    )) +
    geom_line() +
    labs(title = "Cumulative Variance", x = "Components", y = "Cumulative Variance") +
    theme_minimal()
  
  # Add horizontal and vertical lines to Plot 5
  plt5 <- plt5 + geom_hline(yintercept = 0.8, linetype = "dashed") +
    geom_vline(xintercept = which(floor(
      cumsum(EVD$values) / sum(EVD$values) * 10
    ) == 8)[1], linetype = "dashed")
  
  # Combine the two plots
  plots <- grid.arrange(plt1, plt2, plt3, plt4, plt5)
  ggsave(file.path(out_path, "Eigen.pdf"), plots)
}

# Assign individuals to respective folds based on CV technique
prep_data_for_cv <- function(y_data, folds, cv1, cv2, cv0, cv00) {
  
  data <- y_data$data
  gid <- y_data$gid_col
  eid <- y_data$eid_col
  trait_col <- y_data$trait_col
  
  # Auto-assign as true since these are precursors for CV0 and CV00
  #if (cv0) { cv2 <- TRUE }
  if (cv00) { cv1 <- TRUE; cv0<-TRUE }
  
  gids <- unique(data[, gid])
  eids<-unique(data[,eid])
  col_folds <- ncol(data)
  
  # Replace NA values with 0 for now
  # na_rows <- which(is.na(data[, trait_col]))
  # data[na_rows, trait_col] <- 0
  
  # Prepare data.frame to return
  retdf <- as.data.frame(data[, gid])
  retdf$env <- data[, eid]
  retdf$obs <- data[, trait_col]
  
  # Assign fold by line/variety ID
  if (cv1) {
    fold <- rep(1:folds, each=ceiling(length(gids)/folds))[order(runif(length(gids)))]
    y_cv1 <- NA
    for (f in 1:folds) {
      # get all lines IDs in gids that correspond to specific fold
      fold_lines <- gids[fold == f]
      y_cv1[data[[gid]] %in% fold_lines] <- f
    }
    
    retdf$cv1 <- y_cv1
  }
  
  # Assign fold by individual row
  if (cv2) {
    y_cv2 <- NA
    
    for (id in gids) {
      # Find rows with the same line ID
      line_idx <- which(data[, gid] == id)
      
      num_lines <- length(line_idx)
      fold <- sample(1:folds, size = num_lines, replace = num_lines > folds)
      y_cv2[line_idx] <- fold
    }
    
    retdf$cv2 <- y_cv2
  }
  
  if (cv0) {
    y_cv0 <- NA
    fold <- rep(1:folds, times=ceiling(length(eids)/folds))[order(runif(length(eids)))]

    for (f in 1:length(unique(fold))) {
      # Find the index for all rows with the same fold value
      fold_envs <- eids[fold == f]
      # Replace y_cv0 values with NA for these indexes
      
      y_cv0[data[[eid]] %in% fold_envs] <- f
      #y_cv0[eid %in% fold_envs] <- f
      
    }
    retdf$cv0 <- y_cv0
  }
  
  if (cv00) {
    y_cv00<-matrix(NA,nrow=dim(data),ncol=folds)
    y_cv00[,] <- data[,trait_col]
    
    for(fold in 1:folds)
    {
      y_cv00[retdf$cv0==fold,fold] <- NA
    }
    
    y_cv00 <- do.call(cbind, lapply(1:ncol(y_cv00), function(i) matrix(rep(y_cv00[, i], folds), ncol=folds)))
    colnames(y_cv00) <- paste0(colnames(data)[trait_col],"_cv00_fold_",1:folds^2)
    
    for (i in seq(1,folds^2,by=folds)) {
      
    for (j in 1:folds){
      cl_ind<-i+(j-1)
      
      y_cv00[retdf$cv1==j ,cl_ind] <- NA

    }
    }
    retdf <- cbind(retdf, y_cv00)
  }
  return(retdf)
  
}

# Get predictions based on CV methods
get_predictions <- function(data, tmpdir, model_selected, model_eq, cv, folds, 
                            nIter=NULL, burnIn=NULL, esc=FALSE, model="RKHS") {
  
  # Create output directory
  outdir <- file.path(tmpdir, "output", model_selected)
  if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
  
  # Create base prediciton data frame
  predictions <- data.frame(testing = integer(),
                            fold = integer(),
                            individual = character(),
                            env = character(),
                            observed = numeric(),
                            predicted = numeric())
  
  y <- data[, 3]            # Observed values
  eid <- data[, 2]          # Environment IDs
  gid <- data[, 1]          # Line IDs
  
  # Base ETA
  eta <- list()
  
  if (esc) { y <- scale(y, center = TRUE, scale = TRUE) }
  
  for (i in seq_along(model_eq)) {
    # Get file for the term in model
    file_name <- file.path(tmpdir, model_eq[i], "EVD.rda")
    Z <- get(load(file_name))
    
    # TODO -- Ask Diego about algorithm choice here
    if (model == "RKHS") {
      eta[[i]] <- list(V = EVD$vectors, d=EVD$values, model="RKHS")
      rm(EVD)
    }
    rm(Z)
  }
  
  if(!is.null(cv)) {
    #cv <- tolower(cv)
    
    if (cv %in% c("cv00")) {
      cv_data <- data[, grep(paste0(cv,"_"), colnames(data))]
    } else {
      cv_data <- data[, match(cv, colnames(data))]
    }
    
    # Get predictions

    predictions <- fit_cv(cv, cv_data, data, folds, predictions, eta, nIter, burnIn)
    
    write.csv(predictions, file=file.path(outdir, paste0(cv, ".csv")), row.names = FALSE)
    
  } else {
    # Using full data to get variance components
    fm <- BGLR(y=y,ETA=eta,nIter=nIter,burnIn=burnIn,verbose=TRUE)
    fm$equation <-  model_eq
    
    # Save fm
    save(fm, file=file.path(outdir, 'fm.rda'))
  }
  
  # Clean work space
  unlink("*.dat")
  
  # Empty data.frame if cv is NULL
  return(predictions)
}

fit_cv <- function(cv, cv_data, data, folds, predictions, eta, nIter, burnIn) {
  
  y <- data[, 3]            # Observed values
  eid <- data[, 2]          # Environment IDs
  gid <- data[, 1]          # Line IDs
  eids<- unique(eid)

  if (cv %in% c("cv00")) {

    idcv00<-list()
    k<-1
    for (i in seq(1,folds^2,by=folds)) {
     
      for (j in 1:folds){
        cl_ind<-i+(j-1)
        
        idcv00[[cl_ind]]<-which((data$cv0==k) & (data$cv1==j))
      
      }
    k <- (k %% folds) + 1  # 
    }
    idcv00<-do.call(cbind,idcv00)

            for (fold in seq_along(cv_data)) {
      # Train data
      y_na <- cv_data[, fold]
      # Test data index

      #testing <- which(is.na(y_na))
      testing<-idcv00[,fold]
      
      if(all(is.na(testing))){
        next
      }
      
      # Fit BGLR model with training data
      fm <- BGLR(y=y_na, ETA=eta, nIter=nIter, burnIn=burnIn, verbose=FALSE)
      
      # Predict with the fitted model
      preds <- data.frame(testing = testing, fold = fold, env = eid[testing],
                          individual = gid[testing], observed = y[testing],
                          predicted = round(fm$yHat[testing],3))
      
      predictions <- rbind(predictions, preds)
    }
  }
  
  if (cv %in% c("cv0")) {
    
    for (fold in 1:length(unique(data$cv0))) {
      y_na <- y
      
      if (fold != -999) {
        # Get idx for testing
        testing <- which(cv_data == fold)
        # Set testing data for fold to NA
        y_na[testing] <- NA
      
      # Fit BGLR model with training data
      fm <- BGLR(y=y_na, ETA=eta, nIter=nIter, burnIn=burnIn, verbose=FALSE)
      
      # Predict with the fitted model
      preds <- data.frame(testing = testing, fold = fold, env = eid[testing],
                          individual = gid[testing], observed = y[testing],
                          predicted = round(fm$yHat[testing],3))
      
      predictions <- rbind(predictions, preds)
    }
  }
  }
  
  if (cv %in% c("cv1","cv2")) {
    
    for (fold in 1:folds) {
      y_na <- y
      
      if (fold != -999) {
        # Get idx for testing
        testing <- which(cv_data == fold)
        # Set testing data for fold to NA
        y_na[testing] <- NA
        
        # Fit BGLR model
        fm <- BGLR(y=y_na, ETA = eta, nIter=nIter, burnIn=burnIn, verbose=TRUE)
        
        # Generate predictions
        
        preds <- data.frame(testing = testing, fold = fold, env = eid[testing],
                            individual = gid[testing], observed = y[testing],
                            predicted = round(fm$yHat[testing],3))
        # Add predictions to dataset
        predictions <- rbind(predictions, preds)
      }
    }
  }
  
  # Clean work space
  rm(fm)
  closeAllConnections()
  
  return(predictions)
}



