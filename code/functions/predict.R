extract_omics <- function(formula_str) {
  clean_formula <- gsub("\\s+", "", formula_str)
  terms <- unlist(strsplit(clean_formula, "\\+"))
  vars <- unique(unlist(strsplit(terms, "\\*|\\+")))
  return(setdiff(vars, c("E", "L")))
}

predict_helper <- function(trn_list, prd_list, phen_list, prompt, model){
  hash_id <- c('EID', 'GID', 'UID')
  names(hash_id) <- c("Environment ID", "Genotype / Line ID", "Compound ID / UID")
  
  phen_merged <- phen_list
  merged <- trn_list
  
  omics <- names(prd_list)
  pred_booleans <- list()
  #req(all(setdiff(names(trn), omics) %in% c('E', 'L')))
  
  for (o in omics){
    merged[[o]]$id_col <- 1
    trn_o <- trn_list[[o]]
    prd_o <- prd_list[[o]]
    
    all_ids <- sort(unique(c(prd_o$ID, trn_o$data[,trn_o$id_col])))
    id_type <- hash_id[trn_o$linkage_type]
    pred_booleans[[o]] <- prompt[,id_type] %in% all_ids
    trn_feats <- setdiff(colnames(trn_o$data), id_type)
    prd_feats <- setdiff(colnames(prd_o), 'ID')
    
    if (trn_o$type == 'pedigree data'){
      c11 <- as.matrix(trn_o$data %>% dplyr::select(-id_type))
      rownames(c11) = trn_o$data[,trn_o$id_col]
      c21 <- as.matrix(prd_o %>% .[,trn_feats])
      rownames(c21) = prd_o$ID
      c22 <- as.matrix(prd_o %>% .[,rownames(c21)])
      rownames(c22) = rownames(c21)
      c <- cbind(rbind(c11,c21),rbind(t(c21),c22))[as.character(all_ids), as.character(all_ids)]
      prd_ready <- data.frame('ID' = all_ids, c, check.names = FALSE)
      colnames(prd_ready)[1] <- id_type
      merged[[o]]$data <- prd_ready
    }else{
      common_feats <- intersect(trn_feats, prd_feats)
      if (!length(common_feats)){
        stop(sprintf("New data for omic %s has no common features with the training data", o))
      }
      prd_ready <- (prd_o[, c('ID',common_feats)])
      colnames(prd_ready)[1] <- id_type
      merged[[o]]$data <- rbind(trn_o$data[,c(id_type, common_feats)], prd_ready)
    }
  }
  prompt$Set <- 'PRD'
  prompt$Predictable <- apply(bind_rows(pred_booleans), 1, all)
  prompt$Model <- model
  prompt$Y <- NA
  
  # merge phenotypic training
  ntrn <- nrow(phen_list$data)
  prompt$index <- ntrn + cumsum(as.numeric(prompt$Predictable))
  prompt$index[!prompt$Predictable] <- NA
  newp_prompt <- prompt %>% filter(Predictable) %>% dplyr::select(EID, GID, UID, Y)
  
  phen_merged$data <- rbind(phen_list$data, newp_prompt)
  merged$E$data <- data.frame(phen_merged$data[,phen_merged$eid_col])
  merged$E$uid <- data.frame(phen_merged$data[,phen_merged$uid_col])
  merged$L$data <- data.frame(phen_merged$data[,phen_merged$gid_col])
  merged$L$uid <- data.frame(phen_merged$data[,phen_merged$uid_col])
  
  output = list(phen_merged = phen_merged, omics_merged = merged, prompt = prompt)
  return(output)
}
