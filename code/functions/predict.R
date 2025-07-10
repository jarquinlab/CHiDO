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
    prd_feats <- setdiff(colnames(prd_o), c('ID'))
    
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
  new_prompt <- prompt %>% filter(Predictable) %>% dplyr::select(EID, GID, UID, Y)
  
  phen_merged$data <- rbind(phen_list$data, new_prompt)
  merged$E$data <- data.frame(phen_merged$data[,phen_merged$eid_col])
  merged$E$uid <- data.frame(phen_merged$data[,phen_merged$uid_col])
  merged$L$data <- data.frame(phen_merged$data[,phen_merged$gid_col])
  merged$L$uid <- data.frame(phen_merged$data[,phen_merged$uid_col])
  
  output = list(phen_merged = phen_merged, omics_merged = merged, prompt = prompt)
  return(output)
}


get_predictions2 <- function(data, tmpdir, model_selected, model_eq, cv, folds, 
                            nIter=NULL, burnIn=NULL, esc=FALSE, model="RKHS") {
  
  # Create output directory
  outdir <- file.path(tmpdir, "output", model_selected)
  if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
  
  y <- data[, 3]            # Observed values
  eid <- data[, 2]          # Environment IDs
  gid <- data[, 1]          # Line IDs
  
  # Base ETA
  eta <- list()
  for (i in seq_along(model_eq)) {
    file_name <- file.path(tmpdir, model_eq[i], "EVD.rda")
    Z <- get(load(file_name))
    eta[[i]] <- list(V = EVD$vectors, d=EVD$values, model="RKHS")
    rm(EVD)
    rm(Z)
  }
  
  fm <- BGLR(y=y,ETA=eta,nIter=nIter,burnIn=burnIn,verbose=TRUE)
  fm$equation <-  model_eq
  return(fm)
}

generate_predict_results <- function(predict_ds, fm){
  prompt <- predict_ds$prompt
  prompt$yHat <- fm$yHat[c(prompt$index)]
  
  # file 1 -> original prompt + yHats
  prompt_out <- prompt
  prompt_out$Y <- NULL
  prompt_out$index <- NULL
  
  # file 2 <- TRN + original prompt - useful for plots
  aux_out <- predict_ds$phen_merged$data
  ntrn <- nrow(aux_out) - nrow(prompt)
  aux_out <- aux_out[1:ntrn,]
  aux_out$Predictable <- TRUE
  aux_out$Set <- 'TRN'
  aux_out$Model <- unique(prompt$Model)
  aux_out$index <- 1:ntrn
  aux_out$yHat <- fm$yHat[1:ntrn]
  aux_out <- aux_out[, colnames(prompt)]
  aux_out <- rbind(aux_out, prompt)
  
  return(list(prompt_out = prompt_out, aux_out = aux_out))
}

cv1_pred_plot <- function(aux_out, p, top = TRUE){
  keep_envs <- aux_out %>% filter(Set %in% 'TRN') %>% .$EID
  quantiles <- aux_out %>% 
    filter(EID %in% keep_envs) %>% 
    group_by(EID) %>% 
    summarize(topq = quantile(yHat, 1-p), botq = quantile(yHat, p)) %>% 
    mutate(x = as.numeric(factor(EID)),
           xstart = x - .6, xend = x + .6)
  quantiles$yint <- quantiles$botq
  if (top){quantiles$yint <- quantiles$topq}
  
  df <- aux_out %>% 
    filter(EID %in% keep_envs) %>% 
    right_join(quantiles[,c('EID', 'topq', 'botq')])
  df$hl <- if (top) df$yHat > df$top else df$yHat < df$bot
  
  ggplot(df, aes(x = EID, y = yHat, color = Set)) +
    geom_jitter(size = 2, aes(alpha = hl)) +
    geom_segment(data = quantiles,
                 aes(x = xstart, xend = xend, y = yint, yend = yint),
                 inherit.aes = FALSE,
                 color = "darkgray", size = 0.5, linetype = 'dashed')+
    scale_alpha_manual(values = c(0.05,1))+
    scale_color_manual(values = c("#FA4616","#0021A5"))+
    guides(alpha = 'none')+
    theme_light()+
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18, face = 'bold'),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18, face = 'bold'))
}

cv0_pred_plot <- function(aux_out, k, new_env, top = TRUE, sc = TRUE){

  keep_envs <- unique(aux_out %>% filter(Set %in% 'TRN') %>% .$EID)
  quantiles <- aux_out %>% filter(EID %in% new_env) %>% 
    mutate(topq = rank(-yHat, ties.method = 'random'), botq = rank(yHat, ties.method = 'random'))
  best_gids <- quantiles %>% filter(botq <= k) %>% .$GID
  if (top){best_gids <- quantiles %>% filter(topq <= k) %>% .$GID}
  df <- aux_out %>% mutate(hl = GID %in% best_gids) %>% 
    filter(EID %in% union(new_env, keep_envs)) %>% 
    mutate(EID = factor(EID, levels = union(new_env, keep_envs))) 
  
  if (sc){
    df <- df %>% group_by(EID) %>% mutate(ymean = mean(yHat, na.rm = TRUE), 
                                          ysd = sd(yHat, na.rm = TRUE)) %>% 
      ungroup() %>% mutate(yHat = (yHat-ymean)/ysd)
    }

  #df$yHat[!is.na(df$Y)] <- df$Y[!is.na(df$Y)]
  df %>% ggplot(aes(x = EID, y = yHat, color = Set, alpha = hl))+
    geom_line(data = filter(df, hl), aes(group = GID), color = 'black')+
    geom_point(size = 3)+
    scale_alpha_manual(values = c(0.05,1))+
    scale_color_manual(values = c("#FA4616","#0021A5"))+
    guides(alpha = 'none')+
    labs(y = 'scaled yHat')+
    theme_light()+
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18, face = 'bold'),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18, face = 'bold'))
}
