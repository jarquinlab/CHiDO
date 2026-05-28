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
      combined <- rbind(trn_o$data[,c(id_type, common_feats)], prd_ready)
      merged[[o]]$data <- combined[!duplicated(combined[, id_type]), ]
    }
  }
  prompt$Set <- 'PRD'
  prompt$Predictable <- if (length(pred_booleans) == 0) {
    rep(TRUE, nrow(prompt))
  } else {
    Reduce("&", pred_booleans)
  }
  prompt$Model <- model
  trn_gids <- unique(phen_list$data[, phen_list$gid_col])
  trn_eids <- unique(phen_list$data[, phen_list$eid_col])
  new_gen <- !(prompt$GID %in% trn_gids)
  new_env <- !(prompt$EID %in% trn_eids)
  prompt$PRD_type <- dplyr::case_when(
    !new_env & !new_gen ~ "CV2",
    !new_env &  new_gen ~ "CV1",
     new_env & !new_gen ~ "CV0",
     new_env &  new_gen ~ "CV00"
  )

  # merge phenotypic training
  ntrn <- nrow(phen_list$data)
  prompt$index <- ntrn + cumsum(as.numeric(prompt$Predictable))
  prompt$index[!prompt$Predictable] <- NA
  predictable_rows <- prompt %>% filter(Predictable)
  new_prompt <- as.data.frame(matrix(NA, nrow = nrow(predictable_rows), ncol = ncol(phen_list$data)))
  colnames(new_prompt) <- colnames(phen_list$data)
  new_prompt[, phen_list$eid_col] <- predictable_rows$EID
  new_prompt[, phen_list$gid_col] <- predictable_rows$GID
  if (!is.na(phen_list$uid_col)) {
    new_prompt[, phen_list$uid_col] <- predictable_rows$UID
  }

  phen_merged$data <- rbind(phen_list$data, new_prompt)
  merged$E$data <- data.frame(phen_merged$data[, phen_merged$eid_col])
  merged$E$uid <- if (!is.na(phen_merged$uid_col)) data.frame(phen_merged$data[, phen_merged$uid_col]) else data.frame(rep(NA, nrow(phen_merged$data)))
  merged$L$data <- data.frame(phen_merged$data[, phen_merged$gid_col])
  merged$L$uid <- if (!is.na(phen_merged$uid_col)) data.frame(phen_merged$data[, phen_merged$uid_col]) else data.frame(rep(NA, nrow(phen_merged$data)))
  
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
  
  print('predicting')
  fm <- BGLR(y=y,ETA=eta,nIter=nIter,burnIn=burnIn,verbose=FALSE)
  fm$equation <-  model_eq
  return(fm)
}

generate_predict_results <- function(predict_ds, fm){
  prompt <- predict_ds$prompt
  prompt$yHat <- fm$yHat[c(prompt$index)]
  
  # file 1 -> original prompt + yHats
  prompt_out <- prompt
  prompt_out$index <- NULL
  
  # file 2 <- TRN + original prompt - useful for plots
  aux_out <- predict_ds$phen_merged$data
  ntrn <- nrow(aux_out) - sum(prompt$Predictable)
  aux_out <- aux_out[1:ntrn,]
  aux_out$Predictable <- TRUE
  aux_out$Set <- 'TRN'
  aux_out$Model <- unique(prompt$Model)
  aux_out$index <- 1:ntrn
  aux_out$yHat <- fm$yHat[1:ntrn]
  aux_out <- dplyr::bind_rows(aux_out, prompt)
  
  return(list(prompt_out = prompt_out, aux_out = aux_out))
}

cv1_pred_plot <- function(aux_out, p, top = TRUE){
  keep_envs <- sort(unique(aux_out %>% filter(Set %in% 'TRN') %>% .$EID))
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
    right_join(quantiles[,c('EID', 'topq', 'botq')]) %>%
    mutate(
      EID = factor(EID, levels = keep_envs),
      plot_group = case_when(
      Set == "TRN" ~ "TRN",
      new_gen       ~ "PRD (CV1)",
      TRUE          ~ "PRD (CV2)"
    ))
  df$hl <- if (top) df$yHat > df$top else df$yHat < df$bot
  df$color_group <- ifelse(df$hl, df$plot_group, "other")

  trn_obs <- aux_out %>%
    filter(Set == "TRN") %>%
    group_by(GID) %>%
    summarize(first_env = sort(unique(as.character(EID)))[1],
              n_envs = n_distinct(EID), .groups = "drop")

  df <- df %>%
    left_join(trn_obs, by = "GID") %>%
    mutate(obs_label = case_when(
      new_gen     ~ "Never observed in training",
      n_envs == 1 ~ paste0("Observed in: ", first_env),
      TRUE        ~ paste0("Observed in: ", first_env, " and ", n_envs - 1, " other(s)")
    ))

  p <- ggplot(df, aes(x = EID, y = yHat, color = color_group, alpha = hl,
                      tooltip = paste0("GID: ", GID, " | EID: ", EID, "\nyHat: ", round(yHat, 3), "\n", obs_label),
                      data_id = GID)) +
    geom_point_interactive(data = subset(df, color_group == "other"), size = 1, position = position_quasirandom()) +
    geom_point_interactive(data = subset(df, color_group != "other"), size = 1, position = position_quasirandom()) +
    geom_segment(data = quantiles,
                 aes(x = xstart, xend = xend, y = yint, yend = yint),
                 inherit.aes = FALSE,
                 color = "darkgray", size = 0.5, linetype = 'dashed')+
    scale_alpha_manual(values = c(0.15, 1))+
    scale_color_manual(name = "Set",
                       values = c("TRN" = "#0021A5", "PRD (CV2)" = "#FA4616", "PRD (CV1)" = "#27AE60", "other" = "gray80"),
                       breaks = c("TRN", "PRD (CV2)", "PRD (CV1)"))+
    guides(alpha = 'none')+
    coord_flip()+
    theme_light()+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12, face = 'bold'),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12, face = 'bold'))

  girafe(ggobj = p)
}

cv0_radial_plot <- function(aux_out, k, target_env, top = TRUE, cv0_only = FALSE) {

  trn_envs <- sort(unique(aux_out$EID[aux_out$Set == "TRN"]))
  n_envs   <- length(trn_envs)

  # Target env PRD candidates
  prd_cands <- aux_out %>% filter(EID == target_env, Set == "PRD")
  if (cv0_only) prd_cands <- prd_cands %>% filter(new_gen == FALSE)

  trn_obs <- aux_out %>%
    filter(Set == "TRN") %>%
    group_by(GID) %>%
    summarize(first_env = sort(unique(as.character(EID)))[1],
              n_envs    = n_distinct(EID), .groups = "drop") %>%
    mutate(obs_label = case_when(
      n_envs == 1 ~ paste0("Observed in: ", first_env),
      TRUE        ~ paste0("Observed in: ", first_env, " and ", n_envs - 1, " other(s)")
    ))

  prd_ranked <- prd_cands %>%
    mutate(target_rank = rank(-yHat, ties.method = "random")) %>%
    arrange(target_rank) %>%
    left_join(trn_obs %>% dplyr::select(GID, obs_label), by = "GID") %>%
    mutate(obs_label = if_else(new_gen | is.na(obs_label), "Never observed in training", obs_label))

  n_genos   <- nrow(prd_ranked)
  top_k_gids <- if (top) {
    prd_ranked %>% slice_head(n = k) %>% pull(GID)
  } else {
    prd_ranked %>% slice_tail(n = k) %>% pull(GID)
  }

  r_in       <- 0.4
  r_out      <- 0.8
  label_r    <- 0.93
  sector_gap <- 0.15
  sector_span <- (2 * pi / n_envs) * (1 - sector_gap)

  inner_list <- list()
  outer_list <- list()
  lines_list <- list()

  for (ei in seq_along(trn_envs)) {
    env           <- trn_envs[ei]
    sector_center <- (ei - 1) / n_envs * 2 * pi - pi / 2

    # Inner: target env genotypes ranked by target yHat
    inner_angles <- sector_center - sector_span / 2 +
      (prd_ranked$target_rank / (n_genos + 1)) * sector_span

    inner_df <- prd_ranked %>%
      mutate(
        angle      = inner_angles,
        x          = r_in * cos(angle),
        y          = r_in * sin(angle),
        ring       = "inner",
        sector_env = env,
        is_topk    = GID %in% top_k_gids,
        color_group = case_when(!is_topk ~ "other", new_gen ~ "PRD (CV00)", TRUE ~ "PRD (CV0)"),
        tooltip    = paste0("GID: ", GID, " | EID: ", target_env, "\nyHat: ", round(yHat, 3), "\n", obs_label)
      )

    # Outer: same genotypes ranked by TRN env yHat
    trn_data <- aux_out %>%
      filter(EID == env, Set == "TRN", GID %in% prd_ranked$GID) %>%
      mutate(trn_rank = rank(-yHat, ties.method = "random")) %>%
      arrange(trn_rank)

    if (nrow(trn_data) == 0) next
    n_trn <- nrow(trn_data)

    outer_angles <- sector_center - sector_span / 2 +
      (trn_data$trn_rank / (n_trn + 1)) * sector_span

    outer_df <- trn_data %>%
      left_join(prd_ranked %>% dplyr::select(GID, obs_label), by = "GID") %>%
      mutate(
        angle      = outer_angles,
        x          = r_out * cos(angle),
        y          = r_out * sin(angle),
        ring       = "outer",
        sector_env = env,
        is_topk    = GID %in% top_k_gids,
        color_group = if_else(is_topk, "TRN", "other"),
        tooltip    = paste0("GID: ", GID, " | EID: ", env, "\nyHat: ", round(yHat, 3), "\n", obs_label)
      )

    line_df <- inner_df %>%
      dplyr::select(GID, x, y, is_topk) %>%
      dplyr::rename(x_in = x, y_in = y) %>%
      inner_join(outer_df %>% dplyr::select(GID, x, y) %>% dplyr::rename(x_out = x, y_out = y),
                 by = "GID")

    inner_list[[ei]] <- inner_df
    outer_list[[ei]] <- outer_df
    lines_list[[ei]] <- line_df
  }

  all_inner <- bind_rows(inner_list)
  all_outer <- bind_rows(outer_list)
  all_lines <- bind_rows(lines_list)

  sector_labels <- data.frame(
    env   = trn_envs,
    angle = (seq_along(trn_envs) - 1) / n_envs * 2 * pi - pi / 2
  ) %>% mutate(
    x     = label_r * cos(angle),
    y     = label_r * sin(angle),
    hjust = case_when(cos(angle) >  0.1 ~ 0, cos(angle) < -0.1 ~ 1, TRUE ~ 0.5),
    vjust = case_when(sin(angle) >  0.1 ~ 0, sin(angle) < -0.1 ~ 1, TRUE ~ 0.5)
  )

  color_vals <- c("TRN" = "#0021A5", "PRD (CV0)" = "#FA4616",
                  "PRD (CV00)" = "#27AE60", "other" = "gray80")

  p <- ggplot() +
    annotate("path",
             x = r_in  * cos(seq(0, 2 * pi, length.out = 300)),
             y = r_in  * sin(seq(0, 2 * pi, length.out = 300)),
             color = "gray90", linewidth = 0.3) +
    annotate("path",
             x = r_out * cos(seq(0, 2 * pi, length.out = 300)),
             y = r_out * sin(seq(0, 2 * pi, length.out = 300)),
             color = "gray90", linewidth = 0.3) +
    geom_segment(data = filter(all_lines, !is_topk),
                 aes(x = x_in, xend = x_out, y = y_in, yend = y_out),
                 color = "gray60", alpha = 0.25, linewidth = 0.25) +
    geom_segment(data = filter(all_lines,  is_topk),
                 aes(x = x_in, xend = x_out, y = y_in, yend = y_out),
                 color = "black", linewidth = 0.8) +
    geom_point_interactive(data = filter(all_inner, !is_topk),
                           aes(x = x, y = y, tooltip = tooltip, data_id = GID),
                           color = "gray80", size = 0.8, alpha = 0.5) +
    geom_point_interactive(data = filter(all_outer, !is_topk),
                           aes(x = x, y = y, tooltip = tooltip, data_id = GID),
                           color = "gray80", size = 0.8, alpha = 0.5) +
    geom_point_interactive(data = filter(all_inner, is_topk),
                           aes(x = x, y = y, color = color_group,
                               tooltip = tooltip, data_id = GID), size = 2.5) +
    geom_point_interactive(data = filter(all_outer, is_topk),
                           aes(x = x, y = y, color = color_group,
                               tooltip = tooltip, data_id = GID), size = 2.5) +
    geom_text(data = sector_labels,
              aes(x = x, y = y, label = env, hjust = hjust, vjust = vjust),
              size = 3.5, fontface = "bold") +
    scale_color_manual(name = "Set", values = color_vals,
                       breaks = c("TRN", "PRD (CV0)", "PRD (CV00)")) +
    coord_equal(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
    theme_void() +
    theme(legend.text  = element_text(size = 12),
          legend.title = element_text(size = 14, face = "bold"))

  girafe(ggobj = p)
}

cv0_pred_plot <- function(aux_out, k, target_env, top = TRUE, sc = TRUE, cv0_only = FALSE){

  keep_envs <- sort(unique(aux_out %>% filter(Set %in% 'TRN') %>% .$EID))
  prd_candidates <- aux_out %>% filter(EID %in% target_env)
  if (cv0_only) prd_candidates <- prd_candidates %>% filter(new_gen == FALSE)
  quantiles <- prd_candidates %>%
    mutate(topq = rank(-yHat, ties.method = 'random'), botq = rank(yHat, ties.method = 'random'))
  best_gids <- quantiles %>% filter(botq <= k) %>% .$GID
  if (top){best_gids <- quantiles %>% filter(topq <= k) %>% .$GID}

  env_levels <- c(target_env, setdiff(keep_envs, target_env))
  df <- aux_out %>%
    mutate(hl = GID %in% best_gids) %>%
    filter(EID %in% env_levels) %>%
    mutate(
      EID = factor(EID, levels = env_levels),
      plot_group = case_when(
        Set == "TRN" ~ "TRN",
        new_gen       ~ "PRD (CV00)",
        TRUE          ~ "PRD (CV0)"
      ),
      color_group = if_else(hl, plot_group, "other")
    )

  if (sc){
    df <- df %>% group_by(EID) %>% mutate(ymean = mean(yHat, na.rm = TRUE),
                                          ysd = sd(yHat, na.rm = TRUE)) %>%
      ungroup() %>% mutate(yHat = (yHat-ymean)/ysd)
  }

  p <- df %>% ggplot(aes(x = EID, y = yHat, color = color_group, alpha = hl,
                         tooltip = paste0("GID: ", GID, "\nyHat: ", round(yHat, 3)),
                         data_id = GID))+
    geom_point_interactive(data = subset(df, color_group == "other"), size = 3)+
    geom_line(data = filter(df, hl), aes(group = GID), color = 'black')+
    geom_point_interactive(data = subset(df, color_group != "other"), size = 3)+
    scale_alpha_manual(values = c(0.15, 1))+
    scale_color_manual(name = "Set",
                       values = c("TRN" = "#0021A5", "PRD (CV0)" = "#FA4616", "PRD (CV00)" = "#27AE60", "other" = "gray80"),
                       breaks = c("TRN", "PRD (CV0)", "PRD (CV00)"))+
    guides(alpha = 'none')+
    labs(y = 'scaled yHat')+
    theme_light()+
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18, face = 'bold'),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18, face = 'bold'))

  girafe(ggobj = p)
}
