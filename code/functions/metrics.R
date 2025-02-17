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

# Created by: Julian Garcia
# Last updated by: Francisco Gonzalez, Julian Garcia  
# Last updated: 01/10/2025

# This code iteratively gathers all the information in files of the form 'output/[model-name]/[cv].csv'
library(dplyr)
library(ggplot2)
library(scales)
library(viridis)

RMSE = function(x,y){
  keep = intersect(which(!is.na(x)),which(!is.na(y)))
  x = as.numeric(x)[keep]
  y = as.numeric(y)[keep]
  return(sqrt(mean((x - y)^2)))
}

Tiezzi <- function(rs,ns){
  top = 0
  bot = 0
  for (i in 1:length(rs)){
    Vi = (1-rs[i])/(ns[i]-2)
    top = top + rs[i]/Vi
    bot = bot + 1/Vi
  }
  return(top/bot)
}

weighted_mean <- function(x,w){
  return(sum(w*x)/sum(w))
}

create_output_list <- function(eid_col, model_dir, out_dir, metric=NULL, intervals=NULL) {
  models <- setdiff(dir(model_dir), "visualizations")
  output_list <- list()
  
  for (model in models) {
    pre_path = file.path(model_dir, model)
    files = list.files(pre_path, pattern="\\.csv$", recursive = TRUE)
    
    for (file in files) {
      tmp = read.csv(file.path(pre_path, file))
      
      if (!is.null(intervals)) {
        intervals = c('xmin' = min(intervals['xmin'],na.omit(tmp)$predicted),
                      'xmax' = max(intervals['xmax'],na.omit(tmp)$predicted),
                      'ymin' = min(intervals['ymin'],na.omit(tmp)$observed),
                      'ymax' = max(intervals['ymax'],na.omit(tmp)$observed))
      }
      
      if (metric == "accuracy") {
        output = tmp %>% group_by_at(eid_col) %>% 
          summarise(n = n(), cor = round(cor(observed, predicted, use="complete.obs"),2)) %>%
          mutate(model = model, cv = gsub('.csv','',file))
      } else if (metric == "rmse") {
        output = tmp %>% group_by_at(eid_col) %>% 
          summarise(n = n(), rmse = round(RMSE(observed, predicted),2)) %>%
          mutate(model = model, cv = gsub('.csv', '', file))
      }
      
      colnames(output)[1] <- 'env'
      output_list[[paste0(model,"_",file)]] <- output
    }
  }
  return(output_list)
}

get_model_accuracy <- function(eid_col, modeldir, group_by_env=FALSE) {
  
  intervals2 <- c(xmin=Inf, xmax=Inf, ymin=Inf, ymax=Inf)
  out_dir <- file.path(modeldir, "visualizations")
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  
  output_list <- create_output_list(eid_col, modeldir, out_dir, "accuracy", intervals2)
  
  # Accuracy for each environment
  out_df <- do.call("rbind", output_list) %>% 
    mutate(cv = factor(cv, levels = c("cv2", "cv1", "cv0", "cv00")),
           env = factor(env)) %>% dplyr::select(cv, model, env, n, cor)
  
  # Overall accuracy
  tmp2 <- out_df %>% group_by(cv, model) %>% 
    summarise(cor = round(Tiezzi(cor, n), 2), n = sum(n)) %>% 
    mutate(env = "all") %>% dplyr::select(cv, model, env, n, cor)
  
  out_df <- rbind(out_df, tmp2) %>% arrange(cv, model, env)
  
  # Save table and auxiliar objects for Shiny
  rownames(out_df) <- NULL
  write.csv(out_df, file = file.path(out_dir, "accuracy.csv"))
  
  if(group_by_env) {
    retplot <- out_df %>%
      ggplot(aes(x=env, y=cor, fill=model)) +
      geom_col(position = 'dodge', aes(alpha=env=="all", color=env=="all")) +
      facet_grid(cv~model) +
      scale_alpha_manual(values = c(0.7, 1)) +
      scale_fill_viridis_d() +
      scale_color_manual(values = c("transparent", "black")) +
      labs(x="Environment", y="Correlation", fill="Model") +
      theme(axis.text.x = element_text(size = 8, angle=40, hjust = 1)) +
      labs(x="Environment ID", y="Prediction Accuracy (out of 1)", fill="Model")+
      guides(alpha = "none", color = "none")
    
    plot_name <- "accuracy_by_env"
    
  } else {
    retplot <- out_df %>% filter(env=="all") %>% 
      ggplot(aes(x=model, y=cor, fill=model)) +
      geom_col(position = "dodge", color = "black") +
      facet_grid(~cv) + 
      scale_fill_viridis_d() +
      labs(x="", y="Prediction Accuracy (out of 1)", fill="Model") +
      theme(axis.text.x = element_text(size = 8, angle=40, hjust = 1)) +
      theme(axis.text.x = element_blank())
    
    plot_name <- "accuracy_overall"
  }
  
  # Save plot
  ggsave(retplot, filename = file.path(out_dir, paste0(plot_name,".png")),width=2100,height=2100,units=c("px"))
  saveRDS(retplot, file.path(out_dir, paste0(plot_name, ".rds")))
  
}

get_model_rmse <- function(eid_col, model_dir, group_by_env=FALSE) {
  out_dir <- file.path(model_dir, "visualizations")
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  
  output_list = create_output_list(eid_col, model_dir, out_dir, "rmse")
  
  out_df = do.call('rbind', output_list) %>% 
    mutate(cv = factor(cv, levels = c('cv2','cv1','cv0','cv00')),
           env = factor(env)) %>% dplyr::select(cv, model, env, n, rmse)
  
  tmp2 = out_df %>% group_by(cv, model) %>%
    summarise(rmse = round(weighted.mean(rmse, n),2), n = sum(n)) %>%
    mutate(env = "all") %>% dplyr::select(cv, model, env, n, rmse)
  
  out_df = rbind(out_df, tmp2) %>% arrange(cv, model, env)
  
  # Save the table and auxiliar objects for shiny
  rownames(out_df) = NULL
  write.csv(out_df, file = file.path(out_dir, "rmse.csv"))
  
  if (group_by_env) {
    retplot <- out_df %>% ggplot(aes(x=env, y=rmse, fill=model)) +
      geom_col(position = "dodge", aes(alpha = env=="all", color = env=="all")) +
      facet_grid(cv~model) + 
      scale_alpha_manual(values = c(0.95, 1)) +
      scale_fill_viridis_d() +
      scale_color_manual(values = c("transparent", "black")) +
      labs(x="Environment ID", y="RMSE", fill="Model") + 
      theme(axis.text.x = element_text(size = 8, angle=40, hjust = 1)) +
      guides(alpha = "none", color = "none")
    
    plot_name <- "rmse_by_env"
  } else {
    retplot <- out_df %>% filter(env == "all") %>%
      ggplot(aes(x=model, y=rmse, fill=model)) +
      geom_col(position = "dodge", color = "black") +
      facet_grid(~cv) +
      scale_fill_viridis_d() +
      labs(x="", y="RMSE", fill="Model") +
      theme(axis.text.x = element_text(size = 8, angle=40, hjust = 1)) +
      theme(axis.text.x = element_blank())
    
    plot_name <- "rmse_overall"
  }
  
  # Save the plot
  ggsave(retplot, filename = file.path(out_dir, paste0(plot_name,".png")),width=2100,height=2100,units=c("px"))
  saveRDS(retplot, file.path(out_dir, paste0(plot_name, ".rds")))
  
}

get_model_varcomps <- function(modeldir) {
  out_dir = file.path(modeldir, "visualizations")
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  
  models = setdiff(dir(modeldir), "visualizations")
  output_list = list()
  curr_levels = c()
  
  for (model in models) {
    fm = get(load(file.path(modeldir, model, "fm.rda")))
    varcomps = c()
    
    for(i in seq_along(fm$ETA)) {
      varcomps = c(varcomps, fm$ETA[[i]]$varU)
    }
    
    output = data.frame(model = model, 
                        term = c(gsub('_',"*",fm$equation), "residual"),
                        value = c(varcomps, fm$varE)
    )
    
    curr_levels = unique(c(curr_levels, output$term))
    output$varcomp = 100*output$value/sum(output$value)
    output_list[[model]] = output
  }
  
  curr_levels = curr_levels[-which(curr_levels == "residual")]
  curr_levels = c(curr_levels, "residual")
  outdf = do.call("rbind", output_list)
  
  n_terms = length(curr_levels) - 1
  viridis_colors <- viridis(n_terms)
  
  retplot <- outdf %>% mutate(term = factor(term, levels = curr_levels),
                       model = factor(model, levels = models)) %>%
    ggplot(aes(x=model, y=value, fill=term)) +
    geom_bar(position = "fill", stat = "identity", color = "black") +
    scale_fill_manual(values = c(viridis_colors, '#c7c9c8')) +
    labs(y="Variance explained (%)", x="Model") + 
    theme(axis.text.x = element_text(size = 8, angle=40, hjust = 1)) +
    theme_classic()
  
  rownames(outdf) = NULL
  colnames(outdf)[4] = c("percentage")
  write.csv(outdf, file = file.path(out_dir, "varcomps.csv"))
  
  ggsave(retplot, file = file.path(out_dir, "varcomps.png"),width=2100,height=2100,units=c("px"))
  saveRDS(retplot, file.path(out_dir, "varcomps.rds"))
  
  return(retplot)
}