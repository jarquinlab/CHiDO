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

library(dplyr)
library(ggplot2)
library(scales)
library(viridis)
library(gridExtra)

average_scores<-function(outdir,k=2,r,cc){
  # k represents the amount of principal components (2)
  alphama <- read.csv(file.path(outdir, 'alpha.csv'), sep = ",")
  col.names <- c("Iteration",paste("alpha_",tratamientos,"_",rep(1:k,rep(r,k)),sep=""))

  alphama <- alphama[,col.names]
  gamma  <- read.csv(file.path(outdir, 'gamma.csv'), sep = ",")
  col.names <- c("Iteration",paste("gamma_",sitios,"_",rep(1:k,rep(cc,k)),sep=""))
  gamma <- gamma[,col.names]
  
  lambda <- read.csv(file.path(outdir, 'lambda.csv'), sep = ",")
  col.names <- c("Iteration","lambda_1","lambda_2")
  lambda <- lambda[,col.names]
  
  svd.mcmc <- merge(alphama,gamma,"Iteration")
  svd.mcmc <- merge(svd.mcmc,lambda,"Iteration")[-1]
  
  row2mat <- function(x,nrow,ncol)matrix(as.numeric(x),nrow=nrow,ncol=ncol,byrow=F)

  #nsim <- nrow(svd.mcmc)
  #s <- (niter-burnin)/thin
  #sample.svd.mcmc <- svd.mcmc[sample(1:nsim,s),]
  sample.svd.mcmc<-svd.mcmc
  scores.mcmc <- list()
  for(i in 1:nrow(sample.svd.mcmc))
  {
    Ut <- row2mat(sample.svd.mcmc[i,1:(2*r)],r,2)
    Vt <- row2mat(sample.svd.mcmc[i,(2*r+1):(2*r+2*cc)],cc,2)
    rownames(Vt) <- paste('E',1:cc,sep='')
    Dt <- diag(as.numeric(sample.svd.mcmc[i,(2*r+2*cc+1):(2*r+2*cc+2)]),ncol=2)
    if(i==1)
    {
      Xt <- Ut%*%Dt%*%t(Vt)
      dimnames(Xt) <- list(tratamientos,sitios)
    }
    res <- bpca.default2(Xt,Ut,Vt,diag(Dt),method="sqrt")
    scores <- scores.bpca.2d(res)
    scores.mcmc[[i]] <- scores
    #  print(i)
  }
  b<-apply(array(unlist(scores.mcmc), c(dim(scores.mcmc[[1]]), length(scores.mcmc))), 2, rowMeans)
  colnames(b)<-c("Component 1","Component 2")
  rownames(b)<-rownames(scores)
  
  if (sum(duplicated(c(tratamientos,sitios)))>0){rownames(b)<-c(paste0("G",1:r),paste0("E",1:cc),1)}
  write.table(b,file=file.path(outdir, 'average_scores.csv'),row.names = TRUE,col.names=TRUE, sep = ",")

  #write.table(b,file='average_scores.csv',row.names = TRUE,col.names=TRUE) #remove
  
  return(b)
}

perc_expl<-function(outdir){
  ev <- as.matrix(read.csv(file.path(outdir, 'lambda.csv'), sep = ","))
  
  ev<-as.matrix(ev[,-1])^2
  ev<-t(apply(ev,1,function (row) row/sum(row)))
  ev<-apply(ev[,1:2],2,mean)*100
  
  names(ev) <- c("CP1","CP2")
  
  return(ev)
}

bayesian_biplot<-function(tmpdir,scores,lambda,r){
  plot_dir <- file.path(tmpdir, "output/visualizations")
  if (!dir.exists(plot_dir)) { dir.create(plot_dir, recursive = TRUE) }
  
  outdf<-as.data.frame(scores[1:(nrow(scores)-1),])
  colnames(outdf)<-c("CP1","CP2")
  
  perc<-perc_expl(file.path(tmpdir,"output/AMMI"))
  
  retplot <- ggplot(outdf, aes(x=CP1, y=CP2)) +
    ggtitle("Biplot - Bayesian AMMI Analysis") + 
    geom_point() + # Show dots
    xlab(paste0(round(perc[1],2),"%",sep="")) +
    ylab(paste0(round(perc[2],2),"%",sep="")) +
    
    geom_hline(aes(yintercept = 0))  +  
    geom_vline (aes(xintercept = 0))  +

    lapply((r+1):nrow(outdf), function(i) {
      geom_segment(
        data = outdf[i, ], 
        aes(x = 0, y = 0, xend = outdf[i, 1], yend = outdf[i,2 ]),
        arrow = arrow(length = unit(0.3, "cm")), 
        color = "darkblue"
      ) 
    })+
    geom_text(label=rownames(outdf),
              nudge_x = 0.02, nudge_y = 0.03,  
              check_overlap = T) +
    theme_minimal()+
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill="white"))
  
  # Save the plot
  ggsave(retplot, filename = file.path(plot_dir, paste0("bayesian_ammi",".png")),width=2100,height=2100,units=c("px"))
  saveRDS(retplot, file.path(plot_dir, paste0("bayesian_ammi", ".rds")))
  
  #ggsave(retplot, filename = "bayesian_ammi.png",width=2100,height=2100,units=c("px")) # remove
  
}

trace_plots<-function(tmpdir,tau,mu,lambda){
  plot_dir <- file.path(tmpdir, "output/visualizations")
  if (!dir.exists(plot_dir)) { dir.create(plot_dir, recursive = TRUE) }
  
  output<-as.data.frame(tau) %>% reframe(Iteration=1:length(tau),tau=tau,mu=mu,lambda1=lambda[,1],lambda2=lambda[,2])
  
  retplot<-output %>% 
    ggplot(aes(x=Iteration,y=tau))+
    geom_line()+
    labs(y=expression(sigma^2))+
    theme_minimal()
  retplot2<-output %>% 
    ggplot(aes(x=tau))+
    geom_histogram(color='black',fill='gray')+
    labs(y="Frequency",x=expression(sigma^2))+
    theme_minimal()
  tau_plot<-grid.arrange(retplot,retplot2,ncol=2,nrow=1)
  
  # Save the plot
  ggsave(tau_plot, filename = file.path(plot_dir, paste0("trace_tau",".png")),width=2100,height=2100,units=c("px"))
  saveRDS(tau_plot, file.path(plot_dir, paste0("trace_tau", ".rds")))
  
  #ggsave(tau_plot, filename = "trace_tau.png",width=2100,height=2100,units=c("px"))
  
  
  retplot<-output %>% 
    ggplot(aes(x=Iteration,y=mu))+
    geom_line()+
    labs(y=expression(mu))+
    theme_minimal()
  retplot2<-output %>% 
    ggplot(aes(x=mu))+
    geom_histogram(color='black',fill='gray')+
    labs(y="Frequency",x=expression(mu))+
    theme_minimal()
  mu_plot<-grid.arrange(retplot,retplot2,ncol=2,nrow=1)
  # Save the plot
  ggsave(mu_plot, filename = file.path(plot_dir, paste0("trace_mu",".png")),width=2100,height=2100,units=c("px"))
  saveRDS(mu_plot, file.path(plot_dir, paste0("trace_mu", ".rds")))
  
  #ggsave(mu_plot, filename = "trace_mu.png",width=2100,height=2100,units=c("px")) # remove
  
  retplot<-output %>% 
    ggplot(aes(x=Iteration,y=lambda1))+
    geom_line()+
    labs(y=expression(lambda[1]))+
    theme_minimal()
  retplot2<-output %>% 
    ggplot(aes(x=lambda1))+
    geom_histogram(color='black',fill='gray')+
    labs(y="Frequency",x=expression(lambda[1]))+
    theme_minimal()
  lambda1_plot<-grid.arrange(retplot,retplot2,ncol=2,nrow=1)
  # Save the plot
  ggsave(lambda1_plot, filename = file.path(plot_dir, paste0("trace_lambda1",".png")),width=2100,height=2100,units=c("px"))
  saveRDS(lambda1_plot, file.path(plot_dir, paste0("trace_lambda1", ".rds")))
  
  #ggsave(lambda1_plot, filename = "trace_lambda1.png",width=2100,height=2100,units=c("px"))
  
  retplot<-output %>% 
    ggplot(aes(x=Iteration,y=lambda2))+
    geom_line()+
    labs(y=expression(lambda[2]))+
    theme_minimal()
  retplot2<-output %>% 
    ggplot(aes(x=lambda2))+
    geom_histogram(color='black',fill='gray')+
    labs(y="Frequency",x=expression(lambda[2]))+
    theme_minimal()
  lambda2_plot<-grid.arrange(retplot,retplot2,ncol=2,nrow=1)
  
  # Save the plot
  ggsave(lambda2_plot, filename = file.path(plot_dir, paste0("trace_lambda2",".png")),width=2100,height=2100,units=c("px"))
  saveRDS(lambda2_plot, file.path(plot_dir, paste0("trace_lambda2", ".rds")))
  
  #ggsave(lambda2_plot, filename = "trace_lambda2.png",width=2100,height=2100,units=c("px"))
  
}

waas_plot<-function(data,tmpdir,r){
  plot_dir <- file.path(tmpdir, "output/visualizations")
  if (!dir.exists(plot_dir)) { dir.create(plot_dir, recursive = TRUE) }
  
data$Genotype <- factor(data$Genotype)
  
  waas <- ggplot(data, aes(x = as.factor(Genotype), y = BWAAS)) +
    geom_point(aes(color = BWAAS), size = 3) +  # Color points based on B_WAAS value
    geom_errorbar(aes(ymin = CI_0.10, ymax = CI_0.90), width = 0.2, color = "black") +  # Black error bars
    labs(title = "",
         x = "Genotype", y = "Bayesian WAAS") +
    scale_x_discrete(breaks = levels(data$Genotype)) +
    theme_minimal() +
    scale_color_viridis_c()+
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major=element_line(color="gray",linewidth=0.5,linetype=2))
  
  # Save the plot
  ggsave(waas, filename = file.path(plot_dir, paste0("waas_plot",".png")),width=2100,height=2100,units=c("px"))
  saveRDS(waas, file.path(plot_dir, paste0("waas_plot", ".rds")))
  
  #ggsave(waas, filename = "waas_plot.png",width=2100,height=2100,units=c("px"))
  
}

waasy_plot<-function(data,tmpdir,r){
  plot_dir <- file.path(tmpdir, "output/visualizations")
  if (!dir.exists(plot_dir)) { dir.create(plot_dir, recursive = TRUE) }
  
  data$Genotype <- factor(data$Genotype)
  
  waasy<-ggplot(data, aes(x = Genotype, y = BWAASY)) +
    geom_point(aes(color = BWAASY), size = 3) +  # Color points based on B_WAASY value
    geom_errorbar(aes(ymin = CI_0.10, ymax = CI_0.90), width = 0.2, color = "black") +  # Black error bars
    labs(title = "",
         x = "Genotype", y = "Bayesian WAASY") +
    scale_x_discrete(breaks = levels(data$Genotype)) +
    theme_minimal() +
    scale_color_viridis_c()+
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major=element_line(color="gray",linewidth=0.5,linetype=2))
  
  # Save the plot
  ggsave(waasy, filename = file.path(plot_dir, paste0("waasy_plot",".png")),width=2100,height=2100,units=c("px"))
  saveRDS(waasy, file.path(plot_dir, paste0("waasy_plot", ".rds")))
  
  #ggsave(waasy, filename = "waasy_plot.png",width=2100,height=2100,units=c("px"))
  
}

stabdist_plot<-function(data,tmpdir,r){
  plot_dir <- file.path(tmpdir, "output/visualizations")
  if (!dir.exists(plot_dir)) { dir.create(plot_dir, recursive = TRUE) }
  
  retplot <- ggplot(data = data, aes(x = B_Geno_Mean, y = B_StabM)) +
    geom_point() +  
    geom_errorbar(aes(ymin = S_CI_0.10, ymax = S_CI_0.90), width = 0.1, color = "#8FD744FF") + 
    geom_errorbar(aes(xmin = CI_0.10, xmax = CI_0.90), color =  "#31688EFF") + 
    geom_label(aes(label = Genotype), nudge_x = 0.03, nudge_y = 0.03) + 
    labs(x = "Trait Posteriori Mean", y = "Bayesian Stability Mahalanobis Distance") +
    ggtitle("") +  
    geom_vline(xintercept = quantile(data[, 3], probs = c(0.5)), linetype = "dashed", color = "gray")+
    geom_hline(yintercept = quantile(data[, 6], probs = c(0.5)), linetype = "dashed", color = "gray") +
    theme_minimal()  +
    geom_text(x = as.vector(mean(data[,3])+sd(data[,3])), y = as.vector(mean(data[,6])-sd(data[,6])), label= "IV", colour= "black", family = "serif")+ 
    geom_text(x = as.vector(mean(data[,3])-sd(data[,3])), y = as.vector(mean(data[,6])-sd(data[,6])), label= "III", colour= "black", family = "serif")+
    geom_text(x = as.vector(mean(data[,3])-sd(data[,3])), y = as.vector(mean(data[,6])+sd(data[,6])), label= "II", colour= "black", family = "serif")+
    geom_text(x = as.vector(mean(data[,3])+sd(data[,3])), y = as.vector(mean(data[,6])+sd(data[,6])), label= "I", colour="black", family = "serif")+
    scale_y_continuous(limits = c(0, ceiling(max(data[,6])+sd(data[,6]))))+
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill="white"))
  
  # Save the plot
  ggsave(retplot, filename = file.path(plot_dir, paste0("stabdist_plot",".png")),width=2100,height=2100,units=c("px"))
  saveRDS(retplot, file.path(plot_dir, paste0("stabdist_plot", ".rds")))
  
  #ggsave(retplot, filename = "stabdist_plot.png",width=2100,height=2100,units=c("px"))
}