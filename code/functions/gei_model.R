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

library(lme4)
library(Matrix)
library(mvtnorm)
library(MASS)
library(rstiefel)
library(tmvtnorm)
library(coda)

bilinear_model<-function(archivo,colPhen,colSite,colEntry,colRep){
  
  cultivo <- archivo
  if(!is.na(colRep)){
  cultivo$rep<-as.factor(cultivo[,colRep])
  }
  
  cultivo$site <- as.factor(cultivo[,colSite])
  cultivo$entry <- as.factor(cultivo[,colEntry])
  cultivo$y <- cultivo[,colPhen]
  
  if(!is.na(colRep)){
  cultivo$rep <- as.factor(cultivo[,colRep])
  modelo<<-lmer(y~site+entry+site:entry+(1|site:rep),cultivo)
  }else{
  modelo<<-lm(y~site+entry,cultivo)
  }
  
  sigma2<<-sigma(modelo)^2
  
  #tratamientos<<-matrix(as.numeric(levels(cultivo$entry)),ncol=1)
  tratamientos<<-matrix((levels(cultivo$entry)),ncol=1)
  #sitios<<-matrix(as.numeric(levels(cultivo$site)),ncol=1)
  sitios<<-matrix((levels(cultivo$site)),ncol=1)
  
  num_sitios<<-length(sitios)
  num_tratam<<-length(tratamientos)
  
  ce<-num_sitios
  r<<-num_tratam
  cc<<-ce
  
  
  jr<<-matrix(rep(1,r),ncol=1)
  jc<<-matrix(rep(1,ce),ncol=1)
  
  if(!is.na(colRep)){
  cultivo<<-cultivo[ order(cultivo$entry,cultivo$site,cultivo$rep),]
  }else{
  cultivo<<-cultivo[ order(cultivo$entry,cultivo$site),]
  }
  
  matriz_de_medias <- matrix(NaN,nrow=r,ncol=cc)
  for(i in 1:r)
  {
    for(j in 1:cc)
    {
      matriz_de_medias[i,j] <- mean(cultivo$y[cultivo$entry==tratamientos[i] & cultivo$site==sitios[j]] )
    }
  }
  matriz_de_medias <<- matriz_de_medias
  colnames(matriz_de_medias)<-sitios
  rownames(matriz_de_medias)<-tratamientos
  
  sigm<<-function(muda){
    (E<-
       matriz_de_medias
     -mu*jr%*%t(jc)
     -matrix(tao,ncol=1)%*%t(jc) 
     -jr%*%matrix(delta,ncol=ce)
     -alphas%*%diag(lambdas)%*%t(gammas)
    )
  }
  
  if(!is.na(colRep)){
  ne<<-length(levels(cultivo$rep))
  }else{
  ne<<-1
  }
  
  k<<-min(dim(matriz_de_medias))-1

  
  gran_media<<-mean(matriz_de_medias)
  mu<<-gran_media
  
  media_de_tratam <<- apply(matriz_de_medias,1,mean) 
  tao <- matrix(media_de_tratam - mean(matriz_de_medias)  , ncol=1)
  
  media_de_sitio <<- apply(matriz_de_medias,2,mean)
  delta <- matrix(media_de_sitio - mean(matriz_de_medias),ncol=1)
  
  
  zeta3 <<- matriz_de_medias - mu - tao%x%t(jc) - t(delta)%x%jr
  
  
  inicio2<<-svd(zeta3)
  lambdas<<-inicio2$d[1:k]
  alphas<<-inicio2$u[,1:k]
  gammas<<-inicio2$v[,1:k]
  
  
  matriz_de_varianzas <- matrix(NaN,nrow=r,ncol=cc)
  for(i in 1:r)
  {
    for(j in 1:cc)
    {
      matriz_de_varianzas[i,j] <- var( cultivo$y[cultivo$entry==tratamientos[i] & cultivo$site==sitios[j]] )
    }
  }
  rownames(matriz_de_varianzas)<-tratamientos
  colnames(matriz_de_varianzas)<-sitios
  ssy<<-sum(matriz_de_varianzas,na.rm=T)
  
  
  ########################

  ######################
  
  matriz_desv_model<<-function(nada){
    iter<-1
    cultivose2<-NULL
    while(iter<=r)
    {
      cultivos<-cultivo[cultivo$entry==tratamientos[iter],]
      iter2<-1
      cultivose1<-NULL
      while(iter2<=cc)
      {

        cultivose<-((mean(cultivos[cultivos$site==sitios[iter2],]$y)-mu-tao[iter]-delta[iter2]-suma_parci3(lambdas)[iter,iter2])^2)
        cultivose1<-cbind(cultivose1,cultivose)
        iter2<-iter2+1
      }
      cultivose2<-rbind(cultivose2,cultivose1)
      iter<-iter+1
    }
    rownames(cultivose2)<-tratamientos
    colnames(cultivose2)<-sitios
    return(cultivose2)
  }
  
  ##### Functions ####
  matriz_k<-function(ne){
    j<-0
    kn22<-NULL
    while(j<ne-1){
      long<-sqrt((ne-j-1)*(ne-j))	
      i<-0
      kn2<-NULL
      while(i<ne){
        if(i==j) k=(ne-j-1) else k=NA
        if(i<j)  k=0 else k=k
        if(i>j)  k=-1 else k=k
        kn<-k/long
        kn2<-cbind(kn2,kn)
        i<-i+1	
      }
      kn22<-rbind(kn22,kn2)
      j<-j+1	
    }
    return(kn22)
  }
  
  kr<<-t(matriz_k(r))
  kc<<-t(matriz_k(ce))
  
  
  lam<<-function(muda){
    
    med_lamb<-function(iter)(ne*sig_lamb_ap[iter]*sum(alphas[,iter]%*%t(gammas[,iter])*matriz_de_medias)+mu_lamb_ap[iter]*sigma2)/(ne*sig_lamb_ap[iter]+sigma2)
    desv_lamb<-function(sigma2)sqrt(sig_lamb_ap[iter]*sigma2/(ne*sig_lamb_ap[iter]+sigma2))
    
    for(iter in 1:k){
      phi<-0
      uni<-1
      u_neg<- -med_lamb(iter)/desv_lamb(sigma2)
      alpha_star <- (u_neg+sqrt(u_neg^2+4))/2
      if(iter!=1){
        while(phi<uni){
          u_pos<-(lambdas[iter-1]-med_lamb(iter))/desv_lamb(sigma2)
          z<-runif(1,u_neg,u_pos)
          phi<-exp((ifelse(u_pos<0,u_pos^2,ifelse(u_neg>0,u_neg^2,0))-z^2)/2)				
          uni<-runif(1) 
        }
      } 
      else{
        while(phi<uni){
          z <- rexp(1,alpha_star)+u_neg
          phi <- exp(-(z-alpha_star)^2/2)
          uni<-runif(1)
        }
      }
      lambdas[iter]<-desv_lamb(sigma2)*z+med_lamb(iter)		
    }
    return(lambdas)
  }
  
  
  lamI<<-function(mean_lambI,var_dr)
  {
    
    med_lamb<-function(iter) mean_lambI[iter]
    desv_lamb<-sqrt(var_dr)
    lambdas_m<- rep(0,k)
    
    for(iter in 1:k){
      phi<-0
      uni<-1
      u_neg<- -med_lamb(iter)/desv_lamb
      alpha_star <- (u_neg+sqrt(u_neg^2+4))/2
      
      if(iter!=1){
        while(phi<uni)
        {                   
          u_pos<-(lambdas_m[iter-1]-med_lamb(iter))/desv_lamb
          z<-runif(1,u_neg,u_pos)
          phi<-exp((ifelse(u_pos<0,u_pos^2,ifelse(u_neg>0,u_neg^2,0))-z^2)/2)				
          uni<-runif(1) 
        }
      } 
      else{
        while(phi<uni)
        {
          z <- rexp(1,alpha_star)+u_neg
          phi <- exp(-(z-alpha_star)^2/2)
          uni<-runif(1)
        }
      }
      
      lambdas_m[iter]<-desv_lamb*z+med_lamb(iter)
      
    }
    return(lambdas_m)
  }
  
  
  orthonormalization2 <<- function (u = NULL, basis = TRUE, norm = TRUE) 
  {
    if (is.null(u)) 
      return(NULL)
    if (!(is.matrix(u))) 
      u <- as.matrix(u)
    p <- nrow(u)
    ne <- ncol(u)
    if (p < ne) {
      warning("too much vectors to orthogonalize.")
      u <- as.matrix(u[, 1:p])
      ne <- p
    }
    if (basis & (p > ne)) {
      base <- diag(p)
      coef.proj <- crossprod(u, base)/diag(crossprod(u))
      base2 <- base - u %*% matrix(coef.proj, nrow = ne, ncol = p)
      norm.base2 <- diag(crossprod(base2))
      base <- as.matrix(base[, order(norm.base2) > ne])
      u <- cbind(u, base)
      ne <- p
    }
    v <- u
    if (ne > 1) {
      for (i in 2:ne) {
        coef.proj <- c(crossprod(u[, i], v[, 1:(i - 1)]))/diag(crossprod(v[, 
                                                                           1:(i - 1)]))
        v[, i] <- u[, i] - matrix(v[, 1:(i - 1)], nrow = p) %*% 
          matrix(coef.proj, nrow = i - 1)
      }
    }
    if (norm) {
      coef.proj <- 1/sqrt(diag(crossprod(v)))
      v <- t(t(v) * coef.proj)
    }
    return(v)
  }
  
  
  rmf.matrix.gibbs4 <<- function (M, X, rscol = NULL) 
  {
    if (is.null(rscol)) {
      rscol <- max(2, min(round(log(dim(M)[1])), dim(M)[2]))
    }
    sM <- svd(M)
    H <- sM$u %*% diag(sM$d)
    Y <- X %*% sM$v
    m <- dim(H)[1]
    R <- dim(H)[2]
    for (iter in 1:round(R/rscol)) {
      r <- 1:R          
      N <- NullC(Y[, -r])
      y <- rmf.matrix(t(N) %*% H[, r])
      Y[, r] <- N %*% y
    }
    
    entrada<-(Y%*%t(sM$v))
    entrada1<-cbind(rep(1,m),entrada)
    salida2<-orthonormalization2(entrada1, basis=TRUE, norm=TRUE)
    salida3<-salida2[,2:(k+1)]
    return(salida3)
  }
  
  
  aa<-NULL
  for(i in 1:r)
  {
    for(j in 1:cc)
    {
      lvd<-0;for(lk in 1:k){lvd <-lvd + lambdas[lk]*alphas[i,lk]%*%t(gammas[j,lk])} 
      a0 <- cultivo[ cultivo$site==sitios[j] & cultivo$entry==tratamientos[i],   ]$y
      
      aa<- rbind (aa, 
                  matrix(cultivo[ cultivo$site==sitios[j] & cultivo$entry==tratamientos[i],   ]$y)-matrix(rep( mu+ tao[i,] +delta[j,]+lvd,length(a0) )) )
    }
  }
  
  Zc<-matrix(0,r,cc)
  (sum(aa^2) )
  salidas<<-
    list(mu,  1/sigma2,ne,r,cc,k,jr,jc,matriz_de_medias,ssy,lambdas,alphas,gammas,sum(sigm(1)^2),tao,delta,kr,kc,zeta3,Zc,media_de_tratam,media_de_sitio)
  return(salidas)
}

suma_parci3<-function(lambdas){
  iter4<-1
  suma_parcial<-lambdas[1]*alphas[,1]%*%t(gammas[,1])
  iter4<-iter4+1
  while(iter4<=k)
  {
    suma_parcial<-suma_parcial+lambdas[iter4]*alphas[,iter4]%*%t(gammas[,iter4])
    iter4<-iter4+1
  }
  return(suma_parcial)
}

gibbs_chain<-function(exit1,nombre1,niter,times,alphas0,gammas0,lambdas0){
  while(times <= niter)
  {
    mu <<-rnorm(1,(ne*r*cc*sigma2_u*gran_media)/(ne*r*cc*sigma2_u+sigma2),sqrt((sigma2*sigma2_u)/(ne*r*cc*sigma2_u+sigma2)))
    tao <<-kr%*%mvrnorm(1,(cc*ne*sig_tao_ap*t(kr)%*%(media_de_tratam)+sigma2*t(kr)%*%mu_tao_ap)/(cc*ne*sig_tao_ap+sigma2),sqrt((sig_tao_ap*sigma2)/(cc*ne*sig_tao_ap+sigma2)*diag(rep(1,r-1))))

    delta <<-kc%*%mvrnorm(1,(cc*ne*sig_del_ap*t(kc)%*%(media_de_sitio)+sigma2*t(kc)%*%mu_del_ap)/(cc*ne*sig_del_ap+sigma2),sqrt((sig_del_ap*sigma2)/(cc*ne*sig_del_ap+sigma2)*diag(rep(1,cc-1))))
    
    mean_lambI <- diag(t(gammas)%*%t(matriz_de_medias)%*%alphas)
    lambdas <- lamI(mean_lambI=mean_lambI,var_dr =(ne*1/sigma2)^(-1)) 
    
    alphas <- rmf.matrix.gibbs4( ne * matriz_de_medias%*%gammas%*%diag(lambdas)/sigma2,alphas )[,1:k]  

    gammas <- rmf.matrix.gibbs4(  ne * t(matriz_de_medias)%*%alphas%*%diag(lambdas)/sigma2,gammas)[,1:k] 

    sigma2 <<-1/rgamma(1,r*ne*cc/2+alfa_sig_ap,.5*(ne*sum(matriz_desv_model(1))+ssy)+beta_sig_ap)

    for(i in 1:ncol(alphas))
    {
      maxal0 <- which.max(abs(alphas0[,i]))
      maxga0 <- which.max(abs(gammas0[,i]))
      sgn1 <- (alphas0[,i]<0)!=(alphas[,i]<0)
      sgn2 <- (gammas0[,i]<0)!=(gammas[,i]<0)
      sgn3 <- (alphas0[maxal0,i]<0)!=(alphas[maxal0,i]<0)&&(abs(alphas[maxal0,i])>19.5/20)
      sgn4 <- (gammas0[maxga0,i]<0)!=(gammas[maxga0,i]<0)&&(abs(alphas[maxal0,i])>19.5/20)
      if(any(all(sgn1),all(sgn2),sgn3,sgn4))
      {
        alphas[,i] <- -1*alphas[,i]
        gammas[,i] <- -1*gammas[,i]
      }
    }

    exit<-c(mu,sigma2,as.vector(tao),as.vector(delta),as.vector(alphas),as.vector(gammas),as.vector(lambdas) )
    exit1<-rbind(exit1,exit)
    
    if(times%%100==0) cat("iter =",times,"\n")
    if(times==2) write.table(exit1,file=nombre1,append = TRUE, quote = TRUE, sep = " , ",row.names = FALSE,col.names = FALSE)
    if(times%%100==0) write.table(exit1,file=nombre1,append = TRUE, quote = TRUE, sep = " , ",row.names = FALSE,col.names = FALSE)
    if(times%%100==0) exit1<-NULL
    
    #     print(times)
    times <- times + 1 
  }
}

burnin_filter<-function(chain_list,burnin,r,cc,k){
  col.names <- c("mu","sigma2-1",paste("tau",tratamientos,sep="_"),  
                 paste("delta",sitios,sep="_"),
                 paste("alpha_",tratamientos,"_",rep(1:k,rep(r,k)),sep=""),
                 paste("gamma_",sitios,"_",rep(1:k,rep(cc,k)),sep=""),
                 paste("lambda",1:k,sep="_")#,paste("tau",1:6,sep="_")
  )
  
  burn_in<-lapply(chain_list,function (df){
    df<-df[-seq_len(burnin),]
    colnames(df)<-col.names
    
    return(df)
  })
  return(burn_in)
}

thin_filter<-function(chain_list,nIter,burnIn,thin){
  ind <- seq(1,nIter-burnIn,by=thin)

  thin.list<-lapply(chain_list,function (df){
    df<-df[ind,]
    
    return(df)
  })
  return(thin.list)
}

hpd_intervals<-function(chain,outdir){
  S1 <- summary(as.mcmc(chain))
  S2 <- HPDinterval(as.mcmc(chain),prob=0.90)
  S3 <- HPDinterval(as.mcmc(chain))
  colnames(S2) <- paste(colnames(S2),.9)
  colnames(S3) <- paste(colnames(S3),.95)
  
  full.S <- cbind(S1[[1]][,1:2],S1[[2]][,c(2:4)],S2,S3)
  
  write.table(full.S,file=file.path(outdir, 'full.S.csv'),row.names = F,col.names=T,sep=",")
  
  return(full.S)
}

saving_parameters<-function(outdir,mu,tau,tao,beta,lambda,alpha,gamma){
  Iteration <- 1:length(mu)
  
  write.table(list(Iteration=Iteration,mu),file=file.path(outdir, 'mu.csv'),row.names = FALSE,col.names=TRUE, sep = ",")
  write.table(list(Iteration=Iteration,tau),file=file.path(outdir, 'tau.csv'),row.names = FALSE,col.names=TRUE, sep = ",")
  write.table(list(Iteration=Iteration,tao),file=file.path(outdir, 'tao.csv'),row.names = FALSE,col.names=TRUE, sep = ",")
  write.table(list(Iteration=Iteration,beta),file=file.path(outdir, 'beta.csv'),row.names = FALSE,col.names=TRUE, sep = ",")
  write.table(list(Iteration=Iteration,lambda),file=file.path(outdir, 'lambda.csv'),row.names = FALSE,col.names=TRUE, sep = ",")
  write.table(list(Iteration=Iteration,alpha),file=file.path(outdir, 'alpha.csv'),row.names = FALSE,col.names=TRUE, sep = ",")
  write.table(list(Iteration=Iteration,gamma),file=file.path(outdir, 'gamma.csv'),row.names = FALSE,col.names=TRUE, sep = ",")
  
  #write.table(list(Iteration=Iteration,tao),file="tao.csv",row.names = FALSE,col.names=TRUE, sep = ",")
  #write.table(list(Iteration=Iteration,tau),file="tau.csv",row.names = FALSE,col.names=TRUE, sep = ",")
  #write.table(list(Iteration=Iteration,mu),file="mu.csv",row.names = FALSE,col.names=TRUE, sep = ",")
  
}

bpca.default2 <-function(x, u, v, d, lambda.ini=1, lambda.end=2,
           var.position=2, var.center=TRUE, var.scale=TRUE,
           method=c('hj', 'sqrt', 'jk', 'gh'),
           var.rb=FALSE, var.rd=FALSE, limit=10, ...){
    stopifnot(is.matrix(x) || is.data.frame(x))
    
    n.lambdas <- (lambda.end - lambda.ini + 1)
    if(n.lambdas < 2 || n.lambdas > 3)
      stop('Please, check the parameters: lambda.ini and lambda.end:\n',
           'It should be 2 (to bpca.2d) or 3 (to bpca.3d).\n\n')
    
    if(!any(var.position == 1:2))
      stop('Please, check the parameters: var.position:\n',
           'It should be 1 (rows) or 2 (columns).\n\n')
    
    x <- as.matrix(x)
    #if(var.position == 1)
    #  x <- as.matrix(t(x))
    
    obj.names <- rownames(x)
    var.names <- colnames(x)
    
    #x <- scale(x,
    #                  center=var.center,
    #                  scale=var.scale)  # scalle variables
    
    #svdx <- svd(x)        # svd of variables
    
    # variables
    
    rownames(u) <- obj.names
    rownames(v) <- var.names
    colnames(v) <- paste('PC',1:length(d),sep='')
    
    
    
    # variables scaled
    s2.var <- diag(d[lambda.ini:lambda.end])
    #match.arg(method)
    switch(method,
           hj = {
             g.var  <- u[ ,lambda.ini:lambda.end] %*%
               s2.var
             h.var  <- s2.var %*%
               t(v[ ,lambda.ini:lambda.end])
             hl.var <- t(h.var)
           },
           sqrt = {
             g.var  <- u[ ,lambda.ini:lambda.end]%*%
               sqrt(s2.var)
             h.var  <- sqrt(s2.var) %*%
               t(v[ ,lambda.ini:lambda.end])
             hl.var <- t(h.var)
           },
           jk = {
             g.var  <- u[ ,lambda.ini:lambda.end] %*%
               s2.var
             h.var  <- t(v[ ,lambda.ini:lambda.end])
             hl.var <- t(h.var)
           },
           gh = {
             g.var  <- sqrt(nrow(x)-1) *
               u[,lambda.ini:lambda.end]
             h.var  <- 1/sqrt(nrow(x)-1) *
               s2.var %*%
               t(v[,lambda.ini:lambda.end])
             hl.var <- t(h.var)
           })
    
    if(is.null(rownames(x)))
      row.names <- 1:nrow(x) 
    else
      row.names <- rownames(x)
    
    if(is.null(colnames(x)))
      col.names <- paste('V',
                         1:ncol(x),
                         sep='') 
    else
      col.names <- colnames(x)
    
    pc.names <- paste('Component',
                      lambda.ini:lambda.end)
    
    rownames(g.var)  <- row.names
    colnames(g.var)  <- pc.names
    rownames(hl.var) <- col.names
    colnames(hl.var) <- pc.names
    
    # variables
    if(var.rb)
      var.rb.res <- var.rbf(hl.var) 
    else
      var.rb.res <- NA
    
    if(var.rb & var.rd)
      var.rd.res <- var.rdf(x, var.rb.res, limit) 
    else
      var.rd.res <- NA
    #call=match.call()
    res <- list(#,
      eigenvalues=d,
      eigenvectors=v,
      number=seq(lambda.ini, lambda.end, 1),
      importance=rbind(general=round(sum(d[lambda.ini:lambda.end]^2) /
                                       sum(d^2), 3),
                       partial=round(sum(d[lambda.ini:lambda.end]^2) /
                                       sum(d[lambda.ini:length(d)]^2), 3)),
      coord=list(objects=g.var,
                 variables=hl.var),
      var.rb=var.rb.res,
      var.rd=var.rd.res)
    
    colnames(res$importance) <- 'explained'
    
    if(ncol(g.var) == 2)
      class(res) <- c('bpca.2d', 'list')
    else if(ncol(g.var) == 3)
      class(res) <- c('bpca.3d', 'list')
    
    invisible(res)
  }

## Calculate scores from objects of the class 'bpca.2d'
scores.bpca.2d <- function(x, var.factor=1){
    if (!inherits(x, 'bpca.2d'))
      stop("Use this function only with 'bpca.2d' class!")
    scores <- rbind(x$coord$objects,
                    x$coord$variables*var.factor,
                    rep(0,
                        ncol(x$coord$objects)))
    max.pc1 <- max(abs(scores[,1]))
    max.pc2 <- max(abs(scores[,2]))
    scores[,1] <- scores[,1]/max.pc1
    scores[,2] <- scores[,2]/max.pc2 
    invisible(scores)
  }

hdr.boxplot.2d <- function(x, y, prob=c(0.01,0.05), h, show.points=FALSE, pc.sig=1, allabels=TRUE,contour=contour, xlab="", ylab="", kde.package=c("ash","ks"),
                           shadecols=gray((9:1)/10), pointcol=1, xlim, ylim,...)
{   
  #limx <-quantile(x,c(0.1,0.9))
  #limy <-quantile(y,c(0.1,0.9))
  limx <-HPDinterval(as.mcmc(x),1-max(prob))
  limy <-HPDinterval(as.mcmc(y),1-max(prob))
  if(limx[1]<0 && limx[2]<0)
  {
    limx <- limx[2]
  }
  else if(limx[1]<0 && limx[2]>0)
  {
    limx <- min(abs(limx))
  }
  else if (limx[1]>0 && limx[2]>0)
  {
    limx <- limx[1]
  }
  
  if(limy[1]<0 && limy[2]<0)
  {
    limy <- limy[2]
  }
  else if(limy[1]<0 && limy[2]>0)
  {
    limx <- min(abs(limy))
  }
  else if (limy[1]>0 && limy[2]>0)
  {
    limy <- limy[1]
  }
  x <- c(0,limx,0,limx,x)
  y <- c(0,0,limy,limy,y)
  
  # Plots bivariate HDRs in form of boxplot.
  kde.package <- match.arg(kde.package)
  if(kde.package=="ash")
  {
    require(ash)
    if(missing(h))
      h <- c(5,5)
    den <- ash2(bin2(cbind(x,y)),h)
  }
  else
  {
    require(ks)
    X <- cbind(x,y)
    if(missing(h))
      h <- Hpi.diag(X,binned=TRUE)
    else
      h <- diag(h)
    den <- kde(x=X,H=h)
    den <- list(x=den$eval.points[[1]],y=den$eval.points[[2]],z=den$estimate)
  }
  hdr <- plothdr2d(x, y, den, alpha=prob, show.points=show.points, pc.sig=pc.sig, allabels=allabels,contour=contour, xlab=xlab, ylab=ylab, shadecols=shadecols, pointcol=pointcol, xlim=xlim, ylim=ylim, ...)
  invisible(hdr)
}
plothdr2d <- function(x, y, den, alpha=c(0.01,0.05), shaded=TRUE, show.points=TRUE,
                      outside.points=TRUE, pc.sig=1, allabels=TRUE, contour=TRUE, pch=19, shadecols, pointcol, xlim, ylim, xlab, ylab, ...)
{
  hdr <- hdr.info.2d(x, y, den, alpha=alpha)
  switch(pc.sig+1,
         plot.sig <- (hdr$falpha[which.max(alpha)] > hdr$fxy[1]),
         plot.sig <- (hdr$fxy[2] >= hdr$fxy[1]),
         plot.sig <- (hdr$fxy[3] >= hdr$fxy[1]),
         plot.sig <- (hdr$fxy[2] >= hdr$fxy[1])||(hdr$fxy[3] >= hdr$fxy[1]),
         plot.sig <- TRUE)
  if(shaded && plot.sig && contour)
    hdrcde.filled.contour(den$x,den$y,den$z,levels=c(hdr$falpha,1e10),col=shadecols, xlim=xlim, ylim=ylim,...)
  else if(plot.sig && contour)
    contour(den,levels=hdr$falpha, xlim=xlim, ylim=ylim, drawlabels = FALSE,...)
  title(main = "", xlab = xlab, ylab = ylab)
  Axis(x, side = 1)
  Axis(y, side = 2)
  box()
  if(show.points)
    points(x,y,pch=pch,col=pointcol)
  else if(outside.points)
  {
    index <- (hdr$fxy < 0.99999*min(hdr$falpha))
    points(x[index], y[index], pch=pch,col=pointcol)
  }
  #points(hdr$mode[1],hdr$mode[2],pch="o",col=pointcol)
  
  if(!plot.sig && !allabels) hdr$mode <- c(0,0)
  invisible(hdr)
}

hdrcde.filled.contour <- function (x,y,z, xlim = range(x, finite = FALSE), 
                                   ylim = range(y, finite = FALSE), zlim = range(z, finite = TRUE), 
                                   levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                                   col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                                   asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                                   axes = TRUE, frame.plot = axes, ...) 
{
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  #    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #    on.exit(par(par.orig))
  #    mar <- mar.orig
  #    mar[4] <- 1
  #    par(mar = mar)
  #plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                  col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}


hdr.info.2d <- function(x, y, den, alpha){
  # Calculates falpha needed to compute HDR of bivariate density den.
  # Also finds approximate mode.
  # Input: den = bivariate density on grid.
  #      (x,y) = indpendent observations on den
  #      alpha = level of HDR
  # Called by plothdr2d
  
  fxy <- interp.2d(den$x,den$y,den$z,x,y)
  falpha <- quantile(sort(fxy), alpha)
  index <- (fxy==max(fxy))
  mode <- c(x[index],y[index])
  return(list(falpha=falpha,mode=mode,fxy=fxy))
}

interp.2d <- function(x, y, z, x0, y0){
  # Bilinear interpolation of function (x,y,z) onto (x0,y0).
  # Taken from Numerical Recipies (second edition) section 3.6.
  # Called by hdr.info.2d
  # Vectorized version of old.interp.2d. 
  # Contributed by Mrigesh Kshatriya (mkshatriya@zoology.up.ac.za)
  
  nx <- length(x)
  ny <- length(y)
  n0 <- length(x0)
  z0 <- numeric(length = n0)
  xr <- diff(range(x))
  yr <- diff(range(y))
  xmin <- min(x)
  ymin <- min(y)
  j <- ceiling(((nx - 1) * (x0 - xmin))/xr)
  k <- ceiling(((ny - 1) * (y0 - ymin))/yr)
  j[j == 0] <- 1
  k[k == 0] <- 1
  j[j == nx] <- nx - 1
  k[k == ny] <- ny - 1
  v <- (x0 - x[j])/(x[j + 1] - x[j])
  u <- (y0 - y[k])/(y[k + 1] - y[k]) 
  AA <- z[cbind(j, k)]
  BB <- z[cbind(j + 1, k)]
  CC <- z[cbind(j + 1, k + 1)]
  DD <- z[cbind(j, k + 1)]
  z0 <- (1 - v) * (1 - u) * AA + v * (1 - u) * BB + v * u * CC + (1 - v) * u * DD
  return(z0)
}

## Plot mcmc scores with the package 'graphics'
plot.bpca.mcmc.2d <- function(scores, r, c, prob=c(0.01,0.05,0.10),ref.lines=TRUE, shaded=TRUE, pc.sig=0,
           allabels=TRUE, contour=TRUE, join=TRUE, dir.obj=TRUE, dir.var=TRUE, ref.color='navy', ref.lty='dotted',
           var.factor=1,var.color='black', var.lty='solid', var.pch=20,
           var.pos=4, var.cex=1, var.offset=.2,
           obj.color='black', obj.pch=20, obj.pos=4, obj.cex=1, obj.offset=.2,
           obj.names=TRUE, obj.labels=dimnames(scores[[1]])[[1]],
           xlim, ylim, xlab, ylab, xaxs = "i", yaxs = "i", asp = NA,...){
    
    
    #  if (missing(xlim) || missing(ylim)) {
    s <- length(scores)
    #absxy <- sapply(as.data.frame(scores),abs)
    #maxabsx <- max(absxy[seq(1,(2*s)-1,2)])
    #maxabsy <- max(absxy[seq(2,(2*s),2)])    
    #  }
    if (missing(xlim))
      xlim <- c(-1.1,1.1)
    if (missing(ylim))
      ylim <- c(-1.1,1.1)
    if (missing(xlab))
      xlab <- dimnames(scores[[1]])[[2]][1]
    if (missing(ylab))
      ylab <- dimnames(scores[[1]])[[2]][2]
    nsim <- length(scores)
    
    
    ## on a device where alpha-transparency is supported,
    ##  use 'alpha = 0.3' transparency with the default palette :
    if(shaded){
      opal <- col2rgb(gray(sort(cumsum(rep(length(prob)/20,length(prob))),decreasing=TRUE)), alpha=TRUE)/255; 
      opal["alpha",] <- 0.3
      mycols <- do.call(rgb, as.list(as.data.frame(t(opal))))
      opal <- palette(mycols)
    } 
    
    # Objects
    
    if(!shaded){
      if(!join){
        xfig("row.fig",width = 7, height = 7, horizontal = FALSE, pointsize = 12, onefile = TRUE)
        plot.new()
      }
      else{
        xfig("row_col.fig",width = 7, height = 7, horizontal = FALSE, pointsize = 12, onefile = TRUE)
        plot.new()
      }
    } else{
      tiff(paste("row_",jj,".tiff",sep=''))
      plot.new()}   
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    pardef <- par(cex=16/12) 
    rowest <- NULL
    for(j in 1:r){
      rowcomp <- NULL
      for(i in 1:nsim){
        rowcomp <- rbind(rowcomp,scores[[i]][j,])
      }
      rowcomp <- rowcomp #/c(maxabsx,maxabsy)
      if(!shaded){
        hdr <- hdr.boxplot.2d(rowcomp[,1],rowcomp[,2],prob=prob,outside.points=FALSE, pc.sig=pc.sig, shaded=shaded, allabels=allabels,contour=contour, xlab = xlab, ylab = ylab, xlim=xlim, ylim=ylim, add=TRUE)
      }
      else{
        hdr <- hdr.boxplot.2d(rowcomp[,1],rowcomp[,2],prob=prob,outside.points=FALSE, pc.sig=pc.sig, shaded=shaded, allabels=allabels,contour=contour, xlab = xlab, ylab = ylab, xlim=xlim, ylim=ylim, shadecols= palette(mycols))
      }
      rowest <- rbind(rowest,hdr$mode)
    }
    
    ind.sig  <- as.logical(apply(rowest,1,sum)) 
    #points(rowest[ind.sig,],pch=20)
    
    text(x=rowest[ind.sig,1],
         y=rowest[ind.sig,2],
         labels=obj.labels[1:r][ind.sig],
         pos=obj.pos,
         offset=obj.offset,
         col=obj.color,
         cex=obj.cex, ...)
    print(obj.labels[1:r][ind.sig])
    if(dir.obj){
      arrows(x0=0,
             y0=0,
             x1=rowest[ind.sig,1]*var.factor,
             y1=rowest[ind.sig,2]*var.factor,
             col=ref.color,
             lty=var.lty, length = 0.1,...)
    }
    if (ref.lines)
      abline(h=0,
             v=0,
             col=ref.color,
             lty=ref.lty)
    if(!shaded){
      if(!join){
        graphics.off()
        system("fig2dev -L emf row.fig > row.emf")
      }
    }
    par(pardef) 
    dev.off()
    # Variables
    if(!shaded){
      if(!join){
        xfig("col.fig",width = 7, height = 7, horizontal = FALSE, pointsize = 12, onefile = TRUE)
        plot.new()
      } #else plot.new()
    } else{
      tiff(paste("col_",jj,".tiff",sep=''))
      plot.new()}   
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    pardef <- par(cex=16/12) 
    
    nsim <- length(scores)
    colest <- NULL
    for(j in (r+1):(r+c)){
      colcomp <- NULL
      for(i in 1:nsim){
        colcomp <- rbind(colcomp,scores[[i]][j,])
      }
      colcomp <- colcomp #/c(maxabsx,maxabsy)
      if(!shaded){
        hdr <- hdr.boxplot.2d(colcomp[,1],colcomp[,2],prob=prob,outside.points=FALSE, pc.sig=pc.sig, shaded=shaded,  allabels=allabels,contour=contour, xlab = xlab, ylab = ylab, xlim=xlim, ylim=ylim, add=TRUE)
      }                                              
      else{
        hdr <- hdr.boxplot.2d(colcomp[,1],colcomp[,2],prob=prob,outside.points=FALSE, pc.sig=pc.sig, shaded=shaded,  allabels=allabels,contour=contour, xlab = xlab, ylab = ylab, xlim=xlim, ylim=ylim, shadecols= palette(mycols))
      }
      colest <- rbind(colest,hdr$mode)
    }
    
    ind.sig  <- as.logical(apply(colest,1,sum))
    #points(colest[ind.sig,],pch=20)
    
    text(x=colest[ind.sig,1],
         y=colest[ind.sig,2],
         labels=obj.labels[(r+1):(r+c)][ind.sig],
         pos=obj.pos,
         offset=obj.offset,
         col=var.color,
         cex=obj.cex, ...)
    print(obj.labels[(r+1):(r+c)][ind.sig])
    if(dir.var){
      arrows(x0=0,
             y0=0,
             x1=colest[ind.sig,1]*var.factor,
             y1=colest[ind.sig,2]*var.factor,
             col=ref.color,
             lty=var.lty, length = 0.1,...)
    }
    if (ref.lines)
      abline(h=0,
             v=0,
             col=ref.color,
             lty=ref.lty)
    if(!shaded){
      graphics.off()
      if(!join && !shaded) system("fig2dev -L emf col.fig > col.emf")
      else system("fig2dev -L emf row_col.fig > row_col.emf")
    }
    par(pardef)
    dev.off() 
    palette("default")
  }

waas_index<-function(outdir,r){
  alphama<-read.csv(file.path(outdir, 'alpha.csv'), sep = ",")
  
  gamma<- read.csv(file.path(outdir, 'gamma.csv'), sep = ",")

  lambda<- read.csv(file.path(outdir, 'lambda.csv'), sep = ",")

  WAAS_Matrix <- matrix(NA, r, dim(alphama)[1])
  idg <- seq(2,dim(alphama)[2],r)
  
  for ( it in 1:dim(alphama)[1])
  {
    dm <- (dim(alphama)[2]-1)/r
    U <-matrix(NA,r,dm)
    D <- as.matrix(diag(lambda[it,-c(1,dim(lambda)[2]+1)]))
    EP<-lambda[it,-c(1,dim(lambda)[2]+1)]/sum(lambda[it,-c(1,dim(lambda)[2]+1)])
    
    
    for( h in 0:(r-1))
    {
      U[h+1,]<-as.numeric(alphama[1,idg+h])
    }
    
    Sg <- as.matrix(U) %*% as.matrix(sqrt(D))
    
    for( jj in 1:r)
    {
      WAAS <-  t(as.matrix(abs(Sg[jj,])))%*%as.matrix(as.numeric(EP))
      WAAS_Matrix[jj,it] <- WAAS
    } 
  }
  
  B_WAAS<- cbind(1:r,apply(WAAS_Matrix,1,mean))

  B_WAAS_Table<- matrix(NA, r, 4)
  
  for ( g in 1:r)
  {
    hpd <-HPDinterval(as.mcmc(WAAS_Matrix[g,]),prob=0.90)
    B_WAAS_Table[g,1] <- g
    B_WAAS_Table[g,2]<-hpd[1]
    B_WAAS_Table[g,4]<- hpd[2]
    B_WAAS_Table[g,3] <- B_WAAS[g,2]
  }
  
  B_WAAS_Table<-as.data.frame(B_WAAS_Table)
  colnames(B_WAAS_Table) <- c("Genotype",  "CI_0.10", "BWAAS","CI_0.90")
  B_WAAS_Table$Genotype<-tratamientos
  return(B_WAAS_Table)
}

waasy_index<-function(outdir,r,weight){
  alphama<-read.csv(file.path(outdir, 'alpha.csv'), sep = ",")
  
  mu<- read.csv(file.path(outdir, 'mu.csv'), sep = ",")
  mu<-mu[,-1]
  tao<- read.csv(file.path(outdir, 'tao.csv'), sep = ",")
  tao<-tao[,-1]
  lambda<- read.csv(file.path(outdir, 'lambda.csv'), sep = ",")
  
  WAAS_Matrix <- matrix(NA, r, dim(alphama)[1])
  idg <- seq(2,dim(alphama)[2],r)
  
  for ( it in 1:dim(alphama)[1])
  {
    dm <- (dim(alphama)[2]-1)/r
    U <-matrix(NA,r,dm)
    D <- as.matrix(diag(lambda[it,-c(1,dim(lambda)[2]+1)]))
    EP<-lambda[it,-c(1,dim(lambda)[2]+1)]/sum(lambda[it,-c(1,dim(lambda)[2]+1)])
    
    for( h in 0:(r-1))
    {
      U[h+1,]<-as.numeric(alphama[1,idg+h])
    }
    
    Sg <- as.matrix(U) %*% as.matrix(sqrt(D))
    
    for( jj in 1:r)
    {
      WAAS <-  t(as.matrix(abs(Sg[jj,])))%*%as.matrix(as.numeric(EP))
      WAAS_Matrix[jj,it] <- WAAS
    } 
  }
  
  BBWAAS <- t(WAAS_Matrix)
  W <- ((0-100)/(apply(BBWAAS,1,max)-apply(BBWAAS, 1, min)))*(BBWAAS-apply(BBWAAS, 1, max)) + 0
  
  MU <- matrix(rep(mu,r),dim(tao)[1],r)
  TAO <- tao + MU
  G <- ((100-0)/(apply(TAO, 1, max)- apply(TAO, 1, min)))*(TAO-apply(TAO, 1, max)) + 100
  
  theta1 <- weight # weight for trait
  theta2 <- 100-weight # weight for waas
  
  WWAASBY <- ((G*theta1) + (W*theta2))/100
  
  B_WAASY <- matrix(NA,r,4)
  for( i in 1:r)
  {
    B_WAASY[i,1] <- i
    B_WAASY[i,2] <- mean(WWAASBY[,i])
    hpd <-HPDinterval(as.mcmc(WWAASBY[,i]),prob=0.90)
    B_WAASY[i,3]<-hpd[1]
    B_WAASY[i,4]<- hpd[2]
  }
  bwaasy <- as.data.frame(B_WAASY)
  colnames(bwaasy) <- c("Genotype", "BWAASY", "CI_0.10", "CI_0.90")
  
  bwaasy$Genotype <- tratamientos
  
  return(bwaasy)
}


stabdist_index<-function(outdir,r){
  alphama<-read.csv(file.path(outdir, 'alpha.csv'), sep = ",")
  
  lambda<- read.csv(file.path(outdir, 'lambda.csv'), sep = ",")
  
  genotypes <- r
  iter <- dim(alphama)[1]
  idg <- seq(2,dim(alphama)[2],r)
  
  StabDist_M <- matrix(NA,iter,genotypes)
  dm <- (dim(alphama)[2]-1)/r
  CP <- dm
  
  for ( it in 1:dim(alphama)[1])
  {
    dm <- (dim(alphama)[2]-1)/r
    U <-matrix(NA,r,dm)
    D <- as.matrix(diag(lambda[it,-c(1,dim(lambda)[2]+1)]))
    EP<-lambda[it,-c(1,dim(lambda)[2]+1)]/sum(lambda[it,-c(1,dim(lambda)[2]+1)])
    
    for( h in 0:(r-1))
    {
      U[h+1,]<-as.numeric(alphama[1,idg+h])
    }
    
    Sg <- as.matrix(U) %*% as.matrix(sqrt(D))
    Sg1 <-rbind(Sg, rep(0,CP))
    distM <- sqrt(mahalanobis(Sg, rep(0,CP), diag(diag(1/D))))
    StabDist_M[it,]<-as.matrix(distM)
    
  }
  
  Mah_mean <- colMeans(StabDist_M)

  B_StabM<- matrix(NA,r,4)
  for( i in 1:r)
  {
    B_StabM[i,1] <- i
    hpd <-HPDinterval(as.mcmc(StabDist_M[,i]),prob=0.90)
    B_StabM[i,2] <- hpd[1]
    B_StabM[i,3]<- Mah_mean[i]
    B_StabM[i,4]<- hpd[2]
  }
  
  
  colnames(B_StabM) <- c("Genotype", "S_CI_0.10", "B_StabM",  "S_CI_0.90")
  
  return(B_StabM)
}

mean_gen<-function(outdir,r){
  mu<- read.csv(file.path(outdir, 'mu.csv'), sep = ",")
  mu<-mu[,-1]
  tao<- read.csv(file.path(outdir, 'tao.csv'), sep = ",")
  tao<-tao[,-1]
  
  MU <- matrix(rep(mu,r),dim(tao)[1],r)
  TAO <- tao + MU
  
  Geno_mean <- colMeans(TAO)
  
  B_Geno_Mean<- matrix(NA,r,4)
  for( i in 1:r)
  {
    B_Geno_Mean[i,1] <- i
    hpd <-HPDinterval(as.mcmc(TAO[,i]),prob=0.90)
    B_Geno_Mean[i,2] <- hpd[1]
    B_Geno_Mean[i,3]<- Geno_mean[i]
    B_Geno_Mean[i,4]<- hpd[2]
  }
  
  colnames(B_Geno_Mean) <- c("Genotype", "CI_0.10", "B_Geno_Mean",  "CI_0.90")
  
  return(B_Geno_Mean)
}