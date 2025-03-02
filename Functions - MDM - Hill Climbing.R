#Functions for MDM analysis
# Lilia Costa, Michel Miler

# Packages
library(fpc)
library(GGally)
library(bnlearn)
library(Rgraphviz)
library(reshape2)
library(tidyverse)
library(av)
library(gifski)
library(gganimate)
library(dplyr)

############################################################################################
# Creating function for DLM with Filtering for unknown observational and state variances
############################################################################################

# Input:
# Yt = the vector of observed time series with length T
# Ft = the matrix of covariates with dimension: number of thetas (p) X sample size (T)  
# delta = discount factor | Wt=Ctx(1-delta)/delta
# Gt = the matrix of state equation with dimension: p X p X T. The default is identity matrix block.
# m0 = the vector of prior mean at time t=0 with length p. The default is non-informative prior, with zero mean.
# CS0 = the squared matrix of prior variance - C*0 | C*0Vt = C0, with length p. The default is non-informative prior, with prior variance equal to 3 times the observed variance.
# n0 and d0 = the prior hypermarameters of precision phi ~ G(n0/2; d0/2). The default is non-informative prior, with value of 0.001. n0 has to be higher than 0.

# output: 
# mt = the matrix of posterior mean with dimension p X T
# Ct = the squared matrix of posterior variance with dimension p X p X T 
# Rt = the squared matrix of prior variance with dimension p X p X T
# nt and dt = the vector of prior hypermarameters of precision phi with length T
# ft = the vector of one-step forecast mean with length T
# Qt = the vector of one-step forecast variance with length T
# ets = the vector of standardised residuals with length T
# lpl = Log Predictive Likelihood with length T

dlm_filt <- function(Yt, Ft, delta, Gt = array(diag(nrow(Ft)), dim=c(nrow(Ft),nrow(Ft),length(Yt))), m0 = rep(0,nrow(Ft)), CS0 = 3*diag(nrow(Ft)), n0 = 0.001, d0 = 0.001) {
  
  # defining objects
  p = nrow(Ft) # the number of thetas
  Nt = length(Yt)+1 # the sample size + t=0
  if (n0 == 0){
    n0 = 0.001
    warning("n0 is 0.001")
  }
  Y = rep(0, Nt)
  Y[2:Nt] = Yt
  F = array(0, dim=c(p,Nt))
  F[,2:Nt] = Ft
  G = array(0, dim=c(p,p,Nt))
  G[,,2:Nt] = Gt
  mt = array(0, dim=c(p,Nt))
  mt[,1] = m0
  Ct = array(0, dim=c(p,p,Nt)) 
  Ct[,,1] = CS0*d0/n0
  Rt = array(0, dim=c(p,p,Nt))
  nt = rep(0, Nt)
  nt[1] = n0
  dt = rep(0, Nt)
  dt[1] = d0
  ft = rep(0, Nt)
  Qt = rep(0, Nt)
  ets = rep(0, Nt)
  lpl = rep(0, Nt)
  
  for (i in 2:Nt){
    
    # Posterior at {t-1}: (theta_{t-1}|y_{t-1}) ~ t_{n_{t-1}}[m_{t-1}, C_{t-1} = C*_{t-1}xd_{t-1}/n_{t-1}]
    # Prior at {t}: (theta_{t}|y_{t-1}) ~ t_{n_{t-1}}[a_{t}, R_{t}]
    at = G[,,i] %*% mt[,(i-1)]
    RSt = G[,,i] %*% (Ct[,,(i-1)]*nt[(i-1)]/dt[(i-1)]) %*% t(G[,,i]) / delta
    Rt[,,i] = RSt * dt[(i-1)] / nt[(i-1)]
    
    # One-step forecast: (Y_{t}|y_{t-1}) ~ t_{n_{t-1}}[f_{t}, Q_{t}]
    ft[i] = t(F[,i]) %*% at
    QSt = t(F[,i]) %*% RSt %*% F[,i] + 1
    Qt[i] = QSt * dt[(i-1)] / nt[(i-1)]
    et = Y[i] - ft[i]
    ets[i] = et / sqrt(Qt[i])
    
    # Posterior at t: (theta_{t}|y_{t}) ~ t_{n_{t}}[m_{t}, C_{t}]
    At = Rt[,,i] %*% F[,i] / Qt[i]
    mt[,i] = at + At * et
    nt[i] = nt[(i-1)] + 1
    dt[i] = dt[(i-1)] + (et^2) / QSt
    CSt = RSt - (At %*% t(At)) * QSt[1]
    Ct[,,i] = CSt * dt[i] / nt[i]
    
    # Log Predictive Likelihood
    lpl[i] <- lgamma((nt[i]+1)/2)-lgamma(nt[i]/2)-0.5*log(pi*nt[i]*Qt[i])-((nt[i]+1)/2)*log(1+(1/nt[i])*et^2/Qt[i])         
    
  }
  
  result <- list(mt=mt[,2:Nt], Ct=Ct[,,2:Nt], Rt=Rt[,,2:Nt], nt=nt[2:Nt], dt=dt[2:Nt], ft=ft[2:Nt], Qt=Qt[2:Nt], ets=ets[2:Nt], lpl=lpl[2:Nt])
  
  return(result)
}


############################################################################################
# Creating function for DLM with Smoothing for unknown observational and state variances
############################################################################################

# Input: all objects are resulted from "dlm_filt", except Gt
# mt = the matrix of posterior mean with dimension p X T
# Ct = the squared matrix of posterior variance with dimension p X p X T 
# Rt = the squared matrix of prior variance with dimension p X p X T
# nt and dt = the vector of prior hypermarameters of precision phi with length T
# Gt = the matrix of state equation with dimension: p X p X T. The default is identity matrix block.

# output: 
# smt = the matrix of smoothing posterior mean with dimension p X T
# sCt = the squared matrix of smoothing posterior variance with dimension p X p X T 

dlm_smoo <- function(mt, Ct, Rt, nt, dt, Gt = 0) {
  
  # defining objects
  if (is.vector(mt)){
    mt = array(mt, dim=c(1,length(mt)))
    Ct = array(Ct, dim=c(1,1,length(mt)))
    Rt = array(Rt, dim=c(1,1,length(Rt)))     
  }
  if (Gt == 0){Gt = array(diag(nrow(mt)), dim=c(nrow(mt),nrow(mt),ncol(mt)))}
  p = nrow(mt) # the number of thetas
  Nt = ncol(mt) # the sample size
  smt = array(0, dim=c(p,Nt))
  sCt = array(0, dim=c(p,p,Nt)) 
  
  # in the last time point
  smt[,Nt] = mt[,Nt]
  sCt[,,Nt] = Ct[,,Nt] 
  
  # for other time points
  for (i in (Nt-1):1){
    RSt = Rt[,,(i+1)]*nt[i]/dt[i]
    CSt = Ct[,,i]*nt[i]/dt[i]
    inv.sR = solvecov(RSt, cmax = 1e+10)$inv
    B = CSt %*% t(Gt[,,(i+1)]) %*% inv.sR
    smt[,i] = mt[, i] + B %*% (smt[,(i+1)] - Gt[,,(i+1)] %*% mt[,i])
    sCS = CSt + B %*% (sCt[,,(i+1)]*nt[Nt]/dt[Nt] - RSt) %*% t(B)     
    sCt[,,i] = sCS * dt[Nt] / nt[Nt]
  }
  
  result <- list(smt=smt, sCt=sCt)
  return(result)
}

###################################################################
### Choosing the delta
###################################################################

# Input:
#  dts = the matrix with dataset; Number of timepoints X Number of nodes
#  m_ad = # Square Matrix Adjacent with dimension = Number of nodes # 1 if edge exists; 0 otherwise
#  nbf => the Log Predictive Likelihood will be calculate from this time point. It has to be a positive integer number. The default is 15.
#  delta = the vector with the sequence of all discount factors. The default is seq(from=0.5, to=1.0, by=0.01).

# Output: 
# lpldet = LPL for each value of delta and for each node; length(delta) X Number of nodes;
# DF_hat = vector with delta that maximizes the LPL for each node with dimension = Number of nodes.

CDELT_cfd <- function(dts,m_ad,nbf=15,delta=seq(from=0.5, to=1.0, by=0.01)) {
  nd = length(delta)
  Nn = ncol(dts)
  Nt = nrow(dts)
  lpldet = array(0,dim=c(nd,Nn))
  for (k in 1:nd){
    for (i in 1:Nn){
      # Initials:
      p = sum(m_ad[,i])
      if (m_ad[i,i] == 0) {p = p + 1}
      Ft = array(1, dim=c(Nt,p))
      aux = c(1:Nn)[m_ad[,i]>0]
      aux2 = aux[aux!=i]
      #SOLUÇÃO: PARA Error in F[, 2:Nt] = Ft :number of items to replace is not a multiple of replacement length
      if (length(aux2)>0){
        if (typeof(dts[,aux2])=='list'){
          Ft_aux = matrix(unlist(dts[,aux2]),nrow=Nt,ncol=length(aux2))
          Ft[,2:(length(aux2)+1)] = Ft_aux
        }else{
          Ft[,2:(length(aux2)+1)] = dts[,aux2]
        }
      }
      Yt = dts[,i]
      # DLM
      
      a=dlm_filt(Yt, t(Ft), delta=delta[k])
      lpldet[k,i]=sum(a$lpl[nbf:Nt])
    }
  }
  #DF_hat=delta[max.col(t(lpldet))] # with some deltas provide NA in lpl, it means that this delta is not good for this particular dataset so we have to use:
  DF_hat = rep(0,Nn)
  for (i in 1:Nn){
    DF_hat[i] = na.omit(delta[lpldet[,i]==max(lpldet[,i],na.rm=TRUE)])[1]
  }
  result <- list(lpldet = lpldet, DF_hat = DF_hat)
  return(result)
}

############################################################################################
# Creating function for MDM with Filtering for unknown observational and state variances
############################################################################################

# Input:
# dts = the matrix with dataset; Number of timepoints X Number of nodes
# m_ad = # Square Matrix Adjacent with dimension = Number of nodes # 1 if edge exists; 0 otherwise
# DF_hat = vector with delta that maximizes the LPL for each node with dimension = Number of nodes.

# output: 
# mt = list with the matrix of posterior mean with dimension p X T
# Ct = list with the squared matrix of posterior variance with dimension p X p X T 
# Rt = list with the squared matrix of prior variance with dimension p X p X T
# nt and dt = list with the vector of prior hypermarameters of precision phi with length T
# ft = list with the vector of one-step forecast mean with length T
# Qt = list with the vector of one-step forecast variance with length T
# ets = list with the vector of standardised residuals with length T
# lpl = list with Log Predictive Likelihood with length T

mdm_filt <- function(dts, m_ad, DF_hat) {
  Nn = ncol(dts)
  Nt = nrow(dts)
  mt = vector(Nn, mode = "list")
  names(mt) <- colnames(m_ad)
  conections <- which(m_ad == 1, arr.ind = T)
  Ct = vector(Nn, mode = "list")
  Rt = vector(Nn, mode = "list")
  nt = vector(Nn, mode = "list")
  dt = vector(Nn, mode = "list")
  ft = vector(Nn, mode = "list")
  Qt = vector(Nn, mode = "list")
  ets = vector(Nn, mode = "list")
  lpl = vector(Nn, mode = "list")
  for (i in 1:Nn){
    # Initials:
    p = sum(m_ad[,i])
    if (m_ad[i,i] == 0) {p = p + 1}   	         
    Ft = array(1, dim=c(Nt,p))
    aux = c(1:Nn)[m_ad[,i]>0]
    aux2 = aux[aux!=i] 
    #SOLUÇÃO: PARA Error in F[, 2:Nt] = Ft :number of items to replace is not a multiple of replacement length
    if (length(aux2)>0){
      if (typeof(dts[,aux2])=='list'){
        Ft_aux = matrix(unlist(dts[,aux2]),nrow=Nt,ncol=length(aux2))
        #colnames(Ft_aux) = 
        Ft[,2:(length(aux2)+1)] = Ft_aux
      }else{
        Ft[,2:(length(aux2)+1)] = dts[,aux2]
      }
    }
    Yt = dts[,i]
    # DLM
    a=dlm_filt(Yt, t(Ft), delta=DF_hat[i])
    mt[[i]] = a$mt
    if(sum(m_ad[,i]) != 0){
      rownames(mt[[i]]) <- c(
        paste0("beta0_", colnames(m_ad)[i]),
        paste0(
          names(which(conections[,2] == i)),
          "->", colnames(m_ad)[i]
        )
      )
    }
    Ct[[i]] = a$Ct
    Rt[[i]] = a$Rt
    nt[[i]] = a$nt
    dt[[i]] = a$dt
    ft[[i]] = a$ft
    Qt[[i]] = a$Qt
    ets[[i]] = a$ets
    lpl[[i]] = a$lpl         
  }
  result <- list(mt=mt, Ct=Ct, Rt=Rt, nt=nt, dt=dt, ft=ft, Qt=Qt, ets=ets, lpl=lpl)
  return(result)
}

############################################################################################
# Creating function for MDM with Smoothing for unknown observational and state variances
############################################################################################

# Input: all objects are resulted from "mdm_filt"
# mt = list with the matrix of posterior mean with dimension p X T
# Ct = list with the squared matrix of posterior variance with dimension p X p X T 
# Rt = list with the squared matrix of prior variance with dimension p X p X T
# nt and dt = list with the vector of prior hypermarameters of precision phi with length T

# output: 
# smt = list with the matrix of smoothing posterior mean with dimension p X T
# sCt = list with the squared matrix of smoothing posterior variance with dimension p X p X T 

mdm_smoo <- function(mt, Ct, Rt, nt, dt) {
  Nn = length(mt) # the number of nodes
  smt = vector(Nn, mode = "list")
  sCt = vector(Nn, mode = "list")
  for (i in 1:Nn){
    a=dlm_smoo(mt=mt[[i]], Ct=Ct[[i]], Rt=Rt[[i]], nt=nt[[i]], dt=dt[[i]])
    rownames(a$smt) <- rownames(mt[[i]])
    smt[[i]] = a$smt
    sCt[[i]] = a$sCt
  }
  result <- list(smt=smt, sCt=sCt)
  return(result)
}




mdm_score <- function(data_input, nbf=15,method='Brent', delta=seq(from=0.5, to=1.0, by=0.01),
                      GOLB_print = FALSE, subjects_length = 1,q_pais=NULL) {
  
  Nt = dim(data_input)[1] #number of time points
  Nn = dim(data_input)[2] #number of nodes
  Ns = subjects_length #number of subjects
  if (is.null(q_pais)){
    q_pais = Nn-1
    i_pais = Nn
  }else{
    i_pais = q_pais
  }
  dts = array(0, dim=c(Ns, Nt, Nn))
  
  if(subjects_length == 1){
    for(i in 1:Nn){
      dts[,,i] <- data_input[,i]
    }
  } else {
    for(j in 1:Nn){
      for(k in 1:Ns){
        for(p in 1:Nt){
          dts[k, p, j] <- t(data_input[ p, j,k])
        }
      }
    }
  }
  
  # Add names to nodes
  if(subjects_length == 1){
    dimnames(dts) <- list(1:Ns, 1:Nt, colnames(data_input))
  } else {
    dimnames(dts) <- list(1:Ns, 1:Nt, colnames(data_input[,,1]))
  }
  
  Ns = dim(dts)[1] # the number of subjects
  Nt = dim(dts)[2] # the number of timepoints
  Nn = dim(dts)[3] # the number of nodes
  Nd=0
  for (i in seq(0,q_pais)){
    Nd = Nd + ncol(utils::combn(Nn-1,i))
  }
  #Nd =  Nn*2^(q_pais-1)# the number of possible parent-child
  Nd = Nd*Nn
  delt_hat = array(0,dim=c(Ns,Nd)) # the discount factor chosen for each possibility
  lpl = array(0,dim=c(Ns,Nd)) # scores
  par_chil = array(-9,dim=c(Ns,Nd,(Nn+2))) #col1=subject; col2=model; col3={child,number of parents,parents}
  # generating all possible combinations
  cc = vector(i_pais, mode = "list")
  for (i in 1:(i_pais)){
    cc[[i]] = utils::combn(c(1:Nn),i)
  }
  # finding the scores ----
  for (s in 1:Ns){
    #cat("loop subject ",s,"/",Ns, "\n")
    Np=0
    for (i in seq(from=1,to=i_pais,by=2)){
      if(subjects_length > 1){
        cat("looping subject ",s,"/",Ns," and node ",i,"/",Nn,"\n")
      } else {
        cat("looping node ",i,"/",i_pais,"\n")
      }
      for (j in 1:ncol(cc[[i]])){
        # Matrix Adjacent: 1 for edge; 0 otherwise
        m_ad = diag(Nn)
        p=cc[[i]][,j]
        m_ad[p,]=1
        # choosing the Discount Factor and finding the scores
        #delt_dag = CDELT(dts[s,,],m_ad, nbf=nbf, delta=delta)
        delt_dag = CDELT(dts[s,,],m_ad,nbf=nbf, method=method)
        delt_hat[s,(Np+1):(Np+Nn)]=delt_dag$DF_hat
        #aux=t(delt_dag$lpldet)
        #lpl[s,(Np+1):(Np+Nn)]=diag(aux[,max.col(aux)])
        #print(delt_dag$lpldet)
        #lpl[s,(Np+1):(Np+Nn)]=apply(delt_dag$lpldet,2,max,na.rm=TRUE) #deleting the NA from lpl
        lpl[s,(Np+1):(Np+Nn)] = delt_dag$lpldet 
        # child
        par_chil[s,(Np+1):(Np+Nn),1] = seq(1:Nn)-1
        # number of parents
        par_chil[s,(Np+1):(Np+Nn),2] = i
        par_chil[s,c((Np+1):(Np+Nn))[p],2] = i-1
        # parents
        par_chil[s,(Np+1):(Np+Nn),3:(2+i)] = t(array((p-1), dim=c(i,Nn)))
        if (i==1){par_chil[s,c((Np+1):(Np+Nn))[p],3:(2+i)] = -9}
        else {diag(par_chil[s,c((Np+1):(Np+Nn))[p],3:(2+i)]) = -9}
        
        Np=Np+Nn
      }
    }
  }
  #print(par_chil)
  DF_hat= array(-9,dim=c(Ns,Nd,(Nn+2)))
  DF_hat[,,1]=delt_hat
  DF_hat[,,2]=par_chil[,,1]
  DF_hat[,,3:(Nn+2)]=par_chil[,,3:(Nn+2)]
  dimnames(DF_hat)=list(c(1:Ns),c(1:Nd),c("DF","node",1:Nn))
  
  # creating the file with structure of James'programm
  
  all.score = array(" ",dim=c((Nd+Nn+1),(Nn+1)))
  all.score[1,1] = Nn
  #par_chil is the same for all subjects
  a = cbind(colSums(lpl),par_chil[1,,])
  b=a[order(a[,2]),]
  more=2
  j=0
  for (i in 1:Nd){
    if (j==b[i,2]){
      all.score[more,1:2] = c(j,Nd/Nn)
      j=j+1
      more=more+1
    }
    d = b[i,c(1,3:(Nn+3))]
    all.score[more,1:length(d[d!=-9])] = d[d!=-9]
    more = more + 1
  }
  #all.score = list()
  #for (s in 1:Ns){
  # all.score[[s]] = array(" ",dim=c((Nd+Nn+1),(Nn+1)))
  # all.score[[s]][1,1] = Nn
  # a = cbind(lpl[s,],par_chil[s,,])
  # b=a[order(a[,2]),]
  # more=2
  # j=0
  # for (i in 1:Nd){
  #   if (j==b[i,2]){
  #     all.score[[s]][more,1:2] = c(j,Nd/Nn)
  #     j=j+1
  #     more=more+1
  #   }
  #   d = b[i,c(1,3:(Nn+3))]
  #   all.score[[s]][more,1:length(d[d!=-9])] = d[d!=-9]
  #   more = more + 1
  # }
  #}
  
  #JUST GETTING THE FIRST SUBJECT
  #all.score_dic <- all.score[[1]]
  all.score_dic<- all.score
  
  
  #if(!is.null(dimnames(dts)[[3]])){
  
  for(i in 1:Nn){
    all.score_dic[all.score_dic == i-1] <- gsub(x = dimnames(dts)[[3]][i],
                                                " ", replacement =  "")
  }
  #}
  if(GOLB_print){
    write.table(x = all.score_dic,
                file = paste0("mdm_score_", format(Sys.time(), "%d_%b_%Y")),
                quote = FALSE, row.names = FALSE,
                col.names = FALSE)
    
  }
  
  result <- list(all.score=all.score_dic, DF_hat=DF_hat)
  return(result)
}



logpl<-function(dts,m_ad,delta,i){
  Nt = nrow(dts)
  Nn = ncol(dts)
  nbf=15
  p = sum(m_ad[,i])
  if (delta>1|delta<0){
    return(Inf)
  }
  if (m_ad[i,i] == 0) {p = p + 1}
  Ft = array(1, dim=c(Nt,p))
  aux = c(1:Nn)[m_ad[,i]>0]
  aux2 = aux[aux!=i]
  if (length(aux2)>0){
    if (typeof(dts[,aux2])=='list'){
      Ft_aux = matrix(unlist(dts[,aux2]),nrow=Nt,ncol=length(aux2))
      Ft[,2:(length(aux2)+1)] = Ft_aux
    }else{
      Ft[,2:(length(aux2)+1)] = dts[,aux2]
    }
  }
  Yt = dts[,i]
  # DLM
  a=dlm_filt(Yt, t(Ft), delta=delta)
  lpldet=sum(a$lpl[nbf:Nt])
  return(-lpldet)
}





CDELT <- function(dts,m_ad,nbf=15,method='Brent',call=FALSE) {
  #nd = length(delta)
  Nn = ncol(dts)
  Nt = nrow(dts)
  lpldet = rep(0,Nn)
  
  if (call){
    data_input = dts
    dts = array(0, dim=c(Nt, Nn))
    for(i in 1:Nn){
      dts[,i] <- data_input[,i]
    }
  }
  res = list()
  
  for (i in 1:Nn){
    # Initials:
    k<-optim(logpl,par=0.5,dts=dts,m_ad=m_ad,i=i,method=method,lower=0,upper=1)
    res[[i]]<-k
    #print(res)
    lpldet[i]=-res[[i]]$value
  }
  
  #DF_hat=delta[max.col(t(lpldet))] # with some deltas provide NA in lpl, it means that this delta is not good for this particular dataset so we have to use:
  DF_hat = rep(0,Nn)
  for (i in 1:Nn){
    DF_hat[i] = res[[i]]$par
  }
  result <- list(lpldet = lpldet, DF_hat = DF_hat)
  return(result)
}

#################################################################################
# Score of MDM for bnlearn package
#################################################################################

mdm_score_bn <- function(node, parents, data, args) {
  n_n = dim(data)[2]
  m_ad = array(0, dim=c(n_n,n_n))
  dimnames(m_ad) = list(colnames(data),colnames(data))
  #m_ad = m_ad[dimnames(m_ad)[[1]][order(nchar(dimnames(m_ad)[[1]]), dimnames(m_ad)[[1]])],dimnames(m_ad)[[2]][order(nchar(dimnames(m_ad)[[2]]), dimnames(m_ad)[[2]])]]
  m_ad[parents,node] = 1
  nbf= args$nbf
  method = args$method
  call = args$call
  #nd = length(delta)
  Nn = ncol(data)
  Nt = nrow(data)
  i = grep(node, colnames(m_ad))[1]
  #lpldet = rep(0,Nn)
  if (call){
    data_input = data
    data = array(0, dim=c(Nt, Nn))
    for(i in 1:Nn){
      data[,i] <- data_input[,i]
    }
  }
  
  res = list()
  k<-optim(logpl,par=0.5,dts=data,m_ad=m_ad,i=i,method=method,lower=0,upper=1)
  res<-k
  lpldet=-res$value
  
  
  #DF_hat=delta[max.col(t(lpldet))] # with some deltas provide NA in lpl, it means that this delta is not good for this particular dataset so we have to use:
  DF_hat = rep(0,Nn)
  for (i in 1:Nn){
    DF_hat[i] = res$par
  }
  #result <- list(lpldet = lpldet, DF_hat = DF_hat)
  #return(result)
  return(lpldet)
}

