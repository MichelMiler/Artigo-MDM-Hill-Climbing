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



##################################################################################
# Estimation of MDM - Hill Climbing for the real data teacher-student
##################################################################################


aluno1<-read.table('Dados/aluno1_001_oxyhb_aula.txt')
professor1<-read.table('Dados/prof1_001_oxyhb_aula.txt')

aluno8<-read.table('Dados/aluno8_001_oxyhb_aula.txt')
professor8<-read.table('Dados/prof8_001_oxyhb_aula.txt')

aluno9<-read.table('Dados/aluno9_001_oxyhb_aula.txt')
professor9<-read.table('Dados/prof9_001_oxyhb_aula.txt')

aluno11<-read.table('Dados/aluno11_001_oxyhb_aula.txt')
professor11<-read.table('Dados/prof11_001_oxyhb_aula.txt')

aluno12<-read.table('Dados/aluno12_001_oxyhb_aula.txt')
professor12<-read.table('Dados/prof12_001_oxyhb_aula.txt')

aluno14<-read.table('Dados/aluno14_001_oxyhb_aula.txt')
professor14<-read.table('Dados/prof14_001_oxyhb_aula.txt')

aluno15<-read.table('Dados/aluno15_001_oxyhb_aula.txt')
professor15<-read.table('Dados/prof15_001_oxyhb_aula.txt')

aluno17<-read.table('Dados/aluno17_001_oxyhb_aula.txt')
professor17<-read.table('Dados/prof17_001_oxyhb_aula.txt')

aluno20<-read.table('Dados/aluno20_001_oxyhb_aula.txt')
professor20<-read.table('Dados/prof20_001_oxyhb_aula.txt')

aluno21<-read.table('Dados/aluno21_001_oxyhb_aula.txt')
professor21<-read.table('Dados/prof21_001_oxyhb_aula.txt')

aluno22<-read.table('Dados/aluno22_001_oxyhb_aula.txt')
professor22<-read.table('Dados/prof22_001_oxyhb_aula.txt')

aluno24<-read.table('Dados/aluno24_001_oxyhb_aula.txt')
professor24<-read.table('Dados/prof24_001_oxyhb_aula.txt')




varnames = c("SV1",  "SV2",  "SV4","SV5",  "SV7",  "SV8",  "SV10", "SV11", "SV13", "SV14", "SV15", "SV17", "SV18",
             "SV20", "SV21", "SV23","TV1",  "TV2",  "TV4",  "TV5",  "TV7",  "TV8",  "TV10", "TV11", "TV13", "TV14", "TV15", "TV17", "TV18",
             "TV20", "TV21", "TV23")

for (i in c(1,8,9,11,12,14,15,17,20,21,22,24)) {
  x <- paste0("duplas",i)
  eval(call("<-", as.name(x), cbind(get(paste0('aluno',i)),get(paste0('professor',i)))))
}

colnames(duplas1) = varnames
colnames(duplas8) = varnames
colnames(duplas9) = varnames
colnames(duplas11) = varnames
colnames(duplas12) = varnames
colnames(duplas14) = varnames
colnames(duplas15) = varnames
colnames(duplas17) = varnames
colnames(duplas20) = varnames
colnames(duplas21) = varnames
colnames(duplas22) = varnames
colnames(duplas24) = varnames


duplas_agg = list(tail(duplas1,3000),
                  tail(duplas8,3000),
                  tail(duplas9,3000),
                  tail(duplas11,3000),
                  tail(duplas12,3000),
                  tail(duplas14,3000),
                  tail(duplas15,3000),
                  tail(duplas17,3000),
                  tail(duplas20,3000),
                  tail(duplas21,3000),
                  tail(duplas22,3000),
                  tail(duplas24,3000))


duplas_agg_df = rbind(duplas_agg[[1]],duplas_agg[[2]],duplas_agg[[3]],duplas_agg[[4]],duplas_agg[[5]],duplas_agg[[6]],
                      duplas_agg[[7]],duplas_agg[[8]],duplas_agg[[9]],duplas_agg[[10]],duplas_agg[[11]],duplas_agg[[12]])  


allglobal <- function() list2env(mget(ls(name = parent.frame()), envir = parent.frame()), envir = .GlobalEnv)


approach <-function(shd_bool=FALSE){
  
  if(shd_bool = TRUE){
    
  }
}

#------------------------------------------------------------------------------
# MDM-Hill Climbing using individual structure
#------------------------------------------------------------------------------

adj_matrix_hills_is = list()
time_hills = c()
for (i in seq(1,12)){
  print(i)
  n_n=32
  start.time = Sys.time()
  hill = hc(x = data.frame(duplas_agg[[i]]),score='custom',fun=mdm_score_bn,args=list(nbf=15,method='Brent',call=FALSE))
  adj_matrix_hill = matrix(0,nrow=n_n,ncol=n_n,dimnames=list(colnames(duplas_agg[[i]]),colnames(duplas_agg[[i]])))
  if (length(hill$arcs)!=0){
    for (z in seq(length(hill$arcs)/2)){
      adj_matrix_hill[hill$arcs[z],hill$arcs[z+length(hill$arcs)/2]]=1
    }
  }
  end.time = Sys.time()
  time_hills[i] = difftime(end.time, start.time, units = "secs")
  adj_matrix_hills_is[[i]] = adj_matrix_hill
}


for(i in seq(1,12)){
  dimnames(adj_matrix_hills_is[[i]]) = list(varnames,varnames)
  
}


total_matrix_hills_is =  adj_matrix_hills_is[[1]]
for(i in seq(2,12)){
  total_matrix_hills_is =  total_matrix_hills_is + adj_matrix_hills_is[[i]]
}

total_hills_is_pf = setNames(melt(total_matrix_hills_is), c('Father', 'Child', 'Frequency'))

con_matrix_hills_is<-0

k = diag(rep(0,32))
lower_ind = which(lower.tri(adj_matrix_hills_is[[1]]),arr.ind=T)
upper_ind = t(combn(ncol(adj_matrix_hills_is[[1]]),2))

for (i in seq(1,12)) {
  con_matrix_hills_is <- con_matrix_hills_is + adj_matrix_hills_is[[i]][lower_ind] + adj_matrix_hills_is[[i]][upper_ind]
}





#-------------------------------------------------------------------------------
# MDM-Hill Climbing using commum structure (CS)
#-------------------------------------------------------------------------------
n_n=32
start.time = Sys.time()
hill = hc(x = data.frame(duplas_agg_df),score='custom',fun=mdm_score_bn,args=list(nbf=15,method='Brent',call=FALSE))
adj_matrix_hill_cs = matrix(0,nrow=n_n,ncol=n_n,dimnames=list(colnames(duplas_agg_df),colnames(duplas_agg_df)))
if (length(hill$arcs)!=0){
  for (z in seq(length(hill$arcs)/2)){
    adj_matrix_hill_cs[hill$arcs[z],hill$arcs[z+length(hill$arcs)/2]]=1
  }
}

end.time = Sys.time()
time_hill = difftime(end.time, start.time, units = "secs")
#save.image("C:/Users/Michel/Documents/Estatística/TCC/temp.RData")

dimnames(adj_matrix_hill_cs) = list(varnames,varnames)

total_hill_cs_pf = setNames(melt(adj_matrix_hill_cs), c('Father', 'Child', 'Frequency'))
#-------------------------------------------------------------------------------
#Analysis inside the brain region individually.
#-------------------------------------------------------------------------------
neuro_matrix_real = matrix(0,nrow=6,ncol=6)                                                                                            

colnames(neuro_matrix_real)=rownames(neuro_matrix_real)=c('Child_PFE','Child_PFD','Child_TP',
                                                          'Professor_PFE','Professor_PFD','Professor_TP')

sum(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                     total_hill_cs_pf$Child %in% c("SV1",  "SV2",  "SV4",  "SV5")), 'Frequency'])                                                                                         

sum(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("TV1",  "TV2",  "TV4",  "TV5") & 
                     total_hill_cs_pf$Child %in% c("TV1",  "TV2",  "TV4",  "TV5")),'Frequency'])                                                                                                                                                                        

sum(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV7",  "SV8",  "SV10",  "SV11") & 
                     total_hill_cs_pf$Child %in% c("SV7",  "SV8",  "SV10",  "SV11")),'Frequency'])                                                                                         

sum(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("TV7",  "TV8",  "TV10",  "TV11") & 
                     total_hill_cs_pf$Child %in% c("TV7",  "TV8",  "TV10",  "TV11")),'Frequency'] )    

sum(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23") & 
                     total_hill_cs_pf$Child %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23")),'Frequency'])  

sum(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23") & 
                     total_hill_cs_pf$Child %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency'])  


#----------------------------------------------------------------------------
# Analysis betwern regions of the same brain 
#----------------------------------------------------------------------------

neuro_matrix_real[1,2] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                                               total_hill_cs_pf$Child %in% c("SV7",  "SV8",  "SV10",  "SV11")),'Frequency']) 

neuro_matrix_real[2,1] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                                               total_hill_cs_pf$Father %in% c("SV7",  "SV8",  "SV10",  "SV11")),'Frequency']) 


neuro_matrix_real[1,3] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                                               total_hill_cs_pf$Child %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23")),'Frequency']) 

neuro_matrix_real[3,1] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                                               total_hill_cs_pf$Father %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23")),'Frequency']) 


neuro_matrix_real[2,3] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV7",  "SV8",  "SV10",  "SV11") & 
                                               total_hill_cs_pf$Child %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23")),'Frequency']) 


neuro_matrix_real[3,2] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV7",  "SV8",  "SV10",  "SV11") & 
                                               total_hill_cs_pf$Father %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23")),'Frequency']) 




neuro_matrix_real[4,5] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("TV1",  "TV2",  "TV4",  "TV5") & 
                                               total_hill_cs_pf$Child %in% c("TV7",  "TV8",  "TV10",  "TV11")),'Frequency']) 

neuro_matrix_real[5,4] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("TV1",  "TV2",  "TV4",  "TV5") & 
                                               total_hill_cs_pf$Father %in% c("TV7",  "TV8",  "TV10",  "TV11")),'Frequency']) 


neuro_matrix_real[4,6] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("TV1",  "TV2",  "TV4",  "TV5") & 
                                               total_hill_cs_pf$Child %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency']) 

neuro_matrix_real[6,4] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("TV1",  "TV2",  "TV4",  "TV5") & 
                                               total_hill_cs_pf$Father %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency']) 


neuro_matrix_real[5,6] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("TV7",  "TV8",  "TV10",  "TV11") & 
                                               total_hill_cs_pf$Child %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency']) 


neuro_matrix_real[6,5] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("TV7",  "TV8",  "TV10",  "TV11") & 
                                               total_hill_cs_pf$Father %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency']) 


#---------------------------------------------------------------------
# Analysis of the connections between teacher and student brains
#---------------------------------------------------------------------



neuro_matrix_real[1,4] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                                               total_hill_cs_pf$Child %in% c("TV1",  "TV2",  "TV4",  "TV5")),'Frequency'])  

neuro_matrix_real[4,1] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                                               total_hill_cs_pf$Father %in% c("TV1",  "TV2",  "TV4",  "TV5")),'Frequency']) 



neuro_matrix_real[1,5] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                                               total_hill_cs_pf$Child %in% c("TV7",  "TV8",  "TV10",  "TV11")),'Frequency'])  

neuro_matrix_real[5,1] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                                               total_hill_cs_pf$Father %in% c("TV7",  "TV8",  "TV10",  "TV11")),'Frequency']) 


neuro_matrix_real[1,6] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                                               total_hill_cs_pf$Child %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency'])  

neuro_matrix_real[6,1] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV1",  "SV2",  "SV4",  "SV5") & 
                                               total_hill_cs_pf$Father %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency']) 






neuro_matrix_real[2,4] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV7",  "SV8",  "SV10",  "SV11") & 
                                               total_hill_cs_pf$Child %in% c("TV1",  "TV2",  "TV4",  "TV5")),'Frequency'])  

neuro_matrix_real[4,2] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV7",  "SV8",  "SV10",  "SV11") & 
                                               total_hill_cs_pf$Father %in% c("TV1",  "TV2",  "TV4",  "TV5")),'Frequency']) 



neuro_matrix_real[2,5] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV7",  "SV8",  "SV10",  "SV11") & 
                                               total_hill_cs_pf$Child %in% c("TV7",  "TV8",  "TV10",  "TV11")),'Frequency'])  

neuro_matrix_real[5,2] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV7",  "SV8",  "SV10",  "SV11") & 
                                               total_hill_cs_pf$Father %in% c("TV7",  "TV8",  "TV10",  "TV11")),'Frequency']) 


neuro_matrix_real[2,6] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV7",  "SV8",  "SV10",  "SV11") & 
                                               total_hill_cs_pf$Child %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency'])  

neuro_matrix_real[6,2] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV7",  "SV8",  "SV10",  "SV11") & 
                                               total_hill_cs_pf$Father %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency']) 





neuro_matrix_real[3,4] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23") & 
                                               total_hill_cs_pf$Child %in% c("TV1",  "TV2",  "TV4",  "TV5")),'Frequency'])  

neuro_matrix_real[4,3] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23") & 
                                               total_hill_cs_pf$Father %in% c("TV1",  "TV2",  "TV4",  "TV5")),'Frequency']) 



neuro_matrix_real[3,5] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23") & 
                                               total_hill_cs_pf$Child %in% c("TV7",  "TV8",  "TV10",  "TV11")),'Frequency'])  

neuro_matrix_real[5,3] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23") & 
                                               total_hill_cs_pf$Father %in% c("TV7",  "TV8",  "TV10",  "TV11")),'Frequency']) 


neuro_matrix_real[3,6] = mean(total_hill_cs_pf[(total_hill_cs_pf$Father %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23") & 
                                               total_hill_cs_pf$Child %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency'])  

neuro_matrix_real[6,3] = mean(total_hill_cs_pf[(total_hill_cs_pf$Child %in% c("SV13",  "SV14",  "SV15",  "SV17","SV18","SV20","SV21","SV23") & 
                                               total_hill_cs_pf$Father %in% c("TV13",  "TV14",  "TV15",  "TV17","TV18","TV20","TV21","TV23")),'Frequency']) 
#-----------------------------------------------------------------------------------------
# Connectivity and strength graph
#-----------------------------------------------------------------------------------------

my.lines<-data.frame(x=c(.5,16.5,16.5,16.5), y=c(16.5,0.5,16.5,16.5), 
                     xend=c(16.5,16.5,32.5,16.5), yend=c(16.5,16.5,16.5,32.5))


my.vlines <- data.frame(
  y = c(0.5,0.5, 0.5, 0.5,0.5),
  x = c(4.5,8.5, 20.5,24.5,32.5),
  yend = c(32.5,32.5, 32.5, 32.5, 32.5),
  xend = c(4.5,8.5, 20.5,24.5,32.5)
)

my.hlines <- data.frame(
  x = c(0.5,0.5, 0.5, 0.5,0.5),
  y = c(4.5,8.5, 20.5,24.5,32.5),
  xend = c(32.5,32.5, 32.5, 32.5, 32.5),
  yend = c(4.5,8.5, 20.5,24.5,32.5)
)



# Your modified ggplot code
total_hill_cs_pf |>
  mutate(
    Frequency = if_else(round(Frequency) == 1, 1, 0)
  ) |>
  mutate(Father = factor(Father, levels = varnames),
         Child = factor(Child, levels = varnames)) |>
  ggplot(aes(y = Father, x = Child, fill = factor(Frequency))) +
  geom_tile(color = "gray") +
  theme_bw() +
  scale_fill_manual(values = c("0" = "white", "1" = "red")) +
  labs(y = "Parent node", x = "Child node", fill = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  geom_segment(data = my.vlines, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1, inherit.aes = FALSE, linetype = 'dashed', color = 'orange') +
  geom_segment(data = my.hlines, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1, inherit.aes = FALSE, linetype = 'dashed', color = 'orange')+
  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), linewidth=1, inherit.aes=F,linetype='solid',color='black')




#------------------------------------------------------------------------------------
## SHD calculation

lower_ind = which(lower.tri(adj_matrix_hill_cs),arr.ind=T)
k = diag(rep(0,32))
upper_ind = t(combn(ncol(k),2))
con_matrix_hill_cs <- adj_matrix_hill_cs[lower_ind] + adj_matrix_hill_cs[upper_ind]
k[lower_ind] = k[upper_ind] = con_matrix_hill_cs == 1

tot_pos = length(con_matrix_hill_cs)

shd = c()
total_con = c()
for (i in seq(1,12)){
  con_matrix_hills_is <- adj_matrix_hills_is[[i]][lower.tri(adj_matrix_hills_is[[i]])] + t(adj_matrix_hills_is[[i]])[lower.tri(t(adj_matrix_hills_is[[i]]))]
  acc_inverse = 1-mean(con_matrix_hills_is == con_matrix_hill_cs)
  VP_ =  sum((t(adj_matrix_hill_cs)[lower.tri(t(adj_matrix_hill_cs))&(k==1)] == t(adj_matrix_hills_is[[i]])[lower.tri(t(adj_matrix_hills_is[[i]]))&(k==1)])&(adj_matrix_hill_cs[lower.tri(adj_matrix_hill_cs)&(k==1)] == adj_matrix_hills_is[[i]][lower.tri(adj_matrix_hills_is[[i]])&(k==1)]))
  VP = sum(con_matrix_hills_is[con_matrix_hill_cs==1]==1)
  shd[i] = acc_inverse + (VP-VP_)/tot_pos
  total_con[i] = sum(con_matrix_hills_is)
}




#------------------------------------------------------------------------------------ 
#IS-approach graph (12 dyads)
#------------------------------------------------------------------------------------
all_betas_list_is = list()
for (i in seq(1,12)){
  df= CDELT(dts = duplas_agg[[i]], m_ad = adj_matrix_hills_is[[i]])
  filt = mdm_filt(
    # Data input
    dts = as.matrix(duplas_agg[[i]]),
    # Adjacency matrix
    # It's important that the adjacency matrix follows the order of the data input
    m_ad = adj_matrix_hills_is[[i]],
    # Estimated discount factor
    df$DF_hat)
  
  smoo = mdm_smoo(filt$mt, filt$Ct, filt$Rt, filt$nt, filt$dt)
  
  # Bind all parameters togheter
  all_betas <- t(Reduce(x = smoo$smt, rbind))
  all_betas_list_is[[i]] = all_betas
  
}


# IS approach with smoothed parameter


animation=list()
for (i in seq(1,12)){
  all_betas = all_betas_list_is[[i]]
  fix_plot <- all_betas |> as_tibble() |> 
    select(contains("->")) |> mutate(id = 1:dim(all_betas)[1]) |>
    pivot_longer(cols = - id) |>
    tidyr::separate(col = name, into = c("Father", "Child"), sep = "->") |>
    mutate(
      Father = factor(Father, levels = unique(total_hills_is_pf$Father)),
      Child = factor(Child, levels = unique(total_hills_is_pf$Child))
    ) |>
    ggplot(aes(y = Father, x = Child, fill = value)) +
    geom_tile() +
    scale_y_discrete(drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    scale_fill_gradient2(limits = c(-3, 3))+
    #scale_fill_gradient2()+
    #scale_fill_manual(values = c("0" = "gray90")) +
    #scale_fill_viridis_c() +  # You can change the color scale if you prefer
    labs(title = "Time Points: {frame_time}", x = "Child", y = "Father",
         fill = "Connectivity")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_segment(data = my.vlines, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1, inherit.aes = FALSE, linetype = 'dashed', color = 'orange') +
    geom_segment(data = my.hlines, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1, inherit.aes = FALSE, linetype = 'dashed', color = 'orange')+
    geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), linewidth=1, inherit.aes=F,linetype='solid',color='black')
  #theme_bw()
  
  animation[[i]] <- fix_plot +
    transition_time(id) +
    ease_aes('linear'#, interval = 0.5
    )
  print(animation[[i]],renderer = gifski_renderer())
}


#----------------------------------------------------------------------------------
# CS approach
#---------------------------------------------------------------------------------
all_betas_list_cs = list()
for (i in seq(1,12)){
  df= CDELT(dts = duplas_agg[[i]], m_ad = adj_matrix_hill_cs)
  filt = mdm_filt(
    # Data input
    dts = as.matrix(duplas_agg[[i]]),
    # Adjacency matrix
    # It's important that the adjacency matrix follows the order of the data input
    m_ad = adj_matrix_hill_cs,
    # Estimated discount factor
    df$DF_hat)
  
  smoo=mdm_smoo(filt$mt, filt$Ct, filt$Rt, filt$nt, filt$dt)
  
  # Bind all parameters togheter
  all_betas <- t(Reduce(x = smoo$smt, rbind))
  all_betas_list_cs[[i]] = all_betas
  
}


all_betas_geral_cs = all_betas_list_cs[[1]]

for (i in seq(2,12)){
  all_betas_geral_cs = all_betas_geral_cs + all_betas_list_cs[[i]]
}

all_betas_geral_cs = all_betas_geral_cs/12


fix_plot <- all_betas_geral_cs |> as_tibble() |> 
  select(contains("->")) |> mutate(id = 1:dim(all_betas)[1]) |>
  pivot_longer(cols = - id) |>
  tidyr::separate(col = name, into = c("Father", "Child"), sep = "->") |>
  mutate(
    Father = factor(Father, levels = unique(total_hill_cs_pf$Father)),
    Child = factor(Child, levels = unique(total_hill_cs_pf$Child))
  ) |>
  ggplot(aes(y = Father, x = Child, fill = value)) +
  geom_tile() +
  scale_y_discrete(drop=FALSE) +
  scale_x_discrete(drop=FALSE) +
  scale_fill_gradient2(limits = c(-7, 7))+
  #scale_fill_gradient2()+
  #scale_fill_viridis_c() +  # You can change the color scale if you prefer
  labs(title = "Time Points: {frame_time}", x = "Child", y = "Father",
       fill = "Connectivity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_segment(data = my.vlines, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1, inherit.aes = FALSE, linetype = 'dashed', color = 'orange') +
  geom_segment(data = my.hlines, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1, inherit.aes = FALSE, linetype = 'dashed', color = 'orange')+
  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), linewidth=1, inherit.aes=F,linetype='solid',color='black')
#theme_bw()

animation_geral_cs <- fix_plot +
  transition_time(id) +
  ease_aes('linear'#, interval = 0.5
  )
print(animation_geral_cs,renderer = gifski_renderer())


