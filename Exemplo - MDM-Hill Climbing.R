source('Functions - MDM - Hill Climbing.R')

# Setting the seed
set.seed(1564)
# Simulating the data
  ## Defining the true model
  ## Number o variables
  n_n = 4
  ## Creating a square matrix n_n x n_n. Represents the true adjascent matrix (DAG) 
  m_ad = array(0, dim=c(n_n,n_n)) 
  m_ad[3,1] = 1
  m_ad[3,4] = 1
  m_ad[1,2] = 1
  m_ad[4,2] = 1
  print(m_ad)
  ## k is a symmetric matrix with all the existence connections
  k = diag(rep(0,n_n))
  lower_ind = which(lower.tri(k),arr.ind=T)
  upper_ind = t(combn(ncol(k),2))
  m_con = m_ad[lower_ind] + m_ad[upper_ind]
  k[lower_ind] = k[upper_ind] = m_con == 1
  ## Observational Variance
  V = 100
  ## System Variance
  W = 0.1
  ## Sample size
  n=200
  ## Number of theta parameters (Number of connections + Number of variables)
  p = 8
  ## inital theta 
  theta_ant = rep(0,p)
  ## Data
  y=matrix(0,nrow=n,ncol=n_n)
  for (z in seq(1,n)){
    ## System equation
    theta_i = theta_ant + rnorm(p,0,sqrt(W))
    ## Observation equations
    y_3i = theta_i[1] + rnorm(1,0,sqrt(V))
    y_1i =  theta_i[2] + theta_i[3]*y_3i + rnorm(1,0,sqrt(V))
    y_4i = theta_i[4] + theta_i[5] * y_3i + rnorm(1,0,sqrt(V))
    y_2i = theta_i[6] + theta_i[7]*y_1i + theta_i[8]*y_4i + rnorm(1,0,sqrt(V))
    ## Each row of the data. The order is important.
    y[z,] = c(y_1i,y_2i,y_3i,y_4i)
    ##  New theta updated by the system eqaution.
    theta_ant = theta_i 
  }
  ##Data was tranformed into dataframe
  y = data.frame(y)
  colnames(y) = c('Y1','Y2','Y3','Y4')

# Estimation process
  ## Starting the time to measure the estimation time
  start.time = Sys.time()
  ## Here is where the MDM-Hill Climbing starts. The hc function comes from bnlearn package
  ## The default score is BIC, but notice that I chose the custom and the score selected was
  ## the "mdm_score_bn" aka LPL from MDM. You chan look for the Functions - MDM - Hill Climbing file
  ## to have more information about this function.
  hill = hc(x = data.frame(y),score='custom',fun=mdm_score_bn,args=list(nbf=15,method='Brent',call=FALSE))
  ## Building the estimated adjascent matrix (estimated DAG)
  adj_matrix_sim_hill = matrix(0,nrow=n_n,ncol=n_n,dimnames=list(colnames(y),colnames(y)))
  if (length(hill$arcs)!=0){
    for (z in seq(length(hill$arcs)/2)){
      adj_matrix_sim_hill[hill$arcs[z],hill$arcs[z+length(hill$arcs)/2]]=1
    }
  }
  print(adj_matrix_sim_hill)
  ## Ending the time
  end.time = Sys.time()
  ## Time taken to estimate the DAG
  time_hill = end.time - start.time
  print(time_hill)
  ## The estimated connections
  con_matrix_sim_hill <- adj_matrix_sim_hill[lower_ind] + adj_matrix_sim_hill[upper_ind]

# Evaluation process
  ## Accuracy of the connections
  accuracy_hill = mean(con_matrix_sim_hill == m_con)
  print(accuracy_hill)
  ## Sensitivity of the connections
  sensitivity_hill = mean(con_matrix_sim_hill[m_con==1]==1)
  print(sensitivity_hill)
  ## Specificity of the connections
  specificity_hill = mean(con_matrix_sim_hill[m_con==0]==0)
  print(specificity_hill)
  ## Positive Predicted Value of the connections
  ppv_hill = mean(m_con[con_matrix_sim_hill==1]==1)
  print(ppv_hill)
  ## Negative Predicted Value of the connections
  npv_hill = mean(m_con[con_matrix_sim_hill==0]==0)
  print(npv_hill)
  ## Accuracy of the directions of the connections
  d_accuracy_hill =   mean((m_ad[upper_ind][k[upper_ind]==1] == adj_matrix_sim_hill[upper_ind][k[upper_ind]==1])&(m_ad[lower_ind][k[lower_ind]==1] == adj_matrix_sim_hill[lower_ind][k[lower_ind]==1]))
  print(d_accuracy_hill)
  ## Estimation of discount factor
  DF = CDELT_cfd(y,adj_matrix_sim_hill)
  delta=seq(from=0.5, to=1.0, by=0.01)
  df_hill = delta[rowSums(DF$lpl)==max(rowSums(DF$lpl))]
  print(df_hill)

