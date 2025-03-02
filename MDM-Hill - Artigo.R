source('Functions - MDM - Hill Climbing.R')


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


# Time Series Graph from Student 11
aux = tail(aluno11,3000)
aux$id = seq(1,nrow(aux))
ggplot(melt(aux,id = 'id'),aes(x=id,y=value,color=variable)) + 
  geom_line()+ labs(x ="Time Points", y = "Oxyhemoglobin",color='Channels')+
  ylim(-0.002,0.002)


# Time Series Graph from Teachert 11

aux = tail(professor11,3000)
aux$id = seq(1,nrow(aux))
ggplot(melt(aux,id = 'id'),aes(x=id,y=value,color=variable)) + 
  geom_line()+ labs(x ="Time Points", y = "Oxyhemoglobin",color='Channels')+
  ylim(-0.002,0.002)


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
                                                          'Teacher_PFE','Teacher_PFD','Teacher_TP')

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

#Proportions of estimated connections among all possible connections for each region.
print(neuro_matrix_real)

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


#----------------------------------------------------------------------------------
# CS approach graph
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
  scale_fill_gradient2(limits = c(-3, 3))+
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


#-------------------------------------------------------------
# ROIs' mean representation (time-invariant global pattern).
#-------------------------------------------------------------

regex = c("SV\\d+->SV\\d+","TV\\d+->TV\\d+","SV\\d+->TV\\d+","TV\\d+->SV\\d+")

father_level = list(unique(total_hill_cs_pf$Father)[1:16],unique(total_hill_cs_pf$Father)[17:32],
                    unique(total_hill_cs_pf$Father)[1:16],unique(total_hill_cs_pf$Father)[17:32])

child_level = list(unique(total_hill_cs_pf$Child)[1:16],unique(total_hill_cs_pf$Child)[17:32],
                   unique(total_hill_cs_pf$Child)[17:32],unique(total_hill_cs_pf$Child)[1:16])

rois_mean = list()

father_regions = list(c("V13",  "V14",  "V15",  "V17","V18","V20","V21","V23"),
                      c("V13",  "V14",  "V15",  "V17","V18","V20","V21","V23"),
                      c("V13",  "V14",  "V15",  "V17","V18","V20","V21","V23"),
                      c("V7",  "V8",  "V10",  "V11"),
                      c("V7",  "V8",  "V10",  "V11"),
                      c("V7",  "V8",  "V10",  "V11"),
                      c("V1",  "V2",  "V4",  "V5"),
                      c("V1",  "V2",  "V4",  "V5"),
                      c("V1",  "V2",  "V4",  "V5"))

child_regions = list(c("V1",  "V2",  "V4",  "V5"),
                     c("V7",  "V8",  "V10",  "V11"),
                     c("V13",  "V14",  "V15",  "V17","V18","V20","V21","V23"),
                     c("V1",  "V2",  "V4",  "V5"),
                     c("V7",  "V8",  "V10",  "V11"),
                     c("V13",  "V14",  "V15",  "V17","V18","V20","V21","V23"),
                     c("V1",  "V2",  "V4",  "V5"),
                     c("V7",  "V8",  "V10",  "V11"),
                     c("V13",  "V14",  "V15",  "V17","V18","V20","V21","V23"))

person_father = c('S','T','S','T')
person_child = c('S','T','T','S')
for (i in seq(1,4)){
  roi_mean=c()
  all_betas_geral_new<-all_betas_geral_cs |> as_tibble() |> 
    select(matches(regex[i])) |> mutate(id = 1:dim(all_betas)[1]) |>
    pivot_longer(cols = - id) |>
    tidyr::separate(col = name, into = c("Father", "Child"), sep = "->") |>
    mutate(
      Father = factor(Father, levels = father_level[[i]]),
      Child = factor(Child, levels = child_level[[i]])
    ) 
  for (j in seq(1,9)){
    father_region = paste(person_father[[i]],father_regions[[j]],sep="")
    child_region = paste(person_child[[i]],child_regions[[j]],sep="")
    roi_mean[j] = mean(abs(all_betas_geral_new[(all_betas_geral_new$Father %in% father_region & 
                                                  all_betas_geral_new$Child %in% child_region), 'value']$value))                                                                                     
  }
  rois_mean[[i]] = roi_mean
}

roi_graph<- function(i){
  rois_data = data.frame(matrix(rois_mean[[i]],nrow=3,ncol=3,byrow = TRUE))
  rois_data[is.na(rois_data)]=0
  colnames(rois_data) = c('LPF','RPF','RTPJ')
  rownames(rois_data) = c('RTPJ','RPF','LPF')
  rois_data|> mutate(Father = rownames(rois_data)) |>
    pivot_longer(cols = - Father) |> mutate(Child=name)|>
    mutate(name=NULL) |> 
    ggplot(aes(y = Father, x = Child, fill = value)) +
    geom_tile() +
    scale_y_discrete(drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    scale_fill_gradient2(limits = c(0, 1.06))+
    #scale_fill_gradient(low = "yellow", high = "red")+
    #scale_color_distiller(palette='YlOrBr')+
    #scale_fill_viridis_c() +  # You can change the color scale if you prefer
    labs(x = "Child region", y = "Parent region",
         fill = "Connectivity")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
}


# Inside student brain
roi_graph(1)

# Inside teacher brain
roi_graph(2)

# Inside student to teacher brain
roi_graph(3)

#Inside teacher to student brain
roi_graph(4)








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


