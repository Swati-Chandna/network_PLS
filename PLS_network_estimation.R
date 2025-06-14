#########################################PROPOSED PLS with covariates:
#REFERENCE: Chandna, S., B. Bagozzi, and S. Chatterjee. "Profile least squares estimation in networks with covariates." arXiv preprint arXiv:2412.16298 (2024).


####################TO RUN: 
####Load the adjacency matrix: Y as an n by n MATRIX
###and edge covariates: ccedge[,,] as a n by n by p ARRAY
####EXAMPLE DATASET (loaded by): Fungal_tree_net.R

####OUTPUTS: 
##Linear coefficient vector: beta_pls[1:p,1], 
##Residual term estimate: resid_pls_sorted an n by n MATRIX (sorted basd one estimated clusters under a GMM)
##final Phat from PLS: Phat_clustB_sorted n by n MATRIX
####estimated clustering of nodes in the residual matrix: tauhat_pls


require(statsr)
require(sigmoid)
require(lattice)
require(MASS)
require(matrixcalc)
require(graphon)
require(grdpg)
require(igraph)
require(reshape2);
require(ggplot2);
require(grid);
require(loe);
require(vegan);
require(mgcv);
require(fcd)
require(cccd);
require(GoFKernel)
require(fields)

library(Matrix)
library(dplyr)
library(mclust)
library(jsonlite)
library(irlba)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(igraph)
library(graphstats)
library(grdpg)
library(R.matlab)
library(Rtsne)
library(ggplot2)



MC=1##irrelevant for data settings, set to 1
########################
expand_resid_k_to_n<- function(B_kbyk,n,xihat,...){
  B_nbyn=matrix(NA,n,n)
  B_nbyn_sorted=matrix(NA,n,n)
  pi_xihat=order(xihat)
  
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      B_nbyn[i,j]=B_kbyk[xihat[i],xihat[j]]
      
    }
  }
  B_nbyn_sorted=B_nbyn[pi_xihat,pi_xihat]
  return(list(B_nbyn=B_nbyn,B_nbyn_sorted=B_nbyn_sorted))
}

evaluate_LScrit_iter<- function(n,A,Phatk,...){
  #   ## Input:
  #   # n : order of the graph to sample
  #  # A  : n x n adjacency matrix
  # Phatk: estimated model at the kth step of maxiter
  #   ## Output:value of the LS criterion given the estimates of the model terms and observed data (adjacency and covariates)
  
  LSobj_mat=matrix(0,n,n)
  LSobj_mat=(A-Phatk)^2
  indices_upp=upper.tri(LSobj_mat, diag=FALSE)
  LSobj=sum(LSobj_mat[indices_upp])
  
  return(LSobj=LSobj)
}

####SCALING of COVARIATES:
 maxccedge=matrix(0,p,1)
 minccedge=matrix(0,p,1)
 for (l in 1:p)
 {
   maxccedge[l]=max(ccedge[,,l])
   minccedge[l]=min(ccedge[,,l])
 }


 for (l in 1:p)
 {
   if(maxccedge[l]>1)
   {
     ccedge1=ccedge[,,l];
     ccedge[,,l]=(ccedge1 - minccedge[l]) / (maxccedge[l] - minccedge[l])
   }
 }

Y_=Y



SBMblock=0
SBMblock_final=1 #clustering of the residual matrix, only in the final step (post convergence)

maxiter=1000 ###maximum number of iterations set to 1000 as default.
betaa=array(NA,dim=c(p,MC,maxiter));

resid_hat_pls=array(0,dim=c(n,n,MC)) 
resid_hat_pls_sorted=array(0,dim=c(n,n,MC))

resid_pls_sorted_all=array(NA,dim=c(n,n,MC))

resid_pls_sorted=matrix(0,n,n)
coveffect_final_sorted=matrix(0,n,n)
tauhat_pi=matrix(0,n,1)
tauhat_pls_pi=matrix(0,n,1)

beta_pls=matrix(NA,p,MC) ###final PLS estimates of the linear coefficient

dhat_ours_all=matrix(NA,maxiter,MC)

LS=matrix(NA, maxiter,MC)
LS_diff=matrix(1000, maxiter,MC)####check stopping criterion


tauhat_pls_by_muhats=matrix(NA,n,1) 
resid_hat_pls_sorted_bymeans=array(0,dim=c(n,n,MC))
tauhat_plsmeans_pi=matrix(NA,n,1)
Phat_clustB_sortedbymeans=matrix(NA,n,n)
coveffect_final_sortedbymeans=matrix(0,n,n)
Lambda_pls=matrix(0,n,10)
#beta_pls=matrix(0,n,10)
#coveffectPLS
d=2

for (mc in 1:MC){  
  dmax=10
  dinit_ours=dmax 

  Xtilde=matrix(NA,n^2,p)
  Ytilde=matrix(NA,n^2,p)
  
  Lambda=array(NA,dim=c(n,dinit_ours,maxiter));###d: true here bec we know
  Phat_k=array(NA,dim=c(n,n,maxiter));
  
  
  eta_k_sorted=array(NA,dim=c(n,n,maxiter));
  
  
  Phat_k_sorted=array(NA,dim=c(n,n,maxiter));
  
  
  
  betedge=array(NA,dim=c(n,n,p));
  betedge_update=array(NA,dim=c(n,n,p));
  betedge_update1=array(NA,dim=c(n,n,p));
 
  coveffect=array(NA,dim=c(n,n,maxiter));
  coveffect_sorted=array(NA,dim=c(n,n,maxiter));
  eta_k=array(NA,dim=c(n,n,maxiter));
  
  
  B_k=array(NA,dim=c(n,n,maxiter));
  LSmat=array(NA,dim=c(n,n,maxiter))
  
  
  Y_k=array(data=NA, dim=c(n,n,maxiter))
  
  ############################
  T1=5
  t1=-0.9
  t2=0.9
  init_vals1 <- seq(from = t1, to = t2, length.out = T1)
  init_vals=c(-2,2,init_vals1)
  T=length(init_vals)
  ###########################
  dhat_ours_init=matrix(0,T,1)
  beta_pls_init=matrix(0,p,T)
  Lambda_pls_init=array(0,dim=c(n,dmax,T))
  modelbic_init=matrix(0,T,1)
  ###########################

  for (t in 1:T)
  {
  betaa[1:p,mc,1]=init_vals[t]*matrix(1,p,1); 
  
  onesnbyn=matrix(1,n,n);
  B_k_scaled_hat=array(NA, dim=c(n,n,maxiter))
  
  
  
  for (k in 1:1){
    mat_ones=matrix(1,n,n);
    for(l in 1:p)
    {
      betedge[,,l]=(betaa[l,mc,k]*mat_ones)*ccedge[,,l];###n by n by p
    }
    
    coveffect[,,k]=rowSums(betedge,dims = 2);
    
  }
  
  bkk_init=Y-coveffect[,,1]
  
  dmax=10
  embed_Adj_init<- SpectralEmbedding(bkk_init, dmax, maxit = 10000, work = 50)
  s_init <-embed_Adj_init$D
  dhat_init = dimselect(s_init)$elbow[1]+1
  Lambda[,1:dhat_init,1]= embed_Adj_init$X[,1:dhat_init] %*% sqrt(diag(s_init[1:dhat_init], nrow=dhat_init, ncol=dhat_init))#
  

  LP_B11=array(NA,dim=c(n,dinit_ours,maxiter))
  LP_B=matrix(0,n,n)
  dhat_ours_all[1,mc]=dhat_init
  
  
  for (k in 1:(maxiter-1)){
    for(l in 1:p)
    {
      betedge[,,l]=(betaa[l,mc,k]*mat_ones)*ccedge[,,l];
    }
    
    
    coveffect[,,k]=rowSums(betedge,dims = 2);
    
    
    
    if(k==1){
      Ipq <- getIpq(bkk_init, dhat_init)#
    }else
    {Ipq <- getIpq(bkk, dhat_ours_all[k,mc])#
    }
    
    
    if(SBMblock==1)
    { 
      model <- Mclust(Lambda[,1:dhat_ours_all[k,mc],k],verbose = FALSE)
      xihat <- model$classification
      muhats <- matrix(model$parameters$mean, nrow = dhat_ours_all[k,mc])
      
      B_kbyk=t(muhats) %*% Ipq %*% (muhats)
      B_both=expand_resid_k_to_n(B_kbyk,n,xihat)
      B_nbyn=B_both$B_nbyn  
      
      eta_k[,,k]=coveffect[,,k] + B_nbyn 
    }else
    { eta_k[,,k]=coveffect[,,k] + (Lambda[,1:dhat_ours_all[k,mc],k]) %*% Ipq %*% t(Lambda[,1:dhat_ours_all[k,mc],k])
    }
    
  
    Phat_k[,,k]=eta_k[,,k];
    
    Phat_k[,,k][Phat_k[,,k]>1]=1
    Phat_k[,,k][Phat_k[,,k]<0]=0
    
    
    LS[k,mc]=evaluate_LScrit_iter(n,Y_,Phat_k[,,k])
   
    if(k>1)
    {
      LS_diff[k-1,mc]=abs(round(LS[k-1,mc],4)-round(LS[k,mc],4))
      if(LS_diff[k-1,mc]==0)
      {
        cat("Condition met at iteration k =", k, "\n")
        break  # Exits the k loop
      }
    }
    
  
    if(k==(maxiter))
    {
      LS_diff[k-1,mc]=abs(round(LS[k-1,mc],4)-round(LS[k,mc],4))
      if(LS_diff[k-1,mc]>0)
      {
        cat("Condition not met in maximum number of iterations set to", maxiter, "\n")
        cat("Re-start with another initial value", "\n")
        break  # Exits the innermost loop
      }
      
    }
    

    
    
  
    
    if(SBMblock==1)
    {
      model <- Mclust(Lambda[,1:dhat_ours_all[k,mc],k],verbose = FALSE)
      xihatrr <- model$classification
      muhats <- matrix(model$parameters$mean, nrow = dhat_ours_all[k,mc])
      #####NEW:
      Bhat_kbyk=t(muhats) %*% Ipq %*% (muhats)
      B_both=expand_resid_k_to_n(Bhat_kbyk,n,xihatrr)
      Bhat_nbyn=B_both$B_nbyn
      ####
      
      Y_k[,,k]=Y_-Bhat_nbyn
    }else{
      Y_k[,,k]=Y_-(Lambda[,1:dhat_ours_all[k,mc],k]) %*% Ipq %*% t(Lambda[,1:dhat_ours_all[k,mc],k]) 
    }
    
    Ytilde=vec(t(Y_k[,,k]));
    
    for (ll in 1:p){
      Xtilde[,ll]=vec(t(ccedge[,,ll]));
    }
    
    
    lmobj=lm((Ytilde)~(Xtilde)-1)  
    
    betaa[,mc,k+1]=lmobj$coefficients 
  
    betaa[is.na(betaa)]=0
    
   
    for(l in 1:p)
    {
      betedge_update[,,l]=(betaa[l,mc,k+1]*mat_ones)*ccedge[,,l];
    }

    
    coveffect[,,k+1]=rowSums(betedge_update,dim=2);
    
   
    if(SBMblock==1)
    { 
      model <- Mclust(Lambda[,1:dhat_ours_all[k,mc],k],verbose = FALSE)
      xihatrrr <- model$classification
      muhats <- matrix(model$parameters$mean, nrow = dhat_ours_all[k,mc])
      
      Bhat_kbyk=t(muhats) %*% Ipq %*% (muhats)
      B_both=expand_resid_k_to_n(Bhat_kbyk,n,xihatrrr)
      Bhat_nbyn=B_both$B_nbyn
      
      eta_k[,,k+1]=coveffect[,,k+1] + Bhat_nbyn 
    }else{eta_k[,,k+1]=coveffect[,,k+1] +(Lambda[,1:dhat_ours_all[k,mc],k]) %*% Ipq %*% t(Lambda[,1:dhat_ours_all[k,mc],k])
    }
    
    Phat_k[,,k+1]=eta_k[,,k+1];
    
    
    B_k[,,k]=Y_-coveffect[,,k+1]  
    bkk=B_k[,,k]
    
   
   
    embed_Bscaled <- SpectralEmbedding(bkk, dmax, maxit = 10000, work = 50)
    sB <- embed_Bscaled$D

    dhat_ours_all[k+1,mc]=dimselect(sB)$elbow[1]+1;
    
    Ipqbkk=getIpq(bkk,dhat_ours_all[k+1,mc])
    LP_B[,1:dhat_ours_all[k+1,mc]] <- embed_Bscaled$X[,1:dhat_ours_all[k+1,mc]] %*% sqrt(diag(sB[1:dhat_ours_all[k+1,mc]], nrow=dhat_ours_all[k+1,mc], ncol=dhat_ours_all[k+1,mc]))

    # LP_B11[,1:dhat_ours_all[k+1,mc],k]=LP_B[,1:dhat_ours_all[k+1,mc]]
   
      Lambda[,1:dhat_ours_all[k+1,mc],k+1]=LP_B[,1:dhat_ours_all[k+1,mc]];
  
   
    if(SBMblock==1)
    { 
      model <- Mclust(Lambda[,1:dhat_ours_all[k+1,mc],k+1],verbose = FALSE)
      xihat <- model$classification
      muhats <- matrix(model$parameters$mean, nrow = dhat_ours_all[k+1,mc])
      
      Bhat_kbyk=t(muhats) %*% Ipqbkk %*% (muhats)
      B_both=expand_resid_k_to_n(Bhat_kbyk,n,xihat)
      Bhat_nbyn=B_both$B_nbyn
      
      eta_k[,,k+1]=coveffect[,,k+1] + Bhat_nbyn
    }else
    {eta_k[,,k+1]=coveffect[,,k+1] +(Lambda[,1:dhat_ours_all[k+1,mc],k+1]) %*% Ipqbkk %*% t(Lambda[,1:dhat_ours_all[k+1,mc],k+1])
    }
    
    
    Phat_k[,,k+1]=(eta_k[,,k+1]);
    
  
  }######where the k over maxiter loop ends
  
  maxiter_stop=k ###convergence iteration number
  
  
  beta_pls[1:p,mc]=betaa[1:p,mc,maxiter_stop]
  
  if(SBMblock_final==1)#
  { 
    model <- Mclust(Lambda[,1:dhat_ours_all[maxiter_stop,mc],maxiter_stop],verbose = FALSE)
    
    #plot(model, what = "classification")####check validity of GMM
    
    modeluncert=mean(model$uncertainty)
    modeloglik=model$loglik
    modelbic=model$bic
    ####################################
    
    
  }else 
  {resid_hat_pls[,,mc]=(Lambda[,1:dhat_ours_all[maxiter,mc],maxiter]) %*% Ipqbkk %*% t(Lambda[,1:dhat_ours_all[maxiter,mc],maxiter])
  tauhat_pls_pi=order(tauhat_pls)
  resid_hat_pls_sorted[,,mc]=resid_hat_pls[tauhat_pls_pi,tauhat_pls_pi,mc]
  }
  
  dhat_ours_init[t]=dhat_ours_all[maxiter_stop,mc]
  
  beta_pls_init[,t]=beta_pls
  Lambda_pls_init[,1:dhat_ours_init[t],t]=Lambda[,1:dhat_ours_all[maxiter_stop,mc],maxiter_stop]
  modelbic_init[t]=modelbic
  #rm(beta_pls)
  }####initialization loop over t ends
  
minBICindex=which(modelbic_init==min(modelbic))
beta_pls=beta_pls_init[,minBICindex]
dhat_pls=dhat_ours_init[minBICindex]
Lambda_pls[,1:dhat_pls]=Lambda_pls_init[,1:dhat_pls,minBICindex]

  #############################ordering and PLOTTING etc:

model_pls <- Mclust(Lambda_pls[,1:dhat_pls],verbose = FALSE)
tauhat_pls<- model_pls$classification
muhats <- matrix(model_pls$parameters$mean, nrow = dhat_pls)



Bhat_kbyk=t(muhats) %*% Ipqbkk %*% (muhats)
B_both=expand_resid_k_to_n(Bhat_kbyk,n,tauhat_pls)
Bhat_nbyn=B_both$B_nbyn

resid_hat_pls[,,mc]= Bhat_nbyn

tauhat_pls_pi=order(tauhat_pls)
resid_hat_pls_sorted[,,mc]=resid_hat_pls[tauhat_pls_pi,tauhat_pls_pi,mc]


clustker=diag(Bhat_kbyk)

muhat_pi=rank(-clustker)

tauhatlabel_bymeans=muhat_pi


for (s in 1:max(tauhat_pls))
{
  inds=which(tauhat_pls==s)
  tauhat_pls_by_muhats[inds]=muhat_pi[s]
}


tauhat_plsmeans_pi=order(tauhat_pls_by_muhats)
newlabelsorder=tauhat_pls_by_muhats[tauhat_plsmeans_pi]

resid_hat_pls_sorted_bymeans[,,mc]=resid_hat_pls[tauhat_plsmeans_pi,tauhat_plsmeans_pi,mc]

##################################
  resid_pls_sorted=resid_hat_pls_sorted[,,mc]
  
  
  coveffect_final=coveffect[,,maxiter_stop];
  Phat_clustB=coveffect_final+Bhat_nbyn 
  
  Phat_clustB[Phat_clustB<0]=0
  Phat_clustB[Phat_clustB>1]=1
  #image.plot(Phat_clustB)
  Phat_clustB_sorted=Phat_clustB[tauhat_pls_pi,tauhat_pls_pi]
  
  Phat_clustB_sortedbymeans=Phat_clustB[tauhat_plsmeans_pi,tauhat_plsmeans_pi]
  #image.plot(Phat_clustB_sortedbymeans)
 # image.plot(Phat_clustB_sorted)
  
}##mc ends


#coveffect_final=coveffect[,,maxiter_stop]
for(l in 1:p)
{
  betedge[,,l]=(beta_pls[l]*mat_ones)*ccedge[,,l];###n by n by p
}

coveffect_pls=rowSums(betedge,dims = 2);
coveffect_final_sorted=coveffect_pls[tauhat_pls_pi,tauhat_pls_pi]
coveffect_final_sortedbymeans=coveffect_pls[tauhat_plsmeans_pi,tauhat_plsmeans_pi]


par(mfrow = c(2, 2))
image.plot(resid_pls_sorted)
title('Estimated residual term')
image.plot(Phat_clustB_sorted)
title('Phat from PLS')###order corresponding to resid_pls_sorted 
#plot.new
image.plot(coveffect_final_sorted)
title('Covariate effects')

#image.plot(coveffect_final_sortedbymeans)

print(beta_pls)
#print(modelbic)
#dhat_ours_all[maxiter_stop]
#max(tauhat_pls)

