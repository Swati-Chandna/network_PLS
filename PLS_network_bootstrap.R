######################
#REFERENCE: Chandna, S., B. Bagozzi, and S. Chatterjee. "Profile least squares estimation in networks with covariates." arXiv preprint arXiv:2412.16298 (2024).
##IMPLEMENTS SECTION 4. INFERENCE VIA BOOTSTRAP. 
#INPUTS: (1) Adjacency Y (n x n), 
#       (2) edge covariate array: ccedge (n x n x p)
#       (3) PLS ESTIMATES (as obtained via PLS_network_estimation.R): Lambda_pls[,1:dhat_pls], beta_pls, dhat_pls
#       (4) B: number of bootstrap samples
#       (5) alp=significance level for confidence interval construction, set to 0.05 (as default)

#OUTPUT: betaa_bt:p x 1 x B - bootstrap values of the linear coefficient parameter
#        Confidence intervals as reported in the paper.
#########################
require(statsr)
require(sigmoid)
require(lattice)
require(MASS)
require(matrixcalc)
require(graphon)
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
require(DescTools)
require(devtools)

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



expand_resid_k_to_n_NEW<- function(B_kbyk,n,xihat,...){
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

B=1999

Y_=Y
nA=dim(Y_)
n=nA[1]

MC=1####irrelevant for data examples, so set to 1

dhat_ours_bt=matrix(0,B,MC)
betaa_final=matrix(data=NA,p,MC)
betaa_bt=array(data=NA, dim=c(p,MC,B))
betaa_CI_perc_l=array(NA,dim=c(2,p,MC))
betaa_CI_basic_l=array(NA,dim=c(2,p,MC))

CI_nor_low=matrix(data=NA,p,MC)
CI_nor_high=matrix(data=NA,p,MC)
CI_nor_low_nocorr=matrix(data=NA,p,MC)
CI_nor_high_nocorr=matrix(data=NA,p,MC)



for (mc in 1:MC)
{
  
  dhat_ours_final=dhat_pls
  Lambda_final=Lambda_pls[,1:dhat_pls]
  betaa_final[,mc]=beta_pls
  coveffect_final=coveffect_pls
  
  
 
  Y_k_bt=array(data=NA, dim=c(n,n,B))
  W_b_ij=matrix(0,n^2,1)
  Xtilde_bt=matrix(0,n^2,p)

  betedge_update_bt=array(data=NA, dim=c(n,n,p))
  coveffect_bt=array(data=NA, dim=c(n,n,B))
  B_k_bt=array(data=NA, dim=c(n,n,B))
  
  
  
  model <- Mclust(Lambda_final,verbose = FALSE)
  tauhat_pls <- model$classification
  muhats <- matrix(model$parameters$mean, nrow = dhat_ours_final)
  
  
  Bhat_kbyk=t(muhats) %*% Ipqbkk %*% (muhats) 
  

  B_both=expand_resid_k_to_n_NEW(Bhat_kbyk,n,tauhat_pls)
  Bhat_nbyn=B_both$B_nbyn 
  Bhat_nbynsorted=B_both$B_nbyn_sorted 
  

  
 
  Khat=length(unique(tauhat_pls))


    alpha=n^(-1/2)
    W=matrix(NA, nrow=B,ncol=n)
    
  
  for (b in 1:B){
    
    W_b_ij_mat=matrix(NA,nrow=n,ncol=n) 
    
      W[b,]=sqrt(rexp(n,alpha)) 
      
  
      for (i in 1:n){
        for (j in i:n){
          
          W_b_ij_mat[i,j]=W[b,i]*W[b,j]
          W_b_ij_mat[j,i]=W[b,i]*W[b,j];
        }
      }
      
      W_b_ij=vec(W_b_ij_mat)
      
    
      Y_k_bt[,,b]=Y_-Bhat_nbyn
    
    Ytilde_bt=(W_b_ij)*vec(t(Y_k_bt[,,b]));
    
    
    for (ll in 1:p){
      Xtilde_bt[,ll]=(W_b_ij)*vec(t(ccedge[,,ll]));
    }
    
    
    lmobj_bt=lm((Ytilde_bt)~(Xtilde_bt)-1) 
    betaa_bt[,mc,b]=lmobj_bt$coefficients
    betaa_bt[is.na(betaa_bt)]=0
    
    
    for(l in 1:p)
    {
      betedge_update_bt[,,l]=(betaa_bt[l,mc,b]*mat_ones)*ccedge[,,l];
    }
    
    coveffect_bt[,,b]=rowSums(betedge_update_bt,dim=2);
   
    
    
    
  alp=0.05
  basic_perc_low=round((B+1)*(alp)/2)
  basic_perc_high=floor((B+1)*(1-alp/2))

  for (l in 1:p){
    betaa_bt_l_sorted=sort(betaa_bt[l,mc,])
    betaa_CI_perc_l[,l,mc]=c(betaa_bt_l_sorted[basic_perc_low],betaa_bt_l_sorted[basic_perc_high])
    betaa_CI_basic_l[,l,mc]=c(2*betaa_final[l,mc]-betaa_bt_l_sorted[basic_perc_low],2*betaa_final[l,mc]-betaa_bt_l_sorted[basic_perc_high])
    
  }
  
  bias_beta_bt=matrix(0,p,1)
  var_beta_bt=matrix(0,p,1)
  sd_beta_bt=matrix(0,p,1)
  
  
  znoralpby2=1.96 ##assuming alp=0.05
  for  (l in 1:p){
    betaa_bt_l_sorted=sort(betaa_bt[l,mc,])
    bias_beta_bt[l]=mean(betaa_bt_l_sorted)-betaa_final[l,mc]
    var_beta_bt[l]=var(betaa_bt[l,mc,])
    sd_beta_bt[l]=sqrt(var_beta_bt[l])
    
    CI_nor_low[l,mc]=betaa_final[l,mc]-bias_beta_bt[l]-sd_beta_bt[l]*znoralpby2
    CI_nor_high[l,mc]=betaa_final[l,mc]-bias_beta_bt[l]+sd_beta_bt[l]*znoralpby2
    
    ###without bias correction:
    CI_nor_low_nocorr[l,mc]=betaa_final[l,mc]-sd_beta_bt[l]*znoralpby2
    CI_nor_high_nocorr[l,mc]=betaa_final[l,mc]+sd_beta_bt[l]*znoralpby2
    
  }
  
  }
  
}

print(beta_pls)

  hist(betaa_bt[1,1,])
  hist(betaa_bt[2,1,])
  hist(betaa_bt[3,1,])
  
 
  
  