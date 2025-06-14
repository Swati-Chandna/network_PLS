
require(sbm)
data(fungusTreeNetwork)


tree_mat=fungusTreeNetwork$tree_tree 
dtree=dim(tree_mat)
n=dtree[1]

tree_mat01=matrix(0,n,n)
for (i in 1:(n-1))
{
  for ( j in (i+1):n)
  {
    if(tree_mat[i,j]>=1)
    {tree_mat01[i,j]=1
    tree_mat01[j,i]=1
    }
  }
}

###############edge covariates:
covar_1=fungusTreeNetwork$covar_tree$genetic_dist
covar_2=fungusTreeNetwork$covar_tree$taxonomic_dist
covar_3=fungusTreeNetwork$covar_tree$geographic_dist
###############

Y=tree_mat01

p=3
ccedge=array(NA,dim=c(n,n,p))
ccedge[,,1]=covar_1
ccedge[,,2]=covar_2
ccedge[,,3]=covar_3
#######################

