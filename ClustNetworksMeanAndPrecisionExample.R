rm(list=ls())

setwd("C:\\NetworkClustering")

# In this program, clustering of means and graphs is done in parallel. For means,
# DPM is used and graphs by mixture of multivariate normal. 


library(MCMCpack)
library(mvtnorm)
library(pscl)
library(psych)
library(invgamma)
library(abind)

SLURM_ID=1
nData<-1
nodes<-10
n1<-n2<-100

# Chain networks
# cluster 1
precision1<- matrix(rep(0,(nodes*nodes)),nrow=nodes)

cov1<-cov2<-matrix(NA,nrow=nodes,ncol=nodes)

a<-0.1
diag(cov1)<-1
s<-rep(0,nrow(precision1))
s[1]<-1
set.seed(12345)
for (i in 2:length(s))
{
  s[i]<-s[i-1]+runif(1,0.5,1)
}

for (i in 1:(nrow(precision1)-1))
{
  for (j in (i+1):nrow(precision1))
  {
    cov1[i,j]<-exp(-a*abs(s[j]-s[i]))
    cov1[j,i]<-cov1[i,j]
  }
}
precision1<-solve(cov1)

lower<-precision1*lower.tri(precision1)

loc<-which(abs(lower)<(10*10^(-4)),arr.ind=TRUE)
loc1<-which(abs(lower)>(10*10^(-4)),arr.ind=TRUE)

for (i in 1:nrow(loc))
{
  lower[loc[i,1],loc[i,2]]<-0
}
for (i in 1:nrow(loc1))
{
  lower[loc1[i,1],loc1[i,2]]<-1
}

G1<-lower+t(lower)


# cluster 2
# Precision 2
precision2<- matrix(rep(0,(nodes*nodes)),nrow=nodes)

cov2<-matrix(NA,nrow=nodes,ncol=nodes)

a<-0.01
diag(cov2)<-1
s<-rep(0,nrow(precision2))
s[1]<-1
for (i in 2:length(s))
{
  s[i]<-s[i-1]+runif(1,0.5,1)
}

for (i in 1:(nrow(precision2)-1))
{
  for (j in (i+1):nrow(precision2))
  {
    cov2[i,j]<-exp(-a*abs(s[j]-s[i]))
    cov2[j,i]<-cov2[i,j]
  }
}
precision2<-solve(cov2)

precision2[2,3]<-precision2[3,2]<-precision2[3,4]<-precision2[4,3]<-0

cov2<-solve(precision2)

lower<-precision2*lower.tri(precision2)

loc<-which(abs(lower)<(10*10^(-4)),arr.ind=TRUE)
loc1<-which(abs(lower)>(10*10^(-4)),arr.ind=TRUE)

for (i in 1:nrow(loc))
{
  lower[loc[i,1],loc[i,2]]<-0
}
for (i in 1:nrow(loc1))
{
  lower[loc1[i,1],loc1[i,2]]<-1
}

G2<-lower+t(lower)


# cluster 3 -- a random graph
# Number of nodes
set.seed(12345)
n_nodes <- nodes
n_edges <- 10

# Function to create a random symmetric positive definite matrix
create_precision_matrix <- function(n) {
  mat <- matrix(rnorm(n * n), n, n)
  mat <- mat + t(mat)  # Make it symmetric
  mat <- mat + n * diag(n)  # Make it diagonally dominant, hence positive definite
  return(mat)
}

# Generate the precision matrix
precision_matrix <- create_precision_matrix(n_nodes)

# Function to threshold the precision matrix to get exactly 20 edges
create_adj_matrix <- function(precision_matrix, edges) {
  # Get the upper triangular part of the matrix (excluding diagonal)
  upper_tri <- upper.tri(precision_matrix, diag = FALSE)
  abs_values <- abs(precision_matrix[upper_tri])

  # Get the threshold value to have exactly 'edges' number of edges
  threshold <- sort(abs_values, decreasing = TRUE)[edges]

  # Create adjacency matrix by thresholding the precision matrix
  adj_matrix <- ifelse(abs(precision_matrix) >= threshold, 1, 0)

  # Ensure it's symmetric and remove self-loops
  adj_matrix <- adj_matrix * upper.tri(adj_matrix, diag = TRUE)
  adj_matrix <- adj_matrix + t(adj_matrix)
  diag(adj_matrix) <- 0

  return(adj_matrix)
}

# Create the adjacency matrix
G3 <- create_adj_matrix(precision_matrix, n_edges)

precision3<-G3*(precision_matrix)

loc<-which(precision3<0,arr.ind=TRUE)
if (nrow(loc)>=1)
{
  for (jj in 1:nrow(loc))
  {
    precision3[loc[jj,1],loc[jj,2]]<-precision3[loc[jj,1],loc[jj,2]]-2
  }
}
loc<-which(precision3>0,arr.ind=TRUE)
if (nrow(loc)>=1)
{
  for (jj in 1:nrow(loc))
  {
    precision3[loc[jj,1],loc[jj,2]]<-precision3[loc[jj,1],loc[jj,2]]+1
  }
}
diag(precision3)<-diag(precision_matrix)

#eigen(precision3)$values

mean11<-1
mean12<-5
mean21<-2
mean22<-3
mean3<-2

for (i in 1:nData)
{
  set.seed((SLURM_ID*100*i))
  x11<-rmvnorm((n1/2), mean = rep(mean11, nodes), sigma = solve(precision1))
  x12<-rmvnorm((n1/2), mean = rep(mean12, nodes), sigma = solve(precision1))
  x21<-rmvnorm((n1/2), mean = rep(mean21, nodes), sigma = solve(precision2))
  x22<-rmvnorm((n1/2), mean = rep(mean22, nodes), sigma = solve(precision2))
  x3<-rmvnorm(n1, mean = rep(mean3, nodes), sigma = solve(precision3))

if (i==1) 
  {
    # 3 clusters
    dataAll<-rbind(x11,x12,x21,x22,x3)
    # 2 clusters
 #    dataAll<-rbind(x21,x22,x3)
  } else
  {
    # 3 clusters
    dataAll<-rbind(dataAll,x11,x12,x21,x22,x3)
    # 2 clusters
 #    dataAll<-rbind(dataAll,x21,x22,x3)
  }
}

G1True<-G1
G2True<-G2
G3True<-G3

# the following is for two clusters, one chain and one regular graph
#G1True<-G2
#G2True<-G3

GTrue<-abind(G1True, G2True,G3True,along=3)
#GTrue<-abind(G1True, G2True,along=3)

O1True<-precision1
O2True<-precision2
O3True<-precision3

# the following is for two clusters, one chain and one regular graph
#O1True<-precision2
#O2True<-precision3

# conditional posterior of zeta
dzeta<-function(par1,par2,par3){ #return in log scale
  # par1: zeta, par2: sum of log Pi, par3: length of Pi
  if(par1>1) tmp<-2*log(par1)
  else tmp<-0
  return(lgamma(par3*par1)+{par1-1}*par2-lgamma(par1)*par3-tmp)
}

# sampling of D (dn)
sampdn<-function(sampsize,data,zeta,C,dn,Omega,m.cur,mu.cur)
{
  pdn<-rep(0,C)
  
  prob<-rdirichlet(1,rep(zeta,C)+colSums(dn))
  
  for(i in 1:sampsize){
    mu<-mu.cur[m.cur[i],]
    
    covi<-(data[i,]-mu)%*%t(data[i,]-mu)
    
    for(j in 1:C){
      trace<-sum(diag(covi%*%Omega[,,j]))
      if (det(Omega[,,j])<=0) 
      {
        pdn[j]=-99999
      } else
      {
        pdn[j]<-log(prob[j])+0.5*(log(det(Omega[,,j]))-trace)
      }
    }
    pdn<-pdn-max(pdn)
    dn[i,]<-rmultinom(1,1,exp(pdn))
  }
  return(list(prob=prob, dn=dn))
}

# O:	precision matrix assuming one cluster
# covX:	sample covariance matrix based on the complete data
lhCal<-function(O,covX,n)
{
  Trace<-(sum(diag(covX%*%O)))

  if (det(O)<0)
  {
    lh<- -999999
  } else
  {
    lh<-n/2*(log(det(O))-Trace)
  }  
  return(lh)
}

# this is the posterior for each component in the precision matrix, O
postOij<-function(O,covX,G,i,j,n,nu0)
{
  part1<-lhCal(O,covX,n)
  if (i!=j)
  {
    part2<-G[i,j]*(-log(nu1)-(O[i,j])^2/(2*nu1^2))+(1-G[i,j])*(-log(nu0)-(O[i,j])^2/(2*nu0^2))
  }
  if (i==j)
  {
    # following Wang 2015's suggestion, the parameter lambda in exponential is 
    # taken as 1.
    part2<--O[i,i]
  }
  return(part1+part2)
}

# this is to update one node in precision matrix O. 
UpdateO<-function(O,covX,G,i,j,n,nu0)
{
  probO<-postOij(O,covX,G,i,j,n,nu0)
  
  proposal<-rnorm(1,O[i,j],jump)
  Otemp<-O
  Otemp[i,j]<-Otemp[j,i]<-proposal
  
  if (det(Otemp)<=0)  probO0proposal<- -9999999
  else probO0proposal<-postOij(Otemp,covX,G,i,j,n,nu0)
  
  ratio<-probO0proposal-probO
  judge<-log(runif(1))
  if (ratio>judge)
  {
    O[i,j]<-O[j,i]<-Otemp[i,j]
  }
  return(O)
}#end update O

# this is the posterior for each component in the graph matrix, G.
postGij<-function(O,covX,G,i,j,n,nu0)
{
  part1<-lhCal(O,covX,n)
  
  part2<-G[i,j]*(-log(nu1)-(O[i,j])^2/(2*nu1^2))+(1-G[i,j])*(-log(nu0)-(O[i,j])^2/(2*nu0^2))
  
  pi<-0.1
  partG<-log(pi/(1-pi))*(G[i,j])+log(1-pi)
  postOG<-part1+part2+partG
  return(postOG)
}

UpdateG<-function(O,covX,G,i,j,n,nu0)
{
  accept<-0
  Gtemp<-G
  if (G[i,j]==0)
  {
    probG0<-postGij(O,covX,G,i,j,n,nu0)
    Gtemp[i,j]<-Gtemp[j,i]<-1
    Otemp<-UpdateO(O,covX,Gtemp,i,j,n,nu0)
    probG1<-postGij(Otemp,covX,Gtemp,i,j,n,nu0)
    
    maxp<-max(probG1,probG0)
    ratio<-exp(probG1-maxp)/(exp(probG0-maxp)+exp(probG1-maxp))
    judge<-runif(1)
    if (ratio>judge)
    {
      G[i,j]<-G[j,i]<-1
      O[i,j]<-O[j,i]<-Otemp[i,j]
    }
  }
  if (G[i,j]==1)
  {
    probG1<-postGij(O,covX,G,i,j,n,nu0)
    Gtemp[i,j]<-Gtemp[j,i]<-0
    Otemp<-UpdateO(O,covX,Gtemp,i,j,n,nu0)
    probG0<-postGij(Otemp,covX,Gtemp,i,j,n,nu0)
    maxp<-max(probG1,probG0)
    ratio<-1-exp(probG1-maxp)/(exp(probG0-maxp)+exp(probG1-maxp))
    judge<-runif(1)
    if (ratio>judge)
    {
      G[i,j]<-G[j,i]<-0
      O[i,j]<-O[j,i]<-Otemp[i,j]
    }
  }
  return(list(G=G,O=O))
}

GraphEst<-function(O,covX,G,n,nodes1,nu0)
{
  # diagonal
  for (i in 1:nodes1)
  {
    O<-UpdateO(O,covX,G,i,i,n,nu0)
  } 
  
  for (i in 1:(nodes1-1))
  {
    for (j in (i+1):nodes1)
    {
      GG<-UpdateG(O,covX,G,i,j,n,nu0)
      G<-GG$G
      O<-GG$O
    }
    O<-UpdateO(O,covX,G,i,j,n,nu0)
  } #for (i in 1:(p-1))
  return(list(O=O,G=G))
}

# sample sigma^2 in the base distribution for mu
sampsig2<-function(n,mu)
{
  a1<-a2<-0.5
  tmp<-a1+n/2
  sum.mu<-0
  for(i in 1:nodes)
    sum.mu<-sum.mu+sum((mu[,i])^2)
  scale=a2+sum.mu/2
  sigma2<-1/rgamma(1,tmp,scale=scale)
  return(sigma2)
}

DPmu<-function(data,nodes,alpha,m.cur,mu.cur,sigma2,dn,OclustArr,GclustArr)
{
  #m.cur is the cluster assignment for mu_i's for the n subjects
  #mu.cur is a matrix consistent of unique mu's for each cluster
  for(i in 1:sampsize){
    m.i<-m.cur[i]
    m.rle<-rle(sort(m.cur[-i]))
    m.uniq.len<-length(m.rle$values)
    last<-m.uniq.len+1
    prob<-1:last
    
    c.i<-which(dn[i,]!=0)
    Oi<-OclustArr[,,c.i]
    Sigma.mui<-solve(Oi+diag(1/sigma2,nrow=nodes))
    yOmega<-data[i,]%*%Oi
    for(j in 1:m.uniq.len){
      prob[j]<-log(m.rle$lengths[j])+yOmega%*%mu.cur[j,]-t(mu.cur[j,])%*%Oi%*%mu.cur[j,]/2 #The log prob in DP
    }
    
    prob[last]<-log(alpha)+{log(det(Sigma.mui))+yOmega%*%Sigma.mui%*%t(yOmega)-nodes*(log(sigma2))}/2
    prob<-exp(prob-max(prob))
    newm<-which(rmultinom(1,1,prob)==1)

    del<-which(m.i==sort(unique(m.cur)))
    change<-F
    if(newm==last){
      change<-TRUE
      if(is.element(m.i,m.rle$values)){
        m.cur[i]<-m.rle$values[m.uniq.len]+1
        mu.cur<-rbind(mu.cur,rmvnorm(1,solve(Oi)%*%t(yOmega),solve(Oi)))
      }
      else mu.cur[del,]<-rmvnorm(1,solve(Oi)%*%t(yOmega),solve(Oi))
    }
    else{
      m.cur[i]<-m.rle$values[newm]
      if(!is.element(m.i,m.rle$values)) mu.cur<-matrix(mu.cur[-del,],ncol=nodes)
    }
    
    #relabel the clusters for the means
    uniquem<-sort(unique(m.cur))
    m.cur.update<-seq(1,length(uniquem))
    tmp.m.cur<-m.cur
    for (k in 1:length(uniquem))
    {
      loc<-which(tmp.m.cur==uniquem[k])
      m.cur[loc]<-m.cur.update[k]
    }
    
    # update O and G if the mean of subject i changes
    if(m.i!=m.cur[i]|change){
      c.i<-which(dn[i,]!=0)
      c.loc<-which(dn[,c.i]==1)
      c.mean<-matrix(rep(0,(nodes*length(c.loc))),nrow=length(c.loc))
      m.c.i<-unique(m.cur[c.loc])
      for (ii in 1: (length(m.c.i)))
      {
        m.loc<-which(m.cur[c.loc]==m.c.i[ii])
        c.mean[m.loc,]<-mu.cur[m.c.i[ii],]
      }
      O<-OclustArr[,,c.i]
      G<-GclustArr[,,c.i]
      m.c.data<-data[c.loc,]
      num<-length(c.loc)
      
      if (num==1) 
      {
        covXc<-t(m.c.data-c.mean)%*%((m.c.data-c.mean))
      } else 
      {
        covXc<-t(m.c.data-c.mean)%*%((m.c.data-c.mean))/(num-1)
      }
      Graphs<-GraphEst(O,covXc,G,num,nodes,nu0)
      eig<-eigen(Graphs$O)
      len<-length(which(eig$values<0))
      if (len==0)
      {
        OclustArr[,,c.i]<-Graphs$O
        GclustArr[,,c.i]<-Graphs$G
      }
    }
  }
  
  m.uniq.len<-length(unique(m.cur))
  m.uniq<-sort(unique(m.cur))
  for(i in 1:m.uniq.len){
    mclu<-which(m.cur==m.uniq[i])
    
    tmp2<-matrix(0,nrow=nodes,ncol=nodes)
    tmp3<-matrix(0,nrow=nodes,ncol=1)
    for(j in mclu){
      c.j<-which(dn[j,]!=0)
      O<-OclustArr[,,c.j]
      
      tmp2<-tmp2+O
      tmp3<-tmp3+O%*%data[j,]
    }
    tmp<-solve(tmp2+diag(1/sigma2,nrow=nodes))
    mu.cur[m.uniq[i],]<-rmvnorm(1,tmp%*%tmp3,tmp)
  }
  return(list(m.cur=m.cur,mu.cur=mu.cur,OclustArr=OclustArr,GclustArr=GclustArr))
}

GenPseudo<-function(data,num)
{
  n<-nrow(data)
  nodes<-ncol(data)
  n0<-floor((n/nodes))
  
  pseudoD<-matrix(NA,nrow=(nrow(data)),ncol=num)
  loc<-sample(seq(1:nodes),nodes,replace=FALSE)
  pseudoD<-data[,loc]
  for (j in 1:num)
  {  
    set.seed(12345*j)
    loc<-sample(1:nrow(pseudoD), (nrow(pseudoD)), replace = FALSE)
    pseudoD[,j]<-pseudoD[loc,j]
  }
  return(pseudoD[,1:num])
}

initalGraph<-function(O)
{
  nrows<-nrow(O)
  iGraph<-matrix(1,nrow=nrows,ncol=nrows)
  diag(iGraph)<-0
  Ovect<-rep(0,(nrows*(nrows-1)/2))
  k<-1
  for (i in 1:(nrows-1))
  {
    for (j in (i+1):nrows)
    {
      Ovect[k]<-O[i,j]
      k<-k+1
    }
  }
  Omedian<-quantile(abs(Ovect),prob=0.05)
  k<-1
  for (i in 1:(nrows-1))
  {
    for (j in (i+1):nrows)
    {
      if (abs(O[i,j])<= Omedian) iGraph[i,j]<-iGraph[j,i]<-0
    }
  }  
  return(iGraph)
}

# Tune nu0
TuneNu0<-function(data,nodes,sampsize,tune)
{
  set.seed(12345)
  nu0<-0.01
  low<-nodes-3
  upp<-ceiling(nodes+sqrt(nodes*(1-2/(nodes-1)))*10)

  pseudo<-GenPseudo(data,NumPseudo)
  dataAugment<-cbind(data,pseudo)
  
  nodesAug<-ncol(dataAugment)
  covAugX<-cov(dataAugment)
  O0<-solve(covAugX)
  G<-matrix(0,nrow=nrow(O0),ncol=nrow(O0))
  G0<-initalGraph(O0)

  stop<-0
  unchange<-0
  record<-1
  checkk<-0
  
  while (stop==0)
  {
    count<-0
    G<-G0
    O<-O0
    for (i in 1:tune)
    {
      Graphs<-GraphEst(O,covAugX,G,sampsize,nodesAug,nu0)
      O<-Graphs$O
      G<-Graphs$G
      WrongCon2<-sum(G[seq((nodesAug-NumPseudo+1):nodesAug),(1:(nodesAug-NumPseudo))])
      if (WrongCon2>(FPRate*(NumPseudo*nodes))) count<-count+1
    }
    rate<-count/tune
    if (rate<=record)
    {
      record<-rate
      recordNu0<-nu0
    } else
    {
      checkk<-1
      nucheck<-recordNu0
    }
    if (rate<=(cut+cut0))
    {
      stop<-1
    } else
    {
      if (checkk==1)
      {
        stop<-1
        nu0<-nucheck
      } else
      {
        nu0<-nu0*1.1
        
        if (nu0>nu1)
        {
          nu0<-nu1-0.1
          stop<-1
        }
      }
    }
  }
  
  unchange<-0
  count<-0
  while (stop==1)
  {
    Graphs<-GraphEst(O,covAugX,G,sampsize,nodesAug,nu0)
    O<-Graphs$O
    G<-Graphs$G
    
    NumEdg<-sum(G[1:nodes,1:nodes])/2
    if (NumEdg>upp) 
    {
      nu0<-nu0*1.1
      if (nu0>nu1) 
      {
        nu0<-nu1-0.1
        stop<-2
      }
      unchange<-0
      count<-0
    }
    else if (NumEdg<low) 
    {
      nu0<-nu0*0.9
      if (nu0<0.0004)
      {
        nu0<-0.0004
        stop<-2
      }
      unchange<-0
      count<-0
    }
    else if (low<=NumEdg & NumEdg<=upp) 
    {
      unchange<-unchange+1 
      WrongCon2<-sum(G[seq((nodesAug-NumPseudo+1):nodesAug),(1:(nodesAug-NumPseudo))])
      if (WrongCon2>(FPRate*(NumPseudo*nodes))) count<-count+1
      
      if (unchange==50) 
      {
        stop=2
      }
    }  
  }
  WrongCon2<-sum(G[seq((nodesAug-NumPseudo+1):nodesAug),(1:(nodesAug-NumPseudo))])
  return(list(O=O[1:nodes,1:nodes], G=G[1:nodes,1:nodes],nu0=nu0))
}

nloops<-500 # number of loops
itr<-200
alpha<-1.5 # value of alpha in DP

TrueC<-3 #the true number of clusters based on graphs
maxC<-5

sampsize<-nrow(dataAll)/nData
# pre-specified nu1
nu1<-4

# the following settings are used in the tunning process
tune<-100
tune0<-100
cut<-0.5
cut0<-0.1
FPRate<-0.1
NumPseudo<-nodes

# to help determine the selection of rows in the dn matrix
step<-floor(sampsize/TrueC)

zeta<-0.75
sigma2<-1
# this is the variance of a log-normal, the jumping function for sampling of zeta. 
v<-1.8
# this is for the jumping function for sampling O
jump<-0.1

GTPMeanC<-GTNMeanC<-rep(0,TrueC)

for (d in 1:nData)
{
  cat("Data ",d,"\n")
  
  ptm<-proc.time()
  LocStart<-1+(d-1)*sampsize
  LocEnd<-d*sampsize
  
  data<-dataAll[LocStart:LocEnd,]
  data<-apply(data,2,function(y) y-mean(y))
  data<-as.matrix(data)
  covX<-cov(data)

  GraphTuned<-TuneNu0(data,nodes,sampsize,tune)
  nu0<-GraphTuned$nu0
 
  Oinit<-GraphTuned$O
  Ginit<-GraphTuned$G
  # End of tune nu0
  
  # Clustering based on graphs
  PeBIC<-counts<-rep(0,(maxC-1))
  for (C in 2: maxC)
  {
    cat("Value of C ",C,"\n")
    mu<-matrix(rep(0,(sampsize*nodes)),nrow=sampsize)
    m.cur<-rep(1,sampsize)
    mu.cur<-t(as.matrix(rep(0,nodes)))
    dn<-matrix(rep(0,(C*sampsize)),nrow=sampsize)
    
    for (c in 1:(C-1))
    {
      if (c==1)
      {
        dn[,c]<-rbinom(nrow(data),1,(1/C))
        loc<-which(dn[,c]==0)
      } else
      {
        dn[loc,c]<-rbinom(length(loc),1,(1/C))
        loc<-which(rowSums(dn[,seq(1,c)])==0)
      }
    }
    dn[loc,C]<-1

    GClust<-Ginit[1:nodes,1:nodes]
    OClust<-Oinit[1:nodes,1:nodes]
    loc<-which(GClust==1, arr.ind = TRUE)
    
    # initial values for the clusters of G and O
    for (c in 1:C)
    {
      if (c==1)
      {
        GclustArr<-GClust
        OclustArr<-OClust
      } else
      {
        GClust1<-GClust
        GClust1[loc[c,1],loc[c,2]]<-GClust1[loc[c,2],loc[c,1]]<-0
        OClust1<-OClust
        OClust1[loc[c,1],loc[c,2]]<-OClust1[loc[c,2],loc[c,1]]<-0
        GclustArr<-abind(GclustArr,GClust1,along=3)
        OclustArr<-abind(OclustArr,OClust1,along=3)
      }
    }
    
    Tp<-Fp<-rep(0,(nloops/2))
    dic<-0
    exclude<-0
    freqPTP<-freqPTN<-matrix(NA,nrow=TrueC,ncol=C)

    for (i in 1:nloops)
    {
      cat("i ",i,"\n")
      for (sub in 1:sampsize)
      {
        mu[sub,]<-mu.cur[m.cur[sub],]
      }
      # update sigma^2 in the base
      sigma2<-sampsig2(sampsize,mu)
      
      # update mu's
      muClust<-DPmu(data,nodes,alpha,m.cur,mu.cur,sigma2,dn,OclustArr,GclustArr)

      m.cur<-muClust$m.cur
      mu.cur<-muClust$mu.cur
      OclustArr<-muClust$OclustArr
      GclustArr<-muClust$GclustArr
      
      # update graphs in each cluster
      for (c in 1:C)
      {
        c.loc<-which(dn[,c]==1)
        if (length(c.loc>0))
        {
          c.mean<-matrix(rep(0,(nodes*length(c.loc))),nrow=length(c.loc))
          m.c.i<-unique(m.cur[c.loc])
          for (ii in 1: (length(m.c.i)))
          {
            m.loc<-which(m.cur[c.loc]==m.c.i[ii])
            c.mean[m.loc,]<-mu.cur[m.c.i[ii],]
          }
          
          O<-OclustArr[,,c]
          G<-GclustArr[,,c]
          m.c.data<-data[c.loc,]
          num<-length(c.loc)
          
          if (num==1) 
          {
            covXc<-t(m.c.data-c.mean)%*%((m.c.data-c.mean))
          } else 
          {
            covXc<-t(m.c.data-c.mean)%*%((m.c.data-c.mean))/(num-1)
          }
          Graphs<-GraphEst(O,covXc,G,num,nodes,nu0)
          eig<-eigen(Graphs$O)
          len<-length(which(eig$values<0))
          if (len==0)
          {
            OclustArr[,,c]<-Graphs$O
            GclustArr[,,c]<-Graphs$G
          }
          
          if (i>(nloops/2))
          {
            tmpO<-OclustArr[,,c]*(GclustArr[,,c]+diag(1,nrow=nodes))
            loglh<-lhCal(tmpO,covXc,num)
            if (loglh==-999999)
            {
              exclude=exclude+1
            } else
            {
              dic<-dic+loglh
            }
          }
        }
      }

      # update cluster assignment and pi (graphs)
      ProbAndSigClust<-sampdn(sampsize,data,zeta,C,dn,OclustArr,m.cur,mu.cur)
      dn<-ProbAndSigClust$dn
      
      # MH to update zeta with jumping dist lognormal
      sum.log.P<-sum(log(ProbAndSigClust$prob))
      zeta.l<-log(zeta)
      can<-rlnorm(1,zeta.l,v)
      ratio<-min(exp(dzeta(can,sum.log.P,C)-dzeta(zeta,sum.log.P,C))/ (dlnorm(can,zeta.l,v)/dlnorm(zeta,zeta.l,v)),1)
      if(is.nan(ratio)){
        ratio<-0
      }
      if(runif(1)<ratio){
        zeta<-can
      }
      
      if (i==(nloops/2+1))
      {
        sumGclustArr<-GclustArr
        sumOclustArr<-OclustArr
        
        m.means<-matrix(rep(0,(nodes*nrow(data))),nrow=nrow(data))
        m.unique<-unique(m.cur)
        for (ii in 1: (length(m.unique)))
        {
          m.loc<-which(m.cur==m.unique[ii])
          m.means[m.loc,]<-mu.cur[m.unique[ii],]
        }
        summu<-m.means
        sumdn<-dn
      }
      if (i>(nloops/2+1))
      {
        for (cc in 1:C)
        {
          sumGclustArr[,,cc]<-sumGclustArr[,,cc]+GclustArr[,,cc]
          sumOclustArr[,,cc]<-sumOclustArr[,,cc]+OclustArr[,,cc]
        }
      }
        
      if (i>(nloops/2))
      {
        for (c in 1:TrueC)
        {
          locT<-which(GTrue[,,c]==1, arr.ind = TRUE)
          locF<-which(GTrue[,,c]==0, arr.ind = TRUE)
          
          for (cc in 1:C)
          {
            Diff<-GTrue[,,c]-GclustArr[,,cc]
            sumTP<-sumFP<-sumTN<-sumFN<-0
            if (nrow(locT)>=1)
            {
              for (jj in 1:nrow(locT))
              {
                if (Diff[locT[jj,1],locT[jj,2]]==0)
                {
                  sumTP<-sumTP+1
                } else
                {
                  sumFN<-sumFN+1
                }
              }
            }
            
            if (nrow(locF)>=1)
            {
              for (jj in 1:nrow(locF))
              {
                if (Diff[locF[jj,1],locF[jj,2]]==0)
                {
                  sumTN<-sumTN+1
                } else
                {
                  sumFP<-sumFP+1
                }
              }
            }
            freqPTP[c,cc]<-sumTP/(sumTP+sumFP)
            if (sumTP==0) freqPTP[c,cc]<-0
            freqPTN[c,cc]<-(sumTN-nodes)/(sumTN-nodes+sumFN)
            if (sumTN==0) freqPTN[c,cc]<-0
          }
        }
      }
      if (i>(nloops/2))
      {
        m.unique<-unique(m.cur)
        for (ii in 1: (length(m.unique)))
        {
          m.loc<-which(m.cur==m.unique[ii])
          m.means[m.loc,]<-mu.cur[m.unique[ii],]
        }
        summu<-m.means+summu
        sumdn<-dn+sumdn
      }
      if (i>(nloops/2))
      {
        # calculate sensitivity and specificity for each C  
        freq<-matrix(NA,nrow=TrueC,ncol=C)
        freqTp<-freqFp<-rep(0,TrueC)
        
        for (c in 1:TrueC)
        {
          for (cc in 1:C)
          {
            freq[c,cc]<-sum(dn[(step*(c-1)+1):(step*c),cc])
          }
          freqTp[c]<-max(freq[c,])
          loc<-which(freq[c,]==freqTp[c])
          freqFp[c]<-sum(dn[-((step*(c-1)+1):(step*c)),loc])
        }
        Tp[(i-(nloops/2))]<-sum(freqTp)/nrow(dn)
        Fp[(i-(nloops/2))]=sum(freqFp)/nrow(dn)
      }
    }# end for (in in 1:nloops)
    
    # summarize sensitivities and specificities for each C. 
    tmp1<-quantile(Tp,probs = c(0,0.025,0.5,0.975,1))
    mTp<-mean(Tp)
    tmp1<-c(tmp1,mTp)
    if (C==2)
    {
      TpC<-as.matrix(t(tmp1))
    } else
      TpC<-rbind(TpC,tmp1)

    tmp1<-quantile(Fp,probs = c(0,0.025,0.5,0.975,1))
    mFp<-mean(Fp)
    tmp1<-c(tmp1,mFp)
    if (C==2)
    {
      FpC<-as.matrix(t(tmp1))
    } else
      FpC<-rbind(FpC,tmp1)

    for (c in 1:C)
    {
      sumGclustArr[,,c]<-sumGclustArr[,,c]/(nloops/2)
      sumOclustArr[,,c]<-sumOclustArr[,,c]/(nloops/2)
    }

    for (c in 1:TrueC)
    {
       locT<-which(GTrue[,,c]==1, arr.ind = TRUE)
       locF<-which(GTrue[,,c]==0, arr.ind = TRUE)
       
       for (cc in 1:C)
       {
         Diff<-GTrue[,,c]-round(sumGclustArr[,,cc],digits=0)
         sumTP<-sumFP<-sumTN<-sumFN<-0
         if (nrow(locT)>=1)
         {
           for (jj in 1:nrow(locT))
           {
             if (Diff[locT[jj,1],locT[jj,2]]==0)
             {
               sumTP<-sumTP+1
             } else
             {
               sumFN<-sumFN+1
             }
           }
         }
            
         if (nrow(locF)>=1)
         {
           for (jj in 1:nrow(locF))
           {
             if (Diff[locF[jj,1],locF[jj,2]]==0)
             {
               sumTN<-sumTN+1
             } else
             {
               sumFP<-sumFP+1
             }
           }
         }
         freqPTP[c,cc]<-sumTP/(sumTP+sumFP)
         if (sumTP==0) freqPTP[c,cc]<-0
         freqPTN[c,cc]<-(sumTN-nodes)/(sumTN-nodes+sumFN)
         if (sumTN==0) freqPTN[c,cc]<-0
       }
    }

    # the following is to finish PeBIC calculation
    tmpd<-sumdn/(nloops/2)
    max_cols <- apply(tmpd, 1, which.max)
    
    # Combine with row indices
    locations <- cbind(row = 1:nrow(tmpd), col = max_cols)
    tmpdn<-matrix(0,nrow=nrow(sumdn),ncol=C)
    tmpdn[locations]<-1
    
    tmpmu<-summu/(nloops/2)
    
    PeBIC[(C-1)]<- 0 
    
    for (c in 1:C)
    {
      c.loc<-which(tmpdn[,c]==1)
      c.mean<-tmpmu[c.loc,]
      c.data<-data[c.loc,]
      if (c==1 && C==2)
      {
        if ((is.vector(c.data)==TRUE))
        {
          tmpMat<-as.matrix(t(c.data))
          Cdata<-cbind(tmpMat,rep(c,nrow(tmpMat)))
        } else
          Cdata<-cbind(c.data,rep(c,nrow(c.data)))
      } else
      {
        if ((is.vector(c.data)==TRUE))
        {
          tmpMat<-as.matrix(t(c.data))
          Cdata<-rbind(Cdata,cbind(tmpMat,rep(c,nrow(tmpMat))))
        } else
          Cdata<-rbind(Cdata,cbind(c.data,rep(c,nrow(c.data))))
      }
      
      num<-length(c.loc)
      
      if (num==1) 
      {
        covXc<-(c.data-c.mean)%*%t((c.data-c.mean))
      } else 
      {
        covXc<-t(c.data-c.mean)%*%((c.data-c.mean))/(num-1)
      }
      tmpO<-sumOclustArr[,,c]*(sumGclustArr[,,c]+diag(1,nrow=nodes))
      tmpG<-round(sumGclustArr[,,c],digits=0)
      
      if (num >10)
      {
        loglh<-lhCal(tmpO,covXc,num)        
      } else {
        loglh<-0
        counts[C-1]<-counts[C-1]+1
      }
      
      # calculate PeBIC; take gamma=1
      ga=1
      df<-sum(tmpG[upper.tri(tmpG)] != 0)
      penalty<-4 * ga * log(nodes) * df
      PeBIC[(C-1)]<-PeBIC[(C-1)]+2*loglh -(C-counts[C-1])*(log(sampsize)*df+penalty)
    }
  } # end "for" for C
  
  maxDICTmp<-max(DICTmp)
  loc0Tmp<-which(DICTmp==maxDICTmp)
  locTmp<-loc0Tmp-counts[loc0Tmp]
  
  dataSelect<-Cdata[(sampsize*(loc0Tmp-1)+1):(sampsize*loc0Tmp),]
  clusters<-unique(dataSelect[,(nodes+1)])
  for (cc in 1:length(clusters))
  {
    locSub<-which(dataSelect[,(nodes+1)]==clusters[cc])
    dataSelectSub<-dataSelect[locSub,1:nodes]
    if ((is.vector(dataSelectSub)==TRUE))
    {
      num<-1
    } else
      num<-nrow(dataSelectSub)
    if (cc==1)
    {
      if (num<(nodes*2+5) || (length(num)==0))
      {
        # if sample size is small in a cluster, then do not infer a graph.
        GclustFinal<-array(matrix(0,nrow=nodes,ncol=nodes),dim=c(nodes,nodes,1))
      } else {
      tmpp<-TuneNu0(dataSelectSub,nodes,num,tune0)
      GclustFinal<-array(tmpp$G,dim=c(nodes,nodes,1))
      }
    } else {
    if (num<(nodes*2+5) || (length(num)==0))
    {
      # if sample size is small in a cluster, then do not infer a graph.
      GclustFinal<-abind(GclustFinal,matrix(0,nrow=nodes,ncol=nodes),along=3)
    } else {
      tmpp<-TuneNu0(dataSelectSub,nodes,num,tune0)
      GclustFinal<-abind(GclustFinal,tmpp$G,along=3)
      }
    }
  }
  
  # Assess the quality of estimated graph
  freqGTP<-freqGTN<-freqGACC<-matrix(NA,nrow=TrueC,ncol=(length(clusters)))
  for (c in 1:TrueC)
  {
    locT<-which(GTrue[,,c]==1, arr.ind = TRUE)
    locF<-which(GTrue[,,c]==0, arr.ind = TRUE)
  
    for (cc in 1:(length(clusters)))
    {
      Diff<-GTrue[,,c]-GclustFinal[,,cc]
      sumTP<-sumFP<-sumTN<-sumFN<-0
      if (nrow(locT)>=1)
      {
        for (jj in 1:nrow(locT))
        {
          if (Diff[locT[jj,1],locT[jj,2]]==0)
          {
            sumTP<-sumTP+1
          } else
          {
            sumFN<-sumFN+1
          }
        }
      }
  
      if (nrow(locF)>=1)
      {
        for (jj in 1:nrow(locF))
        {
          if (Diff[locF[jj,1],locF[jj,2]]==0)
          {
            sumTN<-sumTN+1
          } else
          {
            sumFP<-sumFP+1
          }
        }
      }
      freqGTP[c,cc]<-sumTP/(sumTP+sumFN)
      if (sumTP==0) freqGTP[c,cc]<-0
      freqGTN[c,cc]<-(sumTN-nodes)/(sumTN-nodes+sumFP)
      if (sumTN==0) freqGTN[c,cc]<-0
    }
    GTPMeanC[c]<-max(freqGTP[c,])
    loc<-which(freqGTP[c,]==GTPMeanC[c])
    if (length(loc)>1)
    {
      GTNMeanC[c]<-max(freqGTN[c,loc])
    } else
    {
      GTNMeanC[c]<-freqGTN[c,loc]
    }
  }
  
  if (d==1)
  {
    PeBICData<-PeBIC
    
    # Cluster quality (number of clusters, true positive and false positve rates)
    NumClust<-locTmp+1
    TpData<-TpC[loc0Tmp,]
    FpData<-FpC[loc0Tmp,]
    # Graph quality (true positive and false positive rates)
    GTpData<-GTPMeanC
    GFpData<-1-GTNMeanC
  } else
  {
    PeBICData<-rbind(PeBICData,PeBIC)
    NumClust<-c(NumClustTmp,(locTmp+1))
    TpData<-rbind(TpDataTmp,TpC[loc0Tmp,])
    FpData<-rbind(FpDataTmp,FpC[loc0Tmp,])
    GTpData<-rbind(GTpData,GTPMeanC)
    GFpData<-rbind(GFpData,(1-GTNMeanC))
  }
} #end data
write.table(PeBICData, file=paste("TwoClust",n1,"PeBIC",SLURM_ID,".txt", sep=""),col.names=FALSE,row.names=FALSE)
write.table(NumClust, file=paste("TwoClust",n1,"NumClusters",SLURM_ID,".txt", sep=""),col.names=FALSE,row.names=FALSE)
write.table(TpData, file=paste("TwoClust",n1,"Tps",SLURM_ID,".txt",sep=""),col.names=FALSE,row.names=FALSE)
write.table(FpData, file=paste("TwoClust",n1,"Fps",SLURM_ID,".txt",sep=""),col.names=FALSE,row.names=FALSE)
write.table(GTpData, file=paste(clus,"Clust",n1,"GTpData",SLURM_ID,".txt",sep=""),col.names=FALSE,row.names=FALSE)
write.table(GFpData, file=paste(clus,"Clust",n1,"GFpData",SLURM_ID,".txt",sep=""),col.names=FALSE,row.names=FALSE)
