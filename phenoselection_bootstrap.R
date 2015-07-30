# Phenotypic selection gradients from Lande and Arnold 1985.
# Includes bootstrap evaluation of significance

gradientslinear<- function(data, indices) {
  d1 <- data[indices,] # allows boot to select sample 
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  s<- cov(d1)[2:3,1]
  B<-Pinv%*%s
  return(B) 
} 

gradientsnonlinear<- function(data, indices) {
  d1 <- data[indices,] # allows boot to select sample 
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  
  new<-d1[,2:3]^2
  z12<-d1[,2]*d1[,3]
  new<-data.frame(new,z12)
  c<-matrix(ncol=2,nrow=2)
  rawcov<-cov(d1$w,new)
  c<-diag(cov(d1$w,new)[1:2])
  c[1,2]<-rawcov[3]
  c[2,1]<-rawcov[3]
  gamma = Pinv%*% c %*%Pinv 
  return(c) 
} 

selectionlinear<- function(data, indices) {
  d1 <- data[indices,] # allows boot to select sample 
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  s<- cov(d1)[2:3,1]
  B<-Pinv%*%s
  return(s) 
} 


selectionnonlinear<- function(data, indices) {
  d1 <- data[indices,] # allows boot to select sample 
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  
  new<-d1[,2:3]^2
  z12<-d1[,2]*d1[,3]
  new<-data.frame(new,z12)
  c<-matrix(ncol=2,nrow=2)
  rawcov<-cov(d1$w,new)
  c<-diag(cov(d1$w,new)[1:2])
  c[1,2]<-rawcov[3]
  c[2,1]<-rawcov[3]
  gamma = Pinv%*% c %*%Pinv 
  return(c) 
} 