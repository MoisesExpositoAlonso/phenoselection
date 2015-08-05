######################################################################################################
### phenoselection, R script implementing phenotypic selection gradients from Lande & Arnold 1983  ###
######################################################################################################
#######################           by Moises Exposito-Alonso        ###################################

# Phenotypic selection gradients from Lande and Arnold 1985.
# Includes bootstrap evaluation of significance

#@@@@@@@ (1) INPUT DATA TO BE PROVIDED  @@@@@@@#

Var1="vector of numerical values of phenotype 1"
Var2="vector of numerical values of phenotype 2"  # notice that the results will be reported in the same order as these two variables
Fitness= "vector of numerical values of fitness"
Gmatrix=matrix(c(heritabilty1,correlation,correlation, heritability2,ncol=2)

#@@@@@@@ (2) LITTLE COMMAND TO DO THE ANALYSES, but first run (3)   @@@@@@@#

PHENOSELECTION(Var1=Var1,Var2=Var2,Fitness=Fitness,Gmatrix=Gmatrix)

## 



#@@@@@@@ (3) RUN ALL DOWN HERE BEFORE (2)  ##### 
PHENOSELECTION<-function(Var1,Var2,Fitness,Gmatrix=NULL){
library(boot)
if (Gmatrix==NULL){
print "Heritabilities not provided, only selection analyses carried out" }

d1<-preparedata(Fitness,Var1,Var2)



################################# start inside functions #############################
#### prepare data

preparedata<-function(Fitness,Var1,Var2){
data<-data.frame(Fitness=as.numeric(Fitness),
                 Var1=as.numeric(Var1),
                 Var2=as.numeric(Var2) )
data<-na.omit(data)
d1<-(data.frame(w=data$Fitness/mean(data$Fitness),
                       Var1=scale(data$Var1),
                       Var2=scale(data$Var2)))

d1$w
d1$Var1<-scale(d1$Var1)
d1$Var2<-scale(d1$Var2)
return(d1)
}

#### functions from Lande & Arnold 1983 with bootstrap significance ####


gradientlinear<- function(data, indices) {
  d1 <- data[indices,] 
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  s<- cov(d1)[2:3,1]
  B<-Pinv%*%s
  return(B) 
} 

gradientquadratic<- function(data, indices) {
  d1 <- data[indices,] 
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
  return(gamma) 
} 

coeflinear<- function(data, indices) {
  d1 <- data[indices,]
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  s<- cov(d1)[2:3,1]
  B<-Pinv%*%s
  return(s) 
} 

coefquadratic<- function(data, indices) {
  d1 <- data[indices,]
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

responselinear<- function(data, indices) {
  d1 <- data[indices,] # allows boot to select sample 
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  s<- cov(d1)[2:3,1]
  B<-Pinv%*%s
  
  
  deltaZ<-Gmatrix%*%B
  deltaZ
  newZ<-deltaZ 
  meanchange<-newZ
  
  return(meanchange)
} 

responsequadratic<- function(data, indices) {
  d1 <- data[indices,] 
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  
  d1 <- data[indices,] 
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  s<- cov(d1)[2:3,1]
  B<-Pinv%*%s
  
  
  new<-d1[,2:3]^2
  z12<-d1[,2]*d1[,3]
  new<-data.frame(new,z12)
  c<-matrix(ncol=2,nrow=2)
  rawcov<-cov(d1$w,new)
  c<-diag(cov(d1$w,new)[1:2])
  c[1,2]<-rawcov[3]
  c[2,1]<-rawcov[3]
  gamma = Pinv%*% c %*%Pinv 
  
  
  deltaP<- (P*gamma*P) - s%*%t(s)  


deltaG<-Gmatrix%*%(gamma-(B%*%t(B)))%*%Gmatrix 
  deltaG+Gmatrix 
  
  return(deltaP)
  
} 

responsequadratic_gmatrix<- function(data, indices) {
  
  d1 <- data[indices,]
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  s<- cov(d1)[2:3,1]
  B<-Pinv%*%s
  
  
  new<-d1[,2:3]^2
  z12<-d1[,2]*d1[,3]
  new<-data.frame(new,z12)
  c<-matrix(ncol=2,nrow=2)
  rawcov<-cov(d1$w,new)
  c<-diag(cov(d1$w,new)[1:2])
  c[1,2]<-rawcov[3]
  c[2,1]<-rawcov[3]
  gamma = Pinv%*% c %*%Pinv 
  return(gamma) 
  
  
  deltaP<- (P*gamma*P) - s%*%t(s)  # why variance is not 1??
  deltaP
  newP<-deltaP+P
  covariancechange<-deltaP

  deltaG<-Gmatrix%*%(gamma-(B%*%t(B)))%*%Gmatrix 
  deltaG+Gmatrix # the change in heritability and genetic correlation!
  
  return(deltaG)
  
} 


#### function to produce nice bootstrap output #####


extractbootstrap<-function(bootstrapresults){
treatbootstrap<-function(x){

  five<-quantile(x,p=c(0.05,0.95) )
  one<-quantile(x,p=c(0.01,0.99) )
  zeroone<-quantile(x,p=c(0.001,0.999) )
  tempsign<-(five[1]/five[2]) / abs(five[1]/five[2])
  tempsign<-c(tempsign,(one[1]/one[2]) / abs(one[1]/one[2]))
  tempsign<-c(tempsign,(zeroone[1]/zeroone[2]) / abs(zeroone[1]/zeroone[2]))
  tempsign[tempsign==1]<-"*"
  tempsign[tempsign==-1]<-""
  tempsign[tempsign==NA]<-""
  tempsign<-paste(as.character(tempsign)[1],as.character(tempsign)[2],as.character(tempsign)[3],sep="")
  
  as.character(tempsign)
  se<-round(sd(x),digits = 3)
  media<-round(mean(x),digits = 3)
  pasted<-paste(media," (",se, ")",tempsign,sep="")
  return(pasted)
}
extracted<-apply(bootstrapresults$t,2,treatbootstrap)
return(extracted)
}

extractbootstrap_numeric<-function(bootstrapresults){
  treatbootstrap<-function(x){
    media<-round(mean(x),digits = 3)
    return(media)
  }
  extracted<-apply(bootstrapresults$t,2,treatbootstrap)
  return(extracted)
}

################################# end inside functions #############################

#@@@@@@ start LITTLE BIT THAT ACTUALLY DO ANALYSES  @@@@@@@#

result_gradient_linear<- boot(data=d1, statistic=gradientlinear, R=1000)
result_gradient_quadratic<- boot(data=d1, statistic=gradientquadratic, R=1000)
result_coefficient_linear<- boot(data=d1, statistic=coeflinear, R=1000)
result_coefficient_quadratic<- boot(data=d1, statistic=coefquadratic, R=1000)

resa<-extractbootstrap(result_gradient_linear)
resb<-extractbootstrap(result_gradient_quadratic)
resc<-extractbootstrap(result_coefficient_linear)
resd<-extractbootstrap(result_coefficient_quadratic)



if (Gmatrix!=NULL){
result_response_linear<- boot(data=d1, statistic=responselinear, R=1000)
result_response_quadratic<- boot(data=d1, statistic=responsequadratic, R=1000)
result_response_quadratic_gmatrix<- boot(data=d1, statistic=responsequadratic_gmatrix, R=1000)

rese<-extractbootstrap(result_response_linear)
resf<-extractbootstrap(result_response_quadratic)
resg<-extractbootstrap(result_response_quadratic_gmatrix)

analysislist<-list(gradient_linear=resa,gradient_quadratic=resb,
coefficient_linear=resc,coefficient_quadratic=resd,
response_linear=rese,
response_quadratic_Vpheno=resf,
response_quadratic_Vaddit=resg)

} else{
analysislist<-list(gradient_linear=resa,gradient_quadratic=resb,
coefficient_linear=resc,coefficient_quadratic=resd)

}

return(analysislist)
#@@@@@@ end LITTLE BIT THAT ACTUALLY DO ANALYSES  @@@@@@@#

} # end phenoselection
