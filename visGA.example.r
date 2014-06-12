#Maximize a mixture of multivariate normal distributions
 
library(mvtnorm)
 
mnMix<-function(args){
   
  mean.vec.d1<-rep(0.3,5)
  std.vec.d1<-diag(rep(1,5))
   
  mean.vec.d2<-rep(1,5)
  std.vec.d2<-diag(rep(1.5,5))
   
  mean.vec.d3<-c(1, 5, 2, 1, 0)
  std.vec.d3<-diag(rep(0.5, 5))
   
  if (args[1]<0){
    y<-dmvnorm(args, mean.vec.d1, std.vec.d1)+dmvnorm(args, mean.vec.d2, std.vec.d2)+dmvnorm(args, mean.vec.d3, std.vec.d3)
  }
  else {
    y<-dmvnorm(args, mean.vec.d1, std.vec.d1)*dmvnorm(args, mean.vec.d2, std.vec.d2)*dmvnorm(args, mean.vec.d3, std.vec.d3)
  }
   
  y
   
}
 
visGAPath(mnMix, 5, 8, parallel=TRUE, file.path='B:\\RFiles\\RInOut\\pro')
 
Colville<-function(args){
  x0<-args[1]
  x1<-args[2]
  x2<-args[3]
  x3<-args[4]
  ret<-100*(x0-x1^2)^2+(1-x0)^2+90*(x3-x2^2)^2+(1-x2)^2+10.1*((x1-1)^2+(x3-1)^2)+19.8*(x1-1)*(x3-1)
  ret
}
 
domains<-matrix(rep(c(-10,10),4), nrow=4, byrow=TRUE)
 
visGAPath(Colville, 4, 8, parallel=TRUE, file.path='B:\\RFiles\\RInOut\\pro', domains=domains, max=FALSE)