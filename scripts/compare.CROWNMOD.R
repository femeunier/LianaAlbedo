rm(list = ls())

library(purrr)
library(abind)

Npatches = 49
Nwl = 2100
N = 2

temp <- 
  map(1:Npatches,function(i){
  file <- paste0("/home/carya/output/PEcAn_99000000001/out/ICANRAD1/rad_profile",i,".dat")
  temp <- scan(file,quiet = TRUE)
  Ncohort <- length(temp)/2100
  
  mat <- array(temp,dim = c(2100,Ncohort))
  return(t(mat[,c(Ncohort,1)]))})

ICANRAD1 <- aperm(abind(temp,along = 3),c(3,2,1))



temp <- 
  map(1:Npatches,function(i){
    file <- paste0("/home/carya/output/PEcAn_99000000001/out/ICANRAD2/rad_profile",i,".dat")
    temp <- scan(file,quiet = TRUE)
    Ncohort <- length(temp)/2100
    
    mat <- array(temp,dim = c(2100,Ncohort))
    return(t(mat[,c(Ncohort,1)]))})

ICANRAD2 <- aperm(abind(temp,along = 3),c(3,2,1))

temp <- 
  map(1:Npatches,function(i){
    file <- paste0("/home/carya/output/PEcAn_99000000001/out/CROWNMOD1/rad_profile",i,".dat")
    temp <- scan(file,quiet = TRUE)
    Ncohort <- length(temp)/2100
    
    mat <- array(temp,dim = c(2100,Ncohort))
    return(t(mat[,c(Ncohort,1)]))})

CROWNMOD1 <- aperm(abind(temp,along = 3),c(3,2,1))

plot(as.vector(ICANRAD2[,,1]),as.vector(CROWNMOD1[,,1]))

ipatch = 1
plot(ICANRAD2[ipatch,,2],ylim = c(0,0.5),type='l',lty=2)
lines(CROWNMOD1[ipatch,,2],col = 'red',lty=2)

# keep = temp2
