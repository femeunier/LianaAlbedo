rm(list = ls())

library(PEcAn.ED2)
library(dplyr)
library(tidyr)
library(redr)

Nensemble = 1

WL = seq(400, 2500)

# Filters
par <- WL>=400 & WL <= 700
green <- WL>=520 & WL <= 600
red <- WL>=550 & WL <= 700
VNIR <- WL>=750 & WL <= 900
IR <- WL>=700 & WL <= 2500 
NIR <- WL>=800 & WL <= 1400
SIR <- WL>=1500 & WL <= 2500
SW <- WL>=400 & WL <= 2500

# read data
All.posteriors <- readRDS("/home/carya/data/RTM/All_posteriors.RDS")

All.posteriors %>% group_by(Param,pft) %>% summarise(mean(value))

# Mean
ensemble <- All.posteriors %>% group_by(Param,pft) %>% summarise(value = mean(value))  %>% mutate(run = 1)  %>% pivot_wider(names_from = Param,
                                                                                                                             values_from = value)
# ensemble <- All.posteriors %>% filter(ref == "Foster") %>% group_by(Param,pft) %>% summarise(value = mean(value))  %>% mutate(run = 1)  %>% pivot_wider(names_from = Param,
#                                                                                                                                                         values_from = value)

# Samples
# ensemble <- All.posteriors %>% group_by(Param,pft) %>% slice_sample(n = Nensemble) %>% mutate(run = 1:Nensemble) %>% dplyr::select(-c("ref")) %>% pivot_wider(names_from = Param,
#                                                                                                                                                               values_from = value)



data.BCI <- read_css(file = "/home/carya/data/BCI/BCI_PFTs.lat9.000lon-79.000.css")
data.mod <-
  data.BCI %>% mutate(pft = case_when(pft %in% c(2, 3, 4) ~ 3,
                                      pft == 17 ~ 17))

sim.patch <- list()
Npatches <- length(unique(data.mod$patch))

href <- 61.7;  b1Ht <- 0.035; b2Ht <- 0.69; hmax = 35 # Tree
hrefL <- 61.7;  b1HtL <- 0.11; b2HtL <- 2.5 # Liana

dbh.thrshld <- 2.5
for (ipatch in seq(1,Npatches)){
  
  site_dat <- data.mod %>% filter(patch == ipatch)
  
  dbh <- site_dat[["dbh"]]
  pft <-  as.factor(site_dat[["pft"]])
  nplant <- site_dat[["n"]]
  patch <- site_dat[["patch"]]
  
  is_liana <- (pft == 17)
  npft <- length(levels(pft))
  levels(pft) <- seq(npft,1)
  
  h <- pmin(hmax,(href*(1 -exp(-b1Ht*(dbh**b2Ht)))))
  
  h[is_liana] <- pmin(hmax,(hrefL*(1 -exp(-b1HtL*(dbh[is_liana]**b2HtL)))))
  
  tree.maxH <- max(h[!is_liana])
  h[dbh>dbh.thrshld & is_liana] <- tree.maxH + 0.5
  h[is_liana & dbh <= dbh.thrshld] <- pmin(h[is_liana & dbh <= dbh.thrshld],tree.maxH)
  
  temp <- sort(h, decreasing = FALSE, index.return = TRUE)
  h <- temp$x
  dbh <- dbh[temp$ix]
  pft <- pft[temp$ix]
  nplant <- nplant[temp$ix]
  patch <- patch[temp$ix]
  is_liana <- is_liana[temp$ix]
  
  
  sim.patch[[ipatch]] <- data.frame(DBH = dbh, H = h, PFT = as.numeric(as.vector(pft)), N = nplant, PA = patch,is.liana = is_liana)
}

# ggplot(data = sim.patch[[12]]) +
#   geom_point(aes(x = 1, y = H,color = as.factor(PFT))) +
#   theme_bw()

OP <- data.frame()
leaf_reflect_all <- list()

for (i in seq(1,Nensemble)){
  
  print(paste0("Ensemble: ",i/Nensemble))
  
  cparams <- ensemble %>% filter(run == i) 
  SLA <- 1/(10*cparams$Cm);b1Bl <- cparams$b1Bl; b2Bl <- cparams$b2Bl ; N <- cparams$Nlayers ; Cab <- cparams$Cab ; Car <- cparams$Car ; Cm <- cparams$Cm ; Cw <- cparams$Cw
  orient_factor <- cparams$orient.factor ; clumping_factor <- cparams$clumping.factor
  
  leaf_spectra <- list()
  
  for (i in seq(1,npft)){
    leaf_spectra[[i]] <- list()
    temp <- PEcAnRTM::prospect(c(N[i],Cab[i],Car[i],Cw[i],Cm[i]), version = "5")
    leaf_spectra[[i]][["reflectance"]] <- as.vector(temp[,1])
    leaf_spectra[[i]][["transmittance"]] <- as.vector(temp[,2])}
  
  leaf_reflect <- Reduce(cbind,Map(function(x) x[[1]], leaf_spectra))
  leaf_reflect_all[[i]] <- leaf_reflect
  
  LAI_liana <- LAI <- delta_PAR <- delta_NIR <- N_liana <- delta_PARgnd <- c() 
  
  for (ipatch in seq(1,Npatches)){
    
     print(paste0("patch: ",ipatch/Npatches))
     
     cpatch <- sim.patch[[ipatch]] %>% mutate(lai = SLA[PFT]*N*b1Bl[PFT]*(DBH**b2Bl[PFT])) %>% filter( (is.liana & DBH >= 1.5) | (!is.liana & DBH > 8))
     
     LT <-
       edr_r_v2(
         pft = cpatch$PFT, lai = cpatch$lai, wai = cpatch$lai * 0, cai = cpatch$lai ** 0,
         N = N[c(1,2)], Cab = Cab[c(1,2)], Car = Car[c(1,2)], Cm = Cm[c(1,2)], Cw = Cw[c(1,2)],
         orient_factor = orient_factor[c(1,2)],clumping_factor = clumping_factor[c(1,2)],
         soil_moisture = c(0.5), direct_sky_frac = c(0.5),czen = c(0.5),
         wavelengths = WL)

     TT <-
       edr_r_v2(
         pft = cpatch$PFT,lai = cpatch$lai,wai = cpatch$lai * 0,cai = cpatch$lai ** 0,
         N = N[c(2,2)],Cab = Cab[c(2,2)],Car = Car[c(2,2)],Cm = Cm[c(2,2)],Cw = Cw[c(2,2)],
         orient_factor = orient_factor[c(2,2)],clumping_factor = clumping_factor[c(2,2)],
         soil_moisture = c(0.5),direct_sky_frac = c(0.5),czen = c(0.5),
         wavelengths = WL)

     delta_PAR[ipatch] <- mean(LT$albedo[par] - TT$albedo[par])
     delta_NIR[ipatch] <- mean(LT$albedo[NIR] - TT$albedo[NIR])
     delta_PARgnd[ipatch] <- mean((LT$light_level[par,1] - TT$light_level[par,1])/TT$light_level[par,1])
     
     LAI[ipatch] <- sum(cpatch$lai)
     LAI_liana[ipatch] <- sum(cpatch$lai[cpatch$is.liana])
     N_liana[ipatch] <- sum(cpatch$N[cpatch$is.liana])
  }
  
  OP <- bind_rows(list(OP,
                 data.frame(delta_PAR,delta_NIR,delta_PARgnd,LAI,LAI_liana,N_liana, run = i,patch = 1:length(N_liana))))
}

saveRDS(object = OP,file = "/home/carya/R/edr-da/data/OP_BCI2.RDS")
saveRDS(object = ensemble,file = "/home/carya/R/edr-da/data/ensemble.RDS")
saveRDS(object = leaf_reflect_all,file = "/home/carya/R/edr-da/data/leaf_reflect_all.RDS")
