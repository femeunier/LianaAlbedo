rm(list = ls())

# Assemble all posteriors from step1

library(tidyr)
library(BayesianTools)
outdir <- "~/output/"

load(file = "~/data/RTM/Inverse_leaf_spectrum.Rdata")

# Create a matrix of PDA for now
# PDA_all[["Sanchez_FTS"]] <- PDA_all[["Sanchez_PNM"]]

References <- names(PDA_all)
References<-"Kalacska" 
PFTs <- names(PDA_all[[1]][[1]])
nchains <- length(PDA_all[[1]][[1]][[1]])
NperPFT <- length(References)*nchains

N = 10000

Prospect_param_names_short <- c("Cab","Car","Cw","Cm","Nlayers","ssigma")

PDA_list <- list()
for(ipft in seq_along(PFTs)){
  PDA_pft <- list()
  compt <- 1
  for (ref in seq_along(References)){
    for (ichain in seq(nchains)){
      PDA_pft[[compt]] <- PDA_all[[ref]][[1]][[ipft]][[ichain]]
      PDA_pft[[compt]][["setup"]][["names"]] <- Prospect_param_names_short
      compt <- compt + 1
    }
  }
  PDA_list[[PFTs[ipft]]] <- BayesianTools::createMcmcSamplerList(PDA_pft)
}

saveRDS(object = PDA_list,file = file.path(outdir,"PDA_leaf_kalacska.RDS"))

N <- 100000
df_pft <- data.frame()
for(ipft in seq_along(PFTs)){
  sampling <- getSample(PDA_list[[PFTs[ipft]]],numSamples = N)[1:N,]
  colnames(sampling) <- Prospect_param_names_short
  df_pft <- rbind(df_pft,
                  as.data.frame(sampling) %>% dplyr::select(-c("ssigma")) %>% pivot_longer(cols = c("Cab","Car","Cw","Cm","Nlayers")) %>% mutate(pft = PFTs[ipft]))
}

# Plot distributions

ggplot(df_pft) +
  geom_density(aes(x = value, fill = pft),alpha = 0.2) +
  facet_wrap(~name,scales = "free") + theme_bw()

# Plot spectra from posteriors

Nensemble = 250

samples_liana = getSample(PDA_list[["Liana_optical"]],numSamples=Nensemble,end=10000)[,c(5,1,2,3,4)]

# run Prospect5 Liana
df_spectrum_L = data.frame()
for (i in seq(1,min(length(unique(samples)),Nensemble))){
  current_parameter_set <- samples_liana[i,]
  current_model_output <- PEcAnRTM::prospect(current_parameter_set, version = "5")
  df_spectrum_L <- rbind(df_spectrum_L,
                         data.frame(num = i,wv = 400:2500,reflectance = current_model_output[,1],transmittance = current_model_output[,2]))
}

df_spectrum_L <- df_spectrum_L %>% mutate(pft = "Liana")

samples_tree = getSample(PDA_list[["Tree_optical"]],numSamples=Nensemble,end=10000)[,c(5,1,2,3,4)]
# run Prospect5 Tree
df_spectrum_T = data.frame()
for (i in seq(1,min(length(unique(samples)),Nensemble))){
  current_parameter_set <- samples_tree[i,]
  current_model_output <- PEcAnRTM::prospect(current_parameter_set, version = "5")
  df_spectrum_T <- rbind(df_spectrum_T,
                         data.frame(num = i,wv = 400:2500,reflectance = current_model_output[,1],transmittance = current_model_output[,2]))
}
df_spectrum_T <- df_spectrum_T %>% mutate(pft = "Tree")

spectra <- rbind(df_spectrum_T,df_spectrum_L) %>% mutate(band = case_when(wv <= 700 ~ 1,
                                                                          wv <= 2500 ~ 2))

spectra_sum <- spectra %>% group_by(pft,wv) %>% summarise(r_m = mean(reflectance),
                                                          t_m = mean(transmittance),
                                                          r_min = min(reflectance),
                                                          r_max = max(reflectance))

ggplot(spectra_sum, aes(x = wv,colour = pft)) +
  geom_line(aes(y = r_m)) +
  geom_line(aes(y = r_min),linetype = 2) +
  geom_line(aes(y = r_max),linetype = 2) +
  geom_line(aes(y = 1-t_m),linetype = 1)+
  scale_y_continuous(sec.axis = sec_axis(trans = ~ 1 - 1 * ., name = "Transmittance [-]")) +
  # geom_ribbon(aes(ymin = r_m,ymax=1-t_m,fill = pft),alpha = 0.2,linetype=0)+
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]") +
  scale_color_manual(values = c("red","blue")) +
  scale_fill_manual(values = c("red","blue")) +
  theme_bw()

##############################################################################################
# For testing
rm(list = ls())

crown_mod = 0
use_leaf_PDA <- 1
use_meta.analysis <- use_prior <- 0
outdir <- "~/output/"

Colors <- c("#137300","#1E64C8")
df_PFT <- data.frame(names = c("Liana_optical","Tree_optical"),PFTnum = c(17,3),Col = Colors)
df_PFT <- df_PFT %>% arrange(PFTnum)%>% mutate(Col = Colors)
df_PFT <- df_PFT %>% arrange(names)

npft <- nrow(df_PFT)

if (crown_mod == 0){
  dis2find <- c('b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','Nlayers',
                'Cab','Car','Cw','Cm')
  dis2find_prim <- c('b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','Nlayers',
                     'Cab','Car','Cw','Cm')
} else {
  dis2find <- c('b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','Nlayers',
                'Cab','Car','Cw','Cm','b1Ca','b2Ca')
  dis2find_prim <- c('b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','Nlayers',
                     'Cab','Car','Cw','Cm','b1Ca','b2Ca')
}

pft_lowers <- c(b1Bl_small = 0.001, b2Bl_small = 1,b1Bl_large = 0.001, b2Bl_large = 1 , SLA = 1, orient_factor = -0.5,
                clumping_factor = 0.4, b1Ca = 0.1, b2Ca = 0.5, Nlayers = 1,
                Cab = 0, Car = 0,Cw = 0,Cm = 0.)
pft_uppers <- c(b1Bl_small = 0.08, b2Bl_small = 2.5,b1Bl_large = 0.08, b2Bl_large = 2.5, SLA = 100, orient_factor = 0.5,
                clumping_factor = 0.9, b1Ca = 2, b2Ca = 3, Nlayers = 5,
                Ca = 150, Car = 50,Cw = 0.1,Cm = 0.1)

prospect_PDA_file <- file.path(outdir,"PDA_all_leaf.RDS")
if (file.exists(prospect_PDA_file) & use_leaf_PDA) {PDA_results <- readRDS(prospect_PDA_file)} else {PDA_results <- NULL}

# Import Posterior/Optimized distributions
N <- 10000
Ndist <- length(dis2find)
sampler_all <- matrix(NA,N,(Ndist)*npft)
lower_all <- upper_all <- best_all <- rep(NA,(Ndist)*npft)

for (ipft in seq(npft)){
  current_pft <- as.character(df_PFT$names[ipft])
  
  postfile <- file.path(outdir,"pft",current_pft,"post.distns.Rdata")
  priorfile <- file.path(outdir,"pft",current_pft,"prior.distns.Rdata")
  
  if (file.exists(prospect_PDA_file) & use_leaf_PDA) {
    names_dis_PDA <- names(MAP(PDA_results[[current_pft]])$parametersMAP)
  } else {
    names_dis_PDA <- NULL
  }
  
  
  if (file.exists(postfile) & use_meta.analysis){
    load(postfile)
    distns <- post.distns
    Distributions <- rownames(post.distns)
  } else if (file.exists(priorfile) & use_prior){
    load(priorfile)
    distns <- prior.distns
    Distributions <- rownames(prior.distns)
  } else {
    distns <- NULL
    Distributions <- NULL
  }
  
  sampler <- matrix(NA,N,length(dis2find))
  lower <- upper <- best <- rep(NA,length(dis2find))
  
  for (idis in seq(dis2find)){
    current_dis <- dis2find[idis]
    if (current_dis %in% names_dis_PDA){
      pos <- which(current_dis == names_dis_PDA)
      sampling <- as.vector(getSample(PDA_results[[current_pft]],whichParameters=pos,numSamples = N))[1:N]
    } else {
      pos <- which(Distributions == current_dis)
      
      if (isempty(pos)){
        unif_distribution <- list(distn="unif",parama=pft_lowers[current_dis],paramb=pft_uppers[current_dis])
        sampling <- get.sample(unif_distribution,n=N)
        
      } else {
        sampling <- get.sample(distns[pos,],n=N)
      }
    }
    
    sampler[,idis] <- sampling
    lower[idis] <- min(sampler[,idis])
    upper[idis] <- max(sampler[,idis])
    best[idis]  <- median(sampler[,idis])
  }
  
  pos <- ((ipft-1)*(Ndist)+1):(((ipft)*(Ndist)))
  
  sampler_all[,pos]<-sampler
  lower_all[pos]<-lower
  upper_all[pos]<-upper
  best_all[pos]<-best
  
  names(lower_all)[pos] <- names(upper_all)[pos] <- names(best_all)[pos] <- dis2find_prim
}

colnames(sampler_all) <- rep(dis2find_prim,npft)

ssigma_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)
soil_moisture_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)

lower_all <- c(0.001,0.001,lower_all)
upper_all <- c(1,1,upper_all)
best_all <- c(0.5,0.5,best_all)
sampler_all <- cbind(ssigma_sampling,soil_moisture_sampling,sampler_all)

prior <- createPriorDensity(sampler_all, method = "multivariate", eps = 1e-10,
                            lower = lower_all, upper = upper_all, best = best_all)


