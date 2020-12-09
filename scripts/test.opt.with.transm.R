rm(list=ls())

library(cowplot)
library(BayesianTools)
library(dplyr)
library(redr)
library(PEcAn.all)
library(pracma)
library(reshape2)
library(ggplot2)
library(ggridges)

iter <-   5000
nrChains  <- 1
wl.max <- 2500
wl.min  <- 350
alpha  <- 0.05

Nrun_prospect <- 2000

files <- c("Figures5and6_sanchez2009_FTS") #Figures5and6_sanchez2009_FTS

Colors <- c("#137300","#1E64C8")

use_meta.analysis <- FALSE
direct <- NULL
df_PFT <- data.frame(names = c("Liana_optical","Tree_optical"),PFTnum = c(17,3),Col = Colors)
df_PFT <- df_PFT %>% arrange(PFTnum)%>% mutate(Col = Colors)
df_PFT <- df_PFT %>% arrange(names)

npft <- nrow(df_PFT)

settings_MCMC <- list(iterations = iter, nrChains = nrChains)

PDA <- list()
best_set <- best_run <- list()

N=10000

# Define likelihood
create_likelihood <- function(waves,observed,observed2 = NULL) {
  function(params) {
    
    # PROSPECT parameters
    Cab <- params[1]
    Car <- params[2]
    Cw <- params[3]
    Cm <- params[4]
    Nlayers <- params[5]
    
    ssigma <- params[6]
    
    optical_param <- c(Nlayers,Cab,Car,Cw,Cm)
    names(optical_param) <- c('Nlayers','chlab','carotenoids','Cw','Cm')
    
    # Call RTM
    result <- tryCatch(
      PEcAnRTM::prospect(optical_param, version = "5"),
      error = function(e) NULL)
    if (is.null(result)) return(-1e20)
    reflectance <- result[,"reflectance"]
    transmittance <- result[,"transmittance"]
    if (any(!is.finite(reflectance)) | any(!is.finite(transmittance))) return(-1e20)
    if (any(reflectance < 0) | any(transmittance < 0)) return(-1e20)
    if (any(reflectance > 1) | any(transmittance > 1)) return(-1e20)
    
    reflectance_waves <- interp1(x = PEcAnRTM::wavelengths(reflectance),
                                 y = as.vector(matrix(reflectance)),
                                 waves)
    transmittance_waves <- interp1(x = PEcAnRTM::wavelengths(transmittance),
                                 y = as.vector(matrix(transmittance)),
                                 waves)
    # Calculate likelihood
    if (!is.null(observed2)){
      ll <- sum(dnorm(transmittance_waves, observed2, ssigma, log = TRUE) +
                  dnorm(reflectance_waves, observed, ssigma, log = TRUE))
    }else{
      ll <- sum(dnorm(reflectance_waves, observed, ssigma, log = TRUE))
    }
    return(ll)
  }
}

############################################################################################
# Loop over papers

marginalPlots <- spectra_post <- performance_prospect_plot <- PDA_all <- list()

performance_all <- best_set_all <- posterior_all <- ensemble_posterior_all <-
  ensemble_posterior_all2 <- data.frame()

model_performance <- sensitivities_all <- p.values_all <- data.frame()

for (ifile in seq(files)){
  data <- load_rds(file.path("~/data/RTM/",files[ifile]))
  
  #############################################################################################
  # Loop over PFTs
  
  # Uniform prior
  Prospect_param_names <- c("Ca","Cb","Car","Cw","Cm","Nlayers","ssigma")
  Prospect_param_names_short <- c("Cab","Car","Cw","Cm","Nlayers","ssigma")
  pft_lowers <- c(chlorophyll_a = 0,chlorophyll_b = 0, carotenoids = 0,Cw = 0, SLA = 1, Nlayers = 1, ssigma = 0)
  pft_uppers <-  c(chlorophyll_a = 100,chlorophyll_b = 50, carotenoids = 50,Cw = 0.1, SLA = 100, Nlayers = 5., ssigma = 1)
  
  dis2find <- c('chlorophyll_a','chlorophyll_b','carotenoids','Cw','SLA','Nlayers',"ssigma")
  
  for (ipft in seq(npft)){ #seq(npft)
    current_pft <- as.character(df_PFT$names[ipft])
    
    # Define priors
    # From posterior
    if ((!is.null(direct))) {
      postfile <- file.path(direct,"pft",current_pft,"post.distns.Rdata")
      if (use_meta.analysis){
        load(postfile)
        distns <- post.distns
        Distributions <- rownames(post.distns)
      } else {
        priorfile <- file.path(direct,"pft",current_pft,"prior.distns.Rdata")
        load(priorfile)
        distns <- prior.distns
        Distributions <- rownames(prior.distns)
      }
    } else {
      Distributions <- NULL
    }
    
    sampler <- matrix(NA,N,length(dis2find))
    lower <- upper <- best <- rep(NA,length(dis2find))
    
    for (idis in seq(dis2find)){
      current_dis = dis2find[idis]
      pos <- which(Distributions == current_dis)
      
      if (isempty(pos)){
        unif_distribution <- list(distn="unif",parama=pft_lowers[current_dis],paramb=pft_uppers[current_dis])
        sampling <- get.sample(unif_distribution,n=N)
      } else {
        sampling <- get.sample(distns[pos,],n=N)
      }
      
      if (current_dis == 'SLA'){
        sampling = 1/sampling/10 # unit conversions
      }
      sampler[,idis] <- sampling
      lower[idis] <- min(sampler[,idis])
      upper[idis] <- max(sampler[,idis])
      best[idis]  <- median(sampler[,idis])
    }
    
    colnames(sampler) <- names(lower) <- names(upper) <- names(best) <- Prospect_param_names
    
    # combine Ca and Cb
    sampler[,"Ca"] <- sampler[,"Ca"] + sampler[,"Cb"]
    colnames(sampler)[colnames(sampler)=="Ca"] = "Cab"
    sampler <- sampler[,!(colnames(sampler) %in% "Cb")]
    
    lower["Ca"] <- lower["Ca"] + lower["Cb"]
    lower <- lower[!(names(lower) %in% c("Cb"))]
    names(lower)[(names(lower) %in% c("Ca"))] <- "Cab"
    
    upper["Ca"] <- upper["Ca"] + upper["Cb"]
    upper <- upper[!(names(upper) %in% c("Cb"))]
    names(upper)[(names(upper) %in% c("Ca"))] <- "Cab"
    
    best["Ca"] <- best["Ca"] + best["Cb"]
    best <- best[!(names(best) %in% c("Cb"))]
    names(best)[(names(best) %in% c("Ca"))] <- "Cab"
    
    # 1/as.numeric(config_temp$SLA)/10*2
    prior <- createPriorDensity(sampler, method = "multivariate", eps = 1e-10,
                                lower = lower, upper = upper, best = best)
    
    # Read Leaf spectra 
    
    if ("Transmittance" %in% colnames(data)){
      temp <- data %>% filter(pft==current_pft & wavelength>=wl.min & wavelength<=wl.max) %>% dplyr::select(c('wavelength','Reflectance','Transmittance'))
      observation <- temp[["Reflectance"]]
      observation2 <- temp[["Transmittance"]]  
      waves <- temp[["wavelength"]]
      likelihood <- create_likelihood(waves,observed = observation,observed2 = observation2)
    } else{
      temp <- data %>% filter(pft==current_pft & wavelength>=wl.min & wavelength<=wl.max) %>% dplyr::select(c('wavelength','Reflectance'))
      observation <- temp[["Reflectance"]] 
      waves <- temp[["wavelength"]]
      likelihood <- create_likelihood(waves,observation)
    }

    # Run inversion
    setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
    samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
    samples <- BayesianTools::runMCMC(samples, settings = settings_MCMC)
    
    MAP_samples <- MAP(samples)$parametersMAP
    names(MAP_samples) <- Prospect_param_names_short
    best_set[[current_pft]]<- c(MAP_samples['Nlayers'],MAP_samples['Cab'],MAP_samples['Car'],MAP_samples['Cw'],MAP_samples['Cm'])
    
    best_run <- PEcAnRTM::prospect(best_set[[current_pft]], version = "5")
    
    test.plot <- as.data.frame(best_run) %>% mutate(wl = PEcAnRTM::wavelengths(best_run))
    ggplot() +
      geom_line(data = test.plot,aes(x = wl,y = reflectance)) +
      geom_point(data = data %>% filter(pft == current_pft),aes(x = wavelength,y = Reflectance)) +
      theme_bw()
    
    ggplot() +
      geom_line(data = test.plot,aes(x = wl,y = transmittance)) +
      geom_point(data = data %>% filter(pft == current_pft),aes(x = wavelength,y = Transmittance)) +
      theme_bw()
    
  }
}