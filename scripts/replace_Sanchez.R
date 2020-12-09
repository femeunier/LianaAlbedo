rm(list = ls())

library(dplyr)
library(ggplot2)
library(pracma)

data.all <- readRDS(file = "./Sanchez_allfits2.RDS")

GFs <- names(data.all)

####################################################################################
# Functions
run_prospect <- function(params,waves){
  # PROSPECT parameters
  Nlayers <- params[1]
  Cab <- params[2]
  Car <- params[3]
  Cw <- params[4]
  Cm <- params[5]
  
  optical_param <- c(Nlayers,Cab,Car,Cw,Cm)
  names(optical_param) <- c('Nlayers','chlab','carotenoids','Cw','Cm')
  
  # Call RTM
  result <- tryCatch(
    PEcAnRTM::prospect(optical_param, version = "5"),
    error = function(e) NULL)
  if (is.null(result)) return(-1e20)
  reflectance <- result[,"reflectance"]
  if (any(!is.finite(reflectance))) return(-1e20)
  if (any(reflectance < 0)) return(-1e20)
  if (any(reflectance > 1)) return(-1e20)
  
  reflectance_waves <- interp1(x = PEcAnRTM::wavelengths(reflectance),
                               y = as.vector(matrix(reflectance)),
                               waves)
  return(reflectance_waves)
}


# test optimization
create_likelihood <- function(observed, waves) {
  function(params) {
    
    ssigma <- params[6]
    reflectance_waves <- run_prospect(params,waves)
    # Calculate likelihood
    ll <- sum(dnorm(reflectance_waves, observed, ssigma, log = TRUE))
    return(ll)
  }
}

################################################################################
df_formatted <- data.frame()
compt <- 1
for (pft in seq(names(data.all))){
  for (species in seq(names(data.all[[pft]]))){
    for (ind in seq(length(data.all[[pft]][[species]]))) {
      
      print(compt)
      data.temp <- data.all[[pft]][[species]][[ind]]
      
      if (!is.null(data.temp)){
        df_formatted <- rbind(df_formatted,data.temp[["posterior"]] %>% mutate(pft = GFs[pft],species = species,ind = ind,site = data.temp[["site"]]))
        compt <- 1 + compt
      }
    }
  }
}

ggplot(df_formatted) +
  geom_density(aes(x = sim_m -  obs_m, fill = pft),alpha = 0.2) +
  theme_bw()


selection <- df_formatted %>% filter(waves < 680 | waves > 750) %>% group_by(site,pft) %>% sample_n(size = 10)
ggplot(selection,aes(x=obs_m,y = sim_m,color = pft)) +
  geom_point() +
  geom_errorbarh(aes(xmin = obs_low,xmax=obs_high)) +
  geom_errorbar(aes(ymin = sim_low,ymax=sim_high)) +
  geom_abline(slope = 1,linetype = 2) +
  facet_wrap(~ site,scales = "free") +
  theme_bw()


Prospect_param_names <- c("Nlayers","Cab","Car","Cw","Cm","ssigma")
df_selection <- df_formatted %>% group_by(site,pft) %>% sample_n(size =1) %>% dplyr::select(pft,species,ind,site)
Nprospect <- 1000
waves <- 400:2500

# Uniform prior
Prospect_param_names <- c("Nlayers","Cab","Car","Cw","Cm","ssigma")
pft_lowers <- c(Nlayers = 1, Cab = 0 , Car = 0,Cw = 0, Cm = 0, ssigma = 0)
pft_uppers <-  c(Nlayers = 5, Cab = 100, Car = 50,Cw = 0.1, Cm = 0.1, ssigma = 1)

new_All_leaf_spectra <- new_ensemble_posterior_all <- data.frame()
for (i in seq(nrow(df_selection))){
  print(i/nrow(df_selection))
  cdata <- data.all[[df_selection$pft[i]]][[df_selection$species[i]]][[df_selection$ind[i]]]
  
  if (!is.null(cdata)){
    cspectrum <- cdata[["spectrum"]]
    observed <- as.vector(t(cspectrum[,"value"]))
    waves_obs <- as.vector(t(cspectrum[,"wv"]))
    observations <- data.frame(wavelength = unlist(waves_obs), obs = unlist(observed))
    
    prior <- BayesianTools::createUniformPrior(pft_lowers, pft_uppers)
    likelihood <- create_likelihood(observed=observed[waves_obs<680 | waves_obs>750], waves_obs[waves_obs<680 | waves_obs>750])
    settings_MCMC <- list(iterations = 100000, nrChains = 1)
    
    # Run inversion
    setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
    samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
    samples <- BayesianTools::runMCMC(samples, settings = settings_MCMC)
    
    cdata[["samples"]] <- samples
    
    posteriors <- BayesianTools::getSample(cdata[["samples"]],numSamples = Nprospect)
    colnames(posteriors) <- Prospect_param_names
    
    df_runs <- data.frame()
    for (irun in seq(1,min(nrow(posteriors),Nprospect))){
      cparams <- c(posteriors[irun,"Nlayers"],posteriors[irun,"Cab"],posteriors[irun,"Car"],posteriors[irun,"Cw"],posteriors[irun,"Cm"])
      crun <- run_prospect(cparams,waves)
      df_runs <- rbind(df_runs,
                       data.frame(wavelength = waves,sim = crun,run = irun))
    }
    
    # leaves <- as.numeric(unique(as.vector(t(cspectrum[,"name"]))))
    df_runs_sum <- df_runs %>% group_by(wavelength) %>% summarise(rmin = min(sim,na.rm=TRUE),
                                                             rmax = max(sim,na.rm=TRUE),
                                                             alphamin = quantile(sim,0.025),
                                                             alphamax = quantile(sim,0.975),
                                                             median = quantile(sim,0.5)) %>% mutate(pft= case_when(df_selection$pft[i] == "Liana" ~ "Liana_optical",
                                                                                                                   TRUE ~ "Tree_optical"),
                                                                                                    ref = case_when(df_selection$site[i] == "FTS" ~ "Sanchez_FTS",
                                                                                                                    TRUE ~ "Sanchez_PNM"))
                                                             

    df_obs_sum <- observations %>% group_by(wavelength) %>% summarise(Reflectance_min = min(obs,na.rm=TRUE),
                                                                      Reflectance_max = max(obs,na.rm=TRUE),
                                                                      Reflectance_alphamin = quantile(obs,0.025),
                                                                      Reflectance_alphamax = quantile(obs,0.975),
                                                                      Reflectance_median = median(obs)) %>% mutate(pft= case_when(df_selection$pft[i] == "Liana" ~ "Liana_optical",
                                                                                                                                    TRUE ~ "Tree_optical"),
                                                                                                                     ref = case_when(df_selection$site[i] == "FTS" ~ "Sanchez_FTS",
                                                                                                                                     TRUE ~ "Sanchez_PNM"))
    
    new_ensemble_posterior_all <- rbind(new_ensemble_posterior_all,
                                    df_runs_sum)
    new_All_leaf_spectra <- rbind(new_All_leaf_spectra,
                                  as.data.frame(df_obs_sum))
    

  }
}

ensemble_posterior_all <- rbind(ensemble_posterior_all %>% filter(!(ref == "Sanchez_PNM" | ref == "Sanchez_FTS")),
                                new_ensemble_posterior_all)
All_leaf_spectra <- bind_rows(list(All_leaf_spectra %>% filter(!(ref == "Sanchez_PNM" | ref == "Sanchez_FTS")),
                                   new_All_leaf_spectra))

rm(data.all)

save.image("~/data/RTM/Inverse_leaf_spectrum_sanchez2.Rdata")

