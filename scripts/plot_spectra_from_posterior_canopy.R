rm(list = ls())

best_sets <- c("~/data/RTM/current_samples_Kalacska_priors_all.rds",
               "~/data/RTM/current_samples_Marvin.rds",
               "~/data/RTM/current_samples_Sanchez_all.rds",
               "~/data/RTM/current_samples_Foster_all.rds")

best_sets <- c("~/data/RTM/current_samples_Sanchez_all.rds")
best_set <- readRDS(best_sets)

Nensemble = 250
ensemble <- getSample(best_set,numSamples = Nensemble)

all.spectra <- data.frame()
for (i in seq(1,max(Nensemble,nrow(ensemble)))){
  
  temp.param <- ensemble[i,]
  Liana_param <- temp.param[3:(2+((length(temp.param)-2)/2))]
  best_set_L<- c(Liana_param['Nlayers'],Liana_param['Cab'],Liana_param['Car'],Liana_param['Cw'],Liana_param['Cm'])
  best_run_L <- PEcAnRTM::prospect(best_set_L, version = "5")
  
  Tree_param <- temp.param[(3+((length(temp.param)-2)/2)):length(temp.param)]
  best_set_T<- c(Tree_param['Nlayers'],Tree_param['Cab'],Tree_param['Car'],Tree_param['Cw'],Tree_param['Cm'])
  best_run_T <- PEcAnRTM::prospect(best_set_T, version = "5")
  
  test.plot <- rbind(all.spectra,
    rbind(as.data.frame(best_run_L) %>% mutate(wl = PEcAnRTM::wavelengths(best_run_L),
                                                          pft = "Liana"),
    as.data.frame(best_run_T) %>% mutate(wl = PEcAnRTM::wavelengths(best_run_T),
                                                          pft = "Tree")) %>% mutate(run = i))
  
  
}

group.plot <- test.plot %>% group_by(wl,pft) %>% summarise(r_m = mean(reflectance),
                                             r_min = min(reflectance),
                                             r_max = max(reflectance))

MAP_samples <- MAP(best_set)$parametersMAP

Liana_param <- MAP_samples[3:(2+((length(MAP_samples)-2)/2))]
best_set_L<- c(Liana_param['Nlayers'],Liana_param['Cab'],Liana_param['Car'],Liana_param['Cw'],Liana_param['Cm'])
best_run_L <- PEcAnRTM::prospect(best_set_L, version = "5")

Tree_param <- MAP_samples[(3+((length(MAP_samples)-2)/2)):length(MAP_samples)]
best_set_T<- c(Tree_param['Nlayers'],Tree_param['Cab'],Tree_param['Car'],Tree_param['Cw'],Tree_param['Cm'])
best_run_T <- PEcAnRTM::prospect(best_set_T, version = "5")

best.run <- rbind(as.data.frame(best_run_L) %>% mutate(wl = PEcAnRTM::wavelengths(best_run_L),
                                                        pft = "Liana"),
                   as.data.frame(best_run_T) %>% mutate(wl = PEcAnRTM::wavelengths(best_run_T),
                                                        pft = "Tree"))

ggplot() +
  geom_line(data = best.run,aes(x = wl,y = reflectance,color = pft),linetype = 2) +
  geom_ribbon(data = group.plot,aes(x = wl,ymin = r_min,ymax = r_max,color = pft,fill = pft)) +
  theme_bw()



# ggplot() +
#   geom_line(data = test.plot,aes(x = wl,y = transmittance)) +
#   theme_bw()
