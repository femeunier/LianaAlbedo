rm(list = ls())

######################################################################################################
# Leaf reflectance
######################################################################################################
library(redr)
library(dplyr)
library(Hmisc)

load(file = "~/data/RTM/Inverse_leaf_spectrum.Rdata")

posterior <- ensemble_posterior_all %>% dplyr::select(wavelength,alphamin,alphamax,median,ref,pft) %>% 
  rename(posterior_alphamin=alphamin,
         posterior_median=median,
         posterior_alphamax=alphamax)

# Merge sanchez
Sanchez <- readRDS(file = "~/data/RTM/fit_Sanchez.RDS")
posterior <- rbind(posterior %>% filter(!(ref %in% c("Sanchez_FTS","Sanchez_PNM"))),
                   Sanchez %>% mutate(ref = paste0("Sanchez_",site),
                                      pft = case_when(GF == "Liana" ~ "Liana_optical",
                                                      GF == "Tree" ~ "Tree_optical"),
                                      posterior_alphamin = NA,
                                      posterior_alphamax = NA) %>% rename(wavelength = wv,
                                                                posterior_median = sim) %>% ungroup() %>% dplyr::select(-c(site,GF,obs)))   

All_leaf_spectra <- bind_rows(list(All_leaf_spectra %>% filter(!(ref %in% c("Sanchez_FTS","Sanchez_PNM"))),
                   Sanchez %>% mutate(ref = paste0("Sanchez_",site),
                                      pft = case_when(GF == "Liana" ~ "Liana_optical",
                                                      GF == "Tree" ~ "Tree_optical"),
                                      Reflectance_min = NA,
                                      Reflectance_max = NA,
                                      Reflectance_alphamin = NA,
                                      Reflectance_alphamax = NA) %>% rename(wavelength = wv,
                                                                          Reflectance_median = obs) %>% ungroup() %>% dplyr::select(-c(site,GF,sim))))


ready.for.table <- All_leaf_spectra %>% mutate(wavelength = as.integer(round(wavelength))) %>%left_join(posterior,by=c("ref","pft","wavelength")) %>% filter(wavelength < 680 | wavelength > 720)

####################################################################################################
references <- unique(All_leaf_spectra$ref)

ggplot(data = ready.for.table) +
  geom_line(aes(x = wavelength,Reflectance_median,color = pft,group = interaction(pft,ref))) +
  geom_line(aes(x = wavelength,posterior_median,color = pft,group = interaction(pft,ref)),linetype = 2) +
  facet_wrap(~ ref) +
  # scale_x_continuous(limits = c(400,850)) +
  # scale_y_continuous(limits = c(0,0.2)) +
  theme_bw()

ggplot(data = ready.for.table %>% filter(wavelength < 680 | wavelength > 720)) +
  geom_point(aes(Reflectance_median,y = posterior_median,color = pft,group = interaction(pft,ref))) +
  facet_wrap(~ ref) +
  # scale_x_continuous(limits = c(400,850)) +
  # scale_y_continuous(limits = c(0,0.2)) +
  theme_bw()

# saveRDS(object = ready.for.table,"~/data/RTM/obsvssim.RDS")

all_data <- all_models <- data.frame()
for (reference in references){
  
  ##################################################################
  # Data 
  temp_T <- All_leaf_spectra %>% filter(ref == reference &
                                   pft == "Tree_optical")
  wavelengths_T <- temp_T %>% pull(wavelength)
  R_T <- temp_T %>% pull(Reflectance_median)
  
  
  temp_L <- All_leaf_spectra %>% filter(ref == reference &
                                   pft == "Liana_optical")
  wavelengths_L <- temp_L %>% pull(wavelength)
  wavelengths_L_out <- seq(round(min(wavelengths_L))+1,round(max(wavelengths_L))-1,1)
  R_L <- temp_L %>% pull(Reflectance_median)
  
  R_L_out <- approxExtrap(x = wavelengths_L,
                     y = R_L,
                     xout = wavelengths_L_out)$y
  
  
  R_T_extrap <- 
    approxExtrap(x = wavelengths_T,
               y = R_T,
               xout = wavelengths_L_out)$y
  
  all_data <- rbind(all_data,
                    data.frame(wv = wavelengths_L_out,
                               R_T = R_T_extrap,
                               R_L = R_L_out,
                               ref = reference))
  
  ##################################################################
  # Models 
  temp_T <- posterior %>% filter(ref == reference &
                                          pft == "Tree_optical")
  wavelengths_T <- temp_T %>% pull(wavelength)
  R_T <- temp_T %>% pull(posterior_median)
  
  
  temp_L <- posterior %>% filter(ref == reference &
                                          pft == "Liana_optical")
  wavelengths_L <- temp_L %>% pull(wavelength)
  R_L <- temp_L %>% pull(posterior_median)
  
  wavelengths_L_out <- seq(round(min(wavelengths_L))+1,round(max(wavelengths_L))-1,1)
  R_L_out <- approxExtrap(x = wavelengths_L,
                          y = R_L,
                          xout = wavelengths_L_out)$y
  
  R_T_extrap <- 
    approxExtrap(x = wavelengths_T,
                 y = R_T,
                 xout = wavelengths_L_out)$y
  
  all_models <- rbind(all_models,
                    data.frame(wv = wavelengths_L_out,
                               R_T = R_T_extrap,
                               R_L = R_L_out,
                               ref = reference))
}

ggplot() +
  geom_line(data = all_models,
            aes(x = wv,
                y = R_L - R_T)) +
  geom_line(data = all_data,
            aes(x = wv,
                y = R_L - R_T),linetype = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = c(400,700,1400,1500,2500)) +
  facet_wrap(~ ref) + theme_bw()

# Data
Delta_PAR <- all_data %>% filter(wv <= 700 & wv >= 500)

# Delta_PAR %>% group_by(ref) %>% summarise(rmse = sqrt(sum((R_L-R_T)^2)/(length(R_T))))
# all_models %>% group_by(ref) %>% summarise(rmse = sqrt(sum((R_L-R_T)^2)/(length(R_T))))

metrics_delta_PAR <- redr::calc_metrics(obs=Delta_PAR %>% pull(R_L),
                                  sim=Delta_PAR %>% pull(R_T))

Delta_NIR <- all_data %>% filter(wv >=800 & wv <= 1400)
metrics_delta_NIR <- redr::calc_metrics(Delta_NIR %>% pull(R_L),
                                  Delta_NIR %>% pull(R_T))

Delta_SIR <- all_data %>% filter(wv >=1500)
metrics_delta_SIR <- redr::calc_metrics(Delta_SIR %>% pull(R_L),
                                        Delta_SIR %>% pull(R_T))

# Models
Delta_PAR_m <- all_models %>% filter(wv <= 700 & wv >= 500)
metrics_delta_PAR_m <- redr::calc_metrics(Delta_PAR_m %>% pull(R_L),
                                    Delta_PAR_m %>% pull(R_T))

Delta_NIR_m <- all_models %>% filter(wv >=800 & wv <= 1400)
metrics_delta_NIR_m <- redr::calc_metrics(Delta_NIR_m %>% pull(R_L),
                                    Delta_NIR_m %>% pull(R_T))

Delta_SIR_m <- all_models %>% filter(wv >= 1500)
metrics_delta_SIR_m <- redr::calc_metrics(Delta_SIR_m %>% pull(R_L),
                                          Delta_SIR_m %>% pull(R_T))

#####################################################################################################
# Calculation
PAR_tree <- ready.for.table %>% filter(pft == "Tree_optical" & wavelength<= 700 & wavelength >= 500)
metrics_PAR_tree <- redr::calc_metrics(obs=PAR_tree%>%pull(Reflectance_median),
                                 sim=PAR_tree%>%pull(posterior_median))

PAR_liana <- ready.for.table %>% filter(pft == "Liana_optical" & wavelength<= 700 & wavelength >= 500)
metrics_PAR_liana <- redr::calc_metrics(obs=PAR_liana%>%pull(Reflectance_median),
                                 sim=PAR_liana%>%pull(posterior_median))

PAR <- ready.for.table %>% filter(wavelength<= 700 & wavelength >= 500)
all_data  %>% filter(wv >= 1500 & wv<= 2500) %>% summarise(Tree_m = mean(R_T),
                                                         Liana_m = mean(R_L))

metrics_PAR <- redr::calc_metrics(obs=PAR%>%pull(Reflectance_median),
                                  sim=PAR%>%pull(posterior_median))

NIR_tree <- ready.for.table %>% filter(pft == "Tree_optical" & wavelength>=800 & wavelength <= 1400)
metrics_NIR_tree <- redr::calc_metrics(obs=NIR_tree%>%pull(Reflectance_median),
                                 sim=NIR_tree%>%pull(posterior_median))

NIR_liana <- ready.for.table %>% filter(pft == "Liana_optical" & wavelength>=800 & wavelength <= 1400)
metrics_NIR_liana <- redr::calc_metrics(obs=NIR_liana%>%pull(Reflectance_median),
                                 sim=NIR_liana%>%pull(posterior_median))

NIR <- ready.for.table %>% filter(wavelength>=800 & wavelength <= 1400)
metrics_NIR <- redr::calc_metrics(obs=NIR%>%pull(Reflectance_median),
                                        sim=NIR%>%pull(posterior_median))

SIR_tree <- ready.for.table %>% filter(pft == "Tree_optical" & wavelength>=1500)
metrics_SIR_tree <- redr::calc_metrics(obs=SIR_tree%>%pull(Reflectance_median),
                                       sim=SIR_tree%>%pull(posterior_median))

SIR_liana <- ready.for.table %>% filter(pft == "Liana_optical" & wavelength>=1500)
metrics_SIR_liana <- redr::calc_metrics(obs=SIR_liana%>%pull(Reflectance_median),
                                        sim=SIR_liana%>%pull(posterior_median))

SIR <- ready.for.table %>% filter(wavelength>=1500)
metrics_SIR <- redr::calc_metrics(obs=SIR%>%pull(Reflectance_median),
                                       sim=SIR%>%pull(posterior_median))

############################################################################################################
# Canopy reflectance
############################################################################################################

load(file = "~/data/RTM/Inverse_canopy_spectrum.Rdata")

model.results <- model_ensemble_all %>% dplyr::select(scenar,wavelength,rbest,alphamin,alphamax,ref) %>% rename(scenario = scenar,median = rbest)
data <- All_canopy_spectra %>% dplyr::select(ref,wavelength,scenario,Reflectance_median,Reflectance_alphamin,Reflectance_alphamax)

data_canopy <- data %>% left_join(model.results,by=c("ref","scenario","wavelength")) 

# differences 
references <- unique(data_canopy$ref)

all_data <- all_models <- data.frame()
for (reference in references){
  
  #######################################################
  # Data
  temp_T <- data_canopy %>% filter(ref == reference &
                                     scenario == "low")
  wavelengths_T <- temp_T %>% pull(wavelength)
  R_T <- temp_T %>% pull(Reflectance_median)
  
  
  temp_L <- data_canopy %>% filter(ref == reference &
                                     scenario == "high")
  wavelengths_L <- temp_L %>% pull(wavelength)
  R_L <- temp_L %>% pull(Reflectance_median)
  
  wavelengths_T_out <- seq(round(min(wavelengths_T))+1,round(max(wavelengths_T))-1,1)
  R_T_out <- approxExtrap(x = wavelengths_T,
                          y = R_T,
                          xout = wavelengths_T_out)$y
  
  
  R_L_extrap <- 
    approxExtrap(x = wavelengths_L,
                 y = R_L,
                 xout = wavelengths_T_out)$y
  
  all_data <- rbind(all_data,
                    data.frame(wv = wavelengths_T_out,
                               R_T = R_T_out,
                               R_L = R_L_extrap,
                               ref = reference))
  
  #######################################################
  # Models
  temp_T <- model.results %>% filter(ref == reference &
                                     scenario == "low")
  wavelengths_T <- temp_T %>% pull(wavelength)
  R_T <- temp_T %>% pull(median)
  
  
  temp_L <- model.results %>% filter(ref == reference &
                                     scenario == "high")
  wavelengths_L <- temp_L %>% pull(wavelength)
  R_L <- temp_L %>% pull(median)
  
  wavelengths_L_out <- seq(round(min(wavelengths_L))+1,round(max(wavelengths_L))-1,1)
  R_L_out <- approxExtrap(x = wavelengths_L,
                          y = R_L,
                          xout = wavelengths_L_out)$y
  
  
  R_T_extrap <- 
    approxExtrap(x = wavelengths_T,
                 y = R_T,
                 xout = wavelengths_L_out)$y
  
  all_models <- rbind(all_models,
                    data.frame(wv = wavelengths_L_out,
                               R_T = R_T_extrap,
                               R_L = R_L_out,
                               ref = reference))
}

data_canopy %>% ungroup() %>% summarise(summary(lm(formula = median ~ Reflectance_median))$r.squared) 
data_canopy %>% ungroup() %>% summarise(coef(lm(formula = median ~ Reflectance_median))[2]) 
data_canopy %>% group_by(ref,scenario) %>% summarise(summary(lm(formula = median ~ Reflectance_median))$r.squared) 
data_canopy %>% group_by(ref,scenario) %>% summarise(coef(lm(formula = median ~ Reflectance_median))[2]) 

# Data
Delta_PAR_canopy <- all_data %>% filter(wv <= 700 & wv >= 500)
metrics_delta_PAR_canopy <- redr::calc_metrics(Delta_PAR_canopy %>% pull(R_L),
                                         Delta_PAR_canopy %>% pull(R_T))

Delta_NIR_canopy <- all_data %>% filter(wv >=800 & wv <= 1400) 
metrics_delta_NIR_canopy <- redr::calc_metrics(Delta_NIR_canopy %>% pull(R_L),
                                         Delta_NIR_canopy %>% pull(R_T))

Delta_SIR_canopy <- all_data %>% filter(wv >= 1400) 
metrics_delta_SIR_canopy <- redr::calc_metrics(Delta_SIR_canopy %>% pull(R_L),
                                               Delta_SIR_canopy %>% pull(R_T))

# Models
Delta_PAR_canopy_m <- all_models %>% filter(wv <= 700 & wv >= 500)
metrics_delta_PAR_canopy_m <- redr::calc_metrics(Delta_PAR_canopy_m %>% pull(R_L),
                                           Delta_PAR_canopy_m %>% pull(R_T))

Delta_NIR_canopy_m <- all_models %>% filter(wv >=800 & wv <= 1400) 
metrics_delta_NIR_canopy_m <- redr::calc_metrics(Delta_NIR_canopy_m %>% pull(R_L),
                                           Delta_NIR_canopy_m %>% pull(R_T))

Delta_SIR_canopy_m <- all_models %>% filter(wv >= 1500) 
metrics_delta_SIR_canopy_m <- redr::calc_metrics(Delta_SIR_canopy_m %>% pull(R_L),
                                                 Delta_SIR_canopy_m %>% pull(R_T))

#############################################################################################################
# Calculations

PAR_tree_canopy <- data_canopy %>% filter(scenario == "low" & wavelength<= 700 & wavelength >= 500)
metrics_PAR_tree_canopy <- redr::calc_metrics(obs=PAR_tree_canopy%>%pull(Reflectance_median),
                                 sim=PAR_tree_canopy%>%pull(median))

All <- data_canopy %>% filter(scenario == "high")
metric_all <- redr::calc_metrics(obs=All%>%pull(Reflectance_median),
                                              sim=All%>%pull(median))

PAR_liana_canopy <- data_canopy %>% filter(scenario == "high" & wavelength<= 700 & wavelength >= 500)
metrics_PAR_liana_canopy <- redr::calc_metrics(obs=PAR_liana_canopy%>%pull(Reflectance_median),
                                 sim=PAR_liana_canopy%>%pull(median))

PAR_canopy <- data_canopy %>% filter(wavelength<= 700 & wavelength >= 500)
metrics_PAR_canopy <- redr::calc_metrics(obs=PAR_canopy%>%pull(Reflectance_median),
                                               sim=PAR_canopy%>%pull(median))

NIR_tree_canopy <- data_canopy %>% filter(scenario == "low" & wavelength>=800 & wavelength <= 1400)
metrics_NIR_tree_canopy <- redr::calc_metrics(obs=NIR_tree_canopy%>%pull(Reflectance_median),
                                 sim=NIR_tree_canopy%>%pull(median))

NIR_liana_canopy <- data_canopy %>% filter(scenario == "high" & wavelength>=800 & wavelength <= 1400)
metrics_NIR_liana_canopy <- redr::calc_metrics(obs=NIR_liana_canopy%>%pull(Reflectance_median),
                                 sim=NIR_liana_canopy%>%pull(median))

NIR_canopy <- data_canopy %>% filter(wavelength>=800 & wavelength <= 1400)
metrics_NIR_canopy <- redr::calc_metrics(obs=NIR_canopy%>%pull(Reflectance_median),
                                               sim=NIR_canopy%>%pull(median))

SIR_tree_canopy <- data_canopy %>% filter(scenario == "low" & wavelength>=1500)
metrics_SIR_tree_canopy <- redr::calc_metrics(obs=SIR_tree_canopy%>%pull(Reflectance_median),
                                              sim=SIR_tree_canopy%>%pull(median))

SIR_liana_canopy <- data_canopy %>% filter(scenario == "high" & wavelength>=1500)
metrics_SIR_liana_canopy <- redr::calc_metrics(obs=SIR_liana_canopy%>%pull(Reflectance_median),
                                               sim=SIR_liana_canopy%>%pull(median))

SIR_canopy <- data_canopy %>% filter(wavelength>=1500)
metrics_SIR_canopy <- redr::calc_metrics(obs=SIR_canopy%>%pull(Reflectance_median),
                                              sim=SIR_canopy%>%pull(median))

table3 <- 
  rbind(
    signif(
      rbind(cbind(rbind(metrics_PAR_tree,metrics_PAR_liana),rbind(metrics_NIR_tree,metrics_NIR_liana),rbind(metrics_SIR_tree,metrics_SIR_liana)),
            rbind(c(metrics_delta_PAR,metrics_delta_NIR,metrics_delta_SIR),
                  c(metrics_delta_PAR_m,metrics_delta_NIR_m,metrics_delta_SIR_m))),3),
    signif(
      rbind(cbind(rbind(metrics_PAR_tree_canopy,metrics_PAR_liana_canopy),rbind(metrics_NIR_tree_canopy,metrics_NIR_liana_canopy),rbind(metrics_SIR_tree_canopy,metrics_SIR_liana_canopy)),
            rbind(c(metrics_delta_PAR_canopy,metrics_delta_NIR_canopy,metrics_delta_SIR_canopy),
                  c(metrics_delta_PAR_canopy_m,metrics_delta_NIR_canopy_m,metrics_delta_SIR_canopy_m))),3))

Table3 <- signif(table3,digits = 2)
View(Table3)
write.csv(x= Table3,file = "~/data/RTM/table3.csv")



