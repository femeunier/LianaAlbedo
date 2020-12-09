rm(list = ls())

library(BayesianTools)
library(reshape2)
library(PEcAnRTM)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggridges)

Colors <- c("#137300","#1E64C8")

df_PFT <- data.frame(names = c("Liana_optical","Tree_optical"),PFTnum = c(17,3),Col = Colors)
df_PFT <- df_PFT %>% arrange(PFTnum)%>% mutate(Col = Colors)
df_PFT <- df_PFT %>% arrange(names)

posteriors <- readRDS(file = "~/data/RTM/Marvin_edr_LAI.rds")
posteriors <- readRDS(file = "~/data/RTM/Foster_edr_LAI.rds")
posteriors <- readRDS(file = "~/data/RTM/Kalacska_edr_LAI2.rds")
posteriors <- readRDS(file = "~/data/RTM/Sanchez_edr_LAI2.rds")

Select = "Sanchez"

data_all <- as.data.frame(BayesianTools::getSample(posteriors,numSamples=1000,parametersOnly = FALSE)) %>% dplyr::select(-c("Lprior")) %>% top_n(250) %>% dplyr::select(-c(21,22))
# data_all <- BayesianTools::getSample(posteriors[[3]][[1]],numSamples=250,end=11112)
data_all <- data_all[,-seq(1,4)]
# data_all <- data_all[data_all[,19]>600,1:(ncol(data_all) -3)]

colnames(data_all) <- rep( c('b1Bl_large','b2Bl_large','Nlayers',
                'Cab','Car','Cw','Cm','orient_factor','clumping_factor'),2)
data_tree <- as.data.frame(data_all[,seq(1,9)]) %>% mutate(sample = 1:nrow(data_all)) %>% dplyr::select(c("sample","Nlayers","Cab","Car","Cw","Cm")) %>% pivot_longer(c("Nlayers","Cab","Car","Cw","Cm")) %>%
  mutate(pft = "Tree_optical")
data_liana <- as.data.frame(data_all[,seq(11,18)]) %>% mutate(sample = 1:nrow(data_all)) %>% dplyr::select(c("sample","Nlayers","Cab","Car","Cw","Cm")) %>% pivot_longer(c("Nlayers","Cab","Car","Cw","Cm")) %>%
  mutate(pft = "Liana_optical")

data.study <- rbind(data_liana,data_tree) %>% arrange(name) %>% rename(Param = name)

Nensemble = 250

# ggplot(data = posterior_all,aes(x = value, y = site, fill = pft)) +
#   geom_density_ridges(alpha= 0.5) +
#   theme_bw() +
#   scale_color_manual(values = as.character(df_PFT$Col)) +
#   scale_fill_manual(values = as.character(df_PFT$Col)) +
#   facet_wrap(Param ~ .,scales = "free",nrow = 1) +
#   theme(axis.text = element_text(size=12),
#         panel.spacing.x = unit(1.5, "lines"))

samples <- pmax(1,as.integer(runif(Nensemble, min = 1, max = max(data.study$sample))))
# Liana 
data <- data.study %>% filter(sample %in% samples) %>% dplyr::select(sample,Param,value,pft) %>% filter( pft == "Liana_optical") %>%
  pivot_wider(names_from = Param, values_from = value) %>% dplyr::select(Nlayers,Cab,Car,Cw,Cm)
data_liana <- data %>% mutate(pft = "Liana")

# run Prospect5
df_spectrum_L = data.frame()
for (i in seq(1,min(length(unique(samples)),Nensemble))){
  current_parameter_set <- data[i,]
  current_model_output <- PEcAnRTM::prospect(current_parameter_set, version = "5")
  df_spectrum_L <- rbind(df_spectrum_L,
                         data.frame(num = i,wv = 400:2500,reflectance = current_model_output[,1],transmittance = current_model_output[,2]))
}

df_spectrum_L <- df_spectrum_L %>% mutate(pft = "Liana")

# Tree
data <- data.study %>% filter(sample %in% samples) %>% dplyr::select(sample,Param,value,pft) %>% filter( pft == "Tree_optical") %>%
  pivot_wider(names_from = Param, values_from = value) %>% dplyr::select(Nlayers,Cab,Car,Cw,Cm)
data_tree <- data %>% mutate(pft = "Tree")

# run Prospect5
df_spectrum_T = data.frame()
for (i in seq(1,min(length(unique(samples)),Nensemble))){
  current_parameter_set <- data[i,]
  current_model_output <- PEcAnRTM::prospect(current_parameter_set, version = "5")
  df_spectrum_T <- rbind(df_spectrum_T,
                         data.frame(num = i,wv = 400:2500,reflectance = current_model_output[,1],transmittance = current_model_output[,2]))
}
df_spectrum_T <- df_spectrum_T %>% mutate(pft = "Tree")

spectra <- rbind(df_spectrum_T,df_spectrum_L) %>% mutate(band = case_when(wv <= 700 ~ 1,
                                                                          wv <= 2500 ~ 2))

spectra_sum <- spectra %>% group_by(pft,wv) %>% summarise(r_m = mean(reflectance),
                                                          t_m = mean(transmittance))

ggplot(spectra_sum, aes(x = wv,colour = pft)) +
  geom_line(aes(y = r_m)) +
  geom_line(aes(y = 1-t_m),linetype = 1)+
  scale_y_continuous(sec.axis = sec_axis(trans = ~ 1 - 1 * ., name = "Transmittance [-]")) +
  # geom_ribbon(aes(ymin = r_m,ymax=1-t_m,fill = pft),alpha = 0.2,linetype=0)+
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]") +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  theme_bw()


spectra %>% group_by(pft,band) %>% summarise(r_m = mean(reflectance),
                                             t_m = mean(transmittance))


saveRDS(object = spectra,file = file.path("~/data/RTM/",paste0("spectrum_",Select,".RDS")))

###########################################################################

prospect_parameters <- rbind(data_liana,
                             data_tree) %>% pivot_longer(c("Nlayers","Cab","Car","Cw","Cm"))

# ggplot(data = prospect_parameters,aes(x = value, y = 0, fill = pft)) +
#   geom_density_ridges(alpha= 0.5) +
#   facet_wrap(~ name,scales = "free",nrow = 1) +
#   scale_color_manual(values = c("darkblue","darkgreen")) +
#   scale_fill_manual(values = c("darkblue","darkgreen")) +
#   theme_bw() 
