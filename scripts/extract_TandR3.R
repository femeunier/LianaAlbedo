rm(list = ls())

library(BayesianTools)
library(reshape2)
library(PEcAnRTM)
library(ggplot2)
library(tidyr)
library(dplyr)

Colors <- c("#137300","#1E64C8")

df_PFT <- data.frame(names = c("Liana_optical","Tree_optical"),PFTnum = c(17,3),Col = Colors)
df_PFT <- df_PFT %>% arrange(PFTnum)%>% mutate(Col = Colors)
df_PFT <- df_PFT %>% arrange(names)

posterior_all <- readRDS(file = "~/data/RTM/params_sanchez.RDS")
posterior_all2 <- posterior_all %>% filter(site == "FTS")

data.study <- posterior_all2 %>% mutate(pft = case_when(GF == "Liana" ~ "Liana_optical",
                                                                                 GF == "Tree" ~ "Tree_optical"),
                                                                 sample = sort(rep(seq(1,nrow(posterior_all2)/5),5))) %>% rename(Param = params)

data.study <- data.study %>% group_by(pft) %>% mutate(sample = 1 + sample - min(sample)) %>% ungroup()

Nensemble = 250

# ggplot(data = posterior_all,aes(x = value, y = site, fill = pft)) +
#   geom_density_ridges(alpha= 0.5) +
#   theme_bw() +
#   scale_color_manual(values = as.character(df_PFT$Col)) +
#   scale_fill_manual(values = as.character(df_PFT$Col)) +
#   facet_wrap(Param ~ .,scales = "free",nrow = 1) +
#   theme(axis.text = element_text(size=12),
#         panel.spacing.x = unit(1.5, "lines"))

samples <- 1:Nensemble
# Liana 
data <- data.study %>% filter(sample %in% samples) %>% dplyr::select(sample,Param,value,pft) %>% filter(pft == "Liana_optical") %>%
  pivot_wider(names_from = Param, values_from = value) %>% dplyr::select(Nlayers,Cab,Car,Cw,Cm)

# run Prospect5
df_spectrum_L = data.frame()
for (i in seq(1,min(length(unique(samples)),Nensemble))){
  current_parameter_set <- data[i,]
  if(!any(is.na(current_parameter_set))){
  current_model_output <- PEcAnRTM::prospect(current_parameter_set, version = "5")
  df_spectrum_L <- rbind(df_spectrum_L,
                         data.frame(num = i,wv = 400:2500,reflectance = current_model_output[,1],transmittance = current_model_output[,2]))
  }
}

df_spectrum_L <- df_spectrum_L %>% mutate(pft = "Liana")

# Tree
data <- data.study %>% filter(sample %in% samples) %>% dplyr::select(sample,Param,value,pft) %>% filter( pft == "Tree_optical") %>%
  pivot_wider(names_from = Param, values_from = value) %>% dplyr::select(Nlayers,Cab,Car,Cw,Cm)

# run Prospect5
df_spectrum_T = data.frame()
for (i in seq(1,min(length(unique(samples)),Nensemble))){
  current_parameter_set <- data[i,]
  if(!any(is.na(current_parameter_set))){
  current_model_output <- PEcAnRTM::prospect(current_parameter_set, version = "5")
  df_spectrum_T <- rbind(df_spectrum_T,
                         data.frame(num = i,wv = 400:2500,reflectance = current_model_output[,1],transmittance = current_model_output[,2]))
  }
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
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]") +
  theme_bw()

saveRDS(object = spectra,file = file.path("~/data/RTM/",paste0("spectrum_","Sanchez_FTS",".RDS")))

spectra %>% group_by(pft,band) %>% summarise(r_m = mean(reflectance),
                                             t_m = mean(transmittance))

