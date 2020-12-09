rm(list = ls())

library(BayesianTools)
library(reshape2)
library(PEcAnRTM)
library(ggplot2)
library(tidyr)
library(dplyr)

load(file = "~/data/RTM/Inverse_leaf_spectrum_sanchez.Rdata")
rm(list=setdiff(ls(), "posterior_all"))

Select = "Guzman"     #  "Guzman"      "Kalacska"    "Sanchez_PNM" "Sanchez_FTS" "Castro_FTS" "Castro_PNM"

data.study <- posterior_all %>% filter(site == Select)

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
  scale_color_manual(values = c("red","blue")) +
  scale_fill_manual(values = c("red","blue")) +
  theme_bw()

saveRDS(object = spectra,file = file.path("~/data/RTM/",paste0("spectrum_",Select,".RDS")))


spectra %>% group_by(pft,band) %>% summarise(r_m = mean(reflectance),
                                             t_m = mean(transmittance))

