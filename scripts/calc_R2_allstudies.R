rm(list = ls())

library(pracma)
library(dplyr)
library(ggplot2)

load(file = "~/data/RTM/Inverse_leaf_spectrum_sanchez.Rdata")

modelled <-  ensemble_posterior_all %>% select(wavelength,median,ref,pft) %>% rename(mod = median)
observed <- All_leaf_spectra %>% select(Reflectance_median,ref,pft,wavelength) %>% rename(obs = Reflectance_median) %>% mutate(wavelength = round(wavelength))

observed_extrap <- data.frame()

pfts <- unique(observed$pft)
references <- unique(observed$ref)

for (iref in seq(1,length(references))){
  for (ipft in seq(1,length(pfts))){
    data.temp <- observed %>% filter(ref == references[iref],
                                     pft == pfts[ipft])
    
    all.wls <- seq(min(data.temp$wavelength),max(data.temp$wavelength))
    all.R <- interp1(data.temp$wavelength,data.temp$obs,all.wls)
    
    observed_extrap <- bind_rows(list(observed_extrap,
                                      data.frame(wavelength = all.wls,
                                                 obs = all.R,
                                                 ref = references[iref],
                                                 pft = pfts[ipft])))
  }
}

all <- observed_extrap %>% left_join(modelled)
all2plot <- all %>% filter(wavelength %in% seq(500,2500,50)) %>% filter(!(wavelength %in% c(700)))

all2plot.lm <- bind_rows(list(all2plot %>% group_by(wavelength) %>% summarise(x = min(obs),
                                                                              y = coef(lm(mod~obs))[2]*x + coef(lm(mod~obs))[1]),
                              all2plot %>% group_by(wavelength) %>% summarise(x = max(obs),
                                                                              y = coef(lm(mod~obs))[2]*x + coef(lm(mod~obs))[1])))


all2plot$pft <- as.factor(all2plot$pft)
levels(all2plot$pft) <- c("Liana","Tree")

ggplot(data = all2plot) +
  geom_abline(slope = 1, linetype = 3) +
  geom_point(aes(y = obs,x = mod,color = pft),alpha = 0.4) +
  geom_line(data = all2plot.lm,
            aes(x = y,y = x, group = as.factor(wavelength)),linetype = 1,color = "darkgrey") +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  labs(y = "Observed reflectance [-]",x = "Modelled reflectance [-]",color = "PFT") +
  theme_bw() +
  theme(legend.position = c(0.15,0.9),
        text = element_text(size = 18))
  

ggsave(plot = last_plot(),
       dpi = 300,
       width = 18,
       height = 15,
       units = "cm",
       file = "./Figures/QoF_bands.png")


all_short <- all %>% filter(wavelength %in% seq(500,2400,100))
df_r2 <- all_short %>% group_by(wavelength) %>% summarise(r2 = summary(lm(mod ~ obs))[["r.squared"]]) %>% arrange(desc(r2)) %>% pull(r2) %>% mean()
hist(df_r2$r2)

df_r2 %>% filter(r2 < 0.8)
