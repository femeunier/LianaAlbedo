rm(list = ls())

library(pracma)
library(dplyr)
library(ggplot2)

load(file = "~/data/RTM/Inverse_canopy_spectrum_priors.Rdata")

All_canopy_spectra <- readRDS(file= "~/data/RTM/All_canopy_spectra.rds")
All_canopy_spectra <- All_canopy_spectra %>% group_by(ref,scenario,wavelength) %>% summarise(Reflectance_min = min(reflectance),
                                                                                             Reflectance_max = max(reflectance),
                                                                                             Reflectance_median = median(reflectance),
                                                                                             Reflectance_alphamin = min(reflectance),
                                                                                             Reflectance_alphamax = max(reflectance)) %>% ungroup()

subplotc <- model_ensemble_all %>% ungroup() %>% mutate(scenar = case_when(
  scenar == "low" ~'Low',
  scenar == "high" ~'High')) %>% dplyr::select(scenar,wavelength,median,ref) %>% rename(scenario = scenar,
                                                                                        mod = median) %>% mutate(wavelength = round(wavelength))

data_subplotc <- All_canopy_spectra %>% ungroup() %>% mutate(scenario = case_when(
  scenario == "low" ~'Low',
  scenario == "high" ~'High')) 

scenarios <- unique(data_subplotc$scenario)
references <- unique(data_subplotc$ref)

observed_extrap <- data.frame()

for (iref in seq(1,length(references))){
  for (ipft in seq(1,length(scenarios))){
    
    # data
    data.temp <- data_subplotc %>% filter(ref == references[iref],
                                          scenario == scenarios[ipft])
    
    all.wls <- seq(round(min(data.temp$wavelength))+1,-1+round(max(data.temp$wavelength)))
    all.R <- interp1(data.temp$wavelength,data.temp$Reflectance_median,all.wls)
    
    # model
    mod.temp <- subplotc %>% filter(ref == references[iref],
                                    scenario == scenarios[ipft])
    
    all.R.mod <- interp1(mod.temp$wavelength,mod.temp$mod,all.wls)
    
    
    observed_extrap <- bind_rows(list(observed_extrap,
                                      data.frame(wavelength = all.wls,
                                                 obs = all.R,
                                                 mod = all.R.mod,
                                                 ref = references[iref],
                                                 scenario = scenarios[ipft])))
  }
}

select.wv <- seq(500,2400,50)

all2plot <- observed_extrap %>% filter(wavelength %in% select.wv) %>% filter(! (wavelength %in% c(650,700,2350,2300,2400)))

all2plot.lm <- bind_rows(list(all2plot %>% group_by(wavelength) %>% summarise(x = min(obs),
                                                                              y = coef(lm(mod~obs))[2]*x + coef(lm(mod~obs))[1]),
                              all2plot %>% group_by(wavelength) %>% summarise(x = max(obs),
                                                                              y = coef(lm(mod~obs))[2]*x + coef(lm(mod~obs))[1])))

all2plot %>% group_by(wavelength) %>%summarise(r2 = summary(lm(mod ~ obs))[["r.squared"]]) %>% pull(r2) %>% mean()

all2plot$scenario <- as.factor(all2plot$scenario)
levels(all2plot$scenario) <- c("High","Low")

ggplot(data = all2plot) +
  geom_abline(slope = 1, linetype = 3) +
  geom_point(aes(y = obs,x = mod,color = scenario),alpha = 0.4) +
  geom_line(data = all2plot.lm,
            aes(y = x,x = y, group = as.factor(wavelength)),linetype = 1,color = "darkgrey") +
  scale_x_continuous(limits = c(0,0.6)) +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  labs(y = "Observed reflectance [-]",x = "Modelled reflectance [-]",color = "Liana infestation") +
  theme_bw() +
  theme(legend.position = c(0.2,0.9),
        text = element_text(size = 18))

ggsave(plot = last_plot(),
       dpi = 300,
       width = 18,
       height = 15,
       units = "cm",
       file = "./Figures/QoF_bands2.png")
