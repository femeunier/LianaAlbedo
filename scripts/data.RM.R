rm(list = ls())

library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggridges)

file2load <- "~/R/edr-da/data/Marvin.RDS"
Marvin <- readRDS(file2load)
data.Marvin <- Marvin %>% pivot_longer(cols = c(wl427,wl482,wl547,wl604,wl660,wl723,wl824,wl914),
                                names_to = "wl",
                                values_to = "value") %>% mutate(wavelength = as.numeric(substring(wl,3)),
                                                                type = !type)
ggplot(data = data.Marvin) +
  geom_boxplot(aes(fill = type, x = (wavelength), y = value,group = interaction(type,wl))) +
  theme_bw()



data.Marvin.outliers <- data.Marvin %>% group_by(wavelength,type) %>% filter(value > quantile(value,0.25) - 1.5*(quantile(value,0.75) - quantile(value,0.25)) &
                                                                                                       value < quantile(value,0.25) + 1.5*(quantile(value,0.75) - quantile(value,0.25)))

ggplot(data = data.Marvin) +
  geom_boxplot(aes(fill = type, x = wavelength, y = value,group = interaction(type,wl))) +
  labs(y = "Reflectance",fill = "Liana cover", x = "Wavelength [nm]") +
  scale_fill_manual(values = c("#1E64C8","#137300"),labels = c("50-100%","0-50%")) +
  theme_bw()

ggsave(plot = last_plot(),
       filename = file.path(getwd(),"Figures","Marvin.png"),dpi = 300,
       height = 15,width = 20, units = "cm")

# Stats 
data.Marvin %>%  group_by(wavelength) %>% summarise(p_value = kruskal.test(formula = value ~ type)$p.value)

ggplot(data = data.Marvin.outliers) +
  geom_density_ridges(aes(fill = type, x = (value), y = type,group = interaction(type)),alpha = 0.4) +
  facet_wrap(~wl,scales = "fixed") +
  theme_bw()


# Stefan plot
file2load <- "~/R/edr-da/data/Stefan.RDS"
Stefan <- readRDS(file2load)
data.Stefan <- Stefan %>% pivot_longer(cols = c(wl427,wl482,wl547,wl604,wl660,wl723,wl824,wl914),
                                       names_to = "wl",
                                       values_to = "value") %>% mutate(wavelength = as.numeric(substring(wl,3)))

ggplot(data = data.Stefan) +
  geom_boxplot(aes(fill = type, x = (wavelength), y = value,group = interaction(type,wl))) +
  theme_bw()

data.Stefan.outliers <- data.Stefan %>% group_by(wavelength,type) %>% filter(value > quantile(value,0.25) - 1.5*(quantile(value,0.75) - quantile(value,0.25)) &
                                                                               value < quantile(value,0.25) + 1.5*(quantile(value,0.75) - quantile(value,0.25)))
  

ggplot(data = data.Stefan) +
  geom_boxplot(aes(fill = type, x = (wavelength), y = value,group = interaction(type,wl))) +
  labs(y = "Reflectance",fill = "Plots", x = "Wavelength [nm]") +
  scale_fill_manual(values = c("#1E64C8","#137300"),labels = c("Control","Removal")) +
  theme_bw()

ggsave(plot = last_plot(),
       filename = file.path(getwd(),"Figures","Removal.png"),dpi = 300,
       height = 15,width = 20, units = "cm")

# Stats 
data.Stefan %>%  group_by(wavelength) %>% filter(!is.na(value)) %>% summarise(p_value = kruskal.test(formula = value ~ as.factor(type))$p.value)
data.Stefan %>%  ungroup() %>% filter(!is.na(value)) %>% summarise(p_value = kruskal.test(formula = value ~ as.factor(type))$p.value)

ggplot(data = data.Stefan) +
  geom_density_ridges(aes(fill = type, x = (value), y = type,group = interaction(type)),alpha = 0.4) +
  facet_wrap(~wl,scales = "fixed") +
  theme_bw()




