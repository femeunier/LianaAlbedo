rm(list = ls())

library(ggplot2)
library(ggridges)
library(dplyr)
library(purrr)
library(tidyr)
library(cowplot)
library(patchwork)
library(pracma)

data.all <- readRDS("./scripts/Sanchez_allfits2.RDS")

GFs <- names(data.all)
Prospect_param_names <- c("Nlayers","Cab","Car","Cw","Cm")

# Add the full spectra
df.data.all_sim <- data.frame()
Nsamples <- 100

params.sanchez <- data.frame()

compt <- 1
for (pft in seq(names(data.all))){
  for (species in seq(names(data.all[[pft]]))){
    for (ind in seq(length(data.all[[pft]][[species]]))) {
      
      data.temp <- data.all[[pft]][[species]][[ind]]
      
      if (!is.null(data.temp)){
        print(compt)
        
        site <- data.temp[["site"]]
        samples <- data.temp[["samples"]]
        sampling <- as.data.frame(BayesianTools::getSample(samples,numSamples = 1000, parametersOnly=F)) %>% 
          dplyr::select(-c("Lprior")) %>% arrange(desc(Lposterior)) %>% top_n(n = Nsamples) %>% dplyr::select(c("par 1","par 2","par 3","par 4","par 5"))
        colnames(sampling) <- Prospect_param_names
        
        
        params.sanchez <- rbind(params.sanchez,
                                sampling %>% pivot_longer(all_of(Prospect_param_names)) %>% rename(Param = name) %>%
                                  mutate(pft = case_when(pft == 1 ~ "Liana_optical",
                                                         pft == 2 ~ "Tree_optical"),
                                         ref = paste0("Sanchez_",site))
        )
        
        compt <- compt + 1                          
      }
    }
  }
}

load(file = "~/data/RTM/Inverse_leaf_spectrum.Rdata")

# Combine reference with Sanchez
raw_data <- Pall %>% dplyr::filter(Param!= "ssigma") %>% ungroup() %>% dplyr::select(-c(sample,rel_value))
raw_data_all <- rbind(raw_data %>%  dplyr::filter(!(ref %in% c("Sanchez_PNM","Sanchez_FTS"))),
                      params.sanchez)

# Combine with ED-RTM calibration
N <- 1000
dis2find <- c('b1Bl','b2Bl','Nlayers','Cab','Car',
              'Cw','Cm','orient.factor','clumping.factor')

directory <- "/home/carya/data/RTM/"

files <- file.path(directory, c("Foster_edr_LAI.rds", "Marvin_edr_LAI.rds","Kalacska_edr_LAI2.rds","Sanchez_edr_LAI2.rds"))
refs <- c("Foster", "Marvin","Kalacska_RTM","Sanchez")

posterior_sample_all <- data.frame()
for (ifile in seq(1, length(files))) {
  file <- files[ifile]
  posteriors <- readRDS(file)
  posterior_sample <-
    as.data.frame(BayesianTools::getSample(
      posteriors,
      numSamples = 10000,
      parametersOnly = FALSE
    )) %>% dplyr::select(-c("Lprior")) %>% top_n(N) %>% dplyr::select(seq(5, 5 +
                                                                            2 * length(dis2find) - 1))
  colnames(posterior_sample) <-
    c(paste0(dis2find, "_tree"), paste0(dis2find, "_liana"))
  posteriors_long <-
    posterior_sample %>% pivot_longer(everything()) %>% mutate(pft = sub(".*\\_", "", name),
                                                               name = sub("\\_.*", "", name)) %>% rename(Param = name)
  posterior_sample_all <- rbind(posterior_sample_all,
                                posteriors_long %>% mutate(ref = refs[ifile]))
}



All_parameters <- rbind(raw_data_all,
                        posterior_sample_all %>% mutate(pft = case_when(pft == "liana" ~ "Liana_optical",
                                                                        pft == "tree" ~ "Tree_optical")))


saveRDS(file = "/home/carya/data/RTM/All_posteriors.RDS",All_parameters)


ggplot(data = All_parameters,
       aes(x = value, y = as.factor(pft), fill = as.factor(pft))) +
  stat_density_ridges(alpha= 0.6,quantile_lines = TRUE, quantiles = c(0.025,0.975)) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  facet_wrap(~Param,scales = "free") +
  labs(x = "",y = "") + 
  theme_bw() + theme(text = element_text(size = 18)) + guides(fill = FALSE)
