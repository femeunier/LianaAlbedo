rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggridges)

directory <- "/home/carya/data/RTM/"

files <- file.path(directory, c("Foster_edr_LAI.rds", "Marvin_edr_LAI2.rds","Kalacska_edr_LAI2.rds","Sanchez_edr_LAI2.rds"))
refs <- c("Foster", "Marvin","Kalacska","Sanchez")

N <- 1000

dis2find <- c('b1Bl','b2Bl','Nlayers','Cab','Car','Cw','Cm','orient.factor','clumping.factor')
posterior_sample_all <- data.frame()

for (ifile in seq(1,length(files))){
  file <- files[ifile]
  posteriors <- readRDS(file)
  posterior_sample <- as.data.frame(BayesianTools::getSample(posteriors,numSamples=10000,parametersOnly=FALSE)) %>% dplyr::select(-c("Lprior")) %>% top_n(N) %>% dplyr::select(seq(5,5+2*length(dis2find)-1))
  colnames(posterior_sample) <- c(paste0(dis2find,"_tree"),paste0(dis2find,"_liana"))
  posteriors_long <- posterior_sample %>% pivot_longer(everything()) %>% mutate(pft = sub(".*\\_", "", name),
  name = sub("\\_.*", "", name)) %>% rename(Param = name)
  posterior_sample_all <- rbind(posterior_sample_all,
  posteriors_long %>% mutate(ref = refs[ifile]))
}

ggplot(data = posterior_sample_all,aes(x = value, y = ref, fill = pft)) +
  geom_density_ridges(alpha= 0.5,scale=0.95,rel_min_height=0.01) +
  geom_hline(yintercept = seq(1,length(files)))+
  labs(y = "",
       x = "",
       fill = "Growth form") +
  facet_wrap(Param ~ ., scales = "free_x", nrow = 1) +
  theme_ridges(font_size = 16) +
  theme_bw()
