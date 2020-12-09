rm(list = ls())

library(ggplot2)
library(ggridges)
library(dplyr)
library(purrr)
library(tidyr)
library(cowplot)
library(patchwork)

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
raw_data <- Pall %>% filter(Param!= "ssigma") %>% ungroup() %>% dplyr::select(-c(sample,rel_value))
raw_data_all <- rbind(raw_data %>% filter(!(ref %in% c("Sanchez_PNM","Sanchez_FTS"))),
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

All_parameters_Bl <- All_parameters %>% filter(Param %in% c("b1Bl","b2Bl"))

#################################################################################
# Data Bl
data_Bl <- readRDS("/home/carya/data/Allometry/Bl_data.RDS") %>% filter(dbh>0.1)

N = 100
All_parameters_Bl_select <- All_parameters_Bl %>% group_by(Param,pft) %>% sample_n(size = N)
b1Bl_liana <- All_parameters_Bl_select %>% filter(pft == "Liana_optical",Param == "b1Bl") %>% pull(value)
b2Bl_liana <- All_parameters_Bl_select %>% filter(pft == "Liana_optical",Param == "b2Bl") %>% pull(value)
b1Bl_tree <- All_parameters_Bl_select %>% filter(pft == "Tree_optical",Param == "b1Bl") %>% pull(value)
b2Bl_tree <- All_parameters_Bl_select %>% filter(pft == "Tree_optical",Param == "b2Bl") %>% pull(value)

dbhs <- seq(min(data_Bl %>% filter(pft == "tree") %>% pull(dbh)),
            max(data_Bl %>% filter(pft == "tree") %>% pull(dbh)),
            length.out = 1000)
dbh_liana <- seq(min(data_Bl %>% filter(pft == "liana") %>% pull(dbh)),
                 max(data_Bl %>% filter(pft == "liana") %>% pull(dbh)),
                 length.out = 1000)
df_Bl <- data.frame()
for (i in seq(N)){
  Bl_liana <- b1Bl_liana[i]*(dbh_liana**b2Bl_liana[i])
  Bl_tree <- b1Bl_tree[i]*(dbhs**b2Bl_tree[i])
  
  df_Bl <- rbind(df_Bl,
                 data.frame(dbh = dbhs,Bl_liana = Bl_liana,Bl_tree = Bl_tree,simu = i))
}

df_Bl_sum <- df_Bl %>% pivot_longer(cols = c( Bl_liana,Bl_tree)) %>% group_by(dbh,name) %>% 
  summarise(Bl = mean(value),
            Bl_low = confint(lm(formula = value ~ 1),level = 0.95)[1,1],
            Bl_high = confint(lm(formula = value ~ 1),level = 0.95)[1,2]) %>% mutate(name = case_when(sub(".*\\_", "", name) == "liana" ~ "liana",
                                                                                         sub(".*\\_", "", name) == "tree" ~ "tree"))

df_Bl_sum[df_Bl_sum$name=="liana","dbh"]<- dbh_liana
###########################################################################################

ggplot(data = df_Bl_sum,
       aes(x = dbh,ymin = Bl_low,ymax = Bl_high,y = Bl,fill = name,color = name)) +
  geom_point(data = data_Bl,aes(color = pft,ymin = Bl,ymax = Bl,fill = pft),alpha = 0.2,size = 0.5) +
  geom_ribbon(alpha = 0.1,color = NA) +
  geom_line() +
  labs(x = "DBH [cm]",y = "Leaf biomass [kg]") +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() + guides(color = FALSE,fill = FALSE)

ggsave(plot = last_plot(),dpi = 300, width = 10,height = 10, filename = file.path("./Figures","Bl.png"),units = "cm")


#################################################################################