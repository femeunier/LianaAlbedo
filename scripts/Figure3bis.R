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


##############################################################################################################################

list.plot <- list()
bandwidths <- c(5,2,4,0.0005,0.1)
limits_M <- c(110,50,105,0.021,3.25)
limits_m <- c(0,0,0,0,1)
Max <- c()
All_parameters_Prospect <- All_parameters %>% filter(Param %in% c("Nlayers","Cab","Car","Cw","Cm")) 
All_parameters_Prospect <- All_parameters_Prospect %>% mutate(id = rep(1:(nrow(All_parameters_Prospect)/5),5))

All_parameters_Prospect2 <- All_parameters_Prospect %>% group_by(Param,ref) %>% filter(value <= (mean(value) + 1.5*IQR(value)) & value >= (mean(value) - 1.5*IQR(value)))   %>% 
  arrange(Param)  %>% mutate(id = 1:length(value)) %>% ungroup()

All_parameters_Prospect2 <- All_parameters_Prospect2 %>% pivot_wider(names_from = Param,
                                                                     values_from = value) %>%
  mutate(Cw = (1-(1/(1+Cw/Cm)))*100) %>% pivot_longer(c("Nlayers","Cab","Car","Cw","Cm")) %>% mutate(Param = name) %>%
  dplyr::select(c(pft,ref,value,Param))

All_parameters_Prospect2 %>% filter(ref %in% c("Guzman","Kalacska","Castro_PNM","Castro_FTS","Sanchez_PNM","Sanchez_FTS")) %>% 
  group_by(ref,Param) %>% summarise(m = mean(value[pft == "Liana_optical"],na.rm = TRUE)- mean(value[pft == "Tree_optical"],na.rm = TRUE)) %>% arrange(Param) %>%
  filter(Param == "Car")


All_parameters_Prospect2 %>% filter(ref %in% c("Guzman","Kalacska","Castro_PNM","Castro_FTS","Sanchez_PNM","Sanchez_FTS")) %>% 
  group_by(ref,Param) %>% summarise(m = mean(value[pft == "Liana_optical"],na.rm = TRUE)- mean(value[pft == "Tree_optical"],na.rm = TRUE)) %>% arrange(Param) %>%
  filter(ref %in% c("Sanchez_PNM","Sanchez_FTS"))

All_parameters_Prospect2 %>% filter(ref %in% c("Guzman","Kalacska","Castro_PNM","Castro_FTS","Sanchez_PNM","Sanchez_FTS")) %>% 
  group_by(ref,Param) %>% summarise(m = mean(value[pft == "Liana_optical"],na.rm = TRUE)- mean(value[pft == "Tree_optical"],na.rm = TRUE)) %>% group_by(Param) %>%
  summarise(m = mean(m))

# Method 2
All_parameters_Prospect2 %>% filter(ref %in% c("Guzman","Kalacska","Castro_PNM","Castro_FTS","Sanchez_PNM","Sanchez_FTS")) %>%
  group_by(Param) %>% summarise(m = mean(value[pft == "Liana_optical"],na.rm = TRUE) - mean(value[pft == "Tree_optical"],na.rm = TRUE)) %>% 
  filter(Param == "Car")


Params <- c("Cab","Car","Cw","Cm","Nlayers")

All_parameters_Prospect2 <- All_parameters_Prospect2 %>% filter(ref %in% c("Castro_FTS","Castro_PNM","Guzman","Sanchez_FTS","Sanchez_PNM","Kalacska"))
All_parameters_Prospect2$ref <- as.factor(All_parameters_Prospect2$ref)
All_parameters_Prospect2$Param <- as.factor(All_parameters_Prospect2$Param)
levels(All_parameters_Prospect2$ref) <- c("Castro (FTS)","Castro (PNM)","Gùzman","Kalacska","Sánchez (FTS)","Sánchez (PNM)")
levels(All_parameters_Prospect2$Param)[4] <- c("Water content")

Params <- c("Cab","Car","Water content","Cm","Nlayers")

for(i in seq(2,length(Params))){
   if (i == 4){ 
    list.plot[[i-1]] <- ggplot(data = All_parameters_Prospect2 %>% filter(Param == Params[i]),
                             aes(x = value, y = as.factor(ref), fill = pft)) +
      geom_density_ridges(alpha= 0.5,bandwidth = bandwidths[i],scale=0.95,rel_min_height=0.01) +
      geom_hline(yintercept = seq(1,length(unique(All_parameters_Prospect2$ref)))) +
      scale_color_manual(values = as.character(df_PFT$Col)) +
      scale_fill_manual(values = as.character(df_PFT$Col)) +
      labs(y="",x = "",fill = "Growth form") +
      theme_ridges(font_size = 16) + 
      facet_wrap(Param ~ .,scales = "free_x",nrow = 1) +
      scale_x_continuous(limits=c(limits_m[i],limits_M[i]),breaks = seq(limits_m[i],limits_M[i],length.out = 4)) +
      theme_bw() +
      theme(text = element_text(size=20),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.margin = unit(c(0, -0., 0, 1), "cm")) + guides(fill = FALSE)
   } else {
     
     list.plot[[i-1]] <- ggplot(data = All_parameters_Prospect2 %>% filter(Param == Params[i]),
                                aes(x = value, y = as.factor(ref), fill = pft)) +
       geom_density_ridges(alpha= 0.5,bandwidth = bandwidths[i],scale=0.95,rel_min_height=0.01) +
       geom_hline(yintercept = seq(1,length(unique(All_parameters_Prospect2$ref)))) +
       scale_color_manual(values = as.character(df_PFT$Col)) +
       scale_fill_manual(values = as.character(df_PFT$Col)) +
       labs(y="",x = "",fill = "Growth form") +
       theme_ridges(font_size = 16) + 
       facet_wrap(Param ~ .,scales = "free_x",nrow = 1) +
       scale_x_continuous(limits=c(limits_m[i],limits_M[i])) +
       theme_bw() +
       theme(text = element_text(size=20),
             axis.text.y=element_blank(),
             axis.title.y=element_blank(),
             axis.ticks.y=element_blank(),
             plot.margin = unit(c(0, -0., 0, 1), "cm")) + guides(fill = FALSE)
   }
}

A <- plot_grid(plotlist = list.plot,align = "hv", nrow = 1,rel_widths = c(1,1,1,1,1))
B <-  ggplot(data = All_parameters_Prospect2 %>% filter(Param == Params[1]),
             aes(x = value, y = ref, fill = pft)) +
  geom_density_ridges(alpha= 0.5,bandwidth = bandwidths[1],scale=0.95,rel_min_height=0.01) +
  geom_hline(yintercept = seq(1,length(unique(All_parameters_Prospect$ref)))) +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_x_continuous(limits=c(limits_m[1],limits_M[1])) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  labs(y="",
       x = "",
       fill = "Growth form") +
  facet_wrap(Param ~ .,scales = "free_x",nrow = 1) +
  theme_ridges(font_size = 16) + 
  theme_bw() +
  theme(text = element_text(size=20),
        plot.margin = unit(c(0, -0., 0, 0.5), "cm")) + guides(fill = FALSE)

plot_grid(B,A,rel_widths = c(0.4,1))
ggsave(plot = last_plot(),dpi = 300, width = 40,height = 20, filename = file.path("./Figures","Posterior_distributions.png"),units = "cm")


###############################################################################################

All_parameters_EDRTM <- All_parameters %>% filter(Param %in% c("b1Bl","b2Bl","orient.factor","clumping.factor"))
All_parameters_EDRTM <- All_parameters_EDRTM %>% mutate(Param = as.factor(Param),
                                                        ref = as.factor(ref))

levels(All_parameters_EDRTM$Param)[c(9,10)] <- c("ω","Ω")
levels(All_parameters_EDRTM$ref) <- c("Foster","Kalacska","Marvin","Sanchez")

ggplot(data = All_parameters_EDRTM,
       aes(x = value, y = ref, fill = pft)) +
  geom_density_ridges(alpha= 0.5,scale=0.95,rel_min_height=0.01) +
  geom_hline(yintercept = seq(1,length(unique(All_parameters_EDRTM$ref)))) +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  labs(y="",
       x = "",
       fill = "Growth form") +
  facet_wrap(Param ~ .,scales = "free_x",nrow = 1) +
  theme_ridges(font_size = 16) + 
  theme_bw() +
  theme(text = element_text(size=20)) + guides(fill = FALSE)
          
ggsave(plot = last_plot(),dpi = 300, width = 40,height = 20, filename = file.path("./Figures","Posterior_distributions_EDRTM.png"),units = "cm")

saveRDS(file = file.path("~/data/RTM/All_parameters.RDS"), object = All_parameters)
