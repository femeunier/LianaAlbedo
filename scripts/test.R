rm(list = ls())

library(ggplot2)
library(ggridges)
library(dplyr)

load(file = "~/data/RTM/Inverse_leaf_spectrum.Rdata")
params.sanchez <- readRDS("~/data/RTM/Sanchez_allfits2.RDS")



raw_data <- Pall %>% filter(Param!= "ssigma") %>% ungroup() %>% dplyr::select(-c(sample,rel_value))

raw_data_all <- rbind(raw_data %>% filter(!(ref %in% c("Sanchez_PNM","Sanchez_FTS"))),
                      params.sanchez %>% rename(Param = params,
                                                pft = GF,
                                                ref = site) %>% mutate(pft = case_when(pft == "Liana" ~ "Liana_optical",
                                                                                       pft == "Tree" ~ "Tree_optical"),
                                                                       ref = paste0("Sanchez_",ref)) %>% dplyr::select(-species))

ggplot(data = raw_data_all ,
       aes(x = value, y = ref, fill = pft)) +
  geom_density_ridges(alpha= 0.5) +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  labs(y="",
       x = "",
       fill = "Growth form") +
  theme_ridges(font_size = 13) + 
  facet_wrap(Param ~ .,scales = "free_x",nrow = 1) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        panel.spacing.x = unit(1.5, "lines"))
