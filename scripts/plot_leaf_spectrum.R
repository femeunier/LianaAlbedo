
rm(list=ls())

alpha = 0.05

All_leaf_spectra <-
  rbind(readRDS(file = "~/data/RTM/Figure1_Guzman.rds") %>% mutate(ref = "Guzman"),
        readRDS(file = "~/data/RTM/Figure1_kalacska.rds") %>% mutate(ref = "Kalacska"),
        readRDS(file = "~/data/RTM/Figures4and5_castro_FTS.rds") %>% mutate(ref = "Castro_FTS"),
        readRDS(file = "~/data/RTM/Figures4and5_castro_PNM.rds") %>% mutate(ref = "Castro_PNM"),
        readRDS(file = "~/data/RTM/Figure6_sanchez2009_PNM_all.rds") %>% mutate(ref = "Sanchez_PNM"),
        readRDS(file = "~/data/RTM/Figure6_sanchez2009_FTS_all.rds") %>% mutate(ref = "Sanchez_FTS"))

ref_alls <- unique(All_leaf_spectra$ref)

curves_all <- c()
All_leaf_spectra_interp <- data.frame()
for (cref in (ref_alls)){
  pos_ref <- which(All_leaf_spectra$ref == cref)
  pfts <- as.character(unique(All_leaf_spectra[pos_ref,"pft"]))
  for (cpft in pfts){
    cdata <- All_leaf_spectra %>% filter(ref == cref, pft == cpft)
    x <- cdata$wavelength
    y <- cdata$Reflectance
    xi <- (as.integer(min(x))+1):(as.integer(max(x))-1)
    yi <- interp1(x,y,xi)
    All_leaf_spectra_interp <- rbind(All_leaf_spectra_interp,
                                     data.frame(wavelength = xi,
                                                Reflectance = yi,
                                                pft = cpft,ref = cref))
  }
}


All_leaf_spectra2 <- All_leaf_spectra_interp %>% ungroup() %>% mutate(pft = as.character(pft),
                                                                      ref = as.character(ref)) %>% 
  mutate(pft = case_when(
  pft == "Tree_optical" ~ 'Tree',
  pft == "Liana_optical" ~ 'Liana'
))  %>% mutate(ref = case_when(
  ref == "Castro_FTS" ~ 'Castro (FTS)',
  ref == "Castro_PNM" ~ 'Castro (PNM)',
  ref == "Sanchez_FTS" ~ 'Sanchez (FTS)',
  ref == "Sanchez_PNM" ~ 'Sanchez (PNM)',
  ref == "Guzman" ~ 'Guzmán',
  TRUE ~ ref
)) %>% group_by(ref,pft)


All_leaf_spectra_interp_diff <- All_leaf_spectra2 %>% group_by(ref,wavelength) %>%
  summarise(diff = ifelse(length(which(pft == "Tree"))>0,Reflectance[pft == "Liana"] - Reflectance[pft == "Tree"],NA))
  

ggplot() +
  # geom_ribbon(alpha = 0.5,linetype = 0) +
  geom_line(data=All_leaf_spectra2,
            aes(x = wavelength,
                y = Reflectance,
                fill = pft,
                color = pft)) +
  geom_line(data = All_leaf_spectra_interp_diff,
            aes(x = wavelength,y = diff),color = "black") +
  facet_wrap(ref ~ .,scales = "free_x") +
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]",
       colour = "Growth form") +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
