rm(list = ls())

library(dplyr)
library(ggplot2)
library(pracma)
library(Hmisc)
library(cowplot)

Values <- c("#1E64C8","#137300","#000000")
names(Values) <- c("Liana","Tree","Soil")  
Values["Soil"] <- "dimgray"

load(file = "~/data/RTM/Inverse_leaf_spectrum_SA.Rdata")

UA.step1 <- model_sensitivities_all %>% rename(param = variable) %>% dplyr::select(-par.var) %>% rename(par.var = CV) %>%
 dplyr::select(ref,pft,param,par.var,OP_variable) %>% 
  mutate(ref = case_when(
    ref == "Sanchez_PNM" ~ "Sanchez (PNM)",
    ref == "Sanchez_FTS" ~ "Sanchez (FTS)",
    ref == "Castro_PNM" ~ "Castro (PNM)",
    ref == "Castro_FTS" ~ "Castro (FTS)",
    ref ==  "Guzman" ~  "Guzman",
    ref ==  "Kalacska" ~  "Kalacska")) %>% filter(OP_variable %in% c("PAR","nir","sir"))

load(file = "~/data/RTM/Inverse_canopy_spectrum_SA.Rdata")
df_PFT <- data.frame(names = c("Liana_optical","Tree_optical"),PFTnum = c(17,3),Col = c("#1E64C8","#137300"))
df_PFT <- df_PFT %>% arrange(PFTnum)%>% mutate(Col = c("#1E64C8","#137300"))
df_PFT <- df_PFT %>% arrange(names)

UA.step2 <- model_sensitivities_all  %>% dplyr::select(-par.var) %>% rename(par.var = CV) %>%
  dplyr::select(ref,OP_variable,pft,par.var,variable) %>%
  rename(param = variable) %>% mutate(OP_variable = as.character(OP_variable)) %>% mutate(OP_variable = case_when(
    OP_variable == "nir_liana" ~ "NIR_Liana",
    OP_variable == "nir_tree" ~ "NIR_Tree",
    OP_variable == "sir_liana" ~ "SIR_Liana",
    OP_variable == "sir_tree" ~ "SIR_Tree",
    OP_variable == "PAR_liana" ~ "PAR_Liana",
    OP_variable == "PAR_tree" ~ "PAR_Tree")) %>% filter(OP_variable %in% c("NIR_Liana","NIR_Tree","PAR_Liana","PAR_Tree","SIR_Liana","SIR_Tree"))

# rescale
UA.step2 <- UA.step2 %>% group_by(OP_variable,ref) %>% filter(param != "soil_moisture") %>% mutate(
  par.var = par.var/sum(par.var))

#######################################################################################################
# Subplot A

subplot.a <-
  UA.step1 %>% group_by(OP_variable,param,pft) %>% mutate(par.var = 100*par.var) %>%
  summarise(par.var_m = mean(par.var),
            par.var_min = min(par.var),
            par.var_max = max(par.var),
            par.var_sd = sd(par.var),
            par.var_se = sd(par.var)/sqrt(length(par.var))) %>% filter(param != "ssigma") %>% mutate(
              pft = case_when(
                pft == "Liana_optical" ~ "Liana",
                pft == "Tree_optical"  ~ "Tree")) %>% mutate(OP_var = case_when(
                  OP_variable == "nir" ~ "Near-Infrared \n[801-1400nm]",
                  OP_variable == "sir" ~ "Shortwave Infrared \n[1501-2500nm]",
                  OP_variable == "PAR" ~ "PAR \n[400-700nm]")) %>% mutate(col = 1)

subplot.a$OP_var2 <- factor(subplot.a$OP_var, levels=c("PAR \n[400-700nm]","Near-Infrared \n[801-1400nm]","Shortwave Infrared \n[1501-2500nm]"))


subplotA <-
  ggplot(data = subplot.a %>% filter(OP_var == "PAR \n[400-700nm]"),
         aes(x = param,
             y = par.var_m,
             ymin = par.var_min,
             ymax = par.var_max,
             color = as.factor(pft),
             fill = as.factor(pft))) +
  geom_errorbar(position=position_dodge(width=0.8),width = 0.5) +
  geom_bar(stat = "identity",
           position = position_dodge(), width = 0.7, alpha = 0.5) +
  coord_flip() +
  theme_bw() +     
  scale_y_continuous(limits = c(0,110),
                     expand = c(0.001,0.001)) +
  scale_fill_manual(values =  as.character(df_PFT$Col)) +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  labs(y = "Posterior CV [%]",
       x = "",
       fill = "Growth form",
       colour = "Growth form") +
  theme(text = element_text(size=18),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        strip.background.x = element_rect(colour=NA, fill= NA),
        strip.text.x = element_text(colour = NA)) + guides(colour = FALSE, fill = FALSE)


#######################################################################################################
# Subplot B

subplot.b <- UA.step2 %>% group_by(OP_variable,pft,param) %>% mutate(par.var = 100*par.var) %>%
  summarise(par.var_m = mean(par.var),
            par.var_min = min(par.var),
            par.var_max = max(par.var),
            par.var_sd = sd(par.var),
            par.var_se = sd(par.var)/sqrt(length(par.var))) %>% ungroup() %>% filter(param != "ssigma") %>% mutate(
              pft = case_when(
                pft == "Liana_optical" ~ "Liana",
                pft == "Tree_optical"  ~ "Tree",
                pft == "soil" ~ "Soil")) %>% mutate(OP_var = sub("\\_.*", "", OP_variable),
                                                    scenario = sub(".*\\_", "", OP_variable))



subplot.b2 <- subplot.b %>% mutate(param = as.character(param)) %>%
  mutate(
    scenario = case_when(
      scenario == "Liana" ~ "Liana-rich patches",
      scenario == "Tree" ~ "Liana-free patches"
    )
  ) %>% mutate(
    param = case_when(
      param == "soil_moisture" ~ "θsoil",
      param == "orient_factor" ~ "ω",
      param == "carotenoids" ~ "Car",
      param == "Cab" ~ "Cab",
      param == "Cm" ~ "Cm",
      param == "Cw" ~ "Cw",
      param == "N" ~ "N",
      param == "b1Bl" ~ "b1Bl",
      param == "b2Bl" ~ "b2Bl",
      param == "clumping_factor" ~ "Ω",
      TRUE ~ param)) %>% mutate(
        order = case_when(
          param == "θsoil" ~ 10,
          param == "b2Bl" ~ 9,
          param == "b1Bl" ~ 8,
          param == "Ω" ~ 7,
          param == "ω" ~ 3,
          param == "Car" ~ 2,
          param == "Cab" ~ 4,
          param == "Cm" ~ 6,
          param == "Cw" ~ 5,
          param == "N" ~ 1
        )
      ) %>% mutate(
        OP_var = case_when(
          OP_var == "NIR" ~ "Near-Infrared \n[801-1400nm]",
          OP_var == "SIR" ~ "Shortwave Infrared \n[1501-2500nm]",
          OP_var == "PAR" ~ "PAR \n[400-700nm]"
        )
      )

subplot.b2$OP_var2 <- factor(subplot.b2$OP_var, levels=c("PAR \n[400-700nm]","Near-Infrared \n[801-1400nm]","Shortwave Infrared \n[1501-2500nm]"))

subplotB_all <- 
  ggplot(data = subplot.b2  %>% filter(OP_var == "PAR \n[400-700nm]"),
         aes(x = order,
             y = par.var_m,
             ymin = par.var_min,
             ymax = par.var_max,
             color = as.factor(pft),
             fill = as.factor(pft))) +
  geom_errorbar(position=position_dodge(width=0.8),width = 0.5) +
  geom_bar(stat = "identity",
           position = position_dodge(), width = 0.7,alpha = 0.5) +
  coord_flip()  +
  scale_fill_manual(values = Values) +
  scale_color_manual(values = Values) +
  scale_y_continuous(limits = c(-40,110),
                     expand = c(0.001,0.001)) +
  scale_x_continuous(breaks = seq(1:10),
                     labels=c("Nlayers","Car","ω","Cab","Cw","Cm","Ω","b1Bl","b2Bl","θsoil")) +
  theme_bw() + 
  labs(y = "Posterior CV [%]",
       x = "",
       fill = "Source",
       colour = "Source") +
  theme(text = element_text(size=18),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        panel.spacing.x = unit(1.5, "lines")) + guides(color = FALSE, fill = FALSE)

subplotB_all

plot_grid(subplotA,subplotB_all,
          align = c("h"),
          nrow = 1,rel_widths = c(1,1))

ggsave(plot = last_plot(),
       dpi = 300,
       width = 30,
       height = 20,
       units = "cm",
       file = "~/data/RTM/Figure4_CV.png")
