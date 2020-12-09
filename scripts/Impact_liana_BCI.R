# rm(list = ls())

library(ggplot2)
library(minpack.lm)
library(cowplot)
library(dplyr)

OP_BCI <-
  readRDS(file = "/home/carya/R/edr-da/data/OP_BCI2.RDS") %>% mutate(pa = 1:length(N_liana)) %>% mutate(x = LAI_liana / LAI,
                                                                                                       canopy = case_when(LAI < 3 ~ "Gap (LAI < 3)",
                                                                                                                          TRUE ~ "Dense canopy (LAI \u2265 3)"))


OP_BCI_fit <- OP_BCI %>% filter(canopy == "Dense canopy (LAI \u2265 3)")

delta_PAR_mod <- nlsLM(data = OP_BCI_fit,
                       delta_PAR ~ ymax*(1-exp(-k*x)),
                       start=list(ymax=0.01,k = 0.1), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                                            printEval = TRUE, warnOnly = TRUE))

delta_NIR_mod <- nlsLM(data = OP_BCI_fit,
                       delta_NIR ~ ymax*(1-exp(-k*x)),
                       start=list(ymax=0.08,k = 0.2), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                                            printEval = TRUE, warnOnly = TRUE))

delta_GND_mod <- nlsLM(data = OP_BCI_fit,
                       delta_PARgnd*100 ~ ymax*(1-exp(-k*x)),
                       start=list(ymax=-60,k = 0.1), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                                            printEval = TRUE, warnOnly = TRUE))


X = seq(0,max(OP_BCI_fit$x),length.out = 1000)

A <- ggplot(data = OP_BCI) +
  geom_point(aes(x = x,y = delta_NIR,shape = as.factor(canopy)),color = "black") +
  geom_point(aes(x = x,y = delta_NIR,shape = as.factor(canopy)),color = "#5ab4ac",show.legend = FALSE) +
  geom_point(aes(x = x,y = delta_PAR,shape = as.factor(canopy)),color = "#d8b365",show.legend = FALSE) +
  geom_line(data = data.frame(x = X,y = predict(delta_PAR_mod,data.frame(x = X))),
            aes(x = x, y = y)) +
  geom_line(data = data.frame(x = X,y = predict(delta_NIR_mod,data.frame(x = X))),
            aes(x = x, y = y)) +
  geom_hline(yintercept = 0,linetype = 2) +
  scale_shape_manual(values = c(19,1)) +
  labs(x = "",
       y = "Change in albedo [-]",
       shape = "") +
  theme_bw() +
  theme(text = element_text(size = 24),
        axis.text.x=element_blank(),
        legend.position = c(0.3,0.88),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))


B <- ggplot(data = OP_BCI) +
  geom_point(aes(x = x,y = delta_PARgnd*100,shape = as.factor(canopy)),color = "darkgrey") +
  geom_line(data = data.frame(x = X,y = predict(delta_GND_mod,data.frame(x = X))),
            aes(x = x, y = y)) +
  scale_shape_manual(values = c(19,1)) +
  geom_hline(yintercept = 0,linetype = 2) +
  labs(x = "Fraction of liana LAI [-]",
       y = "Change in \n understorey light [%]",
       shape = "") +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

plot_grid(A,B,nrow = 2,align = "hv")

ggsave(plot = last_plot(),"./Figures/Liana_impact.png",
       dpi = 300,width = 8,height = 12)
