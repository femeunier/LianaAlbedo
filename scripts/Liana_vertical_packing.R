rm(list = ls())

library(ggplot2)
library(minpack.lm)
library(cowplot)
library(dplyr)
library(PEcAn.ED2)
library(tidyr)
library(ggplot2)

data.BCI <- read_css(file = "/home/carya/data/BCI/BCI_PFTs.lat9.000lon-79.000.css")

All.posteriors <- readRDS("/home/carya/data/RTM/All_posteriors.RDS")
ensemble <- All.posteriors %>% group_by(Param,pft) %>% summarise(value = mean(value))  %>% mutate(run = 1)  %>% pivot_wider(names_from = Param,
                                                                                                                            values_from = value)

SLA <- 1/(10*ensemble$Cm);b1Bl <- ensemble$b1Bl; b2Bl <- ensemble$b2Bl


data.mod <-
  data.BCI %>% mutate(pft = case_when(pft %in% c(2, 3, 4) ~ 3,
                                      pft == 17 ~ 17))

sim.patch <- list()
Npatches <- length(unique(data.mod$patch))
df.LAI <- data.frame()

href <- 61.7;  b1Ht <- 0.035; b2Ht <- 0.69 # Tree
hrefL <- 61.7;  b1HtL <- 0.11; b2HtL <- 2.5 # Liana

dbh.thrshld <- 2.5
for (ipatch in seq(1,Npatches)){
  
  site_dat <- data.mod %>% filter(patch == ipatch)
  
  dbh <- site_dat[["dbh"]]
  pft <-  as.factor(site_dat[["pft"]])
  nplant <- site_dat[["n"]]
  patch <- site_dat[["patch"]]
  
  is_liana <- (pft == 17)
  npft <- length(levels(pft))
  levels(pft) <- seq(npft,1)

  h <- pmin(hmax,(href*(1 -exp(-b1Ht*(dbh**b2Ht)))))

  h[is_liana] <- pmin(hmax,(hrefL*(1 -exp(-b1HtL*(dbh[is_liana]**b2HtL)))))
  
  tree.maxH <- max(h[!is_liana])
  h[dbh>dbh.thrshld & is_liana] <- tree.maxH + 0.5
  h[is_liana & dbh <= dbh.thrshld] <- pmin(h[is_liana & dbh <= dbh.thrshld],tree.maxH)
  
  temp <- sort(h, decreasing = TRUE, index.return = TRUE)
  h <- temp$x
  dbh <- dbh[temp$ix]
  pft <- pft[temp$ix]
  nplant <- nplant[temp$ix]
  patch <- patch[temp$ix]
  is_liana <- is_liana[temp$ix]
  
  
  sim.patch[[ipatch]] <- data.frame(DBH = dbh, H = h, PFT = as.numeric(as.vector(pft)), N = nplant, PA = patch,is.liana = is_liana) %>% 
    mutate(lai = SLA[PFT]*N*b1Bl[PFT]*(DBH**b2Bl[PFT])) %>% mutate(cumLAI = cumsum(lai)) 
  
  df.LAI <- bind_rows(list(df.LAI,
                           data.frame(i = ipatch,
                                      lai = max(sim.patch[[ipatch]][["cumLAI"]]),
                                      h = max(sim.patch[[ipatch]][["H"]]))))
}

ipatch = 46

cpatch <- bind_rows(list(data.frame(DBH = NA,H = max(sim.patch[[ipatch]][["H"]]),PFT = NA,N = NA,is_liana = FALSE,lai = 0,cumLAI = 0),
                         sim.patch[[ipatch]]))
cpatch$PFT <- as.factor(cpatch$PFT)
levels(cpatch$PFT) <- c("Liana","Tree")

ggplot(data = cpatch) +
  geom_line(aes(x = cumLAI,y = H),color = "#A9A9A9",linetype = 1) +
  geom_point(aes(x = cumLAI,y = H,color = as.factor(PFT),shape = as.factor(PFT),alpha = as.factor(PFT),size = (lai))) +
  geom_point(data = cpatch %>% filter(PFT == 1),
             aes(x = cumLAI,y = H,color = as.factor(PFT),shape = as.factor(PFT),alpha = as.factor(PFT),size = (lai))) +
  scale_color_manual(values = c("#1E64C8","#137300","white"),
                     breaks = c("Liana", "Tree", "")) +
  scale_x_continuous(limits = c(0,1.1*max(cpatch$cumLAI))) +
  scale_shape_manual(values = c(19,19)) +
  scale_alpha_manual(values = c(0.6,0.4)) +
  # scale_size_manual(values = c(2,4)) +
  labs(x = "Cumulative LAI [m²/m²]",y = "Height [m]", color = "PFT") +
  # scale_y_continuous(limits = c(25,27)) +
  theme_bw() +
  theme(text = element_text(size = 24),
        legend.position = c(0.85,0.85),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  guides(shape = FALSE,alpha = FALSE,size = FALSE)

ggsave(plot = last_plot(),"./Figures/Liana_packing.png",
       dpi = 300,width = 10,height = 7)

