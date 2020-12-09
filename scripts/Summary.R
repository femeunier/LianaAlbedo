rm(list = ls())

library(PEcAn.ED2)
library(dplyr)
library(tidyr)
library(redr)

Nensemble = 1

WL = seq(400, 2500)

# Filters
par <- WL>=400 & WL <= 700
green <- WL>=520 & WL <= 600
red <- WL>=550 & WL <= 700
VNIR <- WL>=750 & WL <= 900
IR <- WL>=700 & WL <= 2500 
NIR <- WL>=800 & WL <= 1400
SIR <- WL>=1500 & WL <= 2500
SW <- WL>=400 & WL <= 2500

# read data
All.posteriors <- readRDS("/home/carya/data/RTM/All_posteriors.RDS")

All.posteriors %>% group_by(Param,pft) %>% summarise(mean(value))

# Mean
ensemble <- All.posteriors %>% group_by(Param,pft) %>% summarise(value_m = mean(value),
                                                                 value_low = confint(lm(value ~ 1))[1],
                                                                 value_high = confint(lm(value ~ 1))[2]) 
