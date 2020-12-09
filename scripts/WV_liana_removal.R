rm(list = ls())

library(purrr)
library(tidyr)
library(ggplot2)
library(dplyr)

LianaRemovalPixelReflectanceWV3 <- readRDS("/home/femeunier/Documents/projects/albedo/data/LianaRemovalPixelReflectanceWV3.rds")

all.data <- LianaRemovalPixelReflectanceWV3 %>% pivot_longer(cols = -c("type","plot"),
                                                             names_to = "wl",
                                                             values_to = "reflectance")


ggplot(data = all.data) +
  geom_boxplot(aes(x = as.factor(wl),y = reflectance,fill = as.factor(type))) +
  theme_bw()

all.data %>% group_by(wl) %>% summarise(p.val = summary(aov(reflectance ~ as.factor(type)))[[1]][1,5])
