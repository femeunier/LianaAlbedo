rm(list = ls())

library(dplyr)
library(ggplot2)
library(pracma)

data.all <- readRDS(file = "./Sanchez_allfits2.RDS")

params.all <- do.call(rbind,map(1:length(GFs),function(iGF) {do.call(rbind,map(1:length(data.all[[GFs[iGF]]]),function(i){
  print(c(iGF,i))
  params <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'param'))
  site <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'site'))
  if (is.null(params)){
    return(NULL)
  } else{
    data.frame(value = params, params = sub(".*\\.", "", names(params)), site = site, species = i, GF = GFs[iGF])
  }
}))
}))

params.all2 <- params.all %>% filter(!(params == "ssigma"))

ggplot(data = params.all2,aes(x = value, y = GF, fill = GF)) +
  geom_density_ridges(alpha= 0.5) +
  facet_grid(~ params,scales = "free") +
  scale_color_manual(values = c("darkblue","darkgreen")) +
  scale_fill_manual(values = c("darkblue","darkgreen")) +
  theme_bw() 

