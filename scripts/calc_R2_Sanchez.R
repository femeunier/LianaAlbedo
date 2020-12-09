rm(list = ls())

library(dplyr)
library(ggplot2)
library(purrr)

data.all <- readRDS(file = "./Sanchez_allfits2.RDS")

temp <- data.all$Liana[["7"]][[1]][["posterior"]]

GFs <- c("Liana","Tree")

wv_select <- seq(500,2500,100)

data.all_formated <- do.call(rbind,map(1:length(GFs),function(iGF) {do.call(rbind,map(1:length(data.all[[GFs[iGF]]]),function(i){
  print(c(iGF,i))
  
  clist <- data.all[[GFs[iGF]]][[i]]
  clist[sapply(clist, is.null)] <- NULL
  
  waves <- unlist(lapply(sapply(clist, '[', 'spectrum'),function(x){x[,1]}))
  pos <- which(!(waves > 680 & waves < 750))
  waves = waves[pos]
  num <- unlist(lapply(sapply(clist, '[', 'spectrum'),function(x){x[,2]}))[pos]
  Reflectance <- unlist(lapply(sapply(clist, '[', 'spectrum'),function(x){x[,3]}))[pos]
  best_run <- unlist(sapply(clist, '[', 'best_run'))
  site <- unlist(sapply(clist, '[', 'site'))
  site <- site[1]
  names(site) <- NULL
  
  select <- which(waves %in% wv_select)
  
  if (is.null(best_run)){
    return(NULL)
  } else{
    data.frame(wv = waves[select], num = num[select], obs = Reflectance[select], site = site,sim = best_run[select], species = i ,GF = GFs[iGF])
  }
}))
}))

# ggplot(data = data.all_formated %>% filter(wv == 1200)) +
#   geom_point(aes(x= obs,y = sim)) +
#   theme_bw()


all2plot.lm <- bind_rows(list(data.all_formated %>% group_by(wv) %>% summarise(x = min(obs),
                                                                              y = coef(lm(sim~obs))[2]*x + coef(lm(sim~obs))[1]),
                              data.all_formated %>% group_by(wv) %>% summarise(x = max(obs),
                                                                              y = coef(lm(sim~obs))[2]*x + coef(lm(sim~obs))[1])))

WL = 1500
ggplot(data = data.all_formated) +
  geom_abline(slope = 1, linetype = 3) +
  geom_point(aes(x = obs,y = sim,color = GF),alpha = 0.1) +
  geom_line(data = all2plot.lm,
            aes(x = x,y = y, group = as.factor(wv))) +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  theme_bw() + 
  theme(legend.position = c(0.1,0.9))

df_r2 <- data.all_formated %>% group_by(wv) %>% summarise(r2 = summary(lm(sim ~ obs))[["r.squared"]])
hist(df_r2$r2)

df_r2 %>% filter(r2 < 0.8)

