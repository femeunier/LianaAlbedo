rm(list=ls())

library(dplyr)
library(purrr)

data2read <- "~/data/RTM/data/Figure3_marvin.csv"# rm(list=ls())

data <- read.csv(data2read)
scenarios <- c("0.75","0.25")
C <- c("#1E64C8","#137300")
plot(NA,NA,xlim=c(400,2500),ylim=c(0,0.5),xlab="Wavelength",ylab="Reflectance")
data_R <- list()
for (iscenario in seq(scenarios)){
  scenario <- scenarios[iscenario]
  columns <- ((iscenario-1)*2+1):((iscenario-1)*2+2)
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  pos <- sort(temp[,1],index.return = TRUE)
  data_R[[scenario]] <- temp[pos$ix,]
  lines(data_R[[scenario]] [,1],data_R[[scenario]] [,2],col = C[iscenario],type='p',pch=19,lwd=0.01)
  names(data_R[[scenario]]) <- c("wavelength","reflectance")
}

saveRDS(data_R,file = "~/data/RTM/Figure3_marvin.rds")


###################################################################################################

rm(list=ls())

library(dplyr)
library(purrr)

data2read <- "~/data/RTM/data/Forster2008_Figure2.csv"# rm(list=ls())

data <- read.csv(data2read)
scenarios <- c("1","0.")
C <- c("#1E64C8","#137300")
plot(NA,NA,xlim=c(400,2500),ylim=c(0,0.5),xlab="Wavelength",ylab="Reflectance")
data_R <- list()
for (iscenario in seq(scenarios)){
  scenario <- scenarios[iscenario]
  columns <- ((iscenario-1)*2+1):((iscenario-1)*2+2)
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  pos <- sort(temp[,1],index.return = TRUE)
  data_R[[scenario]] <- temp[pos$ix,]
  lines(data_R[[scenario]] [,1],data_R[[scenario]] [,2],col = C[iscenario],type='p',pch=19,lwd=0.01)
  names(data_R[[scenario]]) <- c("wavelength","reflectance")
}

saveRDS(data_R,file = "~/data/RTM/Figure2_foster.rds")




##################################################################
rm(list=ls())

data2read <- "~/data/RTM/data/Figure1_kalacslab.csv"

data <- read.csv(data2read)
scenarios <- c("0.4","0.")
C <- c("#1E64C8","#137300")

plot(NA,NA,xlim=c(400,2500),ylim=c(0,1),xlab="Wavelength",ylab="Reflectance",xaxt="n")
axis(side=1, at=c(450,518,614,757,950,1172,1390,1593,1777,1946,2101,2244), labels = c(450,518,614,757,950,1172,1390,1593,1777,1946,2101,2244))

data_R <- list()

for (iscenario in seq(scenarios)){
  scenario <- scenarios[iscenario]
  columns <- ((iscenario-1)*2+1):((iscenario-1)*2+2)
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  pos <- sort(temp[,1],index.return = TRUE)
  data_R[[scenario]] <- temp[pos$ix,]
  names(data_R[[scenario]]) <- c("wavelength","reflectance")
  
  th <- 1390
  fac <- 1.3
  data_R[[scenario]]["wavelength"][data_R[[scenario]]["wavelength"]<th ] <-  data_R[[scenario]]["wavelength"][data_R[[scenario]]["wavelength"]< th]/fac
  
  th1 <- 757
  th2 <- 1390
  fac <- 1.2
  pos <- which(data_R[[scenario]]["wavelength"]> th1 &
    data_R[[scenario]]["wavelength"]<th2)
  data_R[[scenario]][pos,"wavelength"] <-  
    data_R[[scenario]][pos,"wavelength"]*seq(1,fac,length.out = length(pos))
  
  th1 <- 0
  th2 <- 600
  fac <- 1.25
  pos <- which(data_R[[scenario]]["wavelength"]> th1 &
                 data_R[[scenario]]["wavelength"]<th2)
  data_R[[scenario]][pos,"wavelength"] <-  
    data_R[[scenario]][pos,"wavelength"]*seq(fac,1,length.out = length(pos))

  lines(data_R[[scenario]] [,1],data_R[[scenario]] [,2]/max(data_R[[scenario]] [,2],na.rm=TRUE),col = C[iscenario])

}

Marvin <- readRDS(file = "~/data/RTM/Figure3_marvin.rds")
lines(Marvin[[1]][,1],Marvin[[1]][,2]/max(Marvin[[1]][,2],na.rm=TRUE))

saveRDS(data_R,file = "~/data/RTM/Figure1_kalacslab.rds")

####################################################################

rm(list=ls())

data2read <- "~/data/RTM/data/Figure1_kalacslab_mod.csv"
data2read2 <- "~/data/RTM/data/Figure1_kalacslab_mod2.csv"
axis_t <- "~/data/RTM/data/kalacska_axis.csv"

data <- read.csv(data2read)
data2 <- read.csv(data2read2)

fac_reduc = 0.09405567/0.08865586 # y-axis

axis_data <- read.csv(axis_t,header = FALSE)
axis_data <- rbind(axis_data,c(2300,2300))

scenarios <- c("0.4","0.")
C <- c("#1E64C8","#137300")

par(mfrow=c(1,1))
plot(NA,NA,xlim=c(400,2500),ylim=c(0,1),xlab="Wavelength",ylab="Reflectance",xaxt="n")
axis(side=1, 
     at     = c(450,518,614,757,950,1172,1390,1593,1777,1946,2101,2244),
     labels = c(450,518,614,757,950,1172,1390,1593,1777,1946,2101,2244))

data_R <- data_R2 <- list()

for (iscenario in seq(scenarios)){
  scenario <- scenarios[iscenario]
  columns <- ((iscenario-1)*2+1):((iscenario-1)*2+2)
  
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  temp2 <- data2[!is.na(data2[,columns[1]]),c(columns)]
  
  pos <- sort(temp[,1],index.return = TRUE)
  pos2 <- sort(temp2[,1],index.return = TRUE)
  
  data_R[[scenario]] <- temp[pos$ix,]
  data_R2[[scenario]] <- temp2[pos2$ix,]
  
  names(data_R[[scenario]]) <- c("wavelength","reflectance")
  names(data_R2[[scenario]]) <- c("wavelength","reflectance")
  
  data_R[[scenario]][,1] <- interp1(axis_data[,1],axis_data[,2],data_R[[scenario]][,1])
  data_R[[scenario]][,2] <- data_R[[scenario]][,2]/fac_reduc
  data_R[[scenario]] <- rbind(data_R2[[scenario]],
                              data_R[[scenario]][data_R[[scenario]][,1]>700,])
  
  
  data_R[[scenario]] [,1] = data_R[[scenario]] [,1]-40
  
  lines(data_R[[scenario]] [,1],data_R[[scenario]] [,2],col = C[iscenario],type='p',pch=19,lwd=0.01)
  
}


Marvin <- readRDS(file = "~/data/RTM/Figure3_sanchez2006.rds")
lines(Marvin[[1]][,1],Marvin[[1]][,2])

saveRDS(data_R,file = "~/data/RTM/Figure1_kalacslab_mod.rds")




##################################################################
rm(list=ls())

data2read <- "~/data/RTM/data/Figure3_sanchez2006_small.csv"
data2read <- "~/data/RTM/data/Figure3_sanchez2006.csv"

data <- read.csv(data2read)
scenarios <- c("high","low")
C <- c("#1E64C8","#137300")

plot(NA,NA,xlim=c(400,2500),ylim=c(0,0.7),xlab="Wavelength",ylab="Reflectance")
data_R <- list()

for (iscenario in seq(scenarios)){
  scenario <- scenarios[iscenario]
  columns <- ((iscenario-1)*3+1):((iscenario-1)*3+3)
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  pos <- sort(temp[,1],index.return = TRUE)
  data_R[[scenario]] <- temp[pos$ix,]
  names(data_R[[scenario]]) <- c("wavelength","reflectance","group")
  
  lines(data_R[[scenario]] [,1],data_R[[scenario]] [,2],col = C[iscenario],type='p',pch=19,lwd=0.01)
  
}

saveRDS(data_R,file = "~/data/RTM/Figure3_sanchez2006.rds")

###################################################################################
rm(list = ls())
        
alpha = 0.05

files <- c("~/data/RTM/Figure1_kalacslab_mod.rds",
           "~/data/RTM/Figure3_marvin.rds",
           "~/data/RTM/Figure2_foster.rds",
           "~/data/RTM/Figure3_sanchez2006.rds")
Names <- c("Kalacska","Marvin","Foster","Sanchez")


data <- readRDS("~/data/RTM/Figure3_sanchez2006.rds")
data_sanchez <- rbind(cbind(as.data.frame(data[[2]]),scenario='low'),
      cbind(as.data.frame(data[[1]]),scenario='high')) %>% mutate(ref = "Sanchez")

All_canopy_spectrum <- rbind(do.call("rbind",map2(files,Names,function(file,name){
  data <- readRDS(file)
  return(rbind(cbind(as.data.frame(data[[2]]),scenario='low') %>% mutate(group = 1),
               cbind(as.data.frame(data[[1]]),scenario='high') %>% mutate(group = 1)) %>% mutate(
                 ref = name
               ))})),data_sanchez)
  
 
All_canopy_spectrum <- All_canopy_spectrum %>% group_by(ref,scenario,wavelength,group) %>%
  summarise(Reflectance_min = min(reflectance),
            Reflectance_max = max(reflectance),
            reflectance = median(reflectance),
            Reflectance_median = median(reflectance),
            Reflectance_alphamin = quantile(reflectance,alpha/2),
            Reflectance_alphamax = quantile(reflectance,1-alpha/2))

All_canopy_spectrum <- All_canopy_spectrum %>% filter(! (ref == 'Marvin' & wavelength > 1350 & wavelength < 1520))


All_canopy_spectrum2 <- data.frame()
refs <- unique(All_canopy_spectrum$ref)
scenarios <- unique(All_canopy_spectrum$scenario)

for (cref in refs){
  for (cscenario in scenarios){
    temp_data <- All_canopy_spectrum %>% filter(ref ==cref & scenario == cscenario)
    for (cgroup in seq(1,length(unique(temp_data$group)))) {
    temp <- temp_data %>% filter(group == cgroup)
    x <- temp$wavelength
    y <- temp$Reflectance_median
    
    wavelengths <- x
    for (iwl in seq(1,length(wavelengths)-1)){
      if ((wavelengths[iwl + 1]-wavelengths[iwl]) > 100){
        x <- c(x,c(wavelengths[iwl] + 1,wavelengths[iwl + 1] -1))
        y <- c(y,c(NA,NA))
      }
    }

    if (cref == "Sanchez"){
      xi <- 505:2300
    } else {
      xi <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE))
    }

    yi <- interp1(x,y,xi,"cubic") 
    All_canopy_spectrum2 <- rbind(All_canopy_spectrum2,
                                  data.frame(wavelength = xi,
                                             Reflectance_median = yi,
                                             ref = cref,
                                             group = cgroup,
                                             scenario = cscenario))
    }
  }
}

All_canopy_spectrum2 <- All_canopy_spectrum2 %>% group_by(ref,scenario,wavelength) %>%
  summarise(Reflectance_median = median(Reflectance_median))

ggplot(data=All_canopy_spectrum2,
       aes(x = wavelength,
           group = interaction(scenario),
           y = Reflectance_median,
           fill = scenario,
           color = scenario)) +
  # geom_ribbon(alpha = 0.5,linetype = 0) +
  geom_line(size=0.5) +
  facet_wrap(ref ~ .,scales = "fixed") +
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]",
       colour = "Liana infestation") +
  scale_color_manual(values = c("#137300","#1E64C8")) +
  scale_fill_manual(values = c("#137300","#1E64C8")) +
  scale_x_continuous(limits = c(400,2500)) +
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

ggsave(filename = "~/data/RTM/All_canopy_spectra.png",
       plot = last_plot(),
       width = 20, height = 10,
       dpi = 300,units = "cm")

saveRDS(All_canopy_spectrum,file= "~/data/RTM/All_canopy_spectra.rds")

############################################################################
# Data table

bands_m <- c(550,600,400,800,1500,2000,800,680)
bands_M <- c(600,800,700,1400,1800,2500,2500,700)
band_names <- c('green','Red edge','visible','NIR','SWIR','SWIR2',"SWIR_all",'Rededge2')
studies <- unique(All_canopy_spectrum$ref)

differences <- matrix(NA,length(studies),length(band_names))
rownames(differences) <- studies
colnames(differences) <- band_names

for (istudy in seq(studies)){
  for (iband in seq(band_names)){
    
    Lianas_temp <- All_canopy_spectrum %>% filter(ref == studies[istudy] & scenario == "high" & wavelength >= bands_m[iband] &  wavelength <= bands_M[iband]) %>% ungroup() %>% dplyr::select(wavelength,reflectance)
    
    if (nrow(Lianas_temp)>0){
      Lianas <- rbind(unlist(c(bands_m[iband],Lianas_temp[1,2])),
                      Lianas_temp,
                      unlist(c(bands_M[iband],Lianas_temp[nrow(Lianas_temp),2])))
      L_wl <- Lianas[,1]
      L_R <- Lianas[,2]
      Liana_interp <- interp1(unlist(L_wl),unlist(L_R),bands_m[iband]:bands_M[iband])
      
      Trees_temp <- All_canopy_spectrum %>% filter(ref == studies[istudy] & scenario == "low" & wavelength >= bands_m[iband] &  wavelength <= bands_M[iband]) %>% ungroup() %>% dplyr::select(wavelength,reflectance)
      Trees <- rbind(unlist(c(bands_m[iband],Trees_temp[1,2])),
                     Trees_temp,
                     unlist(c(bands_M[iband],Trees_temp[nrow(Trees_temp),2])))
      T_wl <- Trees[,1]
      T_R <- Trees[,2]
      Trees_interp <- interp1(unlist(T_wl),unlist(T_R),bands_m[iband]:bands_M[iband])
      
      if (band_names[iband] == "Red edge"){
        if (nrow(Lianas_temp)>0 & nrow(Trees_temp)>0 
            # & min(max(L_wl[2:(length(L_wl)-1)]),max(T_wl[2:(length(T_wl)-1)])) - bands_M[iband] > -100 &
            # max(min(L_wl[2:(length(L_wl)-1)]),min(T_wl[2:(length(T_wl)-1)])) - bands_m[iband] < 100
        ){
          red_edge_L <- mean(Liana_interp[(length(Liana_interp)-10):length(Liana_interp)]) - mean(Liana_interp[1:10])
          red_edge_T <- mean(Trees_interp[(length(Trees_interp)-10):length(Trees_interp)]) - mean(Trees_interp[1:10])
          
          if (mean(red_edge_T)>mean(red_edge_L)){
            differences[istudy,iband] <- "T"
          } else{
            differences[istudy,iband] <- "L"      
          }
        } 
        
      } else {
        if (nrow(Lianas_temp)>0 & nrow(Trees_temp)>0 
            # & min(max(L_wl[2:(length(L_wl)-1)]),max(T_wl[2:(length(T_wl)-1)])) - bands_M[iband] > -100 &
            # max(min(L_wl[2:(length(L_wl)-1)]),min(T_wl[2:(length(T_wl)-1)])) - bands_m[iband] < 100
        ){
          if (mean(Trees_interp)>mean(Liana_interp)){
            differences[istudy,iband] <- "T"
          } else{
            differences[istudy,iband] <- "L"      
          }
        } 
      }
    }
  }
}

# sanchez green = L, red edge = L, NIR = L, SWIR = L

