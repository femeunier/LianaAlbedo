rm(list = ls())

library(dplyr)
library(pracma)
library(reshape2)


file2read <- '/home/carya/data/smarts295_test.ext.txt'

BB_min_all <-  c(400,800,1500,800,400)
BB_max_all <- c(700,1400,2500,2500,2500)
names_all <- c('PAR','VNIR','SIR','IR','SW')

select_all <- c("Marvin","Kalacska","Sanchez","Foster")

df_albedo <- data.frame()

for (iBB in seq(1,length(BB_min_all))){
  
  BB_min <- BB_min_all[iBB]
  BB_max <- BB_max_all[iBB]
  ########################################################################
  # Spectrum
  
    for (iselect in seq(1,length(select_all))){
          
      # Canopy
      Select <- select_all[iselect]
      load(file = "~/data/RTM/Inverse_canopy_spectrum.Rdata")
      scenars <- c('low','high')
      
        for (itype in seq(1,2)){
          if (itype == 1){
          spectra <- All_canopy_spectra %>% mutate(median = Reflectance_median,
                                                 scenar = scenario)%>% filter(ref == Select)  # Data
          } else {
            spectra <- model_ensemble_all %>% filter(ref == Select)                               # Model    
          }
        
        
        ########################################################################
        # Extrema
        
        spectra_min <- floor(min(spectra$wavelength))
        spectra_max <- ceil(max(spectra$wavelength))
        
        BB_min <- max(BB_min,spectra_min)
        BB_max <- min(BB_max,spectra_max)
        BB <- seq(BB_min,BB_max)
        
        wl_min <- max(400,BB_min)
        wl_max <- min(2500,BB_max)
        all_wl <- seq(wl_min,wl_max)
        
        ########################################################################
        # Irradiance (sun)
        
        data <- read.table(file2read, header = TRUE, sep = "", dec = ".")
        colnames(data) <- c('wl','I') #Wavelength, Irradiance
        data_interp <- data.frame(wl = all_wl,
                                  I = interp1(data$wl,data$I,all_wl))
        
        
        
        for (iscenar in seq(1,length(scenars))){
          
          temp <-  spectra %>% filter(ref == Select,
                                      scenar == scenars[iscenar]) %>% ungroup() %>% dplyr::select(wavelength,median)
          wavelengths <- temp$wavelength
          
          for (iwl in seq(1,length(wavelengths)-1)){
            if ((wavelengths[iwl + 1]-wavelengths[iwl]) > 100){
              temp <- rbind(temp,
                            data.frame(wavelength = c(wavelengths[iwl] + 1,wavelengths[iwl + 1] -1),
                                       median = NA))
            }
          }
          
          if (wl_min < min(temp$wavelength)){
            temp <- rbind(c(wl_min,temp$median[which.min(temp$wavelength)]),
                          temp)
          }
          
          if (wl_max > max(temp$wavelength)){
            temp <- rbind(temp,
                          c(wl_max,temp$median[which.max(temp$wavelength)]))
          }
          
          Model <- temp %>% arrange(wavelength)
          
          Model_interp <- data.frame(wl = all_wl,
                                     alpha_l = interp1(Model$wavelength,Model$median,all_wl,"cubic"))
          
          X = Model_interp %>% filter(wl >= BB_min & wl <= BB_max) %>% pull(alpha_l)
          w = data_interp %>% filter(wl >= BB_min & wl <= BB_max) %>% pull(I)
          
          # Broadband albedo
          alpha_BB <- 
            weighted.mean(x = X,
                          w = w,
                          na.rm = TRUE)
          
          
          df_albedo <- rbind(df_albedo,
                             data.frame(ref = Select,
                                        scenar = scenars[iscenar],
                                        BB_albedo = alpha_BB,
                                        alpha_m = mean(X,na.rm=TRUE),
                                        type = itype,
                                        BB = names_all[iBB]
                                        ))
          
      }
    }
  }
}


# Data
Data <- df_albedo %>% filter(type == 1) %>% group_by(ref,BB) %>% summarise(alpha_d = (BB_albedo[scenar == 'high'] - BB_albedo[scenar == 'low'])/BB_albedo[scenar == 'low'])
# Model
Model <- df_albedo %>% filter(type == 2) %>% group_by(ref,BB) %>% summarise(alpha_d = (BB_albedo[scenar == 'high'] - BB_albedo[scenar == 'low'])/BB_albedo[scenar == 'low'])


dcast(Data, ref ~ BB)
dcast(Model, ref ~ BB)

