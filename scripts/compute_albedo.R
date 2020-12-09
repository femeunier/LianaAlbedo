rm(list = ls())

library(dplyr)
library(pracma)


file2read <- '/home/carya/data/smarts295_test.ext.txt'

BB_min <-  400
BB_max <- 2500

########################################################################
# Spectrum

# Canopy
Select = "Marvin"  # "Kalacska" "Marvin"   "Sanchez" 
load(file = "~/data/RTM/Inverse_canopy_spectrum.Rdata")
scenars <- c('low','high')
spectra <- All_canopy_spectra %>% mutate(median = Reflectance_median,
                                         scenar = scenario)%>% filter(ref == Select)  # Data
# spectra <- model_ensemble_all %>% filter(ref == Select)                               # Model

# Leaf
# Select = "Guzman"           #  "Guzman"      "Kalacska"    "Sanchez_PNM" "Sanchez_FTS" "Castro_FTS" "Castro_PNM"
# load(file = "~/data/RTM/Inverse_leaf_spectrum2.Rdata")
# scenars <- c('Tree_optical','Liana_optical')
# spectra <- All_leaf_spectra %>% mutate(median = Reflectance_median,
#                                        scenar = pft) %>% filter(ref == Select)       # Data
# spectra <- ensemble_posterior_all %>% mutate(scenar = pft) %>% filter(ref == Select) # Model

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


par(mfrow = c(2,1))
plot(c(BB_min,BB_max),c(NA,NA),ylim = c(0.9*min(spectra %>% filter(wavelength >= BB_min & wavelength <= BB_max) %>% pull(median)),
                                        1.1*max(spectra %>% filter(wavelength >= BB_min & wavelength <= BB_max) %>% pull(median))))

df_albedo <- data.frame()

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
  
  lines(BB,X,type='l',lty = iscenar)
  
  abline(h = alpha_BB,col='red',lty = iscenar)
  abline(h = mean(X,na.rm=TRUE),col='blue',lty = iscenar)
  
  df_albedo <- rbind(df_albedo,
                     data.frame(scenar = scenars[iscenar],
                                BB_albedo = alpha_BB,
                                alpha_m = mean(X,na.rm=TRUE),
                                ref = Select))
  
}

(df_albedo[2,2] - df_albedo[1,2])/df_albedo[1,2]
plot(BB,w,type='l')