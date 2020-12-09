rm(list = ls())

library(purrr)
library(abind)
library(redr)
library(dplyr)
library(pracma)

Npatches = 49
par.wl = c(400:2499)
Nwl = length(par.wl)
N = 2

radprof_up <- read_radprof(Npatches = Npatches,dir =  "/home/carya/output/PEcAn_99000000001/out/PDA_EDRTM/TT_prospect/patches/",
             file.name = "radprof_up")
radprof_down <- read_radprof(Npatches = Npatches,dir =  "/home/carya/output/PEcAn_99000000001/out/PDA_EDRTM/TT_prospect/patches/",
                           file.name = "radprof_down")
LAI <- read_LAI(Npatches = Npatches,dir =  "/home/carya/output/PEcAn_99000000001/out/PDA_EDRTM/TT_prospect/patches/",
                file.name = "LAI")

par <- par.wl>=400 & par.wl <= 700

fapar <- (rowSums(radprof_down[,par,1]) - rowSums(radprof_down[,par,2]))/rowSums(radprof_down[,par,1])

ipatch = 29
plot(par.wl,radprof_down[ipatch,,1],ylim = c(0,0.5),type='l',lty=2) # Top
lines(par.wl,radprof_down[ipatch,,2],ylim = c(0,0.5),col='red',lty=2) # Bottom

Albedo <- radprof_up[ipatch,,1]/radprof_down[ipatch,,1]
plot(par.wl,Albedo,type='l')

####################################################################################
SMARTS_dir <- "~/SMARTS_295_Linux"
file2write <- file.path(SMARTS_dir,"Albedo","ALBEDO.DAT")

wl_min <- 400
wl_max <- 2500

temp_spectrum <- as.data.frame(cbind(par.wl,Albedo)) %>% rename(wavelength = par.wl,
                                                                Reflectance_median = Albedo)

spectrum <- rbind(temp_spectrum,
                  c(wl_max,temp_spectrum$Reflectance_median[which.max(temp_spectrum$wavelength)])) %>% mutate (wavelength = wavelength/1000)

write.table(x = spectrum,file = file2write, col.names = FALSE,row.names = FALSE, sep = " ", dec = ".")

NL <- default_NL.SMARTS()
NL_BCI <- NL
NL_BCI$Date_loc <- c(2004,1,1,12,9.2,-79.85,-5)
dumb <- write_NL.SMARTS(file =  file.path(SMARTS_dir,"smarts295.inp.txt"),NL_BCI)

# Run
CD <- getwd()
setwd(SMARTS_dir)
if (file.exists(file.path(SMARTS_dir,"smarts295.out.txt"))) system2("rm",file.path(SMARTS_dir,"smarts295.out.txt"))
if (file.exists(file.path(SMARTS_dir,"smarts295.ext.txt"))) system2("rm",file.path(SMARTS_dir,"smarts295.ext.txt"))
system2("./smarts295.exe",stdout = "", stderr = "")
setwd(CD)

# Plot
output_f <- file.path(SMARTS_dir,'smarts295.ext.txt')

plot(NA,NA,xlim=c(0,4000),ylim=c(0,2))

OP <- read.table(output_f, header = T, sep = "", dec = ".")
Radiation <- interp1(OP[,1],OP[,2],par.wl)
lines(par.wl,Radiation,type='l') # W/m2



