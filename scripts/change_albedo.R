rm(list = ls())

library(dplyr)
library(hdf5r)
library(redr)
library(pracma)
library(ggplot2)
library(Hmisc)
library(cowplot)
library(BayesianTools)
library(zoo)
library(purrr)
library(abind)
library(PEcAn.ED2)

# Outputs of inversion
load(file = "~/data/RTM/Inverse_canopy_spectrum.Rdata")

# We choose one of the studies
Select = "Sanchez"  # "Kalacska" "Marvin"   "Sanchez" "Foster"

Nensemble = 1
alpha = 0.05
crown_mod = 0
par.wl = c(400:2499)
nir.wl = c(2500)
PFTselect = 1
alpha_frac = 0.8

############################################################################################################
# Radiations

wl_min <- 400
wl_max <- 2500

Nwl = length(par.wl)

# Filters
par <- par.wl>=400 & par.wl <= 700
red <- par.wl>=550 & par.wl <= 700
VNIR <- par.wl>=750 & par.wl <= 900
NIR <- par.wl>=800 & par.wl <= 1400
SIR <- par.wl>=1500 & par.wl <= 2500
SW <- par.wl>=400 & par.wl <= 2500

#############################################################################################################
# Prepare BCI runs

h5file <- "~/data/BCI/bci-S-2004-01-01-000000-g01.h5"

Colors <- c("#137300","#1E64C8")
df_PFT <- data.frame(names = c("Liana_optical","Tree_optical"),PFTnum = c(17,3),Col = Colors)
df_PFT <- df_PFT %>% arrange(PFTnum)%>% mutate(Col = Colors)
df_PFT <- df_PFT %>% arrange(names)  

df_PFT_tree <- df_PFT[2,]

rundir <- "/home/carya/output/PEcAn_99000000001/run/SA-median"
modeloutdir <- "/home/carya/output/PEcAn_99000000001/out/PDA_EDRTM"

# edr_exe_path <- file.path("/home/carya/ED2/EDR/run","ed_2.1-opt")
edr_exe_path <- "/home/carya/ED2/EDR/build/ed_2.1-opt-master-5e6fb784"

# Modify ED2IN
ed2in.file <- "/home/carya/output/PEcAn_99000000001/out/PDA_EDRTM/ED2IN"
ed2in <- read_ed2in(ed2in.file)

ed2in$SFILIN <- file.path(dirname(h5file), sub("\\-.*", "",basename(h5file)))
ed2in$RUNTYPE = "HISTORY"
ed2in$POI_LAT = 9.2
ed2in$POI_LON = -79.8

PEcAn.ED2::write_ed2in(ed2in, file.path(modeloutdir,"ED2IN"))
ED2IN_file <- file.path(modeloutdir,"ED2IN")

hfile <- hdf5r::H5File$new(h5file)

dbh <- readDataSet(hfile[["DBH"]])
nplant <- readDataSet(hfile[["NPLANT"]])
Npatch <- readDataSet(hfile[["NPATCHES_GLOBAL"]])
hite <- readDataSet(hfile[["HITE"]])
pft <- match(readDataSet(hfile[["PFT"]]),df_PFT$PFTnum) # Liana = 1, Tree = 2
PACO_N <- readDataSet(hfile[["PACO_N"]])
hfile$close_all()

Ncohort <- length(dbh)
npft <- max(as.numeric(pft))

Npatch <- length(PACO_N)
PAN <- rep(1:Npatch,PACO_N)

inventory <- list(dbh = dbh,pft = pft,hite = hite,nplant = nplant, PAN = PAN,PFTselect = PFTselect,
                  PACO_N = PACO_N, Ncohort = Ncohort, Npatch = Npatch)

OP_dir <- modeloutdir
h5file_paste <- file.path(modeloutdir,"history-S-2004-01-01-000000-g01.h5")
h5file_paste2 <- file.path(modeloutdir,"history-S-2004-01-01-120000-g01.h5")
system2("cp",args = list(h5file,h5file_paste))
system2("cp",args = list(h5file,h5file_paste2))

# Let's read some info that won't change (no fusion)
h5 <- ncdf4::nc_open(h5file_paste2)
area <- ncdf4::ncvar_get(h5,"AREA")
Npatches <- length(area)
ncdf4::nc_close(h5)

#############################################################################################################
current_ref <- Select

Spectrum_canopy_data <- All_canopy_spectra %>% dplyr::select(c(ref,scenario,wavelength,Reflectance_median)) %>%
  rename(reflectance = Reflectance_median) %>% filter (ref == current_ref)

observed <- observation <- Spectrum_canopy_data

ensemble <- getSample(best_set,numSamples = Nensemble,start = 3500)
ensemble <- rbind(MAP(best_set)$parametersMAP,
                  ensemble)

df_results <- df_patches <- data.frame()

#############################################################################
# SMARTS
NL <- default_NL.SMARTS()
SMARTS_dir <- "~/SMARTS_295_Linux"
file2write <- file.path(SMARTS_dir,"Albedo","ALBEDO.DAT")

print('... Ensemble runs ...')

# nrow(ensemble)
for (iensemble in seq(1,nrow(ensemble))){ 
  
  print(iensemble/nrow(ensemble))
  
  #############################################################################################
  # Liana and tree parameters as calibrated for the select site
  
  print("Liana-Tree")
  
  LT_EDRTM <- file.path(modeloutdir,"LT_EDRTM")
  if(!dir.exists(LT_EDRTM)) dir.create(LT_EDRTM)
  if(!dir.exists(file.path(LT_EDRTM,"patches"))) dir.create(file.path(LT_EDRTM,"patches"))
  
  setwd(LT_EDRTM)
  
  params <- ensemble[iensemble,]

  pft_params <- matrix(params[-c(1,2)], ncol = npft)
  rownames(pft_params) <- head(names(params[-c(1,2)]),length(params[-c(1,2)])/2)
  
  outputs_LT <- run_ED_RTM(rundir,
                        outdir = modeloutdir,params,crown_mod,inventory,par.wl,nir.wl,temp = FALSE,outputdir = LT_EDRTM)
  
  radprof_up <- read_radprof(Npatches = Npatches,dir =  file.path(LT_EDRTM,"patches"), file.name = "radprof_up")
  radprof_down <- read_radprof(Npatches = Npatches,dir =  file.path(LT_EDRTM,"patches"), file.name = "radprof_down")
  LAI <- unlist(lapply(read_LAI(Npatches = 49, file.path(LT_EDRTM,"patches"),"LAI"),sum))
  LAI_PY <- weighted.mean(LAI,area)
  
  # fapar
  fapar <- (rowSums(radprof_down[,par,1]) - rowSums(radprof_down[,par,2]))/rowSums(radprof_down[,par,1])
  fapar_PY <- weighted.mean(fapar,area)
  
  PAR_bot <- rowSums(radprof_down[,par,2])
  PAR_bot_PY <- weighted.mean(PAR_bot,area)
  
  # Albedo
  Reflectance <- t(radprof_up[,,1]/radprof_down[,,1])
  
  Irradiance <- array(NA,dim = c(Nwl,Npatches))
  for (ipatch in seq(1,Npatches)){
    temp_spectrum <- as.data.frame(cbind(par.wl,Reflectance[,ipatch])) %>% rename(wavelength = par.wl,
                                                                                  Reflectance_median = V2)
    spectrum <- rbind(temp_spectrum,
                      c(wl_max,temp_spectrum$Reflectance_median[which.max(temp_spectrum$wavelength)])) %>% mutate (wavelength = wavelength/1000)
    
    write.table(x = spectrum,file = file2write, col.names = FALSE,row.names = FALSE, sep = " ", dec = ".")
    dumb <- write_NL.SMARTS(file =  file.path(SMARTS_dir,"smarts295.inp.txt"),NL)
    dumb <- run_SMARTS(SMARTS_dir)
    
    OP <- read.table(file.path(SMARTS_dir,'smarts295.ext.txt'), header = T, sep = "", dec = ".")
    if (sum(OP[,2],na.rm=TRUE) > 0) Irradiance[,ipatch] <- interp1(OP[,1],OP[,2],par.wl)
  }
  Irradiance_site <- array(rowMeans(Irradiance,na.rm = TRUE),c(Nwl,Npatches))
  
  Albedo_PAR <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance[par,ipatch],Irradiance_site[par,ipatch])}))
  Albedo_PAR_PY <- weighted.mean(Albedo_PAR,area)
  
  Albedo_red <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance[red,ipatch],Irradiance_site[red,ipatch])}))
  Albedo_red_PY <- weighted.mean(Albedo_red,area)
  
  Albedo_VNIR <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance[VNIR,ipatch],Irradiance_site[VNIR,ipatch])}))
  Albedo_VNIR_PY <- weighted.mean(Albedo_VNIR,area)
  
  Albedo_NIR <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance[NIR,ipatch],Irradiance_site[NIR,ipatch])}))
  Albedo_NIR_PY <- weighted.mean(Albedo_NIR,area)
  
  Albedo_SIR <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance[SIR,ipatch],Irradiance_site[SIR,ipatch])}))
  Albedo_SIR_PY <- weighted.mean(Albedo_SIR,area)
  
  Albedo <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance[SW,ipatch],Irradiance_site[SW,ipatch])}))
  Albedo_PY <- weighted.mean(Albedo,area)
  
  NDVI <- (Albedo_NIR - Albedo_red)/(Albedo_NIR + Albedo_red)
  NDVI_PY <- weighted.mean(NDVI,area)
  
  df_patches <- rbind(df_patches,
                      data.frame(patch_num = 1:Npatch,
                                 scenario = 'Tree_liana',
                                 run = iensemble,
                                 fapar = fapar,
                                 PAR_bot = PAR_bot,
                                 LAI = LAI,
                                 NDVI = NDVI,
                                 Albedo = Albedo,
                                 Albedo_PAR = Albedo_PAR,
                                 Albedo_red = Albedo_red,
                                 Albedo_VNIR = Albedo_VNIR,
                                 Albedo_NIR = Albedo_NIR,
                                 Albedo_SIR = Albedo_SIR))
  
  system2('rm',c('-r',LT_EDRTM))
  
  #############################################################################################
  # Liana prospect parameters = Tree prospect parameters
  
  print("Liana-Tree, prospect")
  
  param_prospect_trees <- pft_params
  param_prospect_trees[c("Cab","Car","Cw","Cm","Nlayers"),1] <- param_prospect_trees[c("Cab","Car","Cw","Cm","Nlayers"),2]
  params_LT <- c(params[1:2],as.vector(param_prospect_trees))
  names(params_LT) <- names(params)
  
  TT_prospect <- file.path(modeloutdir,"TT_prospect")
  if(!dir.exists(TT_prospect)) dir.create(TT_prospect)
  if(!dir.exists(file.path(TT_prospect,"patches"))) dir.create(file.path(TT_prospect,"patches"))
  
  setwd(TT_prospect)
  
  outputs_TT_prospect <- run_ED_RTM(rundir,
                                 outdir = modeloutdir,params_LT,crown_mod,inventory,par.wl,nir.wl,temp = FALSE,outputdir = TT_prospect)
  
  radprof_up <- read_radprof(Npatches = Npatches,dir =  file.path(TT_prospect,"patches"), file.name = "radprof_up")
  radprof_down <- read_radprof(Npatches = Npatches,dir =  file.path(TT_prospect,"patches"), file.name = "radprof_down")
  LAI_prospect <- unlist(lapply(read_LAI(Npatches = 49, file.path(TT_prospect,"patches"),"LAI"),sum))
  LAI_prospect_PY <- weighted.mean(LAI_prospect,area)
  
  # fapar
  fapar_prospect <- (rowSums(radprof_down[,par,1]) - rowSums(radprof_down[,par,2]))/rowSums(radprof_down[,par,1])
  fapar_prospect_PY <- weighted.mean(fapar_prospect,area)
  
  PAR_bot_prospect <- rowSums(radprof_down[,par,2])
  PAR_bot_prospect_PY <- weighted.mean(PAR_bot_prospect,area)
  
  # Albedo
  Reflectance_TT <- t(radprof_up[,,1]/radprof_down[,,1])
  
  Irradiance <- array(NA,dim = c(Nwl,Npatches))
  for (ipatch in seq(1,Npatches)){
    temp_spectrum <- as.data.frame(cbind(par.wl,Reflectance_TT[,ipatch])) %>% rename(wavelength = par.wl,
                                                                                  Reflectance_median = V2)
    spectrum <- rbind(temp_spectrum,
                      c(wl_max,temp_spectrum$Reflectance_median[which.max(temp_spectrum$wavelength)])) %>% mutate (wavelength = wavelength/1000)
    
    write.table(x = spectrum,file = file2write, col.names = FALSE,row.names = FALSE, sep = " ", dec = ".")
    dumb <- write_NL.SMARTS(file =  file.path(SMARTS_dir,"smarts295.inp.txt"),NL)
    dumb <- run_SMARTS(SMARTS_dir)
    
    OP <- read.table(file.path(SMARTS_dir,'smarts295.ext.txt'), header = T, sep = "", dec = ".")
    if (sum(OP[,2],na.rm=TRUE) > 0) Irradiance[,ipatch] <- interp1(OP[,1],OP[,2],par.wl)
  }
  Irradiance_site <- array(rowMeans(Irradiance,na.rm = TRUE),c(Nwl,Npatches))
  
  Albedo_PAR_prospect <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT[par,ipatch],Irradiance_site[par,ipatch])}))
  Albedo_PAR_prospect_PY <- weighted.mean(Albedo_PAR_prospect,area)
  
  Albedo_red_prospect <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT[red,ipatch],Irradiance_site[red,ipatch])}))
  Albedo_red_prospect_PY <- weighted.mean(Albedo_red_prospect,area)
  
  Albedo_VNIR_prospect <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT[VNIR,ipatch],Irradiance_site[VNIR,ipatch])}))
  Albedo_VNIR_prospect_PY <- weighted.mean(Albedo_VNIR_prospect,area)
  
  Albedo_NIR_prospect <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT[NIR,ipatch],Irradiance_site[NIR,ipatch])}))
  Albedo_NIR_prospect_PY <- weighted.mean(Albedo_NIR_prospect,area)
  
  Albedo_SIR_prospect <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT[SIR,ipatch],Irradiance_site[SIR,ipatch])}))
  Albedo_SIR_prospect_PY <- weighted.mean(Albedo_SIR_prospect,area)
  
  Albedo_prospect <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT[SW,ipatch],Irradiance_site[SW,ipatch])}))
  Albedo_prospect_PY <- weighted.mean(Albedo_prospect,area)
  
  NDVI_prospect <- (Albedo_NIR_prospect - Albedo_red_prospect)/(Albedo_NIR_prospect + Albedo_red_prospect)
  NDVI_prospect_PY <- weighted.mean(NDVI_prospect,area)
  
  df_patches <- rbind(df_patches,
                      data.frame(patch_num = 1:Npatch,
                                 scenario = 'prospect',
                                 run = iensemble,
                                 fapar = fapar_prospect,
                                 PAR_bot = PAR_bot_prospect,
                                 LAI = LAI_prospect,
                                 NDVI = NDVI_prospect,
                                 Albedo = Albedo_prospect,
                                 Albedo_PAR = Albedo_PAR_prospect,
                                 Albedo_red = Albedo_red_prospect,
                                 Albedo_VNIR = Albedo_VNIR_prospect,
                                 Albedo_NIR = Albedo_NIR_prospect,
                                 Albedo_SIR = Albedo_SIR_prospect))
  
  system2('rm',c('-r',TT_prospect))
  
  #############################################################################################
  # All liana parameters = Tree parameters
  
  print("Liana-Tree, EDRTM")
  
  param_EDRTM_trees <- pft_params
  param_EDRTM_trees[rownames(pft_params),1] <- param_EDRTM_trees[rownames(pft_params),2]
  params_TT <- c(params[1:2],as.vector(param_EDRTM_trees))
  names(params_TT) <- names(params)
  
  TT_EDRTM <- file.path(modeloutdir,"TT_EDRTM")
  if(!dir.exists(TT_EDRTM)) dir.create(TT_EDRTM)
  if(!dir.exists(file.path(TT_EDRTM,"patches"))) dir.create(file.path(TT_EDRTM,"patches"))
  
  setwd(TT_EDRTM)
  
  outputs_TT_EDRTM <- run_ED_RTM(rundir,
                                    outdir = modeloutdir,params_TT,crown_mod,inventory,par.wl,nir.wl,temp = FALSE,outputdir = TT_EDRTM)
  
  
  radprof_up <- read_radprof(Npatches = Npatches,dir =  file.path(TT_EDRTM,"patches"), file.name = "radprof_up")
  radprof_down <- read_radprof(Npatches = Npatches,dir =  file.path(TT_EDRTM,"patches"), file.name = "radprof_down")
  
  LAI_EDRTM <- unlist(lapply(read_LAI(Npatches = 49, file.path(TT_EDRTM,"patches"),"LAI"),sum))
  LAI_EDRTM_PY <- weighted.mean(LAI_EDRTM,area)
  
  # fapar
  fapar_EDRTM <- (rowSums(radprof_down[,par,1]) - rowSums(radprof_down[,par,2]))/rowSums(radprof_down[,par,1])
  fapar_EDRTM_PY <- weighted.mean(fapar_EDRTM,area)
  
  PAR_bot_EDRTM <- rowSums(radprof_down[,par,2])
  PAR_bot_EDRTM_PY <- weighted.mean(PAR_bot_EDRTM,area)
  
  # Albedo
  Reflectance_TT_EDRTM <- t(radprof_up[,,1]/radprof_down[,,1])
  
  Irradiance <- array(NA,dim = c(Nwl,Npatches))
  for (ipatch in seq(1,Npatches)){
    temp_spectrum <- as.data.frame(cbind(par.wl,Reflectance_TT_EDRTM[,ipatch])) %>% rename(wavelength = par.wl,
                                                                                  Reflectance_median = V2)
    spectrum <- rbind(temp_spectrum,
                      c(wl_max,temp_spectrum$Reflectance_median[which.max(temp_spectrum$wavelength)])) %>% mutate (wavelength = wavelength/1000)
    
    write.table(x = spectrum,file = file2write, col.names = FALSE,row.names = FALSE, sep = " ", dec = ".")
    dumb <- write_NL.SMARTS(file =  file.path(SMARTS_dir,"smarts295.inp.txt"),NL)
    dumb <- run_SMARTS(SMARTS_dir)
    
    OP <- read.table(file.path(SMARTS_dir,'smarts295.ext.txt'), header = T, sep = "", dec = ".")
    if (sum(OP[,2],na.rm=TRUE) > 0) Irradiance[,ipatch] <- interp1(OP[,1],OP[,2],par.wl)
  }
  Irradiance_site <- array(rowMeans(Irradiance,na.rm = TRUE),c(Nwl,Npatches))
  
  Albedo_PAR_EDRTM <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT_EDRTM[par,ipatch],Irradiance_site[par,ipatch])}))
  Albedo_PAR_EDRTM_PY <- weighted.mean(Albedo_PAR_EDRTM,area)
  
  Albedo_red_EDRTM <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT_EDRTM[red,ipatch],Irradiance_site[red,ipatch])}))
  Albedo_red_EDRTM_PY <- weighted.mean(Albedo_red_EDRTM,area)
  
  Albedo_VNIR_EDRTM <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT_EDRTM[VNIR,ipatch],Irradiance_site[VNIR,ipatch])}))
  Albedo_VNIR_EDRTM_PY <- weighted.mean(Albedo_VNIR_EDRTM,area)
  
  Albedo_NIR_EDRTM <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT_EDRTM[NIR,ipatch],Irradiance_site[NIR,ipatch])}))
  Albedo_NIR_EDRTM_PY <- weighted.mean(Albedo_NIR_EDRTM,area)
  
  Albedo_SIR_EDRTM <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT_EDRTM[SIR,ipatch],Irradiance_site[SIR,ipatch])}))
  Albedo_SIR_EDRTM_PY <- weighted.mean(Albedo_SIR_EDRTM,area)
  
  Albedo_EDRTM <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_TT_EDRTM[SW,ipatch],Irradiance_site[SW,ipatch])}))
  Albedo_EDRTM_PY <- weighted.mean(Albedo_EDRTM,area)
  
  NDVI_EDRTM <- (Albedo_NIR_EDRTM - Albedo_red_EDRTM)/(Albedo_NIR_EDRTM + Albedo_red_EDRTM)
  NDVI_EDRTM_PY <- weighted.mean(NDVI_EDRTM,area)
  
  df_patches <- rbind(df_patches,
                      data.frame(patch_num = 1:Npatch,
                                 scenario = 'EDRTM',
                                 run = iensemble,
                                 fapar = fapar_EDRTM,
                                 LAI = LAI_EDRTM,
                                 NDVI = NDVI_EDRTM,
                                 PAR_bot = PAR_bot_EDRTM,
                                 Albedo = Albedo_EDRTM,
                                 Albedo_PAR = Albedo_PAR_EDRTM,
                                 Albedo_red = Albedo_red_EDRTM,
                                 Albedo_VNIR = Albedo_VNIR_EDRTM,
                                 Albedo_NIR = Albedo_NIR_EDRTM,
                                 Albedo_SIR = Albedo_SIR_EDRTM))
  
  
  system2('rm',c('-r',TT_EDRTM))
  
  #############################################################################################
  # All liana parameters = Tree parameters
  # Could be easily improved ...
  
  print("No liana")
  
  param_EDRTM <- pft_params
  param_EDRTM[,1] <- param_EDRTM[,2] 
  rownames(param_EDRTM) <- head(names(params[-c(1,2)]),length(params[-c(1,2)])/2)
  param_EDRTM["b1Bl_large",1] <- 0.0001
  param_EDRTM["clumping_factor",1] <- 0.0001
  params_noliana <- c(params[1:2],as.vector(param_EDRTM))
  
  no_liana <- file.path(modeloutdir,"no_liana")
  if(!dir.exists(no_liana)) dir.create(no_liana)
  if(!dir.exists(file.path(no_liana,"patches"))) dir.create(file.path(no_liana,"patches"))
  
  outputs_no_liana <- run_ED_RTM(rundir = rundir,
                                 outdir = modeloutdir,params_noliana,crown_mod,inventory,par.wl,nir.wl,temp = FALSE,outputdir = no_liana)
  
  radprof_up <- read_radprof(Npatches = Npatches,dir =  file.path(no_liana,"patches"), file.name = "radprof_up")
  radprof_down <- read_radprof(Npatches = Npatches,dir =  file.path(no_liana,"patches"), file.name = "radprof_down")
  
  LAI_noliana <- unlist(lapply(read_LAI(Npatches = 49, file.path(no_liana,"patches"),"LAI"),sum))
  LAI_noliana_PY <- weighted.mean(LAI_noliana,area)
  
  # fapar
  fapar_noliana  <- (rowSums(radprof_down[,par,1]) - rowSums(radprof_down[,par,2]))/rowSums(radprof_down[,par,1])
  fapar_noliana_PY <- weighted.mean(fapar_noliana,area)
  
  PAR_bot_noliana <- rowSums(radprof_down[,par,2])
  PAR_bot_noliana_PY <- weighted.mean(PAR_bot_noliana,area)
  
  # Albedo
  Reflectance_noliana <- t(radprof_up[,,1]/radprof_down[,,1])
  
  Irradiance <- array(NA,dim = c(Nwl,Npatches))
  for (ipatch in seq(1,Npatches)){
    temp_spectrum <- as.data.frame(cbind(par.wl,Reflectance_noliana[,ipatch])) %>% rename(wavelength = par.wl,
                                                                                           Reflectance_median = V2)
    spectrum <- rbind(temp_spectrum,
                      c(wl_max,temp_spectrum$Reflectance_median[which.max(temp_spectrum$wavelength)])) %>% mutate (wavelength = wavelength/1000)
    
    write.table(x = spectrum,file = file2write, col.names = FALSE,row.names = FALSE, sep = " ", dec = ".")
    dumb <- write_NL.SMARTS(file =  file.path(SMARTS_dir,"smarts295.inp.txt"),NL)
    dumb <- run_SMARTS(SMARTS_dir)
    
    OP <- read.table(file.path(SMARTS_dir,'smarts295.ext.txt'), header = T, sep = "", dec = ".")
    if (sum(OP[,2],na.rm=TRUE) > 0) Irradiance[,ipatch] <- interp1(OP[,1],OP[,2],par.wl)
  }
  Irradiance_site <- array(rowMeans(Irradiance,na.rm = TRUE),c(Nwl,Npatches))
  
  Albedo_PAR_noliana <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_noliana[par,ipatch],Irradiance_site[par,ipatch])}))
  Albedo_PAR_noliana_PY <- weighted.mean(Albedo_PAR_noliana,area)
  
  Albedo_red_noliana <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_noliana[red,ipatch],Irradiance_site[red,ipatch])}))
  Albedo_red_noliana_PY <- weighted.mean(Albedo_red_noliana,area)
  
  Albedo_VNIR_noliana <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_noliana[VNIR,ipatch],Irradiance_site[VNIR,ipatch])}))
  Albedo_VNIR_noliana_PY <- weighted.mean(Albedo_VNIR_noliana,area)
  
  Albedo_NIR_noliana <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_noliana[NIR,ipatch],Irradiance_site[NIR,ipatch])}))
  Albedo_NIR_noliana_PY <- weighted.mean(Albedo_NIR_noliana,area)
  
  Albedo_SIR_noliana <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_noliana[SIR,ipatch],Irradiance_site[SIR,ipatch])}))
  Albedo_SIR_noliana_PY <- weighted.mean(Albedo_SIR_noliana,area)
  
  Albedo_noliana <- unlist(map(1:Npatches,function(ipatch){weighted.mean(Reflectance_noliana[SW,ipatch],Irradiance_site[SW,ipatch])}))
  Albedo_noliana_PY <- weighted.mean(Albedo_noliana,area)
  
  NDVI_noliana <- (Albedo_NIR_noliana - Albedo_red_noliana)/(Albedo_NIR_noliana + Albedo_red_noliana)
  NDVI_noliana_PY <- weighted.mean(NDVI_noliana,area)
  
  df_patches <- rbind(df_patches,
                      data.frame(patch_num = 1:Npatch,
                                 scenario = 'noliana',
                                 run = iensemble,
                                 fapar = fapar_noliana,
                                 LAI = LAI_noliana,
                                 NDVI = NDVI_noliana,
                                 PAR_bot = PAR_bot_noliana,
                                 Albedo = Albedo_noliana,
                                 Albedo_PAR = Albedo_PAR_noliana,
                                 Albedo_red = Albedo_red_noliana,
                                 Albedo_VNIR = Albedo_VNIR_noliana,
                                 Albedo_NIR = Albedo_NIR_noliana,
                                 Albedo_SIR = Albedo_SIR_noliana))
  
  
  system2('rm',c('-r',no_liana))
  
  
  ###############################################################################################################
  # Combine results

  df_results <- rbind(df_results,
                      data.frame(t(params),
                                 run = iensemble,
                                 Albedo_PAR_PY,Albedo_red_PY,Albedo_VNIR_PY,Albedo_NIR_PY,Albedo_SIR_PY,Albedo_PY,LAI_PY,fapar_PY,PAR_bot_PY,NDVI_PY,
                                 Albedo_PAR_prospect_PY,Albedo_red_prospect_PY,Albedo_VNIR_prospect_PY,Albedo_NIR_prospect_PY,Albedo_SIR_prospect_PY,Albedo_prospect_PY,LAI_prospect_PY,fapar_prospect_PY,PAR_bot_prospect_PY,NDVI_prospect_PY,
                                 Albedo_PAR_EDRTM_PY,Albedo_red_EDRTM_PY,Albedo_VNIR_EDRTM_PY,Albedo_NIR_EDRTM_PY,Albedo_SIR_EDRTM_PY,Albedo_EDRTM_PY,LAI_EDRTM_PY,fapar_EDRTM_PY,PAR_bot_EDRTM_PY,NDVI_EDRTM_PY,
                                 Albedo_PAR_noliana_PY,Albedo_red_noliana_PY,Albedo_VNIR_noliana_PY,Albedo_NIR_noliana_PY,Albedo_SIR_noliana_PY,Albedo_noliana_PY,LAI_noliana_PY,fapar_noliana_PY,PAR_bot_noliana_PY,NDVI_noliana_PY))
  
}

df_results_tr <- df_results %>% mutate(Delta_prospect = Albedo_PY - Albedo_prospect_PY,
                                       Delta_EDRTM = Albedo_prospect_PY - Albedo_EDRTM_PY,
                                       Delta = Albedo_PY - Albedo_EDRTM_PY,
                                       Delta_prospect_rel = (Albedo_PY - Albedo_prospect_PY)/Albedo_prospect_PY,
                                       Delta_EDRTM_rel = (Albedo_PY - Albedo_EDRTM_PY)/Albedo_EDRTM_PY,
                                       Delta_noliana_rel = (Albedo_PY - Albedo_noliana_PY)/Albedo_noliana_PY,
                                       Delta_rel = (Albedo_PY - Albedo_EDRTM_PY)/Albedo_EDRTM_PY,
                                       Delta_LAI = LAI_PY - LAI_EDRTM_PY,
                                       Delta_fapar = (fapar_PY - fapar_EDRTM_PY),
                                       Delta_NDVI = NDVI_PY - NDVI_EDRTM_PY,
                                       Delta_Par_bot = PAR_bot_PY - PAR_bot_EDRTM_PY)

hist(df_results_tr$Delta_EDRTM_rel)
hist(df_results_tr$Delta_EDRTM_rel-df_results_tr$Delta_noliana_rel)

boxplot(df_results_tr$Delta_rel)

boxplot(df_results_tr$Delta_LAI)
boxplot(df_results_tr$Delta_fapar)
boxplot(df_results_tr$Delta_Par_bot)
boxplot(df_results_tr$Delta_NDVI)
# save("~/data/RTM/BCI.RData")

saveRDS(object = df_patches,file = "~/data/RTM/Albedo_BCI.RDS")

# Test NDVI

df_patches_LT <- df_patches %>% filter(run < iensemble,scenario == "Tree_liana")
df_patches_TT <- df_patches %>% filter(run < iensemble,scenario == "EDRTM")

hist(100*(df_patches_LT %>% pull(PAR_bot ) - df_patches_TT %>% pull(PAR_bot ))/(df_patches_TT %>% pull(PAR_bot )))
