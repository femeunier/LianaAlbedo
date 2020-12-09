rm(list = ls())

library(dplyr)
library(ggplot2)

Colors <- c("1E64C8","#137300")

directory <- "/home/carya/data/RTM"

studies <- c("Sanchez","Foster",
             "Marvin","Sanchez_FTS",
             "Sanchez_PNM","Guzman",
             "Kalacska","Castro_PNM",
             "Castro_FTS")

cut <- c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE)
cut <- rep(TRUE,length(studies))

files2load <- file.path(directory,paste0("spectrum_",studies,".RDS"))

All.spectra <- All.spectra_uncut <- data.frame()

for (ifile in seq(1,length(files2load))){
  file <- files2load[ifile]
  temp.spectrum <- readRDS(file)
  
  spectra_sum <- temp.spectrum %>% group_by(pft,wv) %>% summarise(r_m = mean(reflectance),
                                                                  t_m = mean(transmittance))
  
  if (cut[ifile]){
    spectra_sum_cut <- spectra_sum %>% filter(wv < 1000)
  }
  
  All.spectra <- bind_rows(list(All.spectra,
                                spectra_sum_cut %>% mutate(ref = studies[ifile])))
  All.spectra_uncut <- bind_rows(list(All.spectra_uncut,
                                      spectra_sum %>% mutate(ref = studies[ifile])))
}

ggplot(All.spectra, aes(x = wv,colour = pft)) +
  geom_line(aes(y = r_m)) +
  geom_line(aes(y = 1-t_m),linetype = 1)+
  scale_y_continuous(sec.axis = sec_axis(trans = ~ 1 - 1 * ., name = "Transmittance [-]")) +
  # geom_ribbon(aes(ymin = r_m,ymax=1-t_m,fill = pft),alpha = 0.2,linetype=0)+
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]") +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  facet_wrap(~ ref) +
  theme_bw()

diff <- All.spectra %>% group_by(ref,wv) %>% summarise(diff_r = r_m[pft == "Liana"] - r_m[pft == "Tree"],
                                                       diff_t = (t_m[pft == "Liana"] - t_m[pft == "Tree"]))

ggplot(diff, aes(x = wv, colour = ref)) +
  geom_line(aes(y = diff_r)) +
  geom_line(aes(y = 1-diff_t),linetype = 1)+
  scale_y_continuous(sec.axis = sec_axis(trans = ~ 1 - 1 * ., name = "Transmittance [-]")) +
  # geom_ribbon(aes(ymin = r_m,ymax=1-t_m,fill = pft),alpha = 0.2,linetype=0)+
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]") +
  geom_hline(yintercept = c(0,1),linetype = 3) +
  theme_bw()

All <- All.spectra %>% group_by(wv,pft) %>% summarise(r_mean = mean(r_m),
                                                      t_mean = mean(t_m),
                                                      r_min = min(r_m),
                                                      r_max = max(r_m),
                                                      t_min = min(t_m),
                                                      t_max = max(t_m))

ggplot() +
  geom_line(data = diff, aes(x = wv,group = ref,y = diff_r),colour = "black",linetype=3) +
  geom_line(data = diff, aes(x = wv,group = ref,y = 1-diff_t),colour = "black",linetype=3) +
  geom_line(data = All, aes(x = wv,colour = pft,y = r_mean)) +
  geom_line(data = All, aes(x = wv,colour = pft,y = 1-t_mean),linetype = 1) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ 1 - 1 * ., name = "Transmittance [-]")) +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]") +
  theme_bw() + theme(text = element_text(size = 14),
                     axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

diff_sum <- diff %>% group_by(wv) %>% summarise(diff_r_min = min(diff_r),
                                                diff_r_max = max(diff_r),
                                                diff_t_min = min(diff_t),
                                                diff_t_max = max(diff_t))

ggplot() +
  geom_ribbon(data = diff_sum, aes(x = wv,ymin = diff_r_min,ymax = diff_r_max),fill = "darkgrey",alpha=0.3) +
  geom_ribbon(data = diff_sum, aes(x = wv,ymax = 1-diff_t_min,ymin = 1-diff_t_max),fill = "darkgrey",alpha=0.3) +
  geom_line(data = All, aes(x = wv,colour = pft,y = r_mean)) +
  geom_line(data = All, aes(x = wv,colour = pft,y = 1-t_mean),linetype = 1) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ 1 - 1 * ., name = "Transmittance [-]")) +
  geom_hline(yintercept = c(0,1),linetype=3)+
  scale_color_manual(values = c("#1E64C8","#137300")) +
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]") +
  theme_bw() + theme(text = element_text(size = 14),
                     axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


All.spectra_uncut %>% mutate(band = case_when(wv <= 700 ~ 1,
                                              wv <= 2500 ~ 2)) %>% group_by(pft,band) %>% summarise(r_m = mean(r_m),
                                                                                                    t_m = mean(t_m))


