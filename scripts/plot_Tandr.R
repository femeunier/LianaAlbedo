rm(list = ls())

library(cowplot)

spectrum <- readRDS("./data/spectrum_example.RDS")

plot1 <- ggplot(spectrum, aes(x = wv,colour = pft)) +
  geom_line(aes(y = r_m)) +
  geom_line(aes(y = 1-t_m),linetype = 1)+
  scale_y_continuous(sec.axis = sec_axis(trans = ~ 1 - 1 * ., name = "Transmittance [-]")) +
  # geom_ribbon(aes(ymin = r_m,ymax=1-t_m,fill = pft),alpha = 0.2,linetype=0)+
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]") +
  scale_color_manual(values = c("red","blue")) +
  scale_fill_manual(values = c("red","blue")) +
  theme_bw()

plot2 <- ggplot(spectrum, aes(x = wv,colour = pft)) +
  geom_line(aes(y = 1-r_m-t_m)) +
  labs(y = "Absorbance [-]",
       x = "Wavelength [nm]") +
  scale_color_manual(values = c("red","blue")) +
  theme_bw()

plot_grid(plot1,plot2,align = c(, "h", "v", "hv"),)