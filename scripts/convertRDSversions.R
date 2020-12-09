rm(list = ls())

file2load <- "/kyukon/data/gent/vo/000/gvo00074/felicien/R/MarvinTreePixelReflectanceWV3.rds"
Marvin <- readRDS(file2load)
saveRDS(Marvin,file = "/kyukon/data/gent/vo/000/gvo00074/felicien/R/Marvin.RDS")

file2load <- "/kyukon/data/gent/vo/000/gvo00074/felicien/R/LianaRemovalPixelReflectanceWV3.rds"
Stefan <- readRDS(file2load)
saveRDS(Stefan,file = "/kyukon/data/gent/vo/000/gvo00074/felicien/R/Stefan.RDS")




