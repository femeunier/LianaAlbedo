rm(list=ls())

file <- "/home/carya/data/Gigante/tempFoster.lat9.000lon-79.000.css"
css   <-   read.table(file,header=TRUE,stringsAsFactors=FALSE)

patches <- unique(css$patch)

css_mod <- data.frame()

for (ipatch in seq(patches)){
  
  select <- which(css$patch == ipatch)
  cpft <- unique(css$pft[select])
  
  if (17 %in% cpft){
    lianas = css[css$patch == ipatch & css$pft == 17,]
    trees = css[css$patch == ipatch & css$pft != 17,]
    css_mod <- rbind(css_mod,
                    rbind(lianas,trees[nrow(trees),]))
  } else {
    css_mod <- rbind(css_mod,
                     css[select,])
  }
  
}


write.table(x = css_mod,file = "/home/carya/data/Gigante/Foster.lat9.000lon-79.000.css",row.names = FALSE)