
library(tidyverse)
library(plyr)
library(deldir)
library(alphahull)
library(xml2)
library(sp)
library(viridis)
library(readxl)

source("./GRANAR/R/granar.R")
source("./GRANAR/R/micro_hydro.R")
`%!in%` <- compose(`!`, `%in%`)
params <- read_param_xml("./GRANAR/www/Zea_mays_CT.xml")
Sampl <- read.csv("./Anat_proc/Sobol_inputGM.csv")
Sampl_vect <- 1:nrow(Sampl)

for(i in 1:50){
  tmp_sampl <- Sampl[i,]
  print(i)
  if(file.exists("./MECHA/cellsetdata/current_root.xml")){
    file.remove("./MECHA/cellsetdata/current_root.xml")
    file.remove("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml")
  }
  sim <- run_granar(params, tmp_sampl)
  file.copy("./MECHA/cellsetdata/current_root.xml", paste0("./MECHA/cellsetdata/root_",i,".xml"), overwrite = T)
  file.copy("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml",
            paste0( "./MECHA/Projects/GRANAR/in/Maize_Geometry_aer_",i,".xml"), overwrite = T)
}


