

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
print("loading params")
params <- read_param_xml("./GRANAR/www/Zea_mays_CT.xml")
source("./Anat_proc/SobolR.R")
Sampl <- read.csv("./Anat_proc/Sobol_inputGM.csv")
fls <- list.files("./MECHA/cellsetdata/")
fls <- fls[grepl("root_", fls)]

for(j in fls[1:45]){
  message("--------------")
  print(j)
  message("--------------")

  if(file.exists("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt")){
    file.remove("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt")
    file.remove("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_2,1.txt")
    file.remove("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_4,2.txt")
    file.remove("./MECHA/cellsetdata/current_root.xml")
    file.remove("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml")
  }

  fc <- file.copy(paste0("./MECHA/cellsetdata/",j), "./MECHA/cellsetdata/current_root.xml", overwrite = T)
  if(fc == FALSE){next()}
  fc <- file.copy(paste0("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer_", parse_number(j), ".xml"),
              paste0("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml"), overwrite = T)
  if(fc == FALSE){next()}

  # MECHA input change
  id <- parse_number(j) 
  microhydro(path = "MECHA/Projects/GRANAR/in/Maize_hydraulics.xml",
             kw = Sampl$kw[parse_number(j)],
             km = Sampl$km[parse_number(j)], 
             kAQP = Sampl$kAQP[parse_number(j)],
             kpl = Sampl$kpl[parse_number(j)])

  wallthick(path = "MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml", Sampl$thickness[parse_number(j)])

  # Run MECHA - - - - - - - 
  system("python3 ./MECHA/MECHAv4_septa.py")
  message("python script has ended")

  # if works well then
  if(file.exists("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt")){
  
    message ("success")
    file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_1,0_",id,".txt"), overwrite = T)
    file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_2,1.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_2,1_",id,".txt"), overwrite = T)
    file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_4,2.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_4,2_",id,".txt"), overwrite = T)
  }else{message ("fail and move to next simulation")}

}



