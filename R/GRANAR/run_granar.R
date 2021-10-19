

run_granar <- function(params, tmp_sampl, i = NULL){
  
  if(file.exists("./MECHA/cellsetdata/current_root.xml")){
    file.remove("./MECHA/cellsetdata/current_root.xml")
    file.remove("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml")
  }
  
  params$value[params$type == "cell_diameter" & params$name == "stele"] <- tmp_sampl$OneC/2.2
  params$value[params$type == "cell_diameter" & params$name == "pericycle"] <- tmp_sampl$OneC/2.2
  params$value[params$type == "cell_diameter" & params$name == "endodermis"] <- tmp_sampl$OneC/1.3
  params$value[params$type == "cell_diameter" & params$name == "cortex"] <- tmp_sampl$OneC
  params$value[params$type == "cell_diameter" & params$name == "exodermis"] <- tmp_sampl$OneC
  params$value[params$type == "cell_diameter" & params$name == "epidermis"] <- tmp_sampl$OneC/3
  params$value[params$type == "n_layers" & params$name == "cortex"] <- round(tmp_sampl$CF-2)
  params$value[params$type == "layer_diameter" & params$name == "stele"] <- tmp_sampl$r_stele*2
  params$value[params$type == "max_size" & params$name == "xylem"] <- tmp_sampl$X_size -2*params$value[params$type == "cell_diameter" & params$name == "stele"]
  params$value[params$type == "n_files" & params$name == "xylem"] <- round(tmp_sampl$nX)
  params$value[params$type == "ratio" & params$name == "xylem"] <- (2+0.07456*tmp_sampl$r_stele*1000)/tmp_sampl$nX
  
  
  params$value[params$type == "proportion" & params$name == "aerenchyma"] <- tmp_sampl$aerenchyma
  params$value[params$type == "n_files" & params$name == "aerenchyma"] <- ceiling(tmp_sampl$aerenchyma*tmp_sampl$TCA/((tmp_sampl$radius - tmp_sampl$r_stele)*3*tmp_sampl$OneC))
  
  if(tmp_sampl$radius > 0.7){
    tmp <- tibble(name = c("pith", "pith"), type = c("cell_diameter", "layer_diameter"), value = c(tmp_sampl$OneC, (tmp_sampl$r_stele-2*params$value[params$type == "max_size" & params$name == "xylem"])*2*0.8))
    params <- rbind(params, tmp)
  }
  
  sim <- create_anatomy(parameters = params)
  write_anatomy_xml(sim, "./MECHA/cellsetdata/current_root.xml")
  if(length(sim$id_aerenchyma) > 0){
    aer_in_geom_xml(sim, "./MECHA/Projects/GRANAR/in/Maize_Geometry.xml")
  }else{
    fc <- file.copy("./MECHA/Projects/GRANAR/in/Maize_Geometry_noAA.xml", "./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml", overwrite = T )
  }
  
  if(!is.null(i)){
    file.copy("./MECHA/cellsetdata/current_root.xml", paste0("./MECHA/cellsetdata/root_",i,".xml"), overwrite = T)
    file.copy("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml",
            paste0( "./MECHA/Projects/GRANAR/in/Maize_Geometry_aer_",i,".xml"), overwrite = T)
  }else{return(sim)}
  
}

run_Mecha <- function(tmp_sampl, i = NULL){
  
  if(file.exists("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt")){
    file.remove("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt")
    file.remove("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_2,1.txt")
    file.remove("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_4,2.txt")
    file.remove("./MECHA/cellsetdata/current_root.xml")
    file.remove("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml")
  }
  
  fc <- file.copy(paste0("./MECHA/cellsetdata/root_",i,".xml"), "./MECHA/cellsetdata/current_root.xml", overwrite = T)
  file.copy(paste0("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer_", i, ".xml"),
            paste0("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml"), overwrite = T)
  
  microhydro(path = "MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml",
             kw = tmp_sampl$kw,
             km = tmp_sampl$km, 
             kAQP = tmp_sampl$kAQP,
             kpl = tmp_sampl$kpl)
  
  wallthick(thickness = tmp_sampl$thickness)
  
  # Run MECHA - - - - - - - 
  
  # system("C:/Users/heymansad/AppData/Local/Continuum/anaconda3/envs/MECHA/python.exe MECHA/MECHAv4_GRANAR.py")
  system("C:/Users/heymansad/AppData/Local/Continuum/anaconda3/envs/MECHA/python.exe MECHA/MECHAv4_septa.py")
  message("python script has ended")
  
  # if works well then
  if(file.exists("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt")){
    
    message ("success")
    file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_1,0_",i,".txt"), overwrite = T)
    file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_2,1.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_2,1_",i,".txt"), overwrite = T)
    file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_4,2.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_4,2_",i,".txt"), overwrite = T)
  }else{message ("fail and move to next simulation")}
  
}
