
set_params <- function(tiss_param, layer_param, params){
  
  # stele cell size
  c_stele <- tiss_param$one[tiss_param$tiss=="stele"]+d*tiss_param$slope[tiss_param$tiss =="stele"]
  params$value[params$name == "stele" & params$type == "cell_diameter"] <- c_stele/1000
  # stele layer size
  stele <- layer_param$one[layer_param$lay =="stele"]+d*layer_param$slope[layer_param$lay =="stele"]
  params$value[params$name == "stele" & params$type == "layer_diameter"] <- stele/1000
  # pericycle cell size
  c_pericycle <- tiss_param$one[tiss_param$tiss =="pericycle"]+d*tiss_param$slope[tiss_param$tiss =="pericycle"]
  params$value[params$name == "pericycle" & params$type == "cell_diameter"] <- c_pericycle/1000
  # endodermis cell size
  c_endodermis <- tiss_param$one[tiss_param$tiss =="endo"]+d*tiss_param$slope[tiss_param$tiss =="endo"]
  params$value[params$name == "endodermis" & params$type == "cell_diameter"] <- c_endodermis /1000
  # cortex cell size
  c_cortex <- tiss_param$one[tiss_param$tiss =="cortex"]+d*tiss_param$slope[tiss_param$tiss =="cortex"]
  params$value[params$name == "cortex" & params$type == "cell_diameter"] <- c_cortex/1000
  # cortex n_layers
  cortex_w <- layer_param$one[layer_param$lay =="cortex"]+d*layer_param$slope[layer_param$lay =="cortex"]
  n_layer <- round(cortex_w/c_cortex)
  params$value[params$name == "cortex" & params$type == "n_layers"] <- n_layer
  # exodermis cell size
  c_exodermis <- tiss_param$one[tiss_param$tiss =="exo"]+d*tiss_param$slope[tiss_param$tiss =="exo"]
  params$value[params$name == "exodermis" & params$type == "cell_diameter"] <- c_exodermis/1000
  # epidermis cell size 
  c_epidermis <- tiss_param$one[tiss_param$tiss =="epi"]+d*tiss_param$slope[tiss_param$tiss =="epi"]
  params$value[params$name == "epidermis" & params$type == "cell_diameter"] <- c_epidermis/1000
  
  # n metaxylem vessels 
  nX <- layer_param$one[layer_param$lay =="xylem"]+d*layer_param$slope[layer_param$lay =="xylem"]
  params$value[params$name == "xylem" & params$type == "n_files"] <- round(nX)
  # xylem size 
  c_xylem <- tiss_param$one[tiss_param$tiss =="xylem"]+d*tiss_param$slope[tiss_param$tiss =="xylem"]
  params$value[params$name == "xylem" & params$type == "max_size"] <- c_xylem/1000-c_stele/1000
  # ratio proto/meta-xylem
  n_proto <- layer_param$one[layer_param$lay =="proto"]+d*layer_param$slope[layer_param$lay =="proto"]
  params$value[params$name == "xylem" & params$type == "ratio"] <- n_proto/nX
  
  if(length(unique(layer_param$lay == "aer"))>1) {
    # n files of aerenchyma
    n_files <- layer_param$one[layer_param$lay == "aer"]+d*layer_param$slope[layer_param$lay == "aer"]
    params$value[params$name == "aerenchyma" & params$type == "n_files"] <- round(n_files)
    # aerenchyma area
    aer_area <- tiss_param$one[tiss_param$tiss =="aer"]+d*tiss_param$slope[tiss_param$tiss =="aer"]
    aer_area <- aer_area*n_files
    
    # cortex area
    stele_r <- stele/2 + c_pericycle + c_endodermis
    R <- stele_r + cortex_w + c_exodermis + c_epidermis
    RXA <- pi*R^2
    SXA <- pi*stele_r^2
    CXA <- RXA-SXA
    # aerenchyma proportion
    prop <- aer_area/CXA
    params$value[params$name == "aerenchyma" & params$type == "proportion"] <- prop
    
  }
  return(params)
}
