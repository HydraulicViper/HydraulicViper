# Author: Adrien Heymans
# input - output function to run hydraulic viper version 1.01
# Updated: 16-09-2022

# Get conductivities based on cell arragement and cell hydraulic prop from rf & svm algorithm
GM <- function(MatSobol, sam, shcut = TRUE, all_roots){

  load("./R/GRANAR/rf/svm_model_kr1.RData")
  model_kr1svm <- fit.svm
  
  load("./R/GRANAR/rf/rf_model_kr1.RData")
  model_kr1rf <- fit.rf


  load("./R/GRANAR/rf/svm_model_kr2.RData")
  model_kr2svm <- fit.svm
  load("./R/GRANAR/rf/rf_model_kr2.RData")
  model_kr2rf <- fit.rf


  load("./R/GRANAR/rf/svm_model_kr3.RData")
  model_kr3svm <- fit.svm 
  load("./R/GRANAR/rf/rf_model_kr3.RData")
  model_kr3rf <- fit.rf
  
  v_stele = MatSobol$stele[sam]/100
  v_xylem = MatSobol$xylem[sam]/100
  trans_time_apo = MatSobol$trans_time_apo[sam]
  trans_time_xyl = MatSobol$trans_time_xyl[sam]
  
  
  CT_root <- all_roots%>%
    dplyr::group_by(K_type)%>%
    dplyr::summarise(radius = mean(radius)*10)%>% # cm to mm
    ungroup()%>%
    mutate(log_RXA = log(pi*radius^2),
           var_stele = v_stele,
           var_xylem = v_xylem,
           RXA = exp(log_RXA),
           log_TSA = -1.421671 + 1.144070 * log_RXA, #Anatomy_proc.Rmd
           TSA = exp(log_TSA)+var_stele*exp(log_TSA),
           log_TSA = log(TSA, exp(1)),
           r_stele = sqrt(TSA/pi),
           log_nX = 2.9116240 + 0.4459025 * log_TSA, #Anatomy_proc.Rmd
           nX = exp(log_nX),
           ratio = (2+0.07456*r_stele*1000)/nX + var_xylem*(2+0.07456*r_stele*1000)/nX, #Anatomy_proc.Rmd
           mod_XVA = 1.395515 - 0.243006 * log_TSA, #Anatomy_proc.Rmd
           MXA = exp(-mod_XVA^2)+var_xylem*exp(-mod_XVA^2),
           log_CW = log(radius-r_stele, exp(1)),
           CF = exp(3.1091221+0.6718735*log_CW),
           OneC = exp(log_CW)/CF,
           oneC = OneC,
           nPX = nX*ratio,
           PXA_1 = 1000^2*(sqrt(radius/35)/10)^2,
           k_protxyl_s = PXA_1^2/(8*pi*200*1E-5/3600/24)*1E-12,
           kx_unM = k_protxyl_s*nPX*200/1E4, # kx when only the proto xylem have their cell wall lignified 
           LMXA = MXA - nPX*PXA_1/1000^2,
           LMXA_1 = LMXA*1000^2/nX,
           k_Mxyl_s = LMXA_1^2/(8*pi*200*1E-5/3600/24)*1E-12,
           kx_M = k_Mxyl_s*nX*200/1E4 + kx_unM, # kx when all xylem elements have their cell wall lignified 
           km = MatSobol$km[sam], 
           kw = MatSobol$kw[sam], 
           kAQP = MatSobol$aquaporine[sam],
           kpl = MatSobol$plasmo[sam], 
           thickness = 1.5,
           aerenchyma = MatSobol$aer[sam],
           TCA = RXA-TSA,
           AA = aerenchyma*TCA,
           vol = pi*radius^2-AA
    )
  
  
  to_analysis <- CT_root%>%
    select(radius, var_stele, var_xylem, aerenchyma, kw, km, kAQP, kpl, thickness, kx_unM, kx_M)%>%
    mutate(sampl_id = CT_root$K_type)
  
  prediction1 <- (predict(model_kr1svm,to_analysis%>%mutate(aerenchyma = 0))+ predict(model_kr1rf,to_analysis%>%mutate(aerenchyma = 0)))/2 
  prediction2 <- (predict(model_kr2svm,to_analysis)+predict(model_kr2rf,to_analysis))/2
  prediction3 <- (predict(model_kr3svm,to_analysis)+predict(model_kr3rf,to_analysis))/2
  data <- to_analysis%>%
    mutate(kr_1 = abs(prediction1),
           kr_2 = abs(prediction2),
           kr_3 = abs(prediction3))

  
  data$K_type <- data$sampl_id
  
  x <- read_xml("CPlantBox/modelparameter/rootsystem/B73.xml")
  elon <- tibble(type = 1:5, r = rep(0,5), la = rep(0,5))
  for(i in 2:6){
    tmp <- xml_children(x)[i]
    temp_r <- xml_find_all(tmp, ".//parameter")[12]
    r <- xml_attr(temp_r, "value")
    elon$r[i-1] = r
    temp_la <- xml_find_all(tmp, ".//parameter")[8]
    la <- xml_attr(temp_la, "value")
    elon$la[i-1] = la
  }
  
  elon$la <- as.numeric(elon$la)
  elon$r <- as.numeric(elon$r)
  
  conductivities <- NULL
  for(i in CT_root$K_type){
    
    la <- elon$la[all_roots$type[all_roots$K_type == i][1]]
    r <- elon$r[all_roots$type[all_roots$K_type == i][1]]
    
    x_kr = c(0, la/r, la/r+trans_time_apo, la/r+2*trans_time_apo, la/r+3*trans_time_apo, la/r+6*trans_time_apo)
    x_kx = sort(c(0, x_kr[3] + MatSobol$gap[sam],  x_kr[3] + MatSobol$gap[sam] + trans_time_xyl, x_kr[3]+MatSobol$gap[sam] + 5*trans_time_xyl))
    
    tmp_kr1 = data$kr_1[data$K_type == i]
    tmp_kr2 = data$kr_2[data$K_type == i]
    tmp_kr3 = data$kr_3[data$K_type == i]
    if(tmp_kr3 > tmp_kr2){tmp_kr3 = tmp_kr2}
    tmp_kxunM = data$kx_unM[data$K_type == i]
    tmp_kxM = data$kx_M[data$K_type == i]
    
    y = c(tmp_kr1, tmp_kr1,
          tmp_kr2, tmp_kr2,
          tmp_kr3, tmp_kr3/10,
          tmp_kxunM, tmp_kxunM, 
          tmp_kxM, tmp_kxM)
    
    tmp <- tibble(order_id = rep(i,10), order = rep(paste0("type_", i),10), type = c(rep("kr",6), rep("kx", 4)), x = c(x_kr,x_kx), y = y)
    
    conductivities <- rbind(conductivities, tmp)
  }
  
  tmp0 <- tibble(order_id = rep(0,4), order = paste0("type_", 0), type = c(rep("kr",2), rep("kx",2)), x = c(0, 10, 0, 10), 
                 y = c(data$kr_3[data$K_type == 1], data$kr_3[data$K_type == 1], data$kx_M[data$K_type == 1], data$kx_M[data$K_type == 1]))
  conductivities <- rbind(conductivities, tmp0)
  conductivities <- cbind(tibble(id = 1:nrow(conductivities)), conductivities)
  conductivities <- left_join(conductivities, CT_root %>%
                                                 mutate(order = paste0("type_", K_type))%>%
                                                 select(order, vol), by = "order")
  return(conductivities)
  
}
# function to give a conductivity type which is a function of the diameter
K_type <- function(all_roots){
  
  radius_latC <- all_roots$radius[all_roots$type == 2]
  radius_latA <- all_roots$radius[all_roots$type == 3]
  radius_B <- all_roots$radius[all_roots$type == 4]
  radius_SBR <- all_roots$radius[all_roots$type == 5]
  breaks_SBRtype = seq(min(radius_SBR), max(radius_SBR), 
                       by=(max(radius_SBR)-min(radius_SBR))/10)
  
  all_roots$K_type = NA
  all_roots$K_type[all_roots$type == 1] = 1
  all_roots$K_type[all_roots$type == 2 & all_roots$radius < quantile(radius_latC)[2]] = 2
  all_roots$K_type[all_roots$type == 2 & all_roots$radius >= quantile(radius_latC)[2] & all_roots$radius < quantile(radius_latC)[3]] = 3
  all_roots$K_type[all_roots$type == 2 & all_roots$radius >= quantile(radius_latC)[3] & all_roots$radius < quantile(radius_latC)[4]] = 4
  all_roots$K_type[all_roots$type == 2 & all_roots$radius >= quantile(radius_latC)[4]] = 5
  
  all_roots$K_type[all_roots$type == 3 & all_roots$radius < quantile(radius_latA)[2]] = 6
  all_roots$K_type[all_roots$type == 3 & all_roots$radius >= quantile(radius_latA)[2] & all_roots$radius < quantile(radius_latA)[3]] = 7
  all_roots$K_type[all_roots$type == 3 & all_roots$radius >= quantile(radius_latA)[3] & all_roots$radius < quantile(radius_latA)[4]] = 8
  all_roots$K_type[all_roots$type == 3 & all_roots$radius >= quantile(radius_latA)[4]] = 9
  
  all_roots$K_type[all_roots$type == 3 & all_roots$radius < quantile(radius_latA)[2]] = 6
  all_roots$K_type[all_roots$type == 3 & all_roots$radius >= quantile(radius_latA)[2] & all_roots$radius < quantile(radius_latA)[3]] = 7
  all_roots$K_type[all_roots$type == 3 & all_roots$radius >= quantile(radius_latA)[3] & all_roots$radius < quantile(radius_latA)[4]] = 8
  all_roots$K_type[all_roots$type == 3 & all_roots$radius >= quantile(radius_latA)[4]] = 9
  
  all_roots$K_type[all_roots$type == 4 & all_roots$radius < mean(radius_B)] = 10
  all_roots$K_type[all_roots$type == 4 & all_roots$radius == mean(radius_B)] = 11
  all_roots$K_type[all_roots$type == 4 & all_roots$radius > mean(radius_B)] = 12
  
  all_roots$K_type[all_roots$type == 5 & all_roots$radius < breaks_SBRtype[2]] = 13
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[2] & all_roots$radius < breaks_SBRtype[3]] = 14
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[3] & all_roots$radius < breaks_SBRtype[4]] = 15
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[4] & all_roots$radius < breaks_SBRtype[5]] = 16
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[5] & all_roots$radius < breaks_SBRtype[6]] = 17
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[6] & all_roots$radius < breaks_SBRtype[7]] = 18
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[7] & all_roots$radius < breaks_SBRtype[8]] = 19
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[8] & all_roots$radius < breaks_SBRtype[9]] = 20
  all_roots$K_type[all_roots$type == 5 & all_roots$radius > breaks_SBRtype[9] ] = 21
  
  
  return(all_roots)
}

# Create CplantBox param
CPlant_param <- function(MatSobol, sam, age_max = 30){
  
  a <- tibble(type = 1:5, radius = c(unlist(MatSobol[sam,"radius_sem"]),
                                   unlist(MatSobol[sam,"radius_lat"]),
                                   unlist(MatSobol[sam,"radius_lat"])*1.25,
                                   unlist(MatSobol[sam,"radius_sem"])*0.75,
                                   unlist(MatSobol[sam,"radius_shoot"])))%>%
        mutate(radius = radius/10) # cm instead of mm 

  x <- read_xml("./CPlantBox/modelparameter/rootsystem/B73.xml")

  for(i in 2:6){
    root <- xml_children(x)[i]
    temp <- xml_find_all(root, ".//parameter")
    if(xml_attr(temp[[3]], "name") != "a"){
      print("encoding miss placement in B73.xml")
    }
    xml_attr(temp[[3]], "value") <- a$radius[i-1]
  }
  pparam <- xml_children(x)[1]
  ptemp <- xml_find_all(pparam, ".//parameter")
  xml_attr(ptemp[[13]], "value") <- age_max

  write_xml(x, file = "./CPlantBox/modelparameter/rootsystem/B73.xml")
  
}

# # Make a CplantBox script
# make_CPscript <- function(age_max){
# 
#   CPscript <- paste0('"""small example in a container"""
# print("Importing libraries")
# import sys
# sys.path.append("./CPlantBox/")
# import plantbox as pb
# print("libraries loaded")
# 
# rootsystem = pb.RootSystem()
# # Open plant and root parameter from a file
# path = "./CPlantBox/modelparameter/rootsystem/"
# output = "',age_max,'_rootsystem"
# name = "B73"
# 
# print("Read parameter xml file")
# rootsystem.readParameters(path + name + ".xml")
# 
# # Create and set geometry
# # Creates a soil container
# print("Create soil container")
# rhizotron = pb.SDF_PlantBox(900, 900, 900)
# 
# # Pick 1, or 2
# print("Set geometry")
# rootsystem.setGeometry(rhizotron)  # soilcore, or rhizotron
# 
# # Initialize
# print("Initialize")
# rootsystem.initialize()
# 
# # Simulate
# print("Run -- CPlantBox -- ")
# rootsystem.simulate(',age_max,')  # days
# print("analyse")
# ana = pb.SegmentAnalyser(rootsystem)
# 
# print("write rootsytem")
# ana.write("{}.txt".format(str(output)))
#   ')
# write(CPscript, "./B73script.py")
# }

# Run the RootSystem model of CPlantBox
CRootBox <- function(sam, MatSobol, age_max = 30){
  
  all_roots <- NULL
  
  if(file.exists(paste0("./",age_max+3,"_rootsystem.txt"))){
    file.remove(paste0("./",age_max+3,"_rootsystem.txt"))
  }
  
  CPlant_param(MatSobol, sam) # change root radius in B37.xml file

  # make_CPscript(age_max+3)
  message("Run CplantBox")
  system("python ./B73script.py")

  fc <- file.copy(from = paste0("./",age_max+3,"_rootsystem.txt"),
                  to = paste0("./CPlantBox/RS/RS_",sam,".txt"),
                  overwrite = TRUE)
  
  all_roots <- data.table::fread(paste0("./",age_max+3,"_rootsystem.txt"), header = TRUE)

  all_roots <- all_roots%>%arrange(time)
  all_roots <- all_roots%>%dplyr::mutate(length = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2))
  all_roots$node2ID <- 1:nrow(all_roots)
  all_roots$node1ID[all_roots$branchID == 1][1] <- 0
  all_roots$node1ID[all_roots$branchID == 1][-1] <- which(all_roots$branchID == 1)[-length(which(all_roots$branchID == 1))] # tap root ordination

  for(i in unique(all_roots$branchID)[-1]){
     all_roots$node1ID[all_roots$branchID == i][-1] <- which(all_roots$branchID == i)[-length(which(all_roots$branchID == 1))]
     if(all_roots$type[all_roots$branchID == i][1] %in% c(4,5)){ # connection with the collar
       all_roots$node1ID[all_roots$branchID == i][1] <- 0
     }
     if(all_roots$type[all_roots$branchID == i][1] %in% c(2,3)){ # connection with the parental root
       x1_child <- all_roots$x1[all_roots$branchID == i][1]
       y1_child <- all_roots$y1[all_roots$branchID == i][1]
       z1_child <- all_roots$z1[all_roots$branchID == i][1]
       
       tmp_time <- all_roots$time[all_roots$branchID == i][1]
       
       nearest <- all_roots%>%filter(branchID != i)%>%
                      mutate(euc = sqrt((x1-x1_child)^2+ (y1 - y1_child)^2 + (z1 - z1_child)^2))
       nearest <- nearest[nearest$euc == min(nearest$euc), ]
       all_roots$node1ID[all_roots$branchID == i][1] <- nearest$node2ID[1] # oldest segments
     }
   }
  oups <- which(all_roots$node1ID == all_roots$node2ID)
  if(length(oups) > 0){
    for(o in oups){
      self_seg_age <- all_roots$time[all_roots$node2ID == o][1]
      self_seg_id <- all_roots$branchID[all_roots$node2ID == o][1]
      
      nearest <- all_roots%>%filter(branchID == self_seg_id, time < self_seg_age)
      all_roots$node1ID[all_roots$node2ID == o][1] <- nearest$node2ID[nearest$time == max(nearest$time)]
    }
  }
  # increasing shoot born root diameter as the numbers of nodes increase
  # first node have 100% of the shoot born radius defined in param
  # after 30 days --> SBR have a radius equals to 300% of the shoot born radius defined in param
  first_SBR <- min(age_max - all_roots$age[all_roots$type == 5])
  
  all_roots$radius[all_roots$type == 5] <- all_roots$radius[all_roots$type == 5]+all_roots$radius[all_roots$type == 5]*((age_max - all_roots$age[all_roots$type == 5] - first_SBR)/22)^1.4 
  # 
  pl = all_roots%>%ggplot()+geom_point(aes(age, radius, colour = type))
  print(pl)
  message("end of RSA formating")
  all_roots <- all_roots[all_roots$time <= age_max,]
  
  return(all_roots) 
}

# Compute the convex hull
CRoot_hull <- function(all_roots, age_filter = seq(10, age_max, by = 10), sam, C_cost){
  # pol_ch = poly_area = NULL
  ALL_ROOTS <- NULL
  
  for(ag in age_filter){
  
    tmp <- all_roots%>%
      filter(time <= ag)%>%
      mutate(age_RS = ag)
    
    ALL_ROOTS <- rbind(ALL_ROOTS, tmp)
  }
  ALL_ROOTS <- left_join(ALL_ROOTS, C_cost, by = "K_type")

  RLD <- ALL_ROOTS%>%
    dplyr::group_by(age_RS, rz1 = round(z1/2)*2, # on every two cm of the profile  
                    K_type)%>%
    dplyr::summarise(root_length = sum(length),
                     root_area = sum(length*radius),
                     root_vol = sum(length*pi*radius^2),
                     root_volumeAA = sum(length*vol),
                     rooting_depth = min(z2))%>%
    ungroup()
  return(RLD) # list(RLD, pol_ch, poly_area) 
}


read.tlevel.out <- function (project.path, out.file = "T_Level.out") {

  col = c("Time", "rTop", "rRoot", "vTop", "vRoot", "vBot", "sum_rTop", "sum_rRoot", "sum_vTop", 
                        "sum_vRoot", "sum_vBot", "hTop", "hRoot", "hBot", 
                        "RunOff", "sum_Runoff", "Volume", "sum_Infil", 
                        "sum_Evap", "TLevel", "Cum_WTrans", "SnowLayer", "PhSeq","PhLeaf")

  tlevel_out = fread(input = file.path(project.path, out.file), skip = 8, fill = TRUE,
                                 header = FALSE)
  colnames(tlevel_out) <- col
  tlevel_out = slice(tlevel_out, 1:(n() - 1)) # remove last row
  tlevel_out[is.nan(tlevel_out)] <- 0
  tlevel_out = as.data.table(tlevel_out)
  tlevel_out[, names(tlevel_out) := lapply(.SD, as.numeric)]
  return(tlevel_out)
}

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

Soil_type <- function(type = "loam"){
  
  
  soilparam <- read_excel("/home/users/a/h/aheymans/GSUA/www/Soil_type.xlsx")
  if(type %!in% soilparam$type){
    warning(paste0("soil type is not defined: tried ",soilparam$type))
  }
  soilparam <- soilparam[soilparam$type == type,]
  soilparam <- soilparam[1,] # only one row if repetition
  
  write.selector.in("./Hydrus/CouvreurV2", mesh = NA, soilparam = soilparam )
}

read.nod_inf <- function (project.path, out.file = "Nod_Inf.out", output = NULL, warn = FALSE, ...) {
  if (is.null(output) | missing(output)) {
    output = c("Head", "Moisture", "K", "C", "Flux", "Sink", 
               "Kappa", "v/KsTop", "Temp")
  }
  options(warn = -1)
  if (warn == TRUE) 
    options(warn = 0)
  nod_inf = data.table::fread(input = file.path(project.path, 
                                                out.file), fill = TRUE, blank.lines.skip = FALSE, skip = 9)
  time_lines = nod_inf[grepl("Time:", nod_inf[["Node"]]), 
                       ]
  times = c(0, as.numeric(time_lines$Depth))
  for (col in colnames(nod_inf)) set(nod_inf, j = col, value = as.numeric(nod_inf[[col]]))
  nod_inf = na.omit(nod_inf)
  nodes = sort(unique(nod_inf[["Node"]]))
  nod_inf[, `:=`(Time, rep(times, each = length(nodes)))]
  nod_split = split(nod_inf, f = nod_inf$Time)
  nrow_split = sapply(nod_split, nrow)
  extra_index = which(nrow_split > length(nodes))
  for (i in extra_index) {
    nod_split[[i]] = nod_split[[i]][1:length(nodes), ]
  }
  nod_inf = rbindlist(nod_split)
  output_names = intersect(output, colnames(nod_inf))
  output_names = c("Time", "Node", "Depth", output_names)
  nod_out = nod_inf[, .SD, .SDcols = output_names]
  options(warn = 0)
  return(nod_out)
}


soil_initial <- function(soil_type = stype, field_capacity = 0.5){
  if(soil_type == "sandy loam"){
      soil <- data.frame(id=1:101,
                         z = sort(seq(-100,0,1), decreasing = TRUE),
                         value = 1,
                         psi = rep(-73.5021, 101))
  }else if (soil_type == "loam") {
      soil <- data.frame(id=1:101,
                         z = sort(seq(-100,0,1), decreasing = TRUE),
                         value = 1,
                         psi = rep(-132.869, 101))

  }else if (soil_type == "silty clay") {
      soil <- data.frame(id=1:101,
                         z = sort(seq(-100,0,1), decreasing = TRUE),
                         value = 1,
                         psi = rep(-354.424, 101))
  }else{
    warning("This soil type is not referenced yet")
  }

  soil$psi <- -200 # soil$psi*(1/field_capacity)

  return(soil)
  
}

write.options.in <- function(project.path = "./hydrus/CouvreurV2", krs, kcomp){
  
  path <- paste0(project.path, "/Options.in")
  # path = "./Hydrus_1D_Couvreur/CouvreurV2/Options.IN"
  op = readLines(con = paste0(project.path,"/Options_1.IN"), n = -1L, encoding = "unknown")
  tmp_krs <- op[10]
  op[10]<- paste0(krs, collapse = " ")
  tmp_kcomp <- op[13]
  op[13] <- paste0(kcomp, collapse = " ")
  
  write(op, path)
}


write.profile.dat <- function(project.path = "./hydrus/CouvreurV2", SSF){
  
  out.file <- "PROFILE.DAT"
  prof = readLines(con = paste0(project.path,"/PROFILE_1.DAT"), n = -1L, encoding = "unknown")
  h_ind = grep("h", prof)
  n <- as.numeric(str_split(prof[h_ind], " ")[[1]][3])
  if(nrow(SSF) != n){message("sum up SSf every two centimeter along the profile")}
  ile = data.table::fread(input = file.path(paste0(project.path,"/PROFILE_1.DAT")),
                          fill = TRUE, 
                          blank.lines.skip = FALSE, 
                          skip = 4)
  nopcol <- rep("na",(length(ile[h_ind+1,])-9))
  for(i in 1:length(nopcol)){
    nopcol[i] <- paste0("nop",i)
  }
  colnames(ile) <- c("id", "z", "h", "Mat", "Lay", "Beta", "Axz", "Bxz", "Dxz", nopcol)
  
  ile <- ile%>%select(-starts_with("nop"))
  ile <- ile[1:101, ]
  ile$h <- SSF$h
  ile$Beta <- SSF$suf

  data_fmt <- cbind(ile,SSF[,3:26])
  data_fmt = apply(data_fmt, MARGIN = 1, FUN = paste0, collapse = " ")
  
  prof_input1 = prof[1:h_ind]
  prof_input2 = data_fmt
  prof_input3 = prof[(1+h_ind+length(prof_input2)):length(prof)]
  
  prof_input_new = c(prof_input1, prof_input2, prof_input3)
  prof_in_file = file.path(project.path, out.file)
  write(prof_input_new, file = prof_in_file, append = F)
}

write.selector.in<- function(project.path, mesh = NA, soilparam = NA){
  
  out.file = "SELECTOR.IN"
  # default.filename = "ATMOSPH.IN"
  selector_data = readLines(con = paste0(project.path,"/SELECTOR_1.IN"), n = -1L, encoding = "unknown")
  
  if(file.exists(file.path(project.path, out.file))){
    file.remove(file.path(project.path, out.file))
  }
  if(!is.na(soilparam)){
    
    selector_data[27]<- paste(as.character(soilparam[,c("Q_r","Q_s","alpha","n","Ks","l")]),collapse = " ")
  }
  if(!is.na(mesh)){
    selector_data[30]
    mesh_1 <- selector_data[30]
    selector_data[30]<- paste(replace(str_split(mesh_1, " ")[[1]], which(str_split(mesh_1, " ")[[1]]== "1"), mesh), collapse = " ")
    selector_data[30]<- paste(replace(str_split(mesh_1, " ")[[1]], which(str_split(mesh_1, " ")[[1]]== "2"), mesh), collapse = " ")
    selector_data[30]<- paste(replace(str_split(mesh_1, " ")[[1]], which(str_split(mesh_1, " ")[[1]]== "3"), mesh), collapse = " ")
  }
  selector_in_file = file.path(project.path, out.file)
  write(selector_data, file = selector_in_file, append = F)
  
}

write.atmosph.in<- function(project.path, maxAL, deltaT, atm.bc.data, hCritS = 0, ...){
  
  out.file = "ATMOSPH.IN"
  # default.filename = "ATMOSPH.IN"
  atm_data = readLines(con = paste0(project.path,"ATMOSPH_1.IN"), n = -1L, encoding = "unknown")
  
  if(file.exists(file.path(project.path, out.file))){
    file.remove(file.path(project.path, out.file))
  }
  extinction_ind = grep("Extinction", atm_data)
  
  
  # write(atm_data, file = "ATMOSPH_IN.BAK", append = F)
  
  hcrits_ind = grep("hCritS", atm_data)
  atm_data[hcrits_ind + 1] = sprintf("%7.0f", hCritS)
  
  maxAL_ind = grep("MaxAL", atm_data)
  tAtm_ind = grep("tAtm", atm_data)
  
  tMax = maxAL*deltaT
  
  atm_data[(maxAL_ind + 1)] = sprintf("%7.0f", maxAL)
  end_line = atm_data[grep("end", atm_data)]
  
  # bc_data = atm_data[(tAtm_ind +1): (end_line - 1)]
  # data_ind = (tMax*(sim_ind-1) + 1):(sim_ind*tMax)
  
  # tAtm = seq(deltaT, tMax, by = deltaT)
  
  bc_data_vars = c("tAtm", "Prec", "rSoil", "rRoot", "hCritA", "rB",
                   "hB", "ht", "RootDepth")
  
  bc_data_new = atm.bc.data[1:maxAL, bc_data_vars]
  # bc_data_new = data.frame(tAtm = seq(deltaT, tMax, deltaT), bc_data_new, row.names = NULL)
  #  bc_data_new = bc_data_new[rep(seq_len(nrow(bc_data_new)), each = 4), ]
  #  bc_data_new$tAtm = seq(deltaT, tMax, by = deltaT)
  row.names(bc_data_new) = NULL
  
  tstep_decimals = hydrusR::get.decimalplaces(deltaT)
  
  fmt_vec = c("%11.0f", "%12.3f", "%12.4f", "%12.4f", "%12.0f", rep("%12.4f",8))
  fmt_vec[1] = sub(pattern = "0", replacement = tstep_decimals, fmt_vec[1])
  
  bc_data_fmt = bc_data_new
  
  for(a in 1:nrow(bc_data_fmt)) {
    bc_data_fmt[a, ] = sprintf(fmt = fmt_vec[1:ncol(bc_data_fmt)], bc_data_new[a, ])
  }
  bc_data_fmt = apply(bc_data_fmt, MARGIN = 1, FUN = paste, collapse = "")
  
  atm_input1 = atm_data[1:tAtm_ind]
  atm_input2 = bc_data_fmt
  atm_input3 = end_line
  
  atmosph_input_new = c(atm_input1, atm_input2, atm_input3)
  atmosph_in_file = file.path(project.path, out.file)
  write(atmosph_input_new, file = atmosph_in_file, append = F)
  
}

create_soil <- function(soil){
  
  Beta <- rep(0, 101)
  tf = matrix(rep(0,101*26),ncol = 26)
  colnames(tf) <- c("suf","h", paste0("SSF_", 1:24) )

  SSF <- as_tibble(tf)%>%
            mutate(h = soil$psi)
  # overwirte the profile boundary condition
  write.profile.dat(project.path = "./hydrus/CouvreurV2", SSF)
  # message("profile.dat is correctly written")
  write.options.in(project.path = "./hydrus/CouvreurV2", rep(0,24), rep(0,24))
  # message("options.in is correctly written")

  atm_bc_data <- data.frame(tAtm = round(seq(1/24,1,1/24),4), Prec = 0, 
                                  rSoil = 0, rRoot = 0, hCritA = 15000, rB = 0, hB = 0, 
                                  ht = 0, RootDepth = 0)
        # overwrite the atmposheric boundary condition of hydrus.
  write.atmosph.in("./hydrus/CouvreurV2/", maxAL = 24,  deltaT = 1/24, atm_bc_data,
                   hCritS = 15000,  input.pet = F)
  # -----------------------------
  # Run HYDRUS
  # -----------------------------
  system("./h1d_calc")

  hydrus <- read.nod_inf(project.path = "./hydrus/CouvreurV2", 
                         out.file = paste0("Nod_Inf.out")) 
#  hydrus$SSF = rep(SSF$suf,25)
  hydrus <- hydrus%>%mutate(Ksrs = 0, id = Node, z = Depth, value = Time, psi = Head, moisture = Moisture)%>%
                    filter(Time %in% round(seq(1/24,1,1/24),4))
  out_data = read.tlevel.out(project.path = "./hydrus/CouvreurV2", out.file = paste0("T_Level.out"))%>%
                mutate(Ksrs = 0, Kcomp= 0, Tact_eq_MA =0,
                Krs = 0)
  return(list(hydrus, out_data))
}


Ccost <- function(MatSobol, sam, all_roots){
  
  v_stele = MatSobol$stele[sam]/100
  v_xylem = MatSobol$xylem[sam]/100
  trans_time_apo = MatSobol$trans_time_apo[sam]
  
  #-----------------------------
  CT_root <- all_roots%>%
    dplyr::group_by(K_type)%>%
    dplyr::summarise(radius = mean(radius)*10)%>% # cm to mm
    ungroup()%>%
    mutate(log_RXA = log(pi*radius^2),
           var_stele = v_stele,
           var_xylem = v_xylem,
           RXA = exp(log_RXA),
           log_TSA = -1.421671 + 1.144070 * log_RXA, #Anatomy_proc.Rmd
           TSA = exp(log_TSA)+var_stele*exp(log_TSA),
           log_TSA = log(TSA, exp(1)),
           r_stele = sqrt(TSA/pi),
           log_nX = 2.9116240 + 0.4459025 * log_TSA, #Anatomy_proc.Rmd
           nX = exp(log_nX),
           ratio = (2+0.07456*r_stele*1000)/nX, #Anatomy_proc.Rmd
           mod_XVA = 1.395515 - 0.243006 * log_TSA, #Anatomy_proc.Rmd
           MXA = exp(-mod_XVA^2)+var_xylem*exp(-mod_XVA^2),
           log_CW = log(radius-r_stele, exp(1)),
           CF = exp(3.1091221+0.6718735*log_CW),
           OneC = exp(log_CW)/CF,
           oneC = OneC,
           nPX = nX*ratio,
           PXA_1 = 1000^2*(OneC/2.2)^2,
           k_protxyl_s = PXA_1^2/(8*pi*200*1E-5/3600/24)*1E-12,
           kx_unM = k_protxyl_s*nPX*200/1E4, # kx when only the proto xylem have their cell wall lignified 
           LMXA = MXA - nPX*PXA_1/1000^2,
           LMXA_1 = LMXA*1000^2/nX,
           k_Mxyl_s = LMXA_1^2/(8*pi*200*1E-5/3600/24)*1E-12,
           kx_M = k_Mxyl_s*nX*200/1E4 + kx_unM, # kx when all xylem elements have their cell wall lignified 
           km = MatSobol$km[sam], 
           kw = MatSobol$kw[sam], 
           kAQP = MatSobol$aquaporine[sam],
           kpl = MatSobol$plasmo[sam], 
           thickness = 1.5,
           aerenchyma = MatSobol$aer[sam],
           TCA = RXA-TSA,
           AA = aerenchyma*TCA,
           vol = pi*radius^2-AA
    )
  #----------------------------
  
  x <- read_xml("CPlantBox/modelparameter/rootsystem/B73.xml")
  elon <- tibble(type = 1:5, r = rep(0,5), la = rep(0,5))
  for(i in 2:6){
    tmp <- xml_children(x)[i]
    temp_r <- xml_find_all(tmp, ".//parameter")[12]
    r <- xml_attr(temp_r, "value")
    elon$r[i-1] = r
    temp_la <- xml_find_all(tmp, ".//parameter")[8]
    la <- xml_attr(temp_la, "value")
    elon$la[i-1] = la
  }
  
  elon$la <- as.numeric(elon$la)
  elon$r <- as.numeric(elon$r)
  #-------------------------
  C_cost <- NULL
  for(i in CT_root$K_type){
    
    la <- elon$la[all_roots$type[all_roots$K_type == i][1]]
    r <- elon$r[all_roots$type[all_roots$K_type == i][1]]
    
    x = c(0, 2*la, 2*la+r*trans_time_apo, 2*la+2*r*trans_time_apo)
    y = c(CT_root$RXA[CT_root$K_type == i], CT_root$RXA[CT_root$K_type == i],
          CT_root$vol[CT_root$K_type == i], CT_root$vol[CT_root$K_type == i])
    
    tmp <- tibble(x = x, y = y, K_type = i)
    C_cost <- rbind(C_cost, tmp)
  }
  # ---------------------------
  
  order_uni=unique(CT_root$K_type)
  
  Carbon <- NULL
  age_max = max(all_roots$time)
  for(t in seq(0.2,age_max,0.1)){
    table_data <- all_roots%>%
      filter(time <= t)
    Nseg = nrow(table_data)
    order <- table_data$K_type
    seg_age <- table_data$time     # segment age
    seg_age <- max(table_data$time) - seg_age
    
    Cc=matrix(0,Nseg,1) # radial conductivity of the segments
    for ( i in 1:length(order_uni)) {
      
      pos = is.element(order,order_uni[i])
      od <- order_uni[i]
      
      x = C_cost$x[C_cost$K_type == od]
      y = C_cost$y[C_cost$K_type == od]
      x <- c(x, 5000)
      y <- c(y, y[length(y)])
      
      xout = data.frame(seg_age[pos])
      temp=data.frame(approx(x,y,xout[,1]))
      Cc[pos]=temp[,2]
      
    }
    
    table_data$Cc <- Cc[,1]
    table_data$Cc <- table_data$Cc*table_data$length*10 # mm3
    
    tmpC <- tibble(Cc = sum(table_data$Cc), age = t, sam = sam)
    Carbon <- rbind(Carbon, tmpC)
  }
  return(Carbon)
  
}

add_hydraulics <- function(temp_roots, hydraulics){
  
  # Merge output of MARSHAL on specific root segment
  temp_roots$suf <- as.vector(hydraulics$suf)
  temp_roots$suf_eq <- as.vector(hydraulics$suf_eq)
  temp_roots$suf1 <- as.vector(hydraulics$suf1)
  temp_roots$kx <- as.vector(hydraulics$kx)
  temp_roots$kr <- as.vector(hydraulics$kr)
  temp_roots$jr <- as.vector(hydraulics$jr)
  temp_roots$jr_eq <- as.vector(hydraulics$jr_eq)
  temp_roots$jxl_eq <- as.vector(hydraulics$jxl_eq)
  temp_roots$psi_eq <- as.vector(hydraulics$psi_eq)
  temp_roots$psi <- as.vector(hydraulics$psi)
  temp_roots$jxl <- as.vector(hydraulics$jxl)
  temp_roots$psi_soil <- as.vector(hydraulics$psi_soil)
  return(temp_roots)
}


HM_RWUM <- function(weather, all_roots, sam, conductivities, Envi, stype){

    age_max <- round(max(all_roots$time))
    # overwrite the LEVEL_01.DIR file with the right directory
    hydrus_dir <- paste0(getwd(),"/hydrus/CouvreurV2")
    write(hydrus_dir, "LEVEL_01.DIR")

    soil_param <- read_excel("./www/Soil_type.xlsx")%>%
      filter(type == stype)%>%
      mutate(Ksat = Ks, lambda = l)

    time_sequence = 3:round(age_max)
    
    dt = 1
    SOIL = OUT = NULL
    
    # Initial soil condition
    soil <- soil_initial(soil_type = stype)
    resting = 5 # waiting time
    for(i in 1:resting){
      ini_soil = create_soil(soil)
      hydrus <- ini_soil[[1]]
      soil <- hydrus%>%
            filter(Time == 1)%>%transmute(id = Node, z = Depth, value = i-resting, psi = Head)
      SOIL <- rbind(SOIL, hydrus%>%
                      filter(Time != 0)%>%
                      mutate(Time = sort(rep(round(seq(1/24,1,1/24),4),101))+i-resting-1))  
      out_data <- ini_soil[[2]]
      out_data$Time = out_data$Time + i -resting-1
      
      if(!is.null(OUT)){
        out_data$sum_rRoot = out_data$sum_rRoot + OUT$sum_rRoot[nrow(OUT)]
        out_data$sum_rTop = out_data$sum_rTop + OUT$sum_rTop[nrow(OUT)]
        out_data$sum_vRoot = out_data$sum_vRoot + OUT$sum_vRoot[nrow(OUT)]
        out_data$sum_vTop = out_data$sum_vTop + OUT$sum_vTop[nrow(OUT)]
        out_data$sum_Runoff = out_data$sum_Runoff + OUT$sum_Runoff[nrow(OUT)]
        out_data$sum_vBot = out_data$sum_vBot + OUT$sum_vBot[nrow(OUT)]
        out_data$sum_Infil = out_data$sum_Infil + OUT$sum_Infil[nrow(OUT)]
        out_data$sum_Evap = out_data$sum_Evap + OUT$sum_Evap[nrow(OUT)]
        out_data$Cum_WTrans = out_data$Cum_WTrans + OUT$Cum_WTrans[nrow(OUT)]
      }
      
      OUT <- rbind(OUT, out_data)
    }
    
    tmp_soil = hydrus%>%transmute(id = Node, z = Depth, value = Time, psi = Head)
    tpots = rep(-15000,24)

    for(t in time_sequence){

      message(paste0("----- ", t, " -----"))
      # Update root system
      temp_roots <- all_roots%>%
        filter(time <= t)
      conti = TRUE
      ntry = -10
      score_old = score = 10^8
      if (t == min(time_sequence)){ntry = -1} # extra trials at the first steps
      if (weather[t]> 0 | t >= 35){ntry = -3} # extra trials during raining day and soil is dry 
      Negat_T = -10

      while(conti){
        HM <- HM_loop(temp_roots, conductivities, tpots, soil_param, soil,tmp_soil,t, Envi, weather)
        hydrus <- HM[[1]]
        tmp_soil <- hydrus%>%transmute(id = Node, z = Depth, value = Time, psi = Head)
        tmp_score = HM[[3]]
        out_data = HM[[2]]
        tpots = out_data$PhLeaf
        tpots = -abs(tpots)
        tpots[tpots > out_data$PhSeq] = out_data$PhSeq[out_data$PhLeaf > out_data$PhSeq]+out_data$PhSeq[out_data$PhLeaf > out_data$PhSeq]*0.01

        message(paste0("score = ", tmp_score))

        if(tmp_score < 1.5){conti = FALSE} # quite good (minize error)
        if(tmp_score < score_old & length(out_data$vRoot[out_data$vRoot <= 0])== 0){ # save the best run
          case_out = out_data
          case_soil = tmp_soil
        }
        ntry = ntry +1
        if(round(tmp_score*10) == score){conti = FALSE} # converge in a local minimum 
        if(round(tmp_score*10) == score_old){conti = FALSE} # converge in a local minimum
        score_old = score
        score = round(tmp_score*10)
        if (ntry > 7){
          conti = FALSE
          out_data = case_out  # use best run instead
          tmp_soil = case_soil
        }

        if(length(out_data$vRoot[out_data$vRoot <= 0])> 0){
          conti = TRUE
          Negat_T = Negat_T+1
          message("WARNING ! NEGATIVE TRANSPIRATION SPOTED")}else{Negat_T = 0}
        if(Negat_T > 5){
          conti = FALSE
          out_data = case_out # use best run instead
          tmp_soil = case_soil
          }
        
      }

      soil = tmp_soil  %>% filter(value == max(tmp_soil$value))
      SOIL <- rbind(SOIL, hydrus%>%
                      filter(Time != 0)%>%
                      mutate(Time = sort(rep(round(seq(1/24,1,1/24),4),101))+t-3))  
      

      out_data$Time = out_data$Time + t-3
      out_data$sum_rRoot = out_data$sum_rRoot + OUT$sum_rRoot[nrow(OUT)]
      out_data$sum_rTop = out_data$sum_rTop + OUT$sum_rTop[nrow(OUT)]
      out_data$sum_vRoot = out_data$sum_vRoot + OUT$sum_vRoot[nrow(OUT)]
      out_data$sum_vTop = out_data$sum_vTop + OUT$sum_vTop[nrow(OUT)]
      out_data$sum_Runoff = out_data$sum_Runoff + OUT$sum_Runoff[nrow(OUT)]
      out_data$sum_vBot = out_data$sum_vBot + OUT$sum_vBot[nrow(OUT)]
      out_data$sum_Infil = out_data$sum_Infil + OUT$sum_Infil[nrow(OUT)]
      out_data$sum_Evap = out_data$sum_Evap + OUT$sum_Evap[nrow(OUT)]
      out_data$Cum_WTrans = out_data$Cum_WTrans + OUT$Cum_WTrans[nrow(OUT)]
      

      OUT <- rbind(OUT, out_data)
    }
    message("End of HM_RWUM")
    return(list(SOIL, OUT))

}


HM_loop <- function(temp_roots, conductivities, tpots, soil_param, soil, tmp_soil, t, Envi, weather){

  KRS = TACT_MA = KSRS = KCOMP = HSEQ = rep(NA,24)
  vec_i = 1
  temp_soil1 <- tmp_soil[tmp_soil$z %in% seq(-100, 0 ,1),c(1,2,4,3)]
  # pb <- txtProgressBar(min = 1, max = 24, style = 3,  width = 24, char = "=")
  SSF_M = matrix(rep(0,2626),ncol = 26)
  colnames(SSF_M) <- c("h","suf", paste0("SSF_", 1:24))
  SSF_M <- as_tibble(SSF_M)

  for(si in unique(temp_soil1$value)){
      # Select specific Soil for the simulation
    temp_soil <- temp_soil1%>%filter(value == si)
    theta_top <- tmp_soil[1,"moisture"]
    if(length(theta_top) == 0){
       theta_top = 0.0866 # initial condition (-300 hPa)
    }
    h_top <- temp_soil[1,"psi"]
    # -----------------------------
    # Run MARSHAL
    # -----------------------------
    mp_roots = temp_roots%>%filter(time <= t-1+si)

    hydraulics <- getSUF(mp_roots, conductivities, temp_soil, hetero =TRUE, 
                         Psi_collar = tpots[vec_i], soil_param)
    # joint the hydraulic prop on the architecture
    mp_roots <- add_hydraulics(mp_roots, hydraulics)
    ksrs <-  hydraulics$ksrs/15/75 # cm4 hPa-1 d-1
    kcomp <- ksrs
  
    Beta <- rep(0, 101)
    RLDWU <- mp_roots%>% # gather information by layer
      mutate(rz2 = round((z1+z2)/2))%>%
      dplyr::group_by(rz2)%>%
      dplyr::summarise(su_eq = sum(suf_eq), jr_eq = sum(jr_eq), jx_eq = sum(jxl_eq))%>%
      ungroup()
    Beta[which(soil$z %in% RLDWU$rz2 )] <- rev(RLDWU$su_eq) # from above to below
    #SSF <- data.frame(suf = Beta, h = soil$psi)
    if (si == 0.5){
      SSF_M$h = soil$psi
      SSF_M$suf = Beta
    }
  
    Q_dou <- rev(RLDWU$jr_eq)/75/15 # 8 plant per sqare meter, then convert from mL*m-2*day-1 to cm*m-2*day-1
    SUF = rev(RLDWU$su_eq)
    Tact_eq = sum(Q_dou)
    Hsr <- temp_soil$psi[temp_soil$z >= min(RLDWU$rz2) & temp_soil$z <= max(RLDWU$rz2)]# - (100 -temp_soil$z[temp_soil$z >= min(RLDWU$rz2) & temp_soil$z <= max(RLDWU$rz2)])
    
    if(length(unique(Hsr))> 1){
      up = tibble(s = SUF,j = Q_dou,h = Hsr)%>%mutate(Hseqi = s*h, as = s*Tact_eq, upper = j - as)
      Hseq <- sum(up$Hseqi,na.omit = TRUE)
      up <- up%>%mutate(bel = h -Hseq, belo = (bel*s)^(-1))
      kcomp = up$upper%*%(up$belo)
    }else{
      Hseq <- Hsr[1]
    }
    kcomp <- kcomp[1]
    if (kcomp <= 0) {kcomp = 1E-8}
    # setTxtProgressBar(pb, vec_i)
    KSRS[vec_i] = ksrs
    KRS[vec_i] = hydraulics$krs/75/15
    KCOMP[vec_i] = kcomp
    TACT_MA[vec_i] = hydraulics$tact_eq/75/15
    HSEQ[vec_i] = Hseq
    SSF_M[,vec_i+2] = rev(Beta) 
    vec_i = vec_i +1
  } 

  # overwirte the profile boundary condition
  write.profile.dat(project.path = "./hydrus/CouvreurV2", SSF_M)
  # message("profile.dat is correctly written")
  write.options.in(project.path = "./hydrus/CouvreurV2", KSRS, KCOMP)
  # message("options.in is correctly written")

 # soil evaporation 
  T_Cel <- 20 # assumption, soil top layer is at 20?C
  rho_vs <- 10^-3*exp(31.3716-6014.79/(T_Cel+273)-7.92*10^-3*(T_Cel+273))
  E_tmp <- Envi[t,]%>%
        mutate(r_s2 = 10*exp(35.63*(0.15-theta_top)),
                Hr_top = exp(h_top*0.018015*9.81/(1.01*0.287*(T_Cel+273))),
                rho_v = rho_vs*Hr_top,
                E1 = (rho_v-rho_a)/(r_a - (r_s2)),
                E2 = ET-Transpi)
  soil_evap = E_tmp$E1[1]
  if(soil_evap < 0){soil_evap = 0}
     # avoid large evaporation rates
  if(soil_evap > 0.5){soil_evap = 0.5 }
  
  Es <-sin(round(seq(1/24,1,1/24),4)*pi)*soil_evap
  Transpiration = Envi$Transpi[t]
  if(T <= 0){T = 1E-9}
  Tj = Transpiration

  # ratio_T = ksrs/(hydraulics$krs/15/75)
  # if no stress 0.0551*exp(1) = 0.15 --> 15% of day transpi
  Tn = Tj*0.15
  # Transpiration*0.0551*exp(ratio_T)  # less night leakage if stress
  cfact = Tj*24/sum(sin(round(seq(1/24,1,1/24),4)*pi)^5*Tj)
  Tpot = sin(round(seq(1/24,1,1/24),4)*pi)^5*Tj*cfact
  Tpot[Tpot <= Tn] = abs(Tn)
  atm_bc_data <- data.frame(tAtm = round(seq(1/24,1,1/24),4), Prec = 0, # rep(weather[t],24)
                                  rSoil = Es, rRoot = Tpot , hCritA = 15000, rB = 0, hB = 0, 
                                  ht = 0, RootDepth = 0)
  if (t >= 40){
    atm_bc_data <- data.frame(tAtm = round(seq(1/24,1,1/24),4), Prec = rep(weather[t],24),
                                  rSoil = Es, rRoot = Tpot , hCritA = 15000, rB = 0, hB = 0, 
                                  ht = 0, RootDepth = 0)
  }
        # overwrite the atmposheric boundary condition of hydrus.
  write.atmosph.in("./hydrus/CouvreurV2/", maxAL = 24,  deltaT = 1/24, atm_bc_data,
                   hCritS = 15000,  input.pet = FALSE)
  # -----------------------------
  # Run HYDRUS
  # -----------------------------
  system("./h1d_calc2401")

  hydrus <- read.nod_inf(project.path = "./hydrus/CouvreurV2", 
                         out.file = paste0("Nod_Inf.out"))
  if(max(hydrus$Time) < 1){ # Happens only if numerical solution has not converged whitin HYDRUS !
    hydrus <- soil
  }

  hydrus <- hydrus%>%mutate(Ksrs = ksrs, id = Node, z = Depth, value = Time, psi = Head, moisture = Moisture)%>%
                    filter(Time %in% round(seq(1/24,1,1/24),4))
  out_data = read.tlevel.out(project.path = "./hydrus/CouvreurV2", out.file = paste0("T_Level.out"))%>%
                filter(Time %in% round(seq(1/24,1,1/24),4))
  exclude = which(duplicated(out_data$Time))
  out_data$id_row = 1:nrow(out_data)
  out_data <- out_data %>% filter(id_row %!in% exclude)%>%select(-id_row)%>%
                          mutate(Ksrs = KSRS, Kcomp= KCOMP, Tact_eq_MA = TACT_MA, Krs = KRS )
  compa = tibble(hy = out_data$PhSeq, ma = HSEQ, 
                 hyl = out_data$PhLeaf, mal = tpots, x = 1:24)%>%
            mutate(square = (hy-ma)^2,
                   squarel = (hyl-mal)^2,
                   square2 = (out_data$Tact_eq_MA-out_data$vRoot)^2)
  if(length(out_data$vRoot[out_data$vRoot <= 0])> 0){
   pl <- compa%>%ggplot()+
               geom_line(aes(x,hy), colour = "blue")+
               geom_line(aes(x,ma), colour = "blue", alpha = 0.5)+
               geom_line(aes(x,hyl),colour = "red")+
               geom_line(aes(x,mal),colour = "red", alpha = 0.5) 
   print(pl)
  }

  score = sqrt(sum(compa$square/24))/abs(mean(c(compa$hy,compa$ma)))
  scorel = sqrt(sum(compa$squarel/24))/abs(mean(c(compa$hyl,compa$mal)))
  score2 = sqrt(sum(compa$square2/24))/abs(mean(c(out_data$Tact_eq_MA, out_data$vRoot)))
  sc = (sum(compa$square/24)+sum(compa$squarel/24))/2
  message(paste0("Hseq err: ",score*100," Hleaf err: ", scorel*100, " Trans err: ", score2*100))
  return(list(hydrus, out_data, sc = 100*(score+scorel+score2)/3))
}

