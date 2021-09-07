
GM <- function(MatSobol, sam, shcut = T, all_roots){
  
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
  
  x <- read_xml("CPlantBox-master/modelparameter/rootsystem/B73.xml")
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
    
    x_kr = c(0, 2*la, 2*la+r*trans_time_apo, 2*la+2*r*trans_time_apo, 2*la+3*r*trans_time_apo, 2*la+6*r*trans_time_apo)/r
    x_kx = sort(c(0, x_kr[3] + x_kr[3]+r*MatSobol$gap[sam],  x_kr[3] + x_kr[3]+r*MatSobol$gap[sam] + r*trans_time_xyl, x_kr[3]+MatSobol$gap[sam] + 5*r*trans_time_xyl))/r
    
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
# function to give a type as a function of the diameter
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

CPlant_param <- function(MatSobol, sam, age_max = 30){
  
  a <- tibble(type = 1:5, radius = c(unlist(MatSobol[sam,"radius_sem"]),
                                   unlist(MatSobol[sam,"radius_lat"]),
                                   unlist(MatSobol[sam,"radius_lat"])*1.25,
                                   unlist(MatSobol[sam,"radius_sem"])*0.75,
                                   unlist(MatSobol[sam,"radius_shoot"])))%>%
        mutate(radius = radius/10) # cm instead of mm 

  x <- read_xml("./CPlantBox-master/modelparameter/rootsystem/B73.xml")

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

  write_xml(x, file = "./CPlantBox-master/modelparameter/rootsystem/B73.xml")
  
}


make_CPscript <- function(age_max){

  CPscript <- paste0('"""small example in a container"""
print("Importing libraries")
import sys
sys.path.append("../../..")
import plantbox as pb
print("libraries loaded")
rootsystem = pb.RootSystem()

# Open plant and root parameter from a file
path = "../../../modelparameter/rootsystem/"
output = "',age_max,'_rootsystem"
name = "B73"
print("Read parameter xml file")
rootsystem.readParameters(path + name + ".xml")

# Create and set geometry
# Creates a soil container
print("Create soil container")
rhizotron = pb.SDF_PlantBox(900, 900, 900)

# Pick 1, or 2
print("Set geometry")
rootsystem.setGeometry(rhizotron)  # soilcore, or rhizotron

# Initialize
print("Initialize")
rootsystem.initialize()

# Simulate
print("Run -- CPlantBox -- ")
rootsystem.simulate(',age_max,')  # days
print("analyse")
ana = pb.SegmentAnalyser(rootsystem)
print("write rootsytem")
ana.write("{}.txt".format(str(output)))
  ')
  write(CPscript, "./CPlantBox-master/tutorial/examples/python/B73script.py")
}

CRootBox <- function(sam, MatSobol, age_max = 30){
  
  all_roots <- NULL
  
  if(file.exists(paste0("./CPlantBox-master/tutorial/examples/python/",age_max+3,"_rootsystem.txt"))){
    file.remove(paste0("./CPlantBox-master/tutorial/examples/python/",age_max+3,"_rootsystem.txt"))
  }
  
  CPlant_param(MatSobol, sam) # change root radius in B37.xml file

  make_CPscript(age_max+3)
  
  system(paste0("cd ./CPlantBox-master/tutorial/examples/python && python3 B73script.py"))
  
  fc <- file.copy(from = paste0("./CPlantBox-master/tutorial/examples/python/",age_max+3,"_rootsystem.txt"),
                  to = paste0("./CPlantBox-master/tutorial/examples/python/RS_",sam,".txt"),
                  overwrite = T)
  
  all_roots <- fread(paste0("./CPlantBox-master/tutorial/examples/python/RS_",sam,".txt"), header = T)
  all_roots <- as.data.table(all_roots)
  all_roots <- all_roots%>%
    transmute(node1ID = node1ID,
              node2ID = node2ID,
              branchID = branchID,
              x1 = x1, y1 = y1, z1 = z1, x2 = x2, y2= y2, z2 = z2,
              radius = radius,
              length = sqrt((all_roots$x2 - all_roots$x1)^2 + (all_roots$y2 - all_roots$y1)^2 + (all_roots$z2 - all_roots$z1)^2),
              R = R,
              G = G,
              B = B,
              time = time,
              type = type,
              age = age)%>%
    arrange(time)

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
  
  # increasing shoot born root diameter as the numbers of nodes increase
  # first node have 100% of the shoot born radius defined in param
  # after 30 days --> SBR have a radius equals to 300% of the shoot born radius defined in param
  first_SBR <- min(age_max - all_roots$age[all_roots$type == 5])
  
  all_roots$radius[all_roots$type == 5] <- all_roots$radius[all_roots$type == 5]+all_roots$radius[all_roots$type == 5]*(age_max - all_roots$age[all_roots$type == 5] - first_SBR)/10 
  # 

  all_roots <- all_roots[all_roots$time <= age_max,]
  
  return(all_roots) 
}

CRoot_hull <- function(all_roots, age_filter = seq(10, age_max, by = 10), sam, C_cost){
  # pol_ch = poly_area = NULL
  ALL_ROOTS <- NULL
  
  for(ag in age_filter){
    # for(xy in c("x", "y")){
    #   tmp <- all_roots%>%
    #     filter(time <= ag)%>%
    #     select(paste0(xy,"2"),z2)
    #   first_col <- tmp[,1]
    #   x <- matrix(c(first_col, tmp$z2), nc = 2)
    #   # Get the indices of the unique points lying on the convex hull
    #   ch <- chull(x= first_col, y = tmp$z2)
    #   ch <- c(ch, ch[1])
    #   
    #   # Area of the convex hull
    #   xy.coords <- cbind(first_col, tmp$z2)
    #   chull.coords <- xy.coords[ch,]
    #   chull.poly <- sp::Polygon(chull.coords, hole = F)
    #   chull.area <- chull.poly@area
    #   
    #   polyg <- cbind(as.data.frame(x[ch,]), sam, rep = paste0(xy), age = ag)
    #   pol_area <- c(chull.area/10000)%>% # Transformation into square meter
    #     as.tibble()%>%
    #     transmute(area = value, rep= paste0(xy),sam, age = as.numeric(ag))
    #   poly_area <- rbind(poly_area, pol_area)
    #   pol_ch <- rbind(pol_ch, polyg)
    # }
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



HM_RWUM <- function(weather, all_roots, sam, conductivities, ET, E_soil, stype, dt = 0.2){
  age_max <- max(all_roots$time)

  # overwrite the LEVEL_01.DIR file with the right directory
  hydrus_dir <- paste0(getwd(),"/Hydrus/CouvreurV2")
  write(hydrus_dir, "LEVEL_01.DIR")
  
  # interpolate ET for the smallest dt 
  ET <- data.frame(approx(0:59,ET, seq(0.05,60,0.05)))[,2]
  
  E_soil <- tibble(rho_a = data.frame(approx(0:59,E_soil[,1], seq(0.05,60,0.05)))[,2],
                   r_a = data.frame(approx(0:59,E_soil[,2], seq(0.05,60,0.05)))[,2],
                   LAI = data.frame(approx(0:59,E_soil[,3], seq(0.05,60,0.05)))[,2],
                   r_s1 = data.frame(approx(0:59,E_soil[,4], seq(0.05,60,0.05)))[,2],
                   rho_a_night = data.frame(approx(0:59,E_soil[,5], seq(0.05,60,0.05)))[,2]

  )
  
  
  soil_global = new_all_roots <- NULL
  # Initial soil condition
  soil <- soil_initial(soil_type = stype, field_capacity = 0.5)

  temp_conduct <- conductivities
  results <- tibble(age = seq(0.05,60,0.05),
                    DN = ifelse(age-floor(age)>= 0.6,"night", "day"), # ifelse(age%/%round(age) == 0, "night", "day")
                    krs = NA,
                    ksrs = NA, 
                    tact = NA,
                    tact_eq = NA, 
                    tpot = NA,
                    tpot_eq = NA,
                    Kcomp = NA,
                    r_soil = NA,
                    psi_sr = NA,
                    E = NA,
                    psi_collar = NA,
                    rep = sam)
  
 # root system is too small, (length(Hsr) == 1) --> T
  psi_sr = soil$psi[1] # at t = 0 # hPa
  T_PM = Tact <- 0.01 # ET[1] # cm3 d-1
  Ksrs = 6E-6 # cm3 d-1 hPa-1
  temp_psi_collar <- psi_sr - Tact / Ksrs # hPa
  psi_night <- NA
  today <- 0
  while(today <= age_max){
    today = today + dt
    message("--------------------------------------------")

    reso_t <- today/0.05 # indice to look for values in table

    if(results$DN[reso_t] == "night"){
      print("night")
        # Night Transpiration of maize is +- 18% of the one during the day # Tamang and Sadok 2018
      tmp_ksrs <- results$ksrs[!is.na(results$tact_eq)]
      ksrs_fut = 2*tmp_ksrs[length(tmp_ksrs)] - tmp_ksrs[length(tmp_ksrs)-1]
      
      # mean transpiration of the day before
      transpi_day <- results$tact_eq[results$DN == "day" & !is.na(results$tact_eq) & floor(results$age) == floor(today)]
      T_day = mean(transpi_day)
      if(is.na(psi_night)){
        temp_psi_collar = psi_sr - ((T_day*0.18)/ksrs_fut)
      }else{
        temp_psi_collar = psi_night
      }
      
    }else{
      print("day")
      temp_psi_collar = -15000 # shoot demand at max
    }
    
    if(temp_psi_collar < -15000){temp_psi_collar = -15000} # threshold value for Psi collar.

    res_day <- try(HM_loop(soil, all_roots, temp_conduct, temp_psi_collar, today, weather, sam, ET = ET/10000/0.75/0.15, 
                   E_soil, results, stype, dt), silent = T)
  
    hydraulics <- res_day[[3]]
    Tact <- hydraulics$tact_eq
    krs <- hydraulics$krs
    ksrs <- hydraulics$ksrs
    
    # control over dt
    if(today > 0.2){
      if(ksrs > krs*0.95 | results$DN[reso_t] == "night"){
        dt = 0.2
      }else{
        if(ksrs > krs*0.90){
          today = today - dt
          dt = 0.1
          today = today + dt
          res_day <- try(HM_loop(soil, all_roots, temp_conduct, temp_psi_collar, today, weather, sam, ET = ET/10000/0.75/0.15, 
                     E_soil, results, stype, dt), silent = T)
          hydraulics <- res_day[[3]]
          krs <- hydraulics$krs
          ksrs <- hydraulics$ksrs
        }else{
          today = today - dt
          dt = 0.05
          today = today + dt
          res_day <- try(HM_loop(soil, all_roots, temp_conduct, temp_psi_collar, today, weather, sam, ET = ET/10000/0.75/0.15, 
                     E_soil, results, stype, dt), silent = T)
          hydraulics <- res_day[[3]]
        }
      }
    }
    reso_t <- today/0.05 # indice to look for values in table


    if(results$DN[reso_t] == "night"){ # night transpiration
      while(Tact > T_day*0.2 | Tact < T_day*0.15){
        message(paste0("low = ",T_day*0.15, " , high = " ,T_day*0.2, " and now Tact = ",Tact))
        if(Tact < T_day*0.15){
          dis_psi = (Tact - T_day*0.175)
          if(abs(dis_psi) < 2 ){dis_psi = -5}
          temp_psi_collar = temp_psi_collar+10*dis_psi
        }
        if(Tact > T_day*0.2){
          dis_psi = (Tact - T_day*0.175)
          if(abs(dis_psi) < 2 ){dis_psi = 5}
          temp_psi_collar = temp_psi_collar+20*dis_psi
        }
        if(temp_psi_collar < -15000){temp_psi_collar = -15000} # threshold value for Psi collar.
        res_day <- try(HM_loop(soil, all_roots, temp_conduct, temp_psi_collar, today, weather, sam, ET = ET/10000/0.75/0.15, 
                     E_soil, results, stype, dt), silent = T)
        hydraulics <- res_day[[3]]
        Tact <- hydraulics$tact_eq
        psi_night <- temp_psi_collar
      }
    }

    if(results$DN[reso_t] == "day" & Tact > ET[reso_t]*10000/8){ # max transpi day
    # The shoot water demand is limiting
    message("the shoot water demand is limiting, setting threshold")
      while(Tact > ET[reso_t]*10000/8| Tact < (ET[reso_t]*10000/8)*0.95){
        if(Tact < (ET[reso_t]*10000/8)*0.95){
          temp_psi_collar = temp_psi_collar-513
        }
        if(Tact > (ET[reso_t]*10000/8)){
          temp_psi_collar = temp_psi_collar+150
        }
        if(temp_psi_collar < -15000){temp_psi_collar = -15000} # threshold value for Psi collar.
         # soil, all_roots, temp_conduct, tpots, today, weather, sam, ET, E_soil, results, stype, dt
        res_day <- try(HM_loop(soil, all_roots, temp_conduct, temp_psi_collar, today, weather, sam, ET = ET, 
                     E_soil, results, soil_param, dt), silent = T)
        hydraulics <- res_day[[3]]
        Tact <- hydraulics$tact_eq
      }
    }

    results$krs[reso_t] <- hydraulics$krs
    results$tact[reso_t] <- hydraulics$tact
    results$tpot[reso_t] <- hydraulics$tpot
    results$ksrs[reso_t] <- hydraulics$ksrs
    results$tact_eq[reso_t] <- hydraulics$tact_eq
    results$tpot_eq[reso_t] <- hydraulics$tpot_eq
    results$Kcomp[reso_t] <- hydraulics$Kcomp
    #results$soilW[reso_t] <- hydraulics$soilW
    results$r_soil[reso_t] <- hydraulics$r_soil
    results$E[reso_t] <- hydraulics$E
    psi_sr <- res_day[[4]]
    results$psi_sr[reso_t] <- psi_sr
    results$psi_collar[reso_t] <- temp_psi_collar

    soil <- res_day[[1]]
    temp_root <- res_day[[2]]

    new_all_roots <- rbind(new_all_roots, temp_root)
    soil_global <- rbind(soil_global, soil)
  }

  dens <- new_all_roots%>%
    dplyr::group_by(type, age, rz1 = round(((z1/2+z2/2)/2)*2))%>%
    dplyr::summarise(root = sum(length*radius^2*pi),
                     radius = mean(radius), 
                     su = sum(suf),
                     sud = sum(suf/length),
                     su1 = sum(suf1),
                     j = sum(jr),
                     jx = sum(jxl),
                     p = sum(psi),
                     psi = mean(psi),
                     kr = mean(kr),
                     kr_eq = mean(kr_eq),
                     l = sum(length))%>%
    ungroup()%>%
    mutate(sam = sam)
  
  H_M <- list(soil_global, results, dens)
  message("End of HM_RWUM")
  return(H_M)
}

HM_loop <- function(soil, all_roots, temp_conduct, tpots, today, weather, sam, ET, E_soil, results, stype, dt){
  reso_t <- today/0.05
  if(length(all_roots$K_type[1]) > 0){ # If root type has been modified as a function of the radius
    print("yes")
    radtype = T
  }else{
    print("no")
    radtype = F
  }
  # --------- MARSHAL ----------------
  # Computing the root system hydraulic architecture
  print(today)
  
  # Select specific Soil for the simulation
  temp_soil <- soil[soil$z %in% seq(-200, 0 ,2),c(1,2,4)]
  theta_top <- soil[1,"moisture"]
  if(length(theta_top) == 0){
    theta_top = 0.0866 # initial condition (-300 hPa)
  }
  h_top <- temp_soil[1,"psi"] 
  # Select specific root system for the simulation
  temp_root <- all_roots%>%
    filter(time <= today)

  soil_param <- read_excel("/home/users/a/h/aheymans/GSUA/www/Soil_type.xlsx")%>%
      filter(type == stype)%>%
      mutate(Ksat = Ks, lambda = l)
  
  setDT(temp_root)
  # Re-arrange the input data
  orders <- unique(temp_conduct$order)
  ids <- unique(temp_conduct$order_id)
  #temp_root$name <- "root"
  for(o in c(1:length(orders))){
    if(radtype){
      temp_root$name[temp_root$K_type == ids[o]] <- orders[o]
    }else{
      temp_root$name[temp_root$type == ids[o]] <- orders[o]
    }
  }
  
  temp_conduct$y[temp_conduct$type == "kx" & temp_conduct$order_id == 0] <-  max(temp_conduct$y[temp_conduct$type == "kx"])
  
  # -----------------------------
  # Run MARSHAL
  # -----------------------------
  hydraulics <- getSUF(temp_root, 
                       temp_conduct, 
                       temp_soil, 
                       hetero = T, 
                       Psi_collar = tpots, # Psi_collar_tmp
                       soil_param)
  
  # Aggregate output from MARSHAL
  results$krs[results$age == today] <- hydraulics$krs
  results$tact[results$age == today] <- hydraulics$tact
  results$tpot[results$age == today] <- hydraulics$tpot
  
  results$ksrs[results$age == today] <- hydraulics$ksrs
  results$tact_eq[results$age == today] <- hydraulics$tact_eq
  results$tpot_eq[results$age == today] <- hydraulics$tpot_eq
  
  # Merge output of MARSHAL on specific root segment
  temp_root$suf <- as.vector(hydraulics$suf)
  temp_root$suf_eq <- as.vector(hydraulics$suf_eq)
  temp_root$suf1 <- as.vector(hydraulics$suf1)
  temp_root$kx <- as.vector(hydraulics$kx)
  temp_root$kr <- as.vector(hydraulics$kr)
  temp_root$kr_eq <- as.vector(hydraulics$kr_eq)
  temp_root$jr <- as.vector(hydraulics$jr)
  temp_root$jr_eq <- as.vector(hydraulics$jr_eq)
  temp_root$psi <- as.vector(hydraulics$psi)
  temp_root$jxl <- as.vector(hydraulics$jxl)
  temp_root$jxl_eq <- as.vector(hydraulics$jxl_eq)
  temp_root$psi_soil <- as.vector(hydraulics$psi_soil)
  # Add the simulation specificity to the dataset
  temp_root$rep <- sam
  temp_root$age <- today
  
  print("marshal did run")
  # Mean properties at every 2 cm deep
  RLDWU <- temp_root%>%
    mutate(rz2 = round((z1/2+z2/2)/2)*2)%>%
    dplyr::group_by(rz2)%>%
    dplyr::summarise(su = sum(suf1),
                     ps = sum(psi),
                     jr = sum(jr),
                     jx = sum(jxl),
                     su_eq = sum(suf_eq),
                     jr_eq = sum(jr_eq),
                     jx_eq = sum(jxl_eq))%>%
    ungroup()%>%
    mutate(age = today)
  
  # t0 <- rbind(t0, marshal = proc.time())
  message(paste0("the root system allows ",hydraulics$tact_eq[1], " ml of water to be absorbed" ))
  # -------- HYDRUS 1D ------------- 
  max_transpi_root <- hydraulics$tact_eq[1]/75/15 # 8 plant per sqare meter, then convert from mL*m-2*day-1 to cm*m-2*day-1
  r_root <- min(c(max_transpi_root, ET[reso_t])) # transpiration due to the plants: T_aq = Kcb*ET0 or max_transpi_root 
  
  T_Cel <- 20 # assumption, soil top layer is at 20?C
  rho_vs <- 10^-3*exp(31.3716-6014.79/(T_Cel+273)-7.92*10^-3*(T_Cel+273))
  E <- E_soil[reso_t, ]%>%
    mutate(r_s2 = 10*exp(35.63*(0.15-theta_top)),
           Hr_top = exp(h_top*0.018015*9.81/(1.01*0.287*(T_Cel+273))),
           rho_v = rho_vs*Hr_top,
           E = (rho_v-rho_a)/(r_a - (r_s2)) )%>%select(E)
  
  r_soil <- ET[reso_t]-r_root # evaporation of the soil, E = Ke*ET0 or ET - Kcb*ET0 = Ke*ET0
  if(E$E[1] < 0){
    E$E[1] = 0
  }
  if(E$E[1] > 0.5){ # avoid large evaporation rates
    E$E[1] = 0.5
  }
  
  message(paste0("the soil evaporate ",E$E[1], " cm d-1" ))
  
  atm_bc_data <- data.frame(tAtm = 1, Prec = weather[ceiling(today)], rSoil = E$E[1], 
                            rRoot = r_root, hCritA = 15000, rB = 0, hB = 0, ht = 0, 
                            RootDepth = 0)
  
  hydraulics$r_soil <- r_soil
  hydraulics$E <- E$E[1]
  # overwrite the atmposheric boundary condition of hydrus.
  write.atmosph.in("Hydrus/CouvreurV2/",
                   maxAL = 1,
                   deltaT = 1,
                   atm_bc_data,
                   hCritS = 15000,
                   input.pet = F)
  
  Beta <- rep(0, 101)
  Beta[which(temp_soil$z %in% RLDWU$rz2 )] <- rev(RLDWU$su_eq) # from above to below
  SSF <- data.frame(suf = Beta, h = temp_soil$psi)
  # overwirte the profile boundary condition
  write.profile.dat(project.path = "./Hydrus/CouvreurV2", SSF)
  ksrs <- hydraulics$ksrs[1]
  
  Q_dou <- rev(RLDWU$jr_eq)/75/15 # 8 plant per sqare meter, then convert from mL*m-2*day-1 to cm*m-2*day-1
  SUF <- rev(RLDWU$su_eq[RLDWU$rz2 >= -200 & RLDWU$rz2 <= 0])
  Tact <- min(c(max_transpi_root, ET[reso_t]))
  message(paste0("a square meter of this plant transpires ", Tact," per day"))
  message(paste0("however, it rain: ", weather[ceiling(today)], " today"))
  Hsr <- temp_soil$psi[which(temp_soil$z %in% RLDWU$rz2 )]
  kcomp <- ksrs
  if(length(unique(Hsr))> 1){
    Hseq <- t(Hsr) %*% t(t(SUF))
    kcomp <- t(Hsr-Hseq) %*% (Q_dou/SUF - Tact) / (t(Hsr-Hseq) %*% t(t(Hsr-Hseq)))
  }else{
    Hseq <- unique(Hsr)
  }
  write.options.in(project.path = "./Hydrus/CouvreurV2", ksrs/75/15, kcomp/75/15)
  
  if(file.exists("./Hydrus/CouvreurV2/Nod_Inf.out")){
    message("The file Nod_Inf.out will be removed before excecuting hydrus1D")
    file.remove("./Hydrus/CouvreurV2/Nod_Inf.out")
  }
  # Run Hydrus1D ---
  system("./h1d_calc")
  message("hydrus has ended")
  #Load hydrus profile information
  if(file.exists("./Hydrus/CouvreurV2/Nod_Inf.out")){
    message("The file Nod_Inf.out is ready to be loaded ...")
  }
  hydrus <- read.nod_inf(project.path = "./Hydrus/CouvreurV2", 
                         out.file = paste0("Nod_Inf.out"))
  
  hydrus$water_g <- hydrus$Moisture*2*10000
  water_kg <- sum(hydrus$water_g)/1000
  
  message("Hydrus output is loaded into the pipeline")
  soil <- hydrus%>%
    filter(Time == dt)%>%
    transmute(id = Node,
              z = Depth,
              value = today,
              psi = Head,
              moisture = Moisture,
              SSF = SSF$suf)
  
  
  # t0 <- rbind(t0, hydrus = proc.time())
  hydraulics$Kcomp <- kcomp
  hmlo <- list(soil,temp_root, hydraulics, Hseq)
  return(hmlo)
}

Soil_type <- function(type = "loam"){
  
  
  soilparam <- read_excel("/home/users/a/h/aheymans/GSUA/www/Soil_type.xlsx")
  if(type %!in% soilparam$type){
    warning(paste0("soil type is not defined: tried ",soilparam$type))
  }
  soilparam <- soilparam[soilparam$type == type,]
  soilparam <- soilparam[1,] # only one row if repetition
  
  write.selector.in("./Hydrus/CouvreurV2", mesh = NA, soilparam = soilparam )
}

read.nod_inf <- function (project.path, out.file = "Nod_Inf.out", output = NULL, 
                          warn = FALSE, ...) 
{
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
                         z = sort(seq(-200,0,2), decreasing = T),
                         value = 1,
                         psi = rep(-73.5021, 101))
  }else if (soil_type == "loam") {
      soil <- data.frame(id=1:101,
                         z = sort(seq(-200,0,2), decreasing = T),
                         value = 1,
                         psi = rep(-132.869, 101))

  }else if (soil_type == "silty clay") {
      soil <- data.frame(id=1:101,
                         z = sort(seq(-200,0,2), decreasing = T),
                         value = 1,
                         psi = rep(-354.424, 101))
  }else{
    warning("This soil type is not referenced yet")
  }

  soil$psi <- soil$psi*(1/field_capacity)

  return(soil)
  
}

write.options.in <- function(project.path = "./Hydrus_1D_Couvreur/CouvreurV2", krs, kcomp){
  
  path <- paste0(project.path, "/Options.in")
  # path = "./Hydrus_1D_Couvreur/CouvreurV2/Options.IN"
  op = readLines(con = paste0(project.path,"/Options_1.IN"), n = -1L, encoding = "unknown")
  tmp_krs <- op[10]
  op[10]<- paste(replace(str_split(tmp_krs, " ")[[1]], which(str_split(tmp_krs, " ")[[1]]== "0.0001"), krs), collapse = " ")
  tmp_kcomp <- op[13]
  op[13] <- paste(replace(str_split(tmp_kcomp, " ")[[1]], which(str_split(tmp_kcomp, " ")[[1]]== "0.0001"), kcomp), collapse = " ")
  
  write(op, path)
}


write.profile.dat <- function(project.path = "./Hydrus/CouvreurV2", SSF){
  
  out.file <- "PROFILE.DAT"
  prof = readLines(con = paste0(project.path,"/PROFILE_1.DAT"), n = -1L, encoding = "unknown")
  h_ind = grep("h", prof)
  n <- as.numeric(str_split(prof[h_ind], " ")[[1]][3])
  if(nrow(SSF) != n){error("sum up SSf every two centimeter along the profile")}
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
  data_fmt <- ile
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
  tAtm_ind = grep(" tAtm", atm_data)
  
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

create_soil <- function(weather){
  time_step = 1
  endTime = 60
  ntimesteps = 60
  profile_depth = 200
  
  RLD <- all_roots%>%
    mutate(age_RS = ceiling(time))%>%
    dplyr::group_by(age_RS)%>%
    dplyr::summarise(rooting_depth = min(z2))%>%
    ungroup()
  
  soil = atmosph <- NULL
  fit_rdepth <- aov(abs(RLD$rooting_depth)~seq(time_step, endTime, by = 1))
  
  atm_bc_data <- data.frame(tAtm = seq(time_step, endTime, by = time_step),
                            # Precipitation (cm/d)
                            Prec = weather,
                            # Soil evaporation (cm/d)
                            rSoil = rep(0, ntimesteps),
                            # Transpiration (cm/d)
                            rRoot = seq(2/ntimesteps, 2, by = 2/ntimesteps),
                            hCritA = rep(15000, ntimesteps),
                            rB = numeric(ntimesteps),
                            hB = numeric(ntimesteps),
                            ht = numeric(ntimesteps),
                            RootDepth = fit_rdepth$coefficients[1]+fit_rdepth$coefficients[2]*c(1:60))
  
  # overwrite the atmposheric boundary condition of hydrus.
  write.atmosph.in("Hydrus/1DRAINAG/",
                   maxAL = ntimesteps,
                   deltaT = time_step,
                   atm_bc_data,
                   hCritS = 15000,
                   input.pet = F)
  # setwd("C:/Program Files (x86)/PC-Progress/Hydrus-1D 4.xx")
  # Run Hydrus1D ---
  
  system("./h1d_calc")
  #Load hydrus profile information
  hydrus <- read.nod_inf(project.path = "Hydrus/1DRAINAG/", 
                         out.file = paste0("Nod_Inf.out"), output = NULL,
                         warning = FALSE)
  #Hydrus dataset is large, then to optimize computation time, a subset is taken.
  soil_tmp <- hydrus%>%
    filter(Time %in% seq(time_step,endTime,time_step), #Time resolution of the overall exemple
           Depth %in% seq(-profile_depth, 0, by = 0.2))
  
  #formating the soil data to be compatible with MARSHAL
  soil <- soil_tmp%>%
    transmute(id = c(1:nrow(soil_tmp)),
              z = Depth,
              value = Time,
              psi = Head)
  
  if(is.na(soil$psi[soil$value == 60][1])){
    message("hydrus crash once")
    write.selector.in("Hydrus/1DRAINAG/", mesh = "2")
    # Run Hydrus1D ---
    system("./h1d_calc")
    #Load hydrus profile information
    hydrus <- read.nod_inf(project.path = "Hydrus/1DRAINAG/", 
                           out.file = paste0("Nod_Inf.out"), output = NULL,
                           warning = FALSE)
    #Hydrus dataset is large, then to optimize computation time, a subset is taken.
    soil_tmp <- hydrus%>%
      filter(Time %in% seq(time_step,endTime,time_step), #Time resolution of the overall exemple
             Depth %in% seq(-profile_depth, 0, by = 0.2))
    
    #formating the soil data to be compatible with MARSHAL
    soil <- soil_tmp%>%
      transmute(id = c(1:nrow(soil_tmp)),
                z = Depth,
                value = Time,
                psi = Head)
    write.selector.in("Hydrus/1DRAINAG/", mesh = "2")
    if(is.na(soil$psi[soil$value == 60][1])){
      message("hydrus crash twice ...")
    }
  }
  return(list(soil, atm_bc_data))
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
  
  x <- read_xml("CPlantBox-master/modelparameter/rootsystem/B73.xml")
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
  for(t in seq(0.2,30,0.1)){
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
