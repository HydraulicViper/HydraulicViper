# Dir -------------
DIR_path <- getwd()
setwd(dir = DIR_path)
T_start = proc.time()

JOB_NUMBER <- 1
# library and io_function
source("/home/users/a/h/aheymans/GSUA/R/MAGIC_function.R")
source("/home/users/a/h/aheymans/GSUA/R/getSUF.R")
source("/home/users/a/h/aheymans/GSUA/R/dependencies.R")
source("/home/users/a/h/aheymans/GSUA/R/GRANAR/granar.R")

# Parameter--------------
params <- read_param_xml("/home/users/a/h/aheymans/GSUA/www/Zea_mays_CT.xml")
MatSobol <- read_excel(paste0("/home/users/a/h/aheymans/GSUA/sampling/MAGIC_Input_01.xlsx"))
Sobol <- read.csv(paste0("/home/users/a/h/aheymans/GSUA/sampling/SamplingSobol_",JOB_NUMBER,".csv"))

Sobol$ID <- rep(MatSobol$ID,5)

MatSobol <- merge(Sobol, MatSobol[,1:27], by = "ID")

WEA <- read.csv("/home/users/a/h/aheymans/GSUA/www/weather_4loc.csv")
WEA <- WEA%>%filter(local == "Indianapolis,United States")
weather <- WEA$PRECTOT/10 # mm to cm of rain
weather[weather < 0.05] <- 0
ET <- WEA$ET/10 # mm to cm of transpiration allow by the shoot part
E_soil <- WEA %>%
  transmute(rho_a = PS/(1.01*(T2M+273)*0.287),
            r_a = 208/WS2M,
            LAI = 5.76 / ( 1 + exp(-0.013 * (GDD - 502.74))),# Campos et al., 2017 LAI as GDD
            r_s1 = 100/(0.5*LAI)) 


# -1.06*T*log10(100/HR)*10000
WEA <- WEA%>%mutate(tpots = ifelse(Wpot >= -15000, Wpot, -15000 + (Wpot + 15000)/2))
tpots <- WEA$tpots 
# tpots[tpots < -15000] <- -15000 # hPa # Threshold value under -15000 hPa, the stomata close.

# Scrpit R ----------------
pol_ch = poly_area = all_rld <- NULL

ANAT = TIME = soil_water = T_res = HYDRO <- NULL
parti<- NULL
sam <- 3

# ----------------------------------------- #
#############################################
# ----------------------------------------- #

for(sam in 1:nrow(MatSobol)){ 
  message(paste0("--------------------------job number :",
                 JOB_NUMBER, "-------  input line #: ", sam,
                 "--------------------------------------------"))
  # Run CRootBox ------------------------------
  all_roots <- NULL
  all_roots <- CRootBox(sam, MatSobol) # With B73 param  
  
  # Get kr and kx for all roots
  # # granar mecha short-cut with the random forest model
  all_roots <- K_type(all_roots) # gather roots in cluster by their root diameter
  conductivities <- try(GM(MatSobol, sam, shcut = T, all_roots, model),silent = T)
  
  
  if(class(conductivities) == "try-error"){
    # Try again
    # Run CRootBox ------------------------------
    all_roots <- NULL
    all_roots <- CRootBox(sam, MatSobol) # With B73 param
    
    # Get kr and kx for all roots
    # # granar mecha short-cut with the random forest model
    all_roots <- K_type(all_roots) # gather roots in cluster by their root diameter
    conductivities <- try(GM(MatSobol, sam, shcut = T, all_roots, model),silent = T)
    if(class(conductivities) == "try-error"){next()}
  }
  
  if(length(which(is.na(conductivities$y)))>0){
    message("error GM function. It generate NAN values")
    # Try again
    # Run CRootBox ------------------------------
    all_roots <- NULL
    all_roots <- CRootBox(sam, MatSobol) # With B73 param
    
    # Get kr and kx for all roots
    # # granar mecha short-cut with the random forest model
    all_roots <- K_type(all_roots) # gather roots in cluster by their root diameter
    conductivities <- try(GM(MatSobol, sam, shcut = T, all_roots, model),silent = T)
    if(class(conductivities) == "try-error"){next()}
  }
  if(length(which(is.na(conductivities$y)))>0){
    message("error GM function. It generate NAN values")
    next()
  }
  # Compute Convex hull info
  root_hull <- CRoot_hull(all_roots, age_filter = seq(10,30,10), sam)
  RLD <- root_hull[[1]]
  pol_ch <- rbind(pol_ch,root_hull[[2]])
  poly_area <- rbind(poly_area, root_hull[[3]])
  RLD$sam <- sam
  all_rld <- rbind(all_rld, RLD)
  
  HYDRO <- rbind(HYDRO, conductivities%>%mutate(sam = sam))
  
  # Soil --------------
  message("Resolving the hydraulic architecture of the RS and the water flow equation in a 1D soil profile")
  # MARSHAL + Hydrus -----------------
  hydrus_marshal <- try(HM_RWUM(weather, all_roots,sam, tpots, conductivities, ET, E_soil),silent = T)
  if(class(hydrus_marshal) == "try-error"){next()}
  # time <- rbind(time, hydrus_marshal[[4]], loop_end = proc.time())
  
  soil <- hydrus_marshal[[1]]
  soil$sam <- sam
  for(d in unique(soil$value)){
    soil$prec[soil$value == d] <- weather[d/1]
  }
  
  tmp_results <- hydrus_marshal[[2]]
  
  dens <- hydrus_marshal[[3]]
  
  T_res <- rbind(T_res, tmp_results)
  message(paste0( length(unique(T_res$rep)), " / ", sam , " did finish successfully"))
  
  soil_water <- rbind(soil_water, soil)
  
  parti <- rbind(parti, dens)
  
  # MECHA outputs
  write.csv(HYDRO, paste0("./res",JOB_NUMBER,"/conductivities_",JOB_NUMBER,".csv"))
  # Save results of Marshal
  write.csv(T_res, paste0("./res",JOB_NUMBER,"/MARSHAL_",JOB_NUMBER,".csv"))
  # root water uptake partitioning
  write.csv(parti, paste0("./res",JOB_NUMBER,"/RWU_partitioning_",JOB_NUMBER,".csv"))
  # Save soil
  write.csv(soil_water, paste0("./res",JOB_NUMBER,"/soil_",JOB_NUMBER,".csv"))
  
}

T_final = proc.time()
print(T_final - T_start)