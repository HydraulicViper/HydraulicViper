# Dir -------------
DIR_path <- getwd()
setwd(dir = DIR_path)

JOB_NUMBER <- 1
# library and io_function
source(".R/B73_function.R")
source("./R/getSUF.R")
source("./R/dependencies.R")
# source("/home/users/a/h/aheymans/GSUA/R/GRANAR/granar.R")

# Parameter--------------
# params <- read_param_xml("/home/users/a/h/aheymans/GSUA/www/Zea_mays_CT.xml")
MatSobol<- read.csv(paste0("./sampling/SamplingB73_Sobol.csv"))
# missing <- read.csv(paste0("/home/users/a/h/aheymans/GSUA/www/missing.csv"))
MatSobol$km <- 3E-5 # hydraulic conductivity of the membrane without aquaporins. fixed because it has no effect
WEATHER <- read.csv("./WeatherData/Weather_2loc30days.csv")
location <- as.character(unique(WEATHER$local))

MatSobol <- MatSobol%>%filter(X %in% 1:5) # 5 different phenotyping states

# tpots[tpots < -15000] <- -15000 # hPa # Threshold value under -15000 hPa, the stomata close.
# Scrpit R ----------------
pol_ch = poly_area = all_rld <- NULL

ANAT = TIME = soil_water = T_res = HYDRO <- NULL
parti<- NULL
w <- location[2]
WEATHER <- WEATHER %>% filter(local == w)
soil_list <- c("loam", "sandy loam", "silty clay")
# Soil type selection /!\ No hysteresis. 

# Scrpit R ----------------
pol_ch = poly_area = all_rld <- NULL

ANAT = TIME = soil_water = results = CC =  HYDRO <- NULL
parti<- NULL
sam <- 1
age_max = 30
 
for(sam in 1:nrow(MatSobol)){ # :nrow(MatSobol)
  message(paste0("------------------------ Starting simulation ",sam," from job ",JOB_NUMBER," -----------------------------------------------------"))

  # Run CRootBox ------------------------------
  all_roots <- CRootBox(sam, MatSobol, age_max) # With B73 param  
  
  # Get kr and kx for all roots
  # # granar mecha short-cut with the random forest model
  all_roots <- K_type(all_roots) # gather roots in cluster by their root diameter
  conductivities <- try(GM(MatSobol, sam, shcut = T, all_roots),silent = T)
  carbon <- Ccost(MatSobol, sam, all_roots)
  CC <- rbind(CC, carbon)

  n_try  <- 0
  while(length(which(is.na(conductivities$y)))>0){
    message("error GM function. It generate NAN values")
    if(n_try > 30){next()}
    all_roots <- NULL
    all_roots <- CRootBox(sam, MatSobol) # With B73 param
    all_roots <- K_type(all_roots) # gather roots in cluster by their root diameter
    conductivities <- try(GM(MatSobol, sam, shcut = T, all_roots), silent = T)
    n_try <- n_try + 1
  }
            # Compute Convex hull info
            
  C_cost = conductivities%>%mutate(K_type = parse_number(order), vol = vol * 0.01)%>%select(K_type, vol)
  C_cost$vol[C_cost$K_type == 0] <- C_cost$vol[C_cost$K_type == 1][1]
  RLD <- CRoot_hull(all_roots, age_filter = seq(5,30,1), sam, C_cost)
  RLD$sam <- MatSobol$X[sam]
  all_rld <- rbind(all_rld, RLD)
  
  HYDRO <- rbind(HYDRO, conductivities%>%mutate(sam =  MatSobol$X[sam]))


    WEA <- WEATHER%>%filter(local == w)
    weather <- WEA$PRECTOT/10 # mm to cm of rain
    weather[weather < 0.05] <- 0
    ET <- WEA$ET*1000*0.75*0.15 # mm to cm3 of transpiration allow by the shoot part
    E_soil <- WEA %>%
       transmute(rho_a = PS/(1.01*(T2M+273)*0.287),
                 r_a = 208/WS2M,
                 LAI = 5.76 / ( 1 + exp(-0.013 * (GDD - 502.74))),# Campos et al., 2017 LAI as GDD
                 r_s1 = 100/(0.5*LAI),
                 rho_a_night = PS/(1.01*(T2M_MIN+273)*0.287))
  # Soil --------------
    message("Resolving the hydraulic architecture of the RS and the water flow equation in a 1D soil profile")
  
      for(stype in soil_list){
        message(paste0("------------------------ Simulation ",sam," with the environment ",w," and soil type ", stype," -----------------------------------------------------"))
        # MARSHAL + Hydrus -----------------

      hydrus_marshal <- try(HM_RWUM(weather, all_roots, sam, conductivities, ET, E_soil, stype),silent = T)
      if(class(hydrus_marshal) == "try-error"){next()}
        # time <- rbind(time, hydrus_marshal[[4]], loop_end = proc.time())
  
      soil <- hydrus_marshal[[1]]
      soil$sam <-  MatSobol$X[sam]
      soil$loc = w
      soil$soil_type <- stype
      for(d in unique(soil$value)){
        soil$prec[soil$value == d] <- weather[d/1]
      }
      tmp_results <- hydrus_marshal[[2]]
            tmp_results <- tmp_results%>%filter(!is.na(krs))

      tmp_results$soil_type = stype
      tmp_results$loc = w
      tmp_results$sam = MatSobol$X[sam]
      dens <- hydrus_marshal[[3]]
      dens$soil_type = stype
      dens$loc = w
      dens$sam =  MatSobol$X[sam]
      results <- rbind(results, tmp_results)
      soil_water <- rbind(soil_water, soil)
      parti <- rbind(parti, dens)  
    }

  # MECHA outputs
    # Save results of Marshal
  write.csv(results, paste0("./res",JOB_NUMBER,"/MARSHAL_",JOB_NUMBER,".csv"))
    # root water uptake partitioning
  write.csv(parti, paste0("./res",JOB_NUMBER,"/RWU_partitioning_",JOB_NUMBER,".csv"))
  
  write.csv(CC, paste0("./res",JOB_NUMBER,"/Carbon_",JOB_NUMBER,".csv"))
  
  write.csv(all_rld, paste0("./res",JOB_NUMBER,"/rld_",JOB_NUMBER,".csv"))
    # Save soil
  write.csv(soil_water, paste0("./res",JOB_NUMBER,"/soil_",JOB_NUMBER,".csv"))

  write.csv(HYDRO, paste0("./res",JOB_NUMBER,"/conductivities_",JOB_NUMBER,".csv"))
}