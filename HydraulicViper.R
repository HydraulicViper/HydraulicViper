
# Dir ------------
DIR_path <- getwd()
setwd(dir = DIR_path)

JOB_NUMBER<-1 # num
# library and io_function
source("./R/io_function.R")
source("./R/getSUF.R")
source("./R/dependencies.R")

MatSobol<- read.csv(paste0("./www/SamplingB73_Sobol.csv"))
MatSobol$km <- 3E-5 # hydraulic conductivity of the membrane without aquaporins. fixed because it has no effect
WEATHER <- read.csv("./www/Weather_2loc30days.csv")

location <- as.character(unique(WEATHER$local))
w <- location[2]
WEATHER <- WEATHER %>% filter(local == w)
soil_list <- c("sandy loam", "silty clay", "loam")

WEA <- WEATHER%>%filter(local == w)
weather <- WEA$PRECTOT/10 # mm to cm of rain 
weather[weather < 0.05] <- 0 
Envi <- WEA %>%
    transmute(rho_a = PS/(1.01*(T2M+273)*0.287),
              r_a = 208/WS2M,
                 LAI = 5.76 / ( 1 + exp(-0.013 * (GDD - 502.74))),# Campos et al., 2017 LAI as GDD
                 r_s1 = 100/(0.5*LAI),
                 rho_a_night = PS/(1.01*(T2M_MIN+273)*0.287),
                 ET = ET/10,
                 #Transpi = (ET*LAI/8)^(0.6),
                 Kt = 1-exp(-0.52*LAI),
                 Kcb = Kt/0.9,
                 Transpi = ET0*Kcb/10)
age_max = 50

sam <- JOB_NUMBER
hydrus_marshal = tibble()
class(hydrus_marshal) = "try-error"

for(sam in c(((JOB_NUMBER-1)*10+1):((JOB_NUMBER)*10))){ # c(((JOB_NUMBER-1)*10+1):((JOB_NUMBER)*10))-40*floor((JOB_NUMBER-1)*1/4
     # allocate table name
   HYDRO = CC = OUT = SOIL = NULL
   id_sam = MatSobol$X[sam]

   for(stype in soil_list){
        serr_rel = 10
        # selection for adapted phenotype     
        
        while(serr_rel > 5){
            if(serr_rel > 5){    
              # Run CRootBox ------------------------------
              all_roots <- try(CRootBox(sam, MatSobol, age_max), silent = T)
              all_roots <- K_type(all_roots) # gather roots in cluster by their root diameter
              conductivities <- try(GM(MatSobol, sam, shcut = T, all_roots),silent = T)
              carbon <- Ccost(MatSobol, sam, all_roots)
            }
        
            hydrus_marshal = try(HM_RWUM(weather, all_roots, sam, conductivities, Envi, stype), silent = T)
            if(class(hydrus_marshal) == "try-error"){next()}
            tmp = hydrus_marshal[[2]]
            tmp$sam = id_sam
            tmp$soil_type = stype
            err = tmp$vRoot[tmp$Time>3] - tmp$Tact_eq_MA[tmp$Time>3]
            serr_rel = abs(sum(err)/max(tmp$sum_vRoot)*100)
            
            print("+-+-+-+-+-+-")
            print(serr_rel)
            print("+-+-+-+-+-+-")

        }
        tmp_soil = hydrus_marshal[[1]]
        tmp_soil$sam = id_sam
        tmp_soil$soil_type = stype
        SOIL = rbind(SOIL, tmp_soil)
        OUT = rbind(OUT, tmp)
        CC <- rbind(CC, carbon)
        
    }
    write.csv(CC, paste0("./res",JOB_NUMBER,"/Carbon_",JOB_NUMBER,"_",id_sam,".csv"))
    write.csv(OUT, paste0("./res",JOB_NUMBER,"/out_",JOB_NUMBER,"_",id_sam,".csv"))
    write.csv(SOIL, paste0("./res",JOB_NUMBER,"/SOIL_",JOB_NUMBER,"_",id_sam,".csv"))
}


