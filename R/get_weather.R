# Adrien Heymans
# december 2020


require(nasapower)
get_weather <- function(lon = 4.63, lat = 50.66, year = 2018, month = 5){
  
  if((month+3)>9){
    dates = c(paste0(year,"-0",month,"-15"), paste0(year,"-",month+3,"-04"))
  }else{
    dates = c(paste0(year,"-0",month,"-15"), paste0(year,"-0",month+3,"-04"))
  }
  daily_ag <- get_power(community = "AG",
                        lonlat = c(lon, lat),
                        pars = c("RH2M", "T2M", "PRECTOTCORR", "T2M_MAX", "T2M_MIN",
                                 "ALLSKY_SFC_SW_DWN", "CLRSKY_SFC_SW_DWN", "PS", "WS2M", "T2MDEW"),
                        dates = dates,
                        temporal_api = "DAILY")
  
  daily_ag <- daily_ag[1:60,]
  daily_ag <- daily_ag%>%
    mutate(Wpot = -1.06*(T2M+273.3)*log10(100/RH2M), # Water potential MPa
           e_t = 0.6108*exp((17.27*T2M)/(T2M+237.3)), # 
           e_tmin = 0.6108*exp((17.27*T2M_MIN)/(T2M_MIN+237.3)),
           e_tmax = 0.6108*exp((17.27*T2M_MAX)/(T2M_MAX+237.3)),
           e_a = 0.6108*exp((17.27*T2MDEW)/(T2MDEW+237.3)), # actual vapor pressure kPa
           e_s = (e_tmax+e_tmin)/2, # saturation vapor pressure kPa
           VPD = e_s-e_a,
           delta = 4098 * e_t / (T2M +237.3)^2,
           gama = (0.665*10^-3)*PS,
           gd = T2M-10)
  daily_ag$GDD <- 0
  daily_ag$GDD = cumsum(daily_ag$gd)

  daily_ag$ALLSKY_SFC_SW_DWN[daily_ag$ALLSKY_SFC_SW_DWN < 1] <- mean(daily_ag$ALLSKY_SFC_SW_DWN[daily_ag$ALLSKY_SFC_SW_DWN > 1])
  
  daily_ag <- daily_ag%>%
    mutate(ratio = ALLSKY_SFC_SW_DWN/CLRSKY_SFC_SW_DWN)
  daily_ag$ratio[daily_ag$ratio < 1] <- 1
  daily_ag$PRECTOT = daily_ag$PRECTOTCORR
  daily_ag$PRECTOT[daily_ag$PRECTOTCORR < 0.2] <- 0
  
  daily_ag <- daily_ag%>%
    mutate(Rns = (1-0.23)*ALLSKY_SFC_SW_DWN,
           Rnl = (4.903*10^-9)*((T2M_MAX+273.16*2+T2M_MIN)/2)*(0.34-0.14*sqrt(e_a))*(1.35*(ratio)-0.35),
           Rn = Rns-Rnl,
           ET0 = (0.408*delta*Rn+gama*(900/(T2M+273))*WS2M*(e_s-e_a))/(delta+gama*(1+0.34*WS2M)),
           Kc = 1.1*(1-exp(-(GDD/620)^3))+0.1, # based on Campos et al 2017 Kc point and FAO chap 6 ET0 - Kc_mid
           ET = Kc*ET0)
  return(daily_ag)
  
}
