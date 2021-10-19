
alometry <- function(data, RXA_TSA_coef, MXN_TSA_coef, XVA_TSA_coef, CF_CW_coef){

data <- data%>%
  mutate(log_RXA = log(pi*radius^2),
         RXA = exp(log_RXA),
         log_TSA = RXA_TSA_coef[1]+RXA_TSA_coef[2]*log_RXA,
         TSA = exp(log_TSA)+var_stele*exp(log_TSA),
         log_TSA = log(TSA, exp(1)),
         r_stele = sqrt(TSA/pi),
         log_nX = MXN_TSA_coef[1]+MXN_TSA_coef[2]*log_TSA,
         nX = exp(log_nX),
         mod_XVA = XVA_TSA_coef[1]+XVA_TSA_coef[2]*log_TSA,
         MXA = exp(-mod_XVA^2)+var_xylem*exp(-mod_XVA^2),
         XVA = MXA,
         X_size = 2*sqrt((MXA/nX)/pi),
         log_CW = log(radius-r_stele, exp(1)),
         CF = exp(CF_CW_coef[1]+CF_CW_coef[2]*log_CW),
         OneC = exp(log_CW)/CF,
         PXA_1 = (OneC/2.2)^2,
         ratio = (2+0.07456*r_stele*1000)/nX,
         nPX = round(nX*ratio),
         PXA = nPX*PXA_1,
         a = PXA_1*1000^2,
         k_protxyl_s = a^2/(8*pi*200*1E-5/3600/24)*1E-12,
         # kx when only the proto xylem have their cell wall lignified 
         kx_unM = k_protxyl_s*nPX*200/1E4,
         LMXA = MXA - PXA,
         LMXA_1 = LMXA/nX,
         b = LMXA_1*1000^2,
         k_Mxyl_s = b^2/(8*pi*200*1E-5/3600/24)*1E-12,
         # kx when all xylem elements have their cell wall lignifiedk
         kx_M = k_Mxyl_s*nX*200/1E4 + kx_unM,
         km = 2.4E-4, 
         kw = 3.0E-5, 
         kAQP = 4.3E-4, 
         kpl = 5.3E-12, 
         thickness = 1.5,
         TCA = RXA-TSA)

  return(data)

}
