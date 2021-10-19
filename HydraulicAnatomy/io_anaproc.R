
library(tidyverse)
library(plyr)
library(deldir)
library(alphahull)
library(xml2)
library(sp)
library(viridis)
library(readxl)
library(caret)
library(sensitivity)
source("./R/granar.R")
source("./R/micro_hydro.R")

makeMCSample <- function(n, vals, p = 1) {
  # Packages to generate quasi-random sequences
  # and rearrange the data
  require(randtoolbox)
  require(plyr)
  
  # Generate a Sobol' sequence
  sob <- sobol(n, length(vals), seed = round(runif(1, 1000,10000)), scrambling = 1)
  
  # Fill a matrix with the values
  # inverted from uniform values to
  # distributions of choice
  samp <- matrix(rep(0,n*(length(vals)+1)), nrow=n)
  samp[,1] <- 1:n
  for (i in 1:length(vals)) {
    # i=1
    l <- vals[[i]]
    dist <- l$dist
    params <- l$params
    fname <- paste("q",dist,sep="")
    samp[,i+1] <- do.call(fname,c(list(p=sob[,i]),params))
  }
  
  # Convert matrix to data frame and add labels
  samp <- as.data.frame(samp)
  names(samp) <- c("n",laply(vals, function(l) l$var))
  return(samp)
}

Rootdata <- function(n, vals){
  require("sensitivity")
  dataS5 <- read.csv("./HydraulicAnatomy/Data for Fig. S5.csv")
  Sum_up <- read_excel("./HydraulicAnatomy/Sum_up_03.xlsx", col_types = c("text", rep("numeric", 22)))
  sct <- read_csv("./HydraulicAnatomy/www/sct_B73.csv")
  
  Sum_up <- Sum_up%>%
    mutate(r_CT = sqrt(RXSA/pi),
           r_stele = sqrt(TSA/pi),
           d_stele = r_stele*2,
           r_cortex = r_CT - r_stele)
  sct <- sct %>% 
    transmute(name = "Heymans et al. 2020",
              RXSA = pi*(CT/1000)^2,
              TSA = pi*(stele/1000)^2,
              XVA = nX*pi*(xylem/2000)^2,
              CF = CF,
              n_X = nX)
  dataS5 <- dataS5%>%
    transmute(name = "Yang et al. 2019",
              RXSA = RXA,
              TSA = SXA,
              XVA = MXA,
              CF = CF,
              n_X = MXN)
  
  Data <- Sum_up%>%select(name, RXSA, TSA, XVA, CF, n_X)
  Data$n_X[Data$name == "Burton et al. 2013"] <- NA
  Data <- rbind(Data, sct, dataS5)
  
  Data <- Data %>% mutate(log_RXA = log(RXSA, exp(1)),
                          log_TSA = log(TSA, exp(1)),
                          log_XVA = log(XVA, exp(1)),
                          log_CF = log(CF, exp(1)),
                          mod_XVA = sqrt(abs(log_XVA)),
                          log_nX = log(n_X, exp(1)),
                          CW = sqrt(RXSA/pi)-sqrt(TSA/pi),
                          log_CW = log(CW))
  
  
  fit <- aov(log_TSA ~ log_RXA, data = Data)
  RXA_TSA_coef <- fit$coefficients
  fit <- aov(mod_XVA ~ log_TSA, data = Data)
  XVA_TSA_coef <- fit$coefficients
  fit <- aov(log_nX ~ log_TSA, data = Data)
  MXN_TSA_coef <- fit$coefficients
  fit <- aov(log_CF ~ log_CW, data = Data)
  CF_CW_coef <- fit$coefficients
  
  detach("package:sensitivity", unload=TRUE)
  X1 <- makeMCSample(n,vals, p = 1)%>%
    select(-n)
  X2 <- makeMCSample(n,vals, p = 2)%>%
    select(-n)
  
  library(sensitivity)
  Sobol <- sobol2007(model = NULL, X1, X2, nboot = 500)
  
  return(list (Sobol, coef <- list(RXA_TSA_coef, XVA_TSA_coef, MXN_TSA_coef, CF_CW_coef)))
  
}

Rootdetails <- function (Sobol_inf){
  
  Sobol <- Sobol_inf[[1]]
  
  RXA_TSA_coef <- Sobol_inf[[2]][[1]]
  XVA_TSA_coef<- Sobol_inf[[2]][[2]]
  MXN_TSA_coef<- Sobol_inf[[2]][[3]]
  CF_CW_coef <- Sobol_inf[[2]][[4]]
  
  dataframe <- Sobol$X
  Sampl <- dataframe%>%
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
           TCA = RXA-TSA,
           aerenchyma = ifelse(radius < 0.2, aerenchyma/10, aerenchyma),
           km = 3E-5)
  return(Sampl)
  
}

do_GranarMECHA <- function(data, params){
  vect <- data$n
  for(i in vect){
    tmp_sampl <- data[data$n == i,]
    check_granar <- try(run_granar(params, tmp_sampl, i), silent = T)
    check_mecha <- try(run_Mecha(tmp_sampl, i), silent = T)
  }
}

revdet <- function(data, Sobol_inf, design){
  
  RXA_TSA_coef <- Sobol_inf[[2]][[1]]
  XVA_TSA_coef<- Sobol_inf[[2]][[2]]
  MXN_TSA_coef<- Sobol_inf[[2]][[3]]
  CF_CW_coef <- Sobol_inf[[2]][[4]]
  
  Sampl <- data%>%
    mutate(radius = sqrt(RXA/pi),
           log_RXA = log(pi*radius^2),
           log_TSA = RXA_TSA_coef[1]+RXA_TSA_coef[2]*log_RXA,
           TSA_expected = exp(log_TSA),
           var_stele = ((TSA_in-TSA_expected)/TSA_expected)*100,
           log_TSA = log(TSA_expected, exp(1)),
           r_stele = sqrt(TSA_in/pi),
           log_nX = MXN_TSA_coef[1]+MXN_TSA_coef[2]*log_TSA,
           nX = exp(log_nX),
           mod_XVA = XVA_TSA_coef[1]+XVA_TSA_coef[2]*log_TSA,
           MXA_expected = exp(-mod_XVA^2),
           var_xylem = ((MXA_in-MXA_expected)/MXA_expected)*100,
           TCA = RXA-TSA,
           aerenchyma = AA/TCA)
  
  hyd <- design%>%
    select("kw", "kAQP", "thickness")%>%
    mutate(sampl_id = 1:nrow(design))
  
  Sampl <- left_join(Sampl, hyd, by = "sampl_id")
  return(Sampl)
  
}

get_GMecha <- function(design, first = 1){
  
  if (first == 1){
  path1 = "~/Thesis/2020-10 HydraulicViper/MECHA/Projects/GRANAR/out/wide/"
  path2 = "~/Thesis/2020-10 HydraulicViper/MECHA/Projects/GRANAR/out/wide_pl/"
  fls <- list.files(path1)
  fls <- fls[grepl(".txt", fls)]
  
  K <- tibble(kr = NULL, kx = NULL, sampl_id = NULL, apo = NULL)
  for (k in fls){
    M <- read_file(paste0(path1,k))
    tmp_M <- strsplit(M, split="\n")[[1]]
    K_xyl_spec <- as.numeric(strsplit(tmp_M[15], " ")[[1]][5])
    kr_M <- as.numeric(strsplit(tmp_M[17], " ")[[1]][4])
    scenario <- round(parse_number(unlist(str_split(k,"_"))[3])/10)
    sampl_id <- parse_number(unlist(str_split(k,"_"))[4])
    if(sampl_id == 1866){
      next
    }
    
    K <- rbind(K, tibble(kr = kr_M, kx = K_xyl_spec, sampl_id = sampl_id, apo = scenario))
  }
  
  K <- K%>%arrange(sampl_id)
  
  missing = K$sampl_id[K$apo == 2]
  
  K <- K %>% filter(sampl_id %in% missing)
  
  K <- tibble(kr_1 = K$kr[K$apo == 1],
              kr_2 = K$kr[K$apo == 2],
              kr_3 = K$kr[K$apo == 4],
              log_kx = log(K$kx[K$apo == 1]),
              sampl_id = K$sampl_id[K$apo == 1],
              kpl = 5.3E-12)
  
  
  fls <- list.files(path2)
  fls <- fls[grepl(".txt", fls)]
    Kplasmo <- tibble(kr = NULL, kx = NULL, sampl_id = NULL, apo = NULL)
  for (k in fls){
    M <- read_file(paste0(path2,k))
    tmp_M <- strsplit(M, split="\n")[[1]]
    K_xyl_spec <- as.numeric(strsplit(tmp_M[15], " ")[[1]][5])
    kr_M <- as.numeric(strsplit(tmp_M[17], " ")[[1]][4])
    scenario <- round(parse_number(unlist(str_split(k,"_"))[3])/10)
    sampl_id <- parse_number(unlist(str_split(k,"_"))[4])
    Kplasmo <- rbind(Kplasmo, tibble(kr = kr_M, kx = K_xyl_spec, sampl_id = sampl_id, apo = scenario))
  }
  
  Kplasmo <- Kplasmo%>%arrange(sampl_id)
  missing = Kplasmo$sampl_id[Kplasmo$apo == 2]
  Kplasmo <- Kplasmo %>% filter(sampl_id %in% missing)
  Kplasmo <- tibble(kr_1 = Kplasmo$kr[Kplasmo$apo == 1],
                    kr_2 = Kplasmo$kr[Kplasmo$apo == 2],
                    kr_3 = Kplasmo$kr[Kplasmo$apo == 4],
                    log_kx = log(Kplasmo$kx[Kplasmo$apo == 1]),
                    sampl_id = Kplasmo$sampl_id[Kplasmo$apo == 1])
  
  hyd <- design%>%
    select("kpl")%>%
    mutate(sampl_id = 1:nrow(design))
  Kplasmo <- left_join(Kplasmo, hyd, by = "sampl_id")
  K <- rbind(K, Kplasmo)
  }
  if(first == 2){
    path = "~/Thesis/2020-10 HydraulicViper/MECHA/Projects/GRANAR/out_plus/"
    fls <- list.files(path)
    fls <- fls[grepl(".txt", fls)]
    K <- tibble(kr = NULL, kx = NULL, sampl_id = NULL, apo = NULL)
    for (k in fls){
      M <- read_file(paste0(path,k))
      tmp_M <- strsplit(M, split="\n")[[1]]
      K_xyl_spec <- as.numeric(strsplit(tmp_M[15], " ")[[1]][5])
      kr_M <- as.numeric(strsplit(tmp_M[17], " ")[[1]][4])
      scenario <- round(parse_number(unlist(str_split(k,"_"))[3])/10)
      sampl_id <- parse_number(unlist(str_split(k,"_"))[4])
      K <- rbind(K, tibble(kr = kr_M, kx = K_xyl_spec, sampl_id = sampl_id, apo = scenario))
    }
    
    K <- K%>%arrange(sampl_id)
    missing = K$sampl_id[K$apo == 2]
    K <- K %>% filter(sampl_id %in% missing)
    K <- tibble(kr_1 = K$kr[K$apo == 1],
                      kr_2 = K$kr[K$apo == 2],
                      kr_3 = K$kr[K$apo == 4],
                      log_kx = log(K$kx[K$apo == 1]),
                      sampl_id = K$sampl_id[K$apo == 1])
    
    hyd <- design%>%
      select("kpl")%>%
      mutate(sampl_id = 1:nrow(design))
    K <- left_join(K, hyd, by = "sampl_id")
    
  }
  
  return(K)
  
}

train_model <- function(data_GM, n_rep = 5){
  
  require(caret)
  # Run algorithms using 10-fold cross validation
  control <- trainControl(method="cv", number=10)
  metric <- "RMSE"
  to_est <- c("kr_1", "kr_2", "kr_3")
  descr = c("radius", "var_stele", "var_xylem", "aerenchyma", "kw", "kAQP", "kpl", "thickness")
  
  ACCU <- NULL
  memo = c(10,10,10)
  op <- options(digits.secs = 6)
  for(i in 1:n_rep){
    for(est in to_est){
      print(est)
      print(i)
      pos <- which(colnames(data_GM) == est)
      seed <- round(parse_number(unlist(str_split(as.character(Sys.time()), ":"))[3])*10000)
      set.seed(seed)
      
      y <- c(t(data_GM[,pos]))
      validation_index <- createDataPartition(y, p=0.80, list=FALSE)
      # select 20% of the data for validation
      validation <- data_GM[-validation_index,]
      # use the remaining 80% of data to training and testing the models
      dataset <- data_GM[validation_index,]
      
      global <- dataset%>%
        mutate(id = 1:nrow(dataset))
      
      descriptor <- global%>%select(descr, "id")#-kr_1,-kr_2,-kr_3,-sampl_id)
      ground_truth <- global%>%select(c(est, "id"))
      colnames(ground_truth)<- c("kr", "id")
      train <- left_join(ground_truth, descriptor, by="id")%>%
        select(-id)
      
      # cart
      fit.cart <- train(kr~., data=train, method="rpart", metric=metric, trControl=control)
      # KNN
      fit.knn <- train(kr~., data=train, method="knn", metric=metric, trControl=control)
      # Advanced algorithms
      fit.svm <- train(kr~., data=train, method="svmRadial", metric=metric, trControl=control)
      # Random Forest
      fit.rf <- train(kr~., data=train, method="rf", metric=metric, trControl=control)
      
      results <- resamples(list(cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf))
      summary(results)
      
      # svm
      predictions_svm <- predict(fit.svm, validation)
      accuracy = tibble(x = validation[,pos], y = predictions_svm)
      t <- t.test(accuracy$x, accuracy$y, alternative = "two.sided" )
      RMSE = sqrt(sum((accuracy$x-accuracy$y)^2)/nrow(accuracy))
      nrrmse = (RMSE/(max(accuracy$x)-min(accuracy$x)))*100
      fit <- lm(unlist(accuracy$x) ~ accuracy$y)
      tmp_svm <- tibble(kr = est, model_type = "svm", rmse = RMSE, nrrmse = nrrmse, 
                        rsquare = summary(fit)$r.squared, t_pval= t$p.value, seed = seed)
      # random forest
      predictions_rf <- predict(fit.rf, validation)
      accuracy = tibble(x = validation[,pos], y = predictions_rf)
      t <- t.test(accuracy$x, accuracy$y, alternative = "two.sided" )
      RMSE = sqrt(sum((accuracy$x-accuracy$y)^2)/nrow(accuracy))
      nrrmse = (RMSE/(max(accuracy$x)-min(accuracy$x)))*100
      fit <- lm(unlist(accuracy$x) ~ accuracy$y)
      tmp_rf <- tibble(kr = est, model_type = "rf", rmse = RMSE, nrrmse = nrrmse, 
                       rsquare = summary(fit)$r.squared,  t_pval= t$p.value,seed = seed)
      
      # mean prediction
      pred <- tibble(svm = predictions_svm, rf = predictions_rf)%>%
        mutate(mm = (svm+rf)/2)
      accuracy = tibble(x = validation[,pos], y = pred$mm)
      t <- t.test(accuracy$x, accuracy$y, alternative = "two.sided" )
      RMSE = sqrt(sum((accuracy$x-accuracy$y)^2)/nrow(accuracy))
      nrrmse = (RMSE/(max(accuracy$x)-min(accuracy$x)))*100
      fit <- lm(unlist(accuracy$x) ~ accuracy$y)
      tmp_mm <- tibble(kr = est, model_type = "mm", rmse = RMSE, nrrmse = nrrmse, 
                       rsquare = summary(fit)$r.squared,  t_pval= t$p.value, seed = seed)
      
      ACCU = rbind(ACCU, tmp_svm, tmp_rf, tmp_mm)
      
      if(est == "kr_1"){
        if(tmp_svm$nrrmse < memo[1]){
          save(fit.svm, file = "./R/GRANAR/rf/svm_model_kr1.RData")
          model_kr1 = fit.svm
          memo[1] <- tmp_svm$nrrmse 
        }
        if(tmp_mm$nrrmse < memo[1]){
          save(fit.rf, file = "./R/GRANAR/rf/rf_model_kr1.RData")
          save(fit.rf, file = "./R/GRANAR/rf/rf_model_kr1.RData")
          model_kr1 = fit.rf
          memo[1] <- tmp_rf$nrrmse 
        }
        
      }
      if(est == "kr_2"){
        if(tmp_svm$nrrmse < memo[2]){
          save(fit.svm, file = "./R/GRANAR/rf/svm_model_kr2.RData")
          model_kr2 = fit.svm
          memo[2] <- tmp_svm$nrrmse 
        }
        if(tmp_mm$nrrmse < memo[2]){
          save(fit.svm, file = "./R/GRANAR/rf/svm_model_kr2.RData")
          save(fit.rf, file = "./R/GRANAR/rf/rf_model_kr2.RData")
          model_kr2 = fit.rf
          memo[2] <- tmp_rf$nrrmse 
        }
      }
      if(est == "kr_3"){
        if(tmp_svm$nrrmse < memo[3]){
          save(fit.svm, file = "./R/GRANAR/rf/svm_model_kr3.RData")
          model_kr3 = fit.svm
          memo[3] <- tmp_svm$nrrmse 
        }
        if(tmp_mm$nrrmse < memo[3]){
          save(fit.svm, file = "./R/GRANAR/rf/svm_model_kr3.RData")
          save(fit.rf, file = "./R/GRANAR/rf/rf_model_kr3.RData")
          model_kr3 = fit.rf
          memo[3] <- tmp_rf$nrrmse 
        }
      }
      
    }
  }
  return(ACCU)
}

Morris_GM <- function(model_list){
  require(sensitivity)
  to_est <- c("kr_1", "kr_2", "kr_3")
  descr = c("radius", "var_stele", "var_xylem", "aerenchyma", "kw", "kAQP", "kpl", "thickness")
  MAX = c(1.5, 40, 30, 0.4, 2.4E-3, 4.3E-3, 5.3E-11, 2) # Kpl 5.3E-10
  MIN = c(0.1,-30,-30, 0 , 2.5E-5 , 4.5E-5, 5.7E-13, 0.8)
  
  
  Morris <- morris(model = NULL, factors = 8, r = 20, design = list(type = "oat", levels = 8, grid.jump = 4),
                   binf= MIN,  bsup = MAX, scale = T)
  
  svm <- model_list[[1]]
  rf <- model_list[[2]]
  da <- momo(Morris, svm, rf)
  da$apo = 1
  
  svm <- model_list[[3]]
  rf <- model_list[[4]]
  db <- momo(Morris, svm, rf)
  db$apo = 2
  
  svm <- model_list[[5]]
  rf <- model_list[[6]]
  dc <- momo(Morris, svm, rf)
  dc$apo = 3
  
  mor = rbind(da, db , dc)
  return (mor)
  
}
momo <- function(Morris, svm, rf){
  descr = c("radius", "var_stele", "var_xylem", "aerenchyma", "kw", "kAQP", "kpl", "thickness")
  design <- as.data.frame(Morris$X)
  colnames(design) <- descr
  
  kr = (predict(svm, design)+ predict(rf, design))/2
  tell(Morris, kr)
  mu.star <- apply(Morris$ee, 2, function(Morris) mean(abs(Morris)))
  sigma <- apply(Morris$ee, 2, sd)
  x <- as.numeric(mu.star)
  y <- as.numeric(sigma)
  da = tibble(mu = x , sigma = y, name = descr)
  return(da)
}

qlunif <- function(p, min, max, base=exp(1)) {
  if (mode(p) != "numeric")
    stop("'p' must be a non-empty numeric vector")
  if (any(missing(min), missing(max)))
    stop("'min' and 'max' not provided, without default.\n")
  return(base ^ (log(min, base) + (log(max, base) - log(min, base)) * p))
}

Fast_GM <- function(model_list){
  descr = c("radius", "var_stele", "var_xylem", "aerenchyma", "kw", "kAQP", "kpl", "thickness")
  MAX = c(1.5, 40, 30, 0.4, 2.4E-3, 4.3E-3, 5.3E-11, 2)
  MIN = c(0.1,-30,-30, 0 , 2.5E-5 , 4.5E-5, 5.7E-13, 0.8)
  FAST <- fast99(model = NULL, descr, n = 2^11, M = 8, omega = NULL, 
                 q = c(rep("qunif", 4), rep("qlunif", 3), "qunif" ),
                 q.arg = list(list(min = MIN[1], max = MAX[1]),
                              list(min = MIN[2], max = MAX[2]),
                              list(min = MIN[3], max = MAX[3]),
                              list(min = MIN[4], max = MAX[4]),
                              list(min = MIN[5], max = MAX[5], base = 10),
                              list(min = MIN[6], max = MAX[6], base = 10),
                              list(min = MIN[7], max = MAX[7], base = 10),
                              list(min = MIN[8], max = MAX[8])))
  
  svm <- model_list[[1]]
  rf <- model_list[[2]]
  Sa <- fast(FAST, svm, rf)
  Sa$apo = 1
  
  svm <- model_list[[3]]
  rf <- model_list[[4]]
  Sb <- fast(FAST, svm, rf)
  Sb$apo = 2
  
  svm <- model_list[[5]]
  rf <- model_list[[6]]
  Sc <- fast(FAST, svm, rf)
  Sc$apo = 3
  
  S = rbind(Sa, Sb , Sc)
  return (S)
  
}
fast <- function (FAST, svm, rf){
  design <- FAST$X
  kr = (predict(svm, design)+ predict(rf, design))/2
  tell(FAST, y = kr)
  S <- cbind(c(FAST$D1 / FAST$V,  1 - FAST$Dt / FAST$V - FAST$D1 / FAST$V),rep(colnames(design),2))
  colnames(S) = c("effect", "name")
  S <- S %>% as_data_frame()%>%
    mutate(main = c(rep("main",8),rep("inter",8)))
  return(S)
}

Sobol_GM <- function(model_list){
  vals <- list(list(var="radius",dist="unif",params=list(min=0.1,max=1.5)),
               list(var="var_stele",dist="unif",params=list(min=-0.3,max=0.4)),
               list(var="var_xylem",dist="unif",params=list(min=-0.3,max=0.3)),
               list(var="aerenchyma",dist="unif",params=list(min=0,max=0.5)),
               list(var="kw",dist="lunif",params=list(min=0.1*2.4E-4,max=10*2.4E-4)),
               list(var="kAQP",dist="lunif",params=list(min=0.1*4.3E-4,max=10*4.3E-4)),
               list(var="kpl",dist="lunif",params=list(min=0.1*5.3E-12,max=10*5.3E-12)),
               list(var="thickness",dist="unif",params=list(min=0.8,max=2)))
  
  Sobol_inf <- Rootdata(n = 2^(round(runif(1,11,15))), vals)
  Sobol<- Sobol_inf[[1]]
  
  svm <- model_list[[1]]
  rf <- model_list[[2]]
  Sa <- sobo(Sobol, svm, rf)
  Sa$apo = 1
  
  svm <- model_list[[3]]
  rf <- model_list[[4]]
  Sb <- sobo(Sobol, svm, rf)
  Sb$apo = 2
  
  svm <- model_list[[5]]
  rf <- model_list[[6]]
  Sc <- sobo(Sobol, svm, rf)
  Sc$apo = 3
  
  S = rbind(Sa, Sb , Sc)
  
}
sobo <- function (Sobol, svm, rf){
  descr = c("radius", "var_stele", "var_xylem", "aerenchyma", "kw", "kAQP", "kpl", "thickness")
  design <- Sobol$X
  kr = (predict(svm, design)+ predict(rf, design))/2
  kr = (kr-mean(kr))/sd(kr)
  if(length(which(is.na(kr)))> 0){print("na values")}
  #kr = (kr - mean(kr))/sd(kr)
  tell(Sobol, y = kr)
  
  vardec <- tibble(name = descr,
                    main = Sobol$S$original,
                    mainlow = Sobol$S$`min. c.i.`,
                    mainhigh = Sobol$S$`max. c.i.`,
                    tot = Sobol$T$original,
                    totlow = Sobol$T$`min. c.i.`,
                    tothigh = Sobol$T$`max. c.i.`)
  
  print(plot(Sobol))
  
  S <- tibble(name = rep(vardec$name,2),
              effect = c(vardec$main, vardec$tot-vardec$main),
              main = c(rep("main", 8),rep("inter", 8)))
  
  return(S)
}


