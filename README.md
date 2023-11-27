# HydraulicViper
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5582750.svg)](https://doi.org/10.5281/zenodo.5582750)
## How to install

### On linux

The needed modules are:

- CMake (tested working version 3.9.6)
- GCC (tested working version 7.3.0-2.30)
- Python (tested working version 3.7.4)
- R (version 3.6.1)

Then, compile:

- CPlantBox
- Hydrus-1D

Once it is done, you can run Hydraulic Viper through the Rscript B73_lite.R

The script load a sampling matrix for 12 root traits (anatomical traits, subcellular hydraulic traits, and maturation rates).

It load 3 soil types (loam, sandy loam and a silty clay soil).

The results are then stored in the "res" folder.

### Folder content:

- CplantBox contains the CplantBox version used to generate the root system structure.
- HydaulicAnatomy has the Procedure for doing the micro hydraulic part. It contains the MECHA code and the script used to generate anatomies with GRANAR, as well as the script used to train and validate meta-models is included.
- Hydrus contains Hydrus-1D version with the Couvreur model implemented.
- R is the folder with most of the actual code for Hydraulic Viper (mainly B73_function.R, and getSUF.R which is the MARSHAL)
- res has some of the output files which can be created when running the "../B73_lite.R" script. 
- sampling contains the sampling matrix used to run the main root water uptake simulations.
- WeatherData is a folder including scripts to gather climatic data from GPS coordinates and Penman-Monteith equations related to the estimation of ET0

### GRANAR - MECHA: Global sensitivity analysis of the hydraulic properties at the cross section level

#### GRANAR
To generate root cross section:

```{r}

params <- read_param_xml("./GRANAR/www/Zea_mays_CT.xml")
Sampl <- read.csv("./Anat_proc/Sobol_inputGM.csv")

for(i in 1:nrow(Sampl)){
    # Get the input parameter from the sampling matrix
  tmp_sampl <- Sampl[i,]

    # Run GRANAR
  sim <- run_granar(params, tmp_sampl)

   # Save the generated root cross section
  file.copy("./MECHA/cellsetdata/current_root.xml", paste0("./MECHA/cellsetdata/root_",i,".xml"), overwrite = T)
   # Save the id of the aerenchyma polygons
  file.copy("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml",
            paste0( "./MECHA/Projects/GRANAR/in/Maize_Geometry_aer_",i,".xml"), overwrite = T)
}

```

#### MECHA
To estimate the root radial hydraulic conductivities of the generate root cross section, the following script was used:

```{r}

fls <- list.files("./MECHA/cellsetdata/")
fls <- fls[grepl("root_", fls)]

for(j in fls){

  fc <- file.copy(paste0("./MECHA/cellsetdata/",j), "./MECHA/cellsetdata/current_root.xml", overwrite = T)
  fc <- file.copy(paste0("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer_", parse_number(j), ".xml"),
              paste0("./MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml"), overwrite = T)

    # MECHA gets the input parameter from the sampling matrix
  id <- parse_number(j) 
  microhydro(path = "MECHA/Projects/GRANAR/in/Maize_hydraulics.xml",
             kw = Sampl$kw[parse_number(j)],
             km = Sampl$km[parse_number(j)], 
             kAQP = Sampl$kAQP[parse_number(j)],
             kpl = Sampl$kpl[parse_number(j)])
  wallthick(path = "MECHA/Projects/GRANAR/in/Maize_Geometry_aer.xml", Sampl$thickness[parse_number(j)])

  # Run MECHA - - - - - - - 
  system("python3 ./MECHA/MECHAv4_septa.py")

  # Save MECHA output
  file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_1,0.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_1,0_",id,".txt"), overwrite = T)
  file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_2,1.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_2,1_",id,".txt"), overwrite = T)
  file.copy("./MECHA/Projects/GRANAR/out/M1v4/Root/Project_Test/results/Macro_prop_4,2.txt",
              paste0("./MECHA/Projects/GRANAR/out/M1v4/Root/Macro_prop_4,2_",id,".txt"), overwrite = T)

}

```

#### Training the meta-model

Once all data are collected from the coupling between GRANAR & MECHA, we train different machine learning algorithms and compared them.
Saved the two best one in each scenarios.

```{r}

All_data <- read.csv("./data_GRANAR_MECHA.csv")

library(caret)
# Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=10)
metric <- "RMSE"

to_est <- c("kr_1", "kr_2", "kr_3")
descr = c("radius", "var_stele", "var_xylem", "aerenchyma", "kw", "km", "kAQP", "kpl", "thickness")

ACCU <- NULL
memo = c(10,10,10)
op <- options(digits.secs = 6)
for(i in 1:10){
  for(est in to_est){
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
    print(dotplot(results))

    predictions <- predict(fit.svm, validation)
    accuracy = tibble(x = validation[,pos], y = predictions)
    RMSE = sqrt(sum((accuracy$x-accuracy$y)^2)/nrow(accuracy))
    nrrmse = (RMSE/(max(accuracy$x)-min(accuracy$x)))*100
    fit <- lm(unlist(accuracy$x) ~ accuracy$y)
    tmp_svm <- tibble(kr = est, model_type = "svm", rmse = RMSE, nrrmse = nrrmse, rsquare = summary(fit)$r.squared, seed = seed)

    predictions <- predict(fit.knn, validation)
    accuracy = tibble(x = validation[,pos], y = predictions)
    RMSE = sqrt(sum((accuracy$x-accuracy$y)^2)/nrow(accuracy))
    nrrmse = (RMSE/(max(accuracy$x)-min(accuracy$x)))*100
    fit <- lm(unlist(accuracy$x) ~ accuracy$y)
    tmp_knn <- tibble(kr = est, model_type = "knn", rmse = RMSE, nrrmse = nrrmse, rsquare = summary(fit)$r.squared, seed = seed)

    predictions <- predict(fit.cart, validation)
    accuracy = tibble(x = validation[,pos], y = predictions)
    RMSE = sqrt(sum((accuracy$x-accuracy$y)^2)/nrow(accuracy))
    nrrmse = (RMSE/(max(accuracy$x)-min(accuracy$x)))*100
    fit <- lm(unlist(accuracy$x) ~ accuracy$y)
    tmp_cart <- tibble(kr = est, model_type = "cart", rmse = RMSE, nrrmse = nrrmse, rsquare = summary(fit)$r.squared, seed = seed)
    
    predictions <- predict(fit.rf, validation)
    accuracy = tibble(x = validation[,pos], y = predictions)
    t <- t.test(accuracy$x, accuracy$y, alternative = "two.sided" )
    RMSE = sqrt(sum((accuracy$x-accuracy$y)^2)/nrow(accuracy))
    nrrmse = (RMSE/(max(accuracy$x)-min(accuracy$x)))*100
    fit <- lm(unlist(accuracy$x) ~ accuracy$y)
    tmp_rf <- tibble(kr = est, model_type = "rf", rmse = RMSE, nrrmse = nrrmse, rsquare = summary(fit)$r.squared, seed = seed)

    ACCU = rbind(ACCU, tmp_svm, tmp_rf, tmp_knn, tmp_cart)
  }
}

```

## Run Hydraulic Viper

Load input

Sampling matrix

| Abbreviation, Unit      | Variable description                                             | Range            |
|-------------------------|-------------------------------------------------------------------|------------------|
| $Cv_{MX}$ (\%)          | Coefficient of variation from the average ratio between            | [-20; 20]        |
|                         | the meta-xylem area over the stele area.                           |                  |
| $Cv_{S}$ (\%)           | Coefficient of variation from the average ratio between the stele  | [-30; 40]        |
|                         | area and the root cross section area.                              |                  |
| $R_{Lat}$ ($\mu m$)     | Lateral root radius                                               | [100; 150]       |
| $R_{Sem}$ ($\mu m$)     | Seminal root radius                                               | [200; 400]       |
| $R_{SBR}$ ($\mu m$)     | Shoot born root radius                                            | [300; 700]       |
| $P_{aer}$ (\%)          | Percentage of aerenchyma area over the root cortex area.           | [0; 0.5]         |
| $k_W$ ($cm^2 hPa^{-1} d^{-1}$) | Cell wall hydraulic conductivity                           | [$1E^{-4}$; $4E^{-4}$] |
| $K_{PL}$ ($cm^3 hPa^{-1} d^{-1}$) | Plasmodesmata hydraulic conductance                    | [$6E^{-13}$; $1E^{-11}$] |
| $k_{AQP}$ ($cm^2 hPa^{-1} d^{-1}$) | Aquaporin contribution to the cell membrane hydraulic conductivity. | [$2E^{-4}$; $7.5E^{-4}$] |
| $dt_{S-MX}$ (day)       | Time delay between the beginning of the endodermis suberization   | [-0.5; 4]        |
|                         | and of the meta-xylem maturation.                                 |                  |
| $dt_{MX}$ (day)         | Duration of the maturation of the meta-xylem vessels               | [1; 4]           |
| $dt_{Hb}$ (day)         | Duration of the deposition of hydrophobic barriers or the delay    | [1; 4]           |
|                         | between two levels of hydrophobic barriers.                        |                  |



WEATHER


soil_list
age_max 

```{r}
MatSobol<- read.csv(paste0("./www/SamplingB73_Sobol.csv"))

```




