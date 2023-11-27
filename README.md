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

To generate root cross section, the following script was used:

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


