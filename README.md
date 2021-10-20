# HydraulicViper
[![DOI](https://zenodo.org/badge/DOI10.5281/zenodo.5582750.svg)](https://doi.org/10.5281/zenodo.5582750)

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
- HydaulicAnatomy has the MECHA code and the script used to generate anatomies with GRANAR. The script used to train and validate meta-models is included.
- Hydrus contains Hydrus-1D version with the Couvreur model implemented.
- R is the folder with most of the actual code for Hydraulic Viper (mainly B73_function.R, and getSUF.R which is the MARSHAL)
- res has some of the output files which can be created when running the "../B73_lite.R" script. 
- sampling contains the sampling matrix used to run the main root water uptake simulations.
- WeatherData is a folder including scripts to gather climatic data from GPS coordinates and Penman-Monteith equations related to the estimation of ET0

