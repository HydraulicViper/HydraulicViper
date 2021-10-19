# HydraulicViper

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
 