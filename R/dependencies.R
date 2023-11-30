# Depedencies ----------------------
library(readxl) #
library(Matrix)
library(data.table)
library(tidyverse) #
library(xml2) #
`%!in%` <- compose(`!`, `%in%`)

library(plyr) 
library(deldir) #
library(alphahull) #
library(sp)
source("./R/GRANAR/granar.R")
library(viridis)

# granar mecha short-cut : meta-model
# load .RData
library(caret)

load("./R/GRANAR/rf/svm_model_kr1.RData")
model_kr1svm <- fit.svm
  
load("./R/GRANAR/rf/rf_model_kr1.RData")
model_kr1rf <- fit.rf

load("./R/GRANAR/rf/svm_model_kr2.RData")
model_kr2svm <- fit.svm

load("./GSUA/R/GRANAR/rf/rf_model_kr2.RData")
model_kr2rf <- fit.rf

load("./R/GRANAR/rf/svm_model_kr3.RData")
model_kr3svm <- fit.svm 

load("./R/GRANAR/rf/rf_model_kr3.RData")
model_kr3rf <- fit.rf
