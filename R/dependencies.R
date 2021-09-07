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

# granar mecha short-cut : random forest model
#source("/home/users/a/h/aheymans/GSUA/R/GRANAR/rf/machine_learning.r")
#PkgTest(c("randomForest"))
# next line is not working ...
# how to load .RData
library(caret)

  load("/home/users/a/h/aheymans/GSUA/R/GRANAR/rf/svm_model_kr1.RData")
  model_kr1svm <- fit.svm
  
  load("/home/users/a/h/aheymans/GSUA/R/GRANAR/rf/rf_model_kr1.RData")
  model_kr1rf <- fit.rf


  load("/home/users/a/h/aheymans/GSUA/R/GRANAR/rf/svm_model_kr2.RData")
  model_kr2svm <- fit.svm
  load("/home/users/a/h/aheymans/GSUA/R/GRANAR/rf/rf_model_kr2.RData")
  model_kr2rf <- fit.rf


  load("/home/users/a/h/aheymans/GSUA/R/GRANAR/rf/svm_model_kr3.RData")
  model_kr3svm <- fit.svm 
  load("/home/users/a/h/aheymans/GSUA/R/GRANAR/rf/rf_model_kr3.RData")
  model_kr3rf <- fit.rf
