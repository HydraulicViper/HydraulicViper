##################################################################################
#                                                                                #
#                                                                                #
#           *********************************************************            #
#           *                                                       *            #
#           *                  INTERNAL FUNCTIONS                   *            #
#           *                                                       *            #
#           *               non-user interface functions            *            #
#           *                                                       *            #
#           *********************************************************            #
#                                                                                #
#                                                                                #
##################################################################################


#----------------------------------------------------------------------
#FUNCTION PkgTest installs required packages if they are not installed
#----------------------------------------------------------------------
#ARGUMENTS:
#str.pkg			the package name given in qoutes, e.g. "xtable"
#
#RETURNS:
#TRUE / stop 	returns TRUE if succesful, otherwise stop message is output
#
#-----------------------------------------------------------------------
PkgTest <- function(vec.pkg){
  for(str.pkg in vec.pkg){
    if (!require(str.pkg,character.only = TRUE)){
      install.packages(str.pkg,dep=TRUE)
      if(!require(str.pkg,character.only = TRUE)) stop("Package not found")
    }
  }
  return(TRUE)
}

#----------------------------------------------------------------------
#FUNCTION ColVar computes the corrected variances of the coloumns,
#i.e. let n be the number of rows, then the function sums each coloumn 
#and divides each sum by n - 1
#----------------------------------------------------------------------
#ARGUMENTS:
#mat.X			   a matrix 
#
#RETURNS:
#A vector of the corrected coloumn variance, 
#where each entry is the variance of one coloumn
#
#----------------------------------------------------------------------
#EXAMPLE:
# > M <- matrix(c(1:6), ncol=2)
# > M
#      [,1] [,2]
# [1,]    1    4
# [2,]    2    5
# [3,]    3    6
# > ColVar(M)
# [1] 1 1
# 
#-----------------------------------------------------------------------
ColVar <- function(mat.X){
  n     <- NROW(mat.X) 
  res   <- (colSums( apply(mat.X,2, function(col) col - mean(col))^2))/( n - 1)
  return(res)
}

#----------------------------------------------------------------------
#FUNCTION ColMax finds the maximum value of each coloumn 
#REMARK the build-in R function maxCol computes confusingly the max of each row
#----------------------------------------------------------------------
#ARGUMENTS: 
#mat.X          a matrix
#
#RETURNS:
#A vector where each entry is the maximum of the corresponding coloumn of mat.X
#
#----------------------------------------------------------------------
#EXAMPLE:
# > M <- matrix(c(1:6), ncol=2)
# > M
#      [,1] [,2]
# [1,]    1    4
# [2,]    2    5
# [3,]    3    6
# > ColMax(M)
# [1] 3 6
# 
#-----------------------------------------------------------------------
ColMax <- function(X){
  return(apply(X,2,max))
}

#-----------------------------------------------------------------------
#FUNCTION DivideSet divides randomly the data set into a test and training set
#corresponding to the given percentage
#-----------------------------------------------------------------------
#ARGUMENTS:
#size			Number of samples in the data
#perc			Percent of data into the test set, by default 30
#rnr      A random number included in seed setting 
#         (helpful for bagging, when the method is called many times)
#
#RETURNS:
#A list of two vectors of indices, 
#where the first vector "test" is the index set for the test set  
#second vector "train" the index set for the training set
#
#-----------------------------------------------------------------------
DivideSet <- function(size, perc=30, rnr=1){
  if ( 0 < perc && perc < 100){
    set.seed(perc * size * rnr)
    if ( 1 <= floor( perc / 100 * size) )
      test <- sample(1 : floor( perc / 100 * size))
    else test <- NULL
    
    if (size >= 1+floor( perc / 100 * size))
      train <- sample((1+floor( perc / 100 * size)): size)
    else train <- NULL
    
    return( list(test=test, train=train) )
  }
  else
    stop("Percentage has to be in between 0 and 100")
}

#----------------------------------------------------------------------
#FUNCTION ValidityTest checks whether the set of predictor variables
#and the set of response variables are valid,
#i.e. containing at least one variable, being within range and
#that the response variables are a subset of the permissible data set
#----------------------------------------------------------------------
#ARGUMENTS:
#vec.f        the coloumn indices of the response variables
#vec.p        the coloumn indices of the permissible variables;
#             have to contain the response variable indices
#n            the total number of coloumns in the data set
#
#RETURNS:
#Stops if a set is not valid, otherwise it simply returns to the toplevel call
#
#-----------------------------------------------------------------------
ValidityTest <- function(vec.f, vec.p, n){
  
  # Check whether 'vec.f' is a valid parameter
  if ( length(vec.f) < 1 )
    stop("Parameter 'f_vec' needs at least one component.")
  if ( length(vec.f) == length(vec.p))
    stop("Number of response variables has to be a subset of the variables of the whole data set.")
  if (max(vec.f) > max(vec.p) || min(vec.f) < min(vec.p) )
    stop("Indices of the parameter 'vec.p' are out of range of the permissible data set.")
  
  
  # Check whether 'vec.p' is a valid parameter
  if ( length(vec.p) < 1 )
    stop("Parameter 'vec.p' needs at least one component.")
  if (max(vec.p) > n || min(vec.p) < 1 )
    stop("Indices of the parameter 'vec.p' are out of range.")
  
  if(length(which(vec.f %*% t(1 / vec.p) ==1)) >= length( vec.p ))
    stop("The response variables have to be a subset of the permissible data set.")
  
  return();
}

#----------------------------------------------------------------------
#FUNCTION PredNAvalues predicts for each coloumn the missing NA values
#with a random forest of 10 trees
#----------------------------------------------------------------------
#ARGUMENTS:
#mat.X			A matrix of data, 
#           where a row is a datum,
#           and coloumn one response variable
#
#RETURNS:
#A matrix where the NA values are replaced by predicted values
#
#-----------------------------------------------------------------------
PredictNAvalues <- function(mat.X){

  nrcol             <- ncol(mat.X)
  tempdata          <- mat.X #matrix(as.numeric(unlist(data)), ncol = nrcol)
  indexset          <- NULL
  
  #the models work only with coloumns where no entry is NA
  #so replace in every correspond coloumn the NA values by its median
  for(i in 1:nrcol){
    indices         <- which(is.na(mat.X[,i]))
    numberOfNAs     <- length(indices)
    if (numberOfNAs > 0){
      tempdata[indices,i] <- median(mat.X[,i], na.rm=T)
      indexset <- c(indexset,i)
    }
  }
  
  #compute a prediction for the corresponding NA indices based on all other values
  for(i in indexset){
    indices         <- which(is.na(mat.X[,i]));
    fit             <- randomForest(tempdata[,i] ~ . -tempdata[,i], data=as.data.frame(tempdata), ntree=10)
    mat.X[indices,i] <- predict(fit, as.data.frame(tempdata[indices,]))
  }
  
  #return the matrix where NAs have been replaced
  return(mat.X);
}

#----------------------------------------------------------------------
#FUNCTION ModelRF generates a random forest model for the response variable
#indexed by 'tofit' for the given data 'mat.data'
#-----------------------------------------------------------------------
#ARGUMENTS:
#mat.data      input data
#vec.colnrs    a vector of indices representing the response variables, i.e.
#              if the response variables are coloumns 20 and 22 of 'mat.data', then
#              vec.colnrs = c(20,22)
#              contains indices of variables that shall not be
#              considered as predictor variables for the random forest;
#REMARK:       'vec.colnrs' must not include the index 'tofit'
#nrtree        number of trees of the random forest
#tofit         the index corresponding to the response variable in mat.data that shall be fitted
#
#RETURNS:
#A random forest model for response variable given by the index 'tofit'
#
#-----------------------------------------------------------------------
ModelRF<-function(mat.data, vec.colnrs, nrtree, tofit){
  #Generate the formula for the RF
  # tofit ~ . - vec.colnrs[1] - vec.colnrs[2] - ... - vec.colnrs[n]
  string          <- NULL
  
  vec.desc        <- c(1:ncol(mat.data))
  
  vec.desc        <- vec.desc[-c(tofit,vec.colnrs)]
  
  
  if (!is.null(vec.desc))
    for (i in (vec.desc)){
      string          <- paste(string, colnames(mat.data)[i], sep=" + ")
    }
  
  formul          <- formula(paste(colnames(mat.data)[tofit], "~  ", string))
  
  resRF <- randomForest(formula    = formul, 
                        data       = mat.data, 
                        importance = TRUE, 
                        ntree      = nrtree)
  return(resRF)
}




#-----------------------------------------------------------------------
#FUNCTION ErrorComputation computes for a set of different model
#configurations their corresponding total as well as single-point error
#-----------------------------------------------------------------------
# The function iterates over two parameters
# 1) the number of trees of a random forest
# 2) the number of random forest for one response variable
# and constructs for each pair of this configuration a mode on a training set;
# then it computes the total and single-point error 
# as given by the GeneralTotalErrorEstimate and GeneralPointErrorEstimate functions
# on a test set.
#-----------------------------------------------------------------------
#ARGUMENTS:
#mat.data      A matrix of data
#vec.models    A vector comprising numbers, 
#              reflecting the number of random forests that shall be used
#              for predicting a response variable
#              
#vec.trees     A vector comprising the number of trees for a random forest
#vec.f         The index vector of response variables
#
#RETURNS:
#A named list of errors (value), 
#where each entry reflects one response variable,
#and consists of a list of two matrices,
#the first for the total error, called 'ErrorTotal'
#and the second for the single-point error, called "ErrorPoint'
#----------------------------------------------------------------------
#COMPATIBLE FUNCTION INPUT:
#DrawLines
#NormError
#SelectModel
#
#-----------------------------------------------------------------------
ErrorComputation <- function(mat.data,vec.models, vec.trees, vec.f){
  result       <- list()
  result.point <- list()
  result.total <- list()
  k            <- 0
  
  # Divide the data into training and test data
  res        <- DivideSet(dim(mat.data)[1])
  data.train <- mat.data[res$train,]
  data.test  <- mat.data[res$test,]
  
  #for each number of models given in the vector vec_models
  for (nrm in vec.models){
    print(paste("Training with nr of models:",nrm))
    k      <- k + 1
    l      <- 0
    
    # CONSTRUCT ERROR MATRICES
    # each row represents the errors for the corresponding number of trees used in the forest
    # each coloumn represents the errors for a responds variable
    # hence row i, col j represents the error of using a number of vec_trees[i] for response variable j
    error.point    <- matrix(0, ncol=length(vec.f), nrow =length(vec.trees))
    error.total    <- matrix(0, ncol=length(vec.f), nrow =length(vec.trees))
    #assign the names of the response variables to the cols of the error matrix 
    colnames(error.point) <- names(data.test)[vec.f]
    colnames(error.total) <- names(data.test)[vec.f]
    #assign for each row the trees that were used for training the model and producting the error
    rownames(error.point) <- vec.trees
    rownames(error.total) <- vec.trees
    
    
    #for each number of trees given in the vector vec_trees
    for (nrt in vec.trees){
      print(paste("Using",nrt,"trees"))
      l <- l + 1
      rfModels   <- list()
      # generate a number of folds corresponding to the number of models used in this iteration
      # each model is trained on a subset, i.e. fold
     
      n <- length(vec.f)
      
      for (i in 1:n){
        if (n==1)
          colnrs         <- NULL
        else if (i==1)
          colnrs         <- vec.f[2:n]
        else if (i==n)
          colnrs         <- vec.f[1:(n-1)]
        else
          colnrs         <- vec.f[c(1:(i-1),(i+1):n)]
        
        #GENERATE RANDOM FOREST MODELS
        lst.response.RFs <- list()
        for (j in 1:nrm){
          agg                <- DivideSet(dim(data.train)[1], rnr=j)
          data.train.model   <- data.train[agg$train,]
          lst.response.RFs[[j]] <- ModelRF(data.train, colnrs, nrt, vec.f[i])
        }
        rfModels[[i]]    <- lst.response.RFs
      }
      
      #CALCULATE THE ERRORS
      
      #select the rooted max point errors for each response variable computed on the test set
      error.point[l,] <- ColMax(GeneralPointErrorEstimate(rfModels, data.test, vec.f))
      #compute the total rooted mean squared error
      error.total[l,] <- GeneralTotalErrorEstimate(rfModels, data.test, vec.f)
    }
    
    #assign the errors for the current number of models to the result list
    result.point[[k]] <- error.point
    result.total[[k]] <- error.total
  }
  #the entries of the list are named by the number of used models
  
  mat.error.point   <- matrix(0, ncol=length(vec.models), nrow=length(vec.trees))
  mat.error.total   <- matrix(0, ncol=length(vec.models), nrow=length(vec.trees))
  
  #assign for each row the nr of random forests that were used for training the model and producing the error
  colnames(mat.error.point) <- vec.models
  colnames(mat.error.total) <- vec.models
  #assign for each row the trees that were used for training the model and producing the error
  rownames(mat.error.point) <- vec.trees
  rownames(mat.error.total) <- vec.trees
  
  for(i in 1:length(vec.f)){
    for(j in 1:length(vec.models)){
      mat.error.total[,j]   <- result.total[[j]][,i]
      mat.error.point[,j]   <- result.point[[j]][,i]
    }
    result[[i]]  <- list(Error.total=mat.error.total, Error.point=mat.error.point)
  }
  
  names(result)  <- names(data.test)[vec.f]
  
  return(result)
}


#----------------------------------------------------------------------------
#FUNCTION NormError normalizes a list of two errors to values between 0 and 1
#----------------------------------------------------------------------------
#ARGUMENTS:
#lst.errors     A list of errors (value), 
#               where each entry is considered to consist of a 2-elemental list,
#               reflecting a response variable error wrt. to total and point error; 
#               i.e. like the error data structure returned by ErrorComputation
#               lst.errors[[i]] corresponds to the errors of response variable i
#RETURNS:
#A named list, where each entry consits of a 2-elemental list
#Let be 'Error.total' and 'Error.point' the components of this 2-elemental list,
#furthermore let them be matrices, 
#then each component/matirx is normalized by its maximal value
#----------------------------------------------------------------------
#COMPATIBLE FUNCTION OUTPUT:
#ErrorComputation
#GenerateModels
#
#COMPATIBLE FUNCTION INPUT:
#WeightError
#----------------------------------------------------------------------------
NormError <- function(lst.errors){
  nrRespVar   <- length(lst.errors)
  normerror   <- lst.errors
  
  for(nrVar in 1:nrRespVar){
    maxerror                 <- max(lst.errors[[nrVar]][[1]])
    normerror[[nrVar]][[1]]  <- normerror[[nrVar]][[1]] / maxerror
    maxerror                 <- max(lst.errors[[nrVar]][[2]])
    normerror[[nrVar]][[2]]  <- normerror[[nrVar]][[2]] / maxerror
    
  }
  return(normerror)
}    

#-----------------------------------------------------------------------------
#FUNCTION WeightError computes a linearly weighted error of two type of errors.
#------------------------------------------------------------------------------
#ARGUMENTS:
#lst.errors     A list of errors (value), 
#               where each entry is considered to consist of a 2-elemental list,
#               reflecting a response variable error wrt. to total and point error; 
#               i.e. like the error data structure returned by ErrorComputation, NormError
#               lst.errors[[i]] corresponds to the errors of response variable i
#REMARK:        lst.errors[[i]][[1]] and lst.errors[[i]][[2]] 
#               have to be of equal dimension in order to be summable
#vec.weights    A vector of weights, one entry for each response variable in 'lst.errors',
#               where each entry has to be between 0 and 1, 
#               reflecting the desired impact on the total error, 
#               e.g. if the i-th weight is 0.3 then 30% of the first error will be taken
#               into account and 70% of second - this correlates to a stronger reduction
#               of the second error
#
#RETURNS:
#A named list of matrices, each matrix representing the linearly weighted sum
#of the error matrices of one response variable
#----------------------------------------------------------------------
#COMPATIBLE FUNCTION OUTPUT:
#GenerateModels
#ErrorComputation
#NormError
#-----------------------------------------------------------------------
#EXAMPLE:
# > height.total <- matrix(c(1:6), ncol=2)
# > height.point <- matrix(c(7:12), ncol=2)
# > height       <- list(total=height.total,point=height.point)
# > width.total  <- matrix(c(13:20), ncol=2)
# > width.point  <- matrix(c(21:28), ncol=2)
# > width        <- list(total=width.total,point=width.point)
# > lst.Error    <- list(height=height, width=width)
# > vec.weigths  <- c(0.3,0.6)
# > WeightError(lst.Error, vec.weigths)
#$height
#[,1] [,2]
#[1,]  5.2  8.2
#[2,]  6.2  9.2
#[3,]  7.2 10.2
#
#$width
#[,1] [,2]
#[1,] 16.2 20.2
#[2,] 17.2 21.2
#[3,] 18.2 22.2
#[4,] 19.2 23.2
#
#-----------------------------------------------------------------------
WeightError <- function(lst.error,vec.weights){
  n                      <-  length(lst.error)
  weighted.models        <-  list()
  for(i in 1:n){
    weighted.models[[i]]   <- vec.weights[i]*lst.error[[i]][[1]] + (1-vec.weights[i])*lst.error[[i]][[2]]
  }
  
  names(weighted.models) <- names(lst.error)
  return(weighted.models)
}

#--------------------------------------------------------------------------------
#FUNCTION FindBestModel finds the model with the smallest for a response variable
#--------------------------------------------------------------------------------
#ARGUMENTS:
#lst.error    A list of error matrices
#col          The response variable whose lowest error is inquired
#RETURNS:
#A two-dimensional vector, 
#where the fst entry is the index of coloumn of the lowest error
#and the snd the index of the according row 
#--------------------------------------------------------------------------------
FindBestModel <- function(lst.error, col){
  
  index       <- which.min(lst.error[[col]])
  poscol      <- (index - 1) %/% nrow(lst.error[[col]]) + 1
  posrow      <- (index - 1) %% nrow(lst.error[[col]]) + 1
  
  return(c(poscol, posrow))
}  


##################################################################################
#                                                                                #
#                                                                                #
#           *********************************************************            #
#           *                                                       *            #
#           *              USER INTERFACE FUNCTIONS                 *            #
#           *                                                       *            #
#           *                                                       *            #
#           *                                                       *            #
#           *********************************************************            #
#                                                                                #
#                                                                                #
##################################################################################
#-----------------------------------------------------------------------
#FUNCTION Plot plots a function and its according errors, 
#either as prediction/error band 'type="l"'
#or as prediction/error bars 'type="b"'
#-----------------------------------------------------------------------
#Whether it is a prediction or error depends on how the 'vec.errors'
#was computed.
#If it is a generalized error as given by the GeneralPointErrorEstimate
#function, then it is a prediction band.
#-----------------------------------------------------------------------
#ARGUMENTS:
#vec.x            values for the x-axis
#vec.y            values for the y-axis
#vec.errors       values of errors, considered as symmetric around the function
#REMARK: all vectors have to be of same length
#type             "l" for bands, "b" for bars
#
#RETURNS:
#A prediction/error band-bar plot
#----------------------------------------------------------------------
#COMPATIBLE FUNCTION OUTPUT:
#GeneralPointErrorEstimate
#-----------------------------------------------------------------------
Plot <- function(vec.x, vec.y, vec.errors, type="l"){
  
  #order the values
  ord      <- order(vec.x)
  xaxis    <- sort(vec.x)
  yaxis    <- vec.y[ord]
  if (type == "l"){
    plot(xaxis,yaxis, 
         ylim=range(c(yaxis - vec.errors[ord], yaxis + vec.errors[ord])),
         type="l", 
         col="red"
    )
    
    polygon(c(xaxis, rev(xaxis)), c(yaxis - vec.errors[ord], rev(yaxis + vec.errors[ord])), col ="grey", border = NA)
    lines(xaxis,yaxis, type="l")
    lines(xaxis, yaxis + vec.errors[ord], 
          lty   = 2, 
          type  = "l", 
          col   = "green"
    )
    lines(xaxis, yaxis - vec.errors[ord], 
          lty   = 2, 
          type  = "l", 
          col   ="green"
    )
  }
  else if (type == "b"){
    plot(xaxis, yaxis, 
         ylim=range(c(yaxis - vec.errors[ord], yaxis + vec.errors[ord])),
         pch=19
    )
    arrows(xaxis, yaxis - vec.errors[ord], xaxis, yaxis + vec.errors[ord], 
           length  = 0.05,
           angle   = 90, 
           code    =3
    )
  }
}

#----------------------------------------------------------------------
#FUNCTION ShowPerformance plots the performance of different model configurations
#for one response variable
#----------------------------------------------------------------------
#ARGUMENTS:
#lst.data     A list (mostly an error list), 
#             where each entry represents a response variable
#             consisting of two error matrices, 
#             as returned by the GenerateModels or the ErrorComputation functions
#col          The index of the response variable
#nrer         The number of the error, i.e. 1 or 2
#
#RETURNS:
#A plot showing the performance of different model parameters
#----------------------------------------------------------------------
#COMPATIBLE FUNCTION OUTPUT:
#ErrorComputation
#GenerateModels
#-----------------------------------------------------------------------
ShowPerformance <- function(lst.data, col, nrer){
  
  n      <- NCOL(lst.data[[col]][[nrer]])
  cols   <- rainbow(n)
  ymin <- min(lst.data[[col]][[nrer]])
  ymax <- max(lst.data[[col]][[nrer]])
  
  xaxis  <- rownames(lst.data[[col]][[nrer]])
  
  for(i in 2:n)
  {
    ymin <- min(c(ymin, lst.data[[col]][[nrer]]))
    ymax <- max(c(ymax, lst.data[[col]][[nrer]]))
  }
  
  plot(xaxis, lst.data[[col]][[nrer]][,1],
       ylim    = c(ymin,ymax),
       lwd     = 3, 
       col     = cols[1], 
       main    = paste(names(lst.data[[col]])[nrer], names(lst.data)[col]), 
       type    = "l", 
       ylab    = "RMSE",
       xlab    = "Nr of trees")
  
  for(i in 2:n)
    lines(xaxis,lst.data[[col]][[nrer]][,i], lwd=3,col=cols[i], type="l")
  
  legend("topright", inset=.05,title="Nr of models",legend=colnames(lst.data[[col]][[nrer]]),fill=cols)    
  
}

#----------------------------------------------------------------------
#FUNCTION GeneralPointErrorEstimate computes an estimate of the
#rooted mean squared error (RSME) generalization error for each point.
#----------------------------------------------------------------------
#That is, for each response variable the point error is computed 
#and returned as a matrix, 
#where each coloumn corresponds to a response variable
#and each entry in the coloumn to a point error 
#
#Due to the Bias-variance decomposition it holds
#      MSE = variance + bias^2 + noise
#where we set the noise to zero.
#-----------------------------------------------------------------------
#ARGUMENTS:
#lst.rfModels         a list of models; each list entry reflects the fit to one
#                     response variable; each list entry is again a list, a list of 
#                     random forests, whose mean is taken as prediction;
#                     i.e. lst.rfModels[[i]] corresponds to response variable i
#                     lst.rfModels[[i]][[j]] corresponds to the j-th random forest model
#                     forresponse variable j
#mat.data             the test data
#vec.f                an index vector of the coloumns of the data that are the 
#                     response variables
#RETURNS:
#A matrix with NROW(mat.data) and NCOL(nr of response variables), 
#where each coloumn is labeled by the name of a response variable
#and containing its point errors
#i.e. in coloumn i, row j, we have the point error of observation j of response variable i
#----------------------------------------------------------------------
#COMPATIBLE FUNCTION OUTPUT:
#GenerateModels
#SelectModel
#
#-----------------------------------------------------------------------
GeneralPointErrorEstimate <- function(lst.rfModels, mat.data, vec.f){
  
  n             <- length(lst.rfModels)
  #the result matrix
  generalError  <- matrix(0,nrow = dim(mat.data)[1],ncol=n)
  
  #for each model/response variable
  for (i in 1:n){
    nrModels  <- length(lst.rfModels[[i]])
    # temporary matrix for storing the predictions of each separate random forest
    modelPred <- matrix(0, ncol = dim(mat.data)[1], nrow = nrModels )
    
    #store the results row-wise, such that a prediction of one and the same point 
    #are in a single coloumn
    for (k in 1:nrModels){
      modelPred[k,]    <- predict(lst.rfModels[[i]][[k]], mat.data) 
    }
    #compute the variance of a point prediction
    variance         <- ColVar(modelPred)
    bias             <- colMeans(modelPred) - mat.data[,vec.f[i]]
    #bias-variance decomposition
    generalError[,i] <- variance + bias^2
  }
  #assign names to the coloumns corresponding to the response variables
  colnames(generalError)  <- names(lst.rfModels)
  return(data.frame(sqrt(generalError)))
}

#----------------------------------------------------------------------
#FUNCTION GeneralTotalErrorEstimate computes an estimate of the
#rooted mean squared error (RSME) generalization error for the whole data set.
#-----------------------------------------------------------------------
#That is, for each response variable the overall error on the whole data set 
#is computed and returned in a vector, 
#where each entry is an error of a response variable.
#
#Due to the Bias-variance decomposition it holds
#      MSE = variance + bias^2 + noise
#where we set the noise to zero.
#-----------------------------------------------------------------------
#ARGUMENTS:
#lst.rfModels         a list of models; each list entry reflects the fit to one
#                     response variable; each list entry is again a list, a list of 
#                     random forests, whose mean is taken as prediction;
#                     i.e. lst.rfModels[[i]] corresponds to response variable i
#                     lst.rfModels[[i]][[j]] corresponds to the j-th random forest model
#                     forresponse variable j
#mat.data             the test data
#vec.f                an index vector of the coloumns of the data that are the 
#                     response variables
#
#RETURNS:
#An entry named vector containing the response variable errors, 
#i.e. each entry is labeled with the name and contains the RSME of one response variable
#----------------------------------------------------------------------
#COMPATIBLE FUNCTION OUPUT:
#GenerateModels
#SelectModel
#
#COMPATIBLE FUNCTION INPUT:
#Plot
#
#-----------------------------------------------------------------------
GeneralTotalErrorEstimate <- function(lst.rfModels, mat.data,vec.f){
  #result vector of length of the number of response variables
  generalError     <- rep(0, length(vec.f))
  #predict the values, using the given models; 
  #a matrix is returned, coloumns represent response variable predictions
  prediction       <- PredictRFs(lst.rfModels, mat.data)
  #vector containing the variance of each coloumn
  variance         <- ColVar(prediction)
  #check whether more than one response variable has to predicted
  #if so, use the colMeans functions
  #otherwise the normal mean
  if (!is.null(dim(mat.data[,vec.f])))
    bias                <- colMeans(prediction) - colMeans(mat.data[,vec.f])
  else
    bias                <- colMeans(prediction, na.rm = T) - mean(mat.data[,vec.f], na.rm = T)
  #bias - variance decomposition 
  generalError        <- variance + bias^2
  #assign names to the entries reflecting the response variables
  names(generalError) <- names(lst.rfModels)
  return(sqrt(generalError))
}

#---------------------------------------------------------------------
#FUNCTION PredictRF computes the prediction of the response variables 
#for a data set given a list of models, 
#where each model corresponds to one response variable
#---------------------------------------------------------------------
#ARGUMENTS:
#lst.rfModels         a list of models; each list entry reflects the fit to one
#                     response variable; each list entry is again a list, a list of 
#                     random forests, whose mean is taken as prediction;
#                     i.e. lst.rfModels[[i]] corresponds to response variable i
#                     lst.rfModels[[i]][[j]] corresponds to the j-th random forest model
#                     forresponse variable j
#mat.data             the data to predict
#
#RETURNS:
#A matrix with NROW(mat.data) and NCOL(nr of response variables), 
#where each coloumn is labeled by the name of a response variable
#and containing its corresponding predictions 
#----------------------------------------------------------------------
#COMPATIBLE FUNCTION OUTPUT:
#GenerateModels 
#
#-----------------------------------------------------------------------
PredictRFs <- function(lst.rfModels, mat.data){

  #lst.rfModels <- model
  #mat.data <- rs1[,c(2:ncol(rs1))]
  
  n             <- length(lst.rfModels)
  # the result matrix
  result        <- matrix(0,nrow = dim(mat.data)[1],ncol=n)
  #for each response variable
  for (i in 1:n){
    res        <- 0
    
    nrForests  <- length(lst.rfModels[[i]])
    #compute the sum of the random forest predictions 
    for(k in 1:nrForests){
      res <- res +  predict(lst.rfModels[[i]][[k]], mat.data) 
    }
    #compute then the mean
    res <- res / nrForests
    #store the result in coloumn i
    result[,i] <- res
  }
  # name the coloumns corresponding to the name of the response variables
  colnames(result) <- names(lst.rfModels)
  return(data.frame(result))
}

#------------------------------------------------------------------------
#FUNCTION GenerateModels prepocesses the data, wiping out NA values,
#and computes then a for each response variable for each model configuration
# - given by vec.trees and vec.models - 
#two types of errors, namely the total and the single point error.
#-------------------------------------------------------------------------
#ARGUMENTS:
#fname              The .csv filename
#mat.data           Alternatively a data matrix can passed; however 'fname' has to be set then explicitly to 'NULL'                   
#vec.models         A vector comprising numbers, 
#                   reflecting the number of random forests that shall be used
#                   for predicting a response variable
#vec.trees          A vector comprising the number of trees for a random forest
#vec.f              The response variables, i.e. an index vector corresponding to the coloumns of the data
#vec.p              The prediction variables, i.e. an index vector corresponding to the permittable data coloumns
#                   permittable means here, those data coloumns that are allowed to be included as predictors
#RETURNS:
#A named list consisting of five components is returned
#'Data'         The data without NA values
#'Errors'       A named list of errors (value), 
#               where each entry reflects one response variable,
#               and consists of a list of two matrices,
#               the first for the total error, called 'ErrorTotal'
#               and the second for the single-point error, called "ErrorPoint'
#'Models'       The 'vec.models' vector
#'Trees'        The 'vec.trees' vector
#'ResponseVar'  The adjusted set of response variables
#----------------------------------------------------------------------
#COMPATIBLE FUNCTION INPUT:
#DrawLines
#SelectModel
#-----------------------------------------------------------------------
GenerateModels <- function(fname, mat.data=NULL, vec.models, vec.trees, vec.f, vec.p){
  
  # vec.f = to_est_ind
  # vec.p = descr_ind
  # fname = NULL
  # mat.data = train
  # vec.models = vec.models
  # vec.trees = vec.trees
  
  # Load the data
  if(! is.null(fname) ){
    data       <- read.csv(get('fname'), header=TRUE)
  }else if (! is.null(mat.data)){
    data=mat.data
  }else{
    stop("No data was given.")
  }
  
  # Check the validity of predictor and response indices
  ValidityTest(vec.f,vec.p,dim(data)[2])
  
  # Update the index vector of response variables
  temp      <- vec.f
  vec.f     <- NULL
  for(i in 1:length(temp))
    vec.f <- c(vec.f, which(vec.p==temp[i]))
  
  # Check existence of and load the packacke for random forests
  PkgTest("randomForest")
  library(randomForest)
  
  
  # PREPROCESSING THE DATA #
  
  cdata      <- data[ ,vec.p]
  
  # Eliminating the NA values
  cdata      <- PredictNAvalues(cdata)
  
  # COMPUTING THE ERRORS FOR DIFFERENT MODELS
  errors     <- ErrorComputation(mat.data   = cdata, 
                                 vec.models = vec.models,
                                 vec.trees  = vec.trees,
                                 vec.f      = vec.f
  )
  
  return(list(Data=cdata, 
              Errors=errors, 
              Models=vec.models, 
              Trees=vec.trees,
              ResponseVar=vec.f))
  
}

#----------------------------------------------------------------------
#FUNCTION SelectModel selects and generates the best models 
#for the response variables according the weighted sum of error
#----------------------------------------------------------------------
#ARGUMENTS:
#lst.GenerateModels A list consisting of five components 
#    [[1]]          The data
#    [[2]]          A list of errors (value), 
#                   where each entry reflects one response variable,
#                   and consists of a list of two matrices,
#                   the first for the total error, called 'ErrorTotal'
#                   and the second for the single-point error, called "ErrorPoint'
#    [[3]]          A comprehensive vector of used models (vec.models) for generating the errors
#    [[4]]          A comprehensive vector of used number of trees (vec.trees) for generating the errors
#    [[5]]          The adjusted set of response variables
#vec.weights        A vector of weights, one entry for each response variable in 'lst.errors',
#                   where each entry has to be between 0 and 1, 
#                   reflecting the desired impact on the total error, 
#                   e.g. if the i-th weight is 0.3 then 30% of the first error will be taken
#                   into account and 70% of second - this correlates to a stronger reduction
#                   of the second error
#RETURNS:
#A list of best models,one for each response variable
#-----------------------------------------------------------------------
##COMPATIBLE FUNCTION INPUT:
#GenerateModels 
#-----------------------------------------------------------------------
SelectModel <- function(lst.GenerateModels, vec.weights){
  mat.data             <-   lst.GenerateModels[[1]]
  lst.errors           <-   lst.GenerateModels[[2]]
  vec.models           <-   lst.GenerateModels[[3]]
  vec.trees            <-   lst.GenerateModels[[4]]
  vec.f                <-   lst.GenerateModels[[5]]
  #lst.errors.normed    <-   NormError(lst.errors)
  lst.errors.normed    <-   NormError(lst.GenerateModels$Errors)
  lst.errors.weighted  <-   WeightError(vec.weights = vec.weights, 
                                        lst.error=lst.errors.normed)
  n                    <-   length(vec.f)
  
  # Divide the data into training and test data
  res        <- DivideSet(dim(mat.data)[1])
  
  rfModels   <- list()
  
  for (i in 1:n){
    posColRow    <- FindBestModel(lst.errors.weighted, i)
    
      if (n==1)
        {colnrs         <- NULL
      }else if (i==1)
        {colnrs         <- vec.f[2:n]
      }else if (i==n)
        {colnrs         <- vec.f[1:(n-1)]
      }else
        {colnrs         <- vec.f[c(1:(i-1),(i+1):n)]
        }
      #GENERATE RANDOM FOREST MODELS
      lst.response.RFs <- list()
      for (j in 1:vec.models[posColRow[1]]){
          # data.train.model   <- data.train[folds[,j],]
          lst.response.RFs[[j]] <- ModelRF( mat.data    = mat.data, 
                                            vec.colnrs  = colnrs, 
                                            nrtree      = vec.trees[posColRow[2]], 
                                            tofit       = vec.f[i]
                                            )
      }
      rfModels[[length(rfModels)+1]]    <- lst.response.RFs
  }
  
  names(rfModels)  <- names(lst.errors)
  
  return(rfModels)
}

#----------------------------------------------------------------------
#FUNCTION ShowWeighted plots the performance of different model configurations
#for one response variable of the weighted error
#----------------------------------------------------------------------
#ARGUMENTS:
#lst.data     A list (mostly an error list), 
#             where each entry represents a response variable
#             consisting of two error matrices, 
#             as returned by the GenerateModels or the ErrorComputation functions
#col          The index of the response variable
#
#RETURNS:
#A plot showing the performance of different model parameters
#----------------------------------------------------------------------
#COMPATIBLE FUNCTION OUTPUT:
#ErrorComputation
#GenerateModels
#WeightedError
#-----------------------------------------------------------------------
ShowWeighted <- function(lst.data, col){
  
  n      <- NCOL(lst.data[[col]])
  cols   <- rainbow(n)
  ymin <- min(lst.data[[col]])
  ymax <- max(lst.data[[col]])
  
  xaxis  <- rownames(lst.data[[col]])
  
  for(i in 2:n)
  {
    ymin <- min(c(ymin, lst.data[[col]]))
    ymax <- max(c(ymax, lst.data[[col]]))
  }
  
  plot(xaxis, lst.data[[col]][,1],
       ylim    = c(ymin,ymax),
       lwd     = 3, 
       col     = cols[1], 
       main    = paste(names(lst.data[[col]]), names(lst.data)[col]), 
       type    = "l", 
       ylab    = "RMSE",
       xlab    = "Nr of trees")
  
  for(i in 2:n)
    lines(xaxis,lst.data[[col]][,i], lwd=3,col=cols[i], type="l")
  
  legend("topright", inset=.05,title="Nr of models",legend=colnames(lst.data[[col]]),fill=cols)    
  
}





