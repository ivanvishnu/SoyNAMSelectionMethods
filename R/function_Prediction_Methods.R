### Prediction Methods


## Fn 9:
##### Function to Predict Phenotype Values given genotable in 012 format

PredictRRBLUPPhenoValues <- function(GenoTable,TrainModel){

  nMarkers <- 4289
  markerNames <-rep(0,nMarkers)
 
  markerNames <- paste("m",c(1:nMarkers),sep="")
  
  givenGenoTable<- GenoTable
  beta_markers <- unlist(TrainModel[1])
  names(beta_markers)<- markerNames


  y.hat <- givenGenoTable %*% beta_markers
  Avg <- unlist(TrainModel[2])

  pred.y <- y.hat+ Avg

  return(pred.y)

}

# testRRBLUPPredict<- PredictRRBLUPPhenoValues(genoTable_New2,model)


## Fn 11:


PredictRRBLUPPhenoValues_WGS_Jannink1 <- function(GenoTable,TrainModel,FreqVec){

  givenGenoTable<- GenoTable
  beta.markers <- unlist(TrainModel[1])
  # names(beta)<- markerNames
  Freq <- FreqVec

  weightsVec <- (1/sqrt(Freq))

  wt.beta <- beta.markers* weightsVec

  Neginf.indices<- which(wt.beta == -Inf)
  inf.indices <-  which(wt.beta == Inf)
  indices.unWtd <- c(Neginf.indices,inf.indices)

  if(length(indices.unWtd)!=0) {
	wt.beta[indices.unWtd] <- beta.markers[indices.unWtd]
  }

  y.hat <- givenGenoTable %*% wt.beta

  Avg <- unlist(TrainModel[2])

  pred.y <- y.hat+ Avg

  return(pred.y)

}


## Fn 3

PredictRRBLUPPhenoValues_WGS_DWModel<- function(GenoTable,TrainModel,FreqVec,alphaShape,betaShape,cycleNumber,timeHorizon){

	  givenGenoTable<- GenoTable
	  beta.markers <- unlist(TrainModel[1])
	  Freq <- FreqVec

	  alphaS <- alphaShape
	  betaS <- betaShape
	  N <- timeHorizon
	  t <- cycleNumber


	  alphaPar <- (alphaS+ (t*(1-alphaS)/N))
	  betaPar <- betaS

	  weightsVec <- (Freq^(alphaS+ (t*((1-alphaS)/N))))/ beta(alphaPar,betaPar)

	  wt.beta <- beta.markers* weightsVec

	  Neginf.indices<- which(wt.beta == -Inf)
	  inf.indices <-  which(wt.beta == Inf)
	  indices.unWtd <- c(Neginf.indices,inf.indices)

	  if(length(indices.unWtd)!=0) {
		wt.beta[indices.unWtd] <- beta.markers[indices.unWtd]
	  }

	  y.hat <- givenGenoTable %*% wt.beta

	  Avg <- unlist(TrainModel[2])

	  pred.y <- y.hat+ Avg

	  return(list(pred.y,alphaPar))

}





## Fn 40: 
##### Function to Predict Phenotype Values given genotable in 012 format

PredictBayesPhenoValues <- function(GenoTable,TrainModel){

  givenGenoTable <- GenoTable

  # Pred.pheno.valid <- (testGeno %*% predBayesA$ETA[[1]]$b)
  # Pred.Accuracy <- cor(Pred.pheno.valid,testPheno)


  # TrainModel in sol <- list(Pred.Accuracy,predBayesA,predBayesA$ETA[[1]]$b,train_MSE,test_MSE)

  beta.markers <- as.vector(TrainModel[[3]])
  # names(beta)<- markerNames


  y.hat <- givenGenoTable %*% beta.markers
  
  pred.y <- y.hat
  return(pred.y)

}

# predictBayesPhenoValues(testGeno,sol)

## Fn 41: 
##### Function to Predict Phenotype Values given genotable in 012 format

PredictBayesPhenoValues_Wght <- function(GenoTable,TrainModel,FreqVec){

  givenGenoTable <- GenoTable
  Freq <- FreqVec
  beta.markers <- as.vector(TrainModel[[3]])

 
  weightsVec <- (1/sqrt(Freq))
  wt.beta <- beta.markers* weightsVec

  Neginf.indices<- which(wt.beta == -Inf)
  inf.indices <-  which(wt.beta == Inf)
  indices.unWtd <- c(Neginf.indices,inf.indices)

  if(length(indices.unWtd)!=0) {
	wt.beta[indices.unWtd] <- beta.markers[indices.unWtd]
  }

  y.hat <- givenGenoTable %*% wt.beta

  pred.y <- y.hat
  return(pred.y)

}
