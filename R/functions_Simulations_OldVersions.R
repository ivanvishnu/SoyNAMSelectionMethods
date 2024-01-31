
## Simulations to perform island model selection

runSimulations20X_IslandModelSelection <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,migrationFrequency,Direction){



	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues
	  nFamilies <- 20
	  nSel_inFamily <- no_selected/20
	  nCrosses_inFamily <- nSel_inFamily
	  migFreq <- migrationFrequency
	  direction <- Direction

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- 100
	  nProgeny <-  nProgeny/nSel_inFamily


	  if(nSel_inFamily ==1){
	      nSel_inFamily <- 2
	   }

       nMarkers <- 4289

######################################################################################################

	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

######################################################################################################


	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList


### Get predicted geno and pheno values for cycle1 ###################################################

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
      newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

########################################################################################################


	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }


	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]

	  }



## Split Geno and Pheno values according to families  ############################################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()



	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }



### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)

	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)

	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()


######## Functions on Cycle 1 Geno Data #################################################################

######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   } else if (selectionOnSimulated==TRUE) {

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	    }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  attainedGenoValues_List[[nFamily]][1] <- max(topNGenoValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


########### Functions on Cycle 1 Pheno Data#############################################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- max(topNPhenoValues)
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
	}

			GSTable_afterMigration <- migrate2X(genoSelectionTable,1,1,direction)
			PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,1,direction)

			genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
			phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
			emigrantGroups <- GSTable_afterMigration[[2]]
			immigrantGroups<- GSTable_afterMigration[[3]]





###########################################################################

# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################

   	for( i in 2:nCycles){

	  	print(i)

		selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }



### Exchange selected geno data among pairs


	   selectedGenoData_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)


### Generate F5 RIL progeny for selected geno data
# Parent_Combn_indices <- getParentCombnIndices(no_selected)


	  Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nSel_inFamily,nMarkers,2,nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nSel_inFamily))

	  for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F2 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))


###################################################################################
		if(nSel_inFamily==2){
				nCrosses_inFamily <- 1
		}else{nCrosses_inFamily <- nSel_inFamily}

	    # nCrosses <- nSel_inFamily

		for( j in 1:(nCrosses_inFamily)){


		  if(j ==1 && nCrosses_inFamily==1){
			  Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=nCrosses_inFamily-1)) {
		       Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][1,,]
		  }

		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)

		  for(m in 1:nProgeny){
			Parent1<- Cycle_Progeny_F1[j,,,m]
			Parent2<- Cycle_Progeny_F1[j,,,m]

			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F2[j,,,m]  <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }



### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){

			for(nSel in 1:nSel_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format


		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nSel_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}




# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

###########################################################################################3

####################################################################################################################################################################################################

#nextGenGenoTable <- generate012GenoFormat(Cycle_Progeny_F5,no_selected)

### This section works with full data table

  		nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,cycle1_nProgeny,nFamilies,no_selected)

		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)



## Split Geno and Pheno values according to families  ########################

	    initIndex <- 1
	    finalIndex <- cycle1_nProgeny
	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()

	    for(nFamily in 1:20){

			genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
			phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
			genoValSimValues_List[[nFamily]] <- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues[initIndex:finalIndex]
			initIndex <- finalIndex+1
			finalIndex <- finalIndex+cycle1_nProgeny

	    }

### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
		}

       if(i==2){
			GSTable_afterMigration <- migrate2X(genoSelectionTable,1,2,direction)
			PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,2,direction)

			genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
			phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
			emigrantGroups <- GSTable_afterMigration[[2]]
			immigrantGroups<- GSTable_afterMigration[[3]]

		} else if(i>2 && i%%migFreq==0){

			genoSelectionTable_afterMigration <- migrate(genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
			phenoSelectionTable_afterMigration <- migrate(phenoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
		}else if(i>2 && i%%migFreq !=0){

			genoSelectionTable_afterMigration <- genoSelectionTable
			phenoSelectionTable_afterMigration <- phenoSelectionTable

		}




##############################################################################################

		for(nFamily in 1:20){



			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable_afterMigration[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable_afterMigration[[nFamily]][,1]

			selectedGenoIndividualIndices_List[[nFamily]]<- genoSelectionTable_afterMigration[[nFamily]][,2]
			selectedGenoSimValues <- genoSimValues_List[[nFamily]][as.vector(selectedGenoIndividualIndices_List[[nFamily]])]

	        attainedGenoValues_List[[nFamily]][i] <- max(selectedGenoSimValues)

			selectedPhenoIndividualIndices_List[[nFamily]] <- phenoSelectionTable_afterMigration[[nFamily]][,2]
			attainedPhenoValues_List[[nFamily]][i] <- max(PhenoVal_NX_N_3c_List[[nFamily]][,i])

		}

	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List)

      return(simResults_List)

}


################################################################################################
## migration frequency , direction - 1/2

runSimulations20X_IslandModelSelection_GM <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,migrationFrequency,Direction,NAM_LinkMap){


### Assign Variables #########################################################################

	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues
	  nFamilies <- 20
	  migFreq <- migrationFrequency
	  direction <- Direction

	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }
	  nCrosses_inFamily <- nSel_inFamily

	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny
	  nProgeny <- nProgeny/nSel_inFamily

	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  nMarkers <- 4289
	  NAM_LinkMap <- NAM_LinkMap_New
##################################################

      # i<-1
	  # k<-1
	  # noCrosses <- 20
	  # no_Progeny<- 100

      # h2 <- h2[i]
	  # noQTLs <- no_QTL[i]

	  # no_selected <- numberSelected
	  # no_selected <-20

	  # cycle1_nProgeny<-100
	  # nCrosses <- noCrosses
	  # nProgeny <- no_Progeny
	  # nCycles <- noCycles

	  # selectionOnGeno <- TRUE
	  # selectionOnSimulated <- FALSE
	  # nFamilies <- 20
	  # F5_Progeny_List<- F5_Progeny_ListReps[[k]][[i]]

# # #############################################################################

	  # PredictionModel_Pheno <- model_Sol_Pheno_List[[i]][[k]]
	  # PredictionModel_Geno <- model_Sol_Geno_List[[i]][[k]]

	  # cycle1_Geno_Data <- F5_Progeny_ListReps[[k]][[i]]


      # genoSimValues <-  genotypicValuesSimListReps[[k]][[i]]
	  # phenoSimValues <- phenotypeSimListReps[[k]][[i]]

	  # varE <- errorVarListReps[[k]][[i]]

# ###############################################################################################

### Get predicted geno and pheno values for cycle1

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
      newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

########  Split F5 Progeny into list based on families


	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }

######## split cycle1_Geno data in to list according to families

	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]

	  }



## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)
	  # percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)


	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)


	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
     for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }else if (selectionOnSimulated==TRUE) {

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  selectedGenoSimValues <- genoSimValues_List[[nFamily]][GenoTableIndices_topN]

	  attainedGenoValues_List[[nFamily]][1] <- max(selectedGenoSimValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


#########################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- max(PhenoVal_NX_N_3c_List[[nFamily]][,1])
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
	}

	GSTable_afterMigration <- migrate2X(genoSelectionTable,1,1,direction)
	PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,1,direction)

	genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
	phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
	emigrantGroups <- GSTable_afterMigration[[2]]
	immigrantGroups<- GSTable_afterMigration[[3]]





###########################################################################

# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################

    for( i in 2:nCycles){

	    # i<-2
		print(i)


		selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		# nSel_inFamily <- no_selected/20

		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }


### Exchange selected geno data among pairs


	   selectedGenoData_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)


### Generate F5 RIL progeny for selected geno data

	  Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nSel_inFamily,nMarkers,2,nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nSel_inFamily))

	  for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F2 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))


	###################################################################################
	    if(nSel_inFamily==2){
				nCrosses_inFamily <- 1
		}else{nCrosses_inFamily <- nSel_inFamily}

		for( j in 1:(nCrosses_inFamily)){

		  if(j ==1 && nCrosses_inFamily ==1){
			  Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=nCrosses_inFamily-1)) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
		  }

		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)

		  for(m in 1:nProgeny){

			  Parent1<- Cycle_Progeny_F1[j,,,m]
			  Parent2<- Cycle_Progeny_F1[j,,,m]

			  progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			  Cycle_Progeny_F2[j,,,m]  <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){

			for(nSel in 1:nSel_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nSel_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}

# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

###########################################################################################3

####################################################################################################################################################################################################

#nextGenGenoTable <- generate012GenoFormat(Cycle_Progeny_F5,no_selected)

### This section works with full data table

  		nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,cycle1_nProgeny,nFamilies,no_selected)

		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)


		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

###################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoValSimValues_List[[nFamily]] <- (genoValSimValues)[initIndex:finalIndex]
		phenoValSimValues_List[[nFamily]] <- (phenoValSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
		}

		if(i==2){
			GSTable_afterMigration <- migrate2X(genoSelectionTable,1,2,direction)
			PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,2,direction)

			genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
			phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
			emigrantGroups <- GSTable_afterMigration[[2]]
			immigrantGroups<- GSTable_afterMigration[[3]]

		}else if(i>2 && i%%migFreq==0){

			genoSelectionTable_afterMigration <- migrate(genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
			phenoSelectionTable_afterMigration <- migrate(phenoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
		}else if(i>2 && i%%migFreq !=0){

			genoSelectionTable_afterMigration <- genoSelectionTable
			phenoSelectionTable_afterMigration <- phenoSelectionTable

		}

#####################################################################################################################



## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

 ####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}

	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List)

      return(simResults_List)

}


#######################################################################################
## Four issues with the current settings
## No_of Progeny (bulk- 100/cross, discrete 10/cross- 10/family)
## Single-round robin design
## (first crossed to the rest)
## (migrate samples ingroup and outgroup every generation)
#######################################################################################
## Fn 26:

runSimulations20X_IslandModelSelection_GM_Mod <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,migrationFrequency,Direction,NAM_LinkMap){


### Assign Variables #########################################################################

  	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues


	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny


	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  nMarkers <- 4289

	  NAM_LinkMap_New <- NAM_LinkMap
############################################################
      # i<-1
	  # k<-1
	  # noCrosses <- 200
	  # no_Progeny<- 100

      # h2 <- h2[i]
	  # noQTLs <- no_QTL[i]

	  # # no_selected <- numberSelected
	  # no_selected <-200

	  # cycle1_nProgeny<-100
	  # nCrosses <- noCrosses
	  # nProgeny <- no_Progeny
	  # nCycles <- noCycles

	  # selectionOnGeno <- TRUE
	  # selectionOnSimulated <- TRUE
	  # nFamilies <- 20
	  # F5_Progeny_List<- F5_Progeny_ListReps[[k]][[i]]

# # #############################################################################

	  # PredictionModel_Pheno <- model_Sol_Pheno_List[[i]][[k]]
	  # PredictionModel_Geno <- model_Sol_Geno_List[[i]][[k]]

	  # cycle1_Geno_Data <- F5_Progeny_ListReps[[k]][[i]]


      # genoSimValues <-  genotypicValuesSimListReps[[k]][[i]]
	  # phenoSimValues <- phenotypeSimListReps[[k]][[i]]

	  # varE <- errorVarListReps[[k]][[i]]

	  # direction <- 1

################

	  nFamilies <- 20
	  migFreq <- migrationFrequency
	  direction <- Direction

	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }

      nCrosses_inFamily <- nSel_inFamily

	  nProgeny <- nProgeny/nSel_inFamily



#### /nSel_inFamily

###############################################################################################
## Get predicted geno and pheno values for cycle1
	  nCyc <-1

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	  newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

########  Split F5 Progeny into list based on families

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }

######## split cycle1_Geno data in to list according to families
######## Get percentHeterozygous & favorableAllele frequency
      if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[1])
	  }else if(selectionOnGeno==FALSE){ 
	    MarkerEffects <- unlist(PredictionModel_Pheno[1])
	  }
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:20){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[
			nFamily]],1,cycle1_nProgeny)

	  }
###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1

		}

## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)

	  #percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)


	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)


	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }else if (selectionOnSimulated==TRUE){

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  selectedGenoSimValues <- genoSimValues_List[[nFamily]][GenoTableIndices_topN]

	  attainedGenoValues_List[[nFamily]][1] <- max(selectedGenoSimValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


#########################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- (PhenoVal_NX_N_3c_List[[nFamily]][,1])
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
	}

	GSTable_afterMigration <- migrate2X(genoSelectionTable,1,1,direction)
	PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,1,direction)



	genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
	phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
	emigrantGroups <- (GSTable_afterMigration[[2]])
	immigrantGroups<- GSTable_afterMigration[[3]]


###########################################################################
# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################

    # nCycles <- 10

	for( i in 2:nCycles){

	    # i<-2
	    nCyc <- i
		print(i)

        selectedGenoData_List <- list()
		selectedGenoData_List <- list()


		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }


### Exchange selected geno data among pairs


	  selectedGenoData_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)


### Generate F5 RIL progeny for selected geno data
    if(nSel_inFamily==2){
			nCrosses_inFamily <- 1
	}else{nCrosses_inFamily <- nSel_inFamily}


	Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))),nFamilies)

	Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nCrosses_inFamily*nProgeny))),nFamilies)

	Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nCrosses_inFamily))

	for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nCrosses_inFamily,nMarkers,2,1))

		Cycle_Progeny_F2 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))


###################################################################################


		for( j in 1:(nCrosses_inFamily)){

		  if(j ==1 && nCrosses_inFamily ==1){
			  Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=(nCrosses_inFamily%/%2))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j>(nCrosses_inFamily%/%2))&&(j<=(nCrosses_inFamily-1))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][2,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][(j-(nCrosses_inFamily%/%2)+1),,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
		  }

		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  # for(m in 1:nProgeny){

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1

		  # }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){

			for(nSel in 1:nCrosses_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nCrosses_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}

# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

####################################################################################################################################################################################################
### This section works with full data table

  		# nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,(nSel_inFamily*cycle1_nProgeny),nFamilies,no_selected)


		nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nFamilies,nCrosses_inFamily*nProgeny)

		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

		for(nFamily in 1:20){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)

		}



		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)


### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	    initIndex <- 1
	    finalIndex <- nProgeny*nCrosses_inFamily


	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()

	    for(nFamily in 1:20){


		    genoValues.Family <- genoValues[initIndex:finalIndex]
			phenoValues.Family <- phenoValues[initIndex:finalIndex]
			genoValSimValues.Family<- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues.Family <- phenoValSimValues[initIndex:finalIndex]
			genoValues_List[[nFamily]] <- genoValues.Family
			phenoValues_List[[nFamily]] <- phenoValues.Family

			genoValSimValues_List[[nFamily]] <- genoValSimValues.Family
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues.Family



			initIndex <- finalIndex+1
			finalIndex <- finalIndex + (nProgeny*nCrosses_inFamily)

		}



### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)

		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)

		}


		if(i>=2 && i%%migFreq==0){

			genoSelectionTable_afterMigration <- migrateN_v1(genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
			phenoSelectionTable_afterMigration <- migrateN_v1(phenoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
		}else if(i>2 && i%%migFreq !=0){

			genoSelectionTable_afterMigration <- genoSelectionTable
			phenoSelectionTable_afterMigration <- phenoSelectionTable

		}

#####################################################################################################################
## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}
		
	

	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List)

	  return(simResults_List)

}

## Fn 27: 



runSimulations20X_IslandModelSelection_GM_Mod_PG <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,migrationFrequency,Direction,NAM_LinkMap){


### Assign Variables #########################################################################

  	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues


	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny


	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  nMarkers <- 4289

	  NAM_LinkMap_New <- NAM_LinkMap

############################################################
      # i<-1
	  # k<-1
	  # noCrosses <- 200
	  # no_Progeny<- 100

      # h2 <- h2[i]
	  # noQTLs <- no_QTL[i]

	  # # no_selected <- numberSelected
	  # no_selected <-200

	  # cycle1_nProgeny<-100
	  # nCrosses <- noCrosses
	  # nProgeny <- no_Progeny
	  # nCycles <- noCycles

	  # selectionOnGeno <- TRUE
	  # selectionOnSimulated <- TRUE
	  # nFamilies <- 20
	  # F5_Progeny_List<- F5_Progeny_ListReps[[k]][[i]]

# # #############################################################################

	  # PredictionModel_Pheno <- model_Sol_Pheno_List[[i]][[k]]
	  # PredictionModel_Geno <- model_Sol_Geno_List[[i]][[k]]

	  # cycle1_Geno_Data <- F5_Progeny_ListReps[[k]][[i]]


      # genoSimValues <-  genotypicValuesSimListReps[[k]][[i]]
	  # phenoSimValues <- phenotypeSimListReps[[k]][[i]]

	  # varE <- errorVarListReps[[k]][[i]]

	  # direction <- 1

################

	  nFamilies <- 20
	  migFreq <- migrationFrequency
	  direction <- Direction


	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }

      nCrosses_inFamily <- nSel_inFamily

	  nProgeny <- nProgeny/nSel_inFamily

	  percentHeterozygous_List<- rep(list(rep(0,20)),nFamilies)
	  percentHeterozygous_QTL_List<- rep(list(rep(0,20)),nFamilies)
	  FavorableAlleleFreq <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableAllele <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableQTLAlleleFreq <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableQTLAllele <- rep(list(rep(list(list()),20)),nFamilies)
	  percentHet <- rep(0,20)

	  QTL_GenoTable_List <- list()


	  Diff_Stats_List_Global <- list()
	  Diff_Stats_List_Local <- list()



#### /nSel_inFamily

###############################################################################################
## Get predicted geno and pheno values for cycle1
	  nCyc <-1

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	  percentHet[nCyc] <- getPercentHeterozygous(nextGenGenoTable)

      newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

########  Split F5 Progeny into list based on families

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }

######## split cycle1_Geno data in to list according to families
######## Get percentHeterozygous & favorableAllele frequency
      if(selectionOnGeno==TRUE){
		 MarkerEffects <- unlist(PredictionModel_Geno[1])
	   }else if (selectionOnGeno==FALSE){
		 MarkerEffects <- unlist(PredictionModel_Pheno[1])
	   }
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:20){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[nFamily]],1,cycle1_nProgeny)
			QTL_GenoTable_List[[nFamily]] <- getQTLTable(nextGenGenoTable_List[[nFamily]],no_QTL)

			Freq_Allele_List <-	getFreq_BasePoln(nextGenGenoTable_List[[nFamily]],MarkerEffects)
			FavorableAlleleFreq[[nFamily]][[nCyc]] <- Freq_Allele_List[[1]]
			FavorableAllele[[nFamily]][[nCyc]]  <- Freq_Allele_List[[2]]
			percentHeterozygous_List[[nFamily]][nCyc] <- getPercentHeterozygous(nextGenGenoTable_List[[nFamily]])


			Freq_QTLAllele_List <-	getFreq_BasePoln(QTL_GenoTable_List[[nFamily]],MarkerEffects)
			FavorableQTLAlleleFreq[[nFamily]][[nCyc]] <- Freq_QTLAllele_List[[1]]
			FavorableQTLAllele[[nFamily]][[nCyc]]  <- Freq_QTLAllele_List[[2]]
			percentHeterozygous_QTL_List[[nFamily]][nCyc] <- getPercentHeterozygous(QTL_GenoTable_List[[nFamily]])
	  }

###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1

		}

### get island PG parameters

		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")

		Diff_Stats <- diff_stats(nextGenGenoTable_GInd)
	    Diff_Stats_List_Local[[nCyc]] <- as.big.matrix(Diff_Stats[[1]])
	    Diff_Stats_List_Global[[nCyc]] <- Diff_Stats[[2]]

		rm(nextGenGenoTable_GInd)
		rm(Diff_Stats)
		gc()


## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)
	  # percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)


	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)


	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }else if (selectionOnSimulated==TRUE){

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  selectedGenoSimValues <- genoSimValues_List[[nFamily]][GenoTableIndices_topN]

	  attainedGenoValues_List[[nFamily]][1] <- max(selectedGenoSimValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


#########################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- (PhenoVal_NX_N_3c_List[[nFamily]][,1])
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
	}

	GSTable_afterMigration <- migrate2X(genoSelectionTable,1,1,direction)
	PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,1,direction)



	genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
	phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
	emigrantGroups <- (GSTable_afterMigration[[2]])
	immigrantGroups<- GSTable_afterMigration[[3]]


###########################################################################
# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################

    # nCycles <- 10

	for( i in 2:nCycles){

	    # i<-2
	    nCyc <- i
		print(i)

        selectedGenoData_List <- list()
		selectedGenoData_List <- list()


		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }


### Exchange selected geno data among pairs


	  selectedGenoData_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)


### Generate F5 RIL progeny for selected geno data
    if(nSel_inFamily==2){
			nCrosses_inFamily <- 1
	}else{nCrosses_inFamily <- nSel_inFamily}


	Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))),nFamilies)

	Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nCrosses_inFamily*nProgeny))),nFamilies)

	Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nCrosses_inFamily))

	for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nCrosses_inFamily,nMarkers,2,1))

		Cycle_Progeny_F2 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))


###################################################################################


		for( j in 1:(nCrosses_inFamily)){

		  if(j ==1 && nCrosses_inFamily ==1){
			  Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=(nCrosses_inFamily%/%2))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j>(nCrosses_inFamily%/%2))&&(j<=(nCrosses_inFamily-1))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][2,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][(j-(nCrosses_inFamily%/%2)+1),,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
		  }

		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  # for(m in 1:nProgeny){

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1

		  # }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){

			for(nSel in 1:nCrosses_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nCrosses_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}

# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

####################################################################################################################################################################################################
### This section works with full data table

  		# nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,(nSel_inFamily*cycle1_nProgeny),nFamilies,no_selected)


		nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nFamilies,nCrosses_inFamily*nProgeny)
		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)


		for(nFamily in 1:20){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)
			QTL_GenoTable_List[[nFamily]] <- getQTLTable(nextGenGenoTable_List[[nFamily]],no_QTL)

			FavorableAlleleFreq[[nFamily]][[nCyc]] <- getFreq(nextGenGenoTable_List[[nFamily]],FavorableAllele[[nFamily]][[1]])
			percentHeterozygous_List[[nFamily]][nCyc] <- getPercentHeterozygous(nextGenGenoTable_List[[nFamily]])

			FavorableQTLAlleleFreq[[nFamily]][[nCyc]] <- getFreq(QTL_GenoTable_List[[nFamily]],FavorableQTLAllele[[nFamily]][[1]])
			percentHeterozygous_QTL_List[[nFamily]][nCyc] <- getPercentHeterozygous(QTL_GenoTable_List[[nFamily]])

		}


		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)


### get island PG parameters
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")

		Diff_Stats <- diff_stats(nextGenGenoTable_GInd)
	    Diff_Stats_List_Local[[nCyc]] <- as.big.matrix(Diff_Stats[[1]])
	    Diff_Stats_List_Global[[nCyc]] <- Diff_Stats[[2]]

		rm(nextGenGenoTable_GInd)
		rm(Diff_Stats)
		gc()

### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	    initIndex <- 1
	    finalIndex <- nProgeny*nCrosses_inFamily


	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()
		# genoValues_List.fam <- list()

	    for(nFamily in 1:20){


		    genoValues.Family <- genoValues[initIndex:finalIndex]
			phenoValues.Family <- phenoValues[initIndex:finalIndex]
			genoValSimValues.Family<- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues.Family <- phenoValSimValues[initIndex:finalIndex]
			genoValues_List[[nFamily]] <- genoValues.Family
			phenoValues_List[[nFamily]] <- phenoValues.Family

			genoValSimValues_List[[nFamily]] <- genoValSimValues.Family
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues.Family



			initIndex <- finalIndex+1
			finalIndex <- finalIndex + (nProgeny*nCrosses_inFamily)

		}



### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)

		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)

		}


		if(i>=2 && i%%migFreq==0){

			genoSelectionTable_afterMigration <- migrateN_v1(genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
			phenoSelectionTable_afterMigration <- migrateN_v1(phenoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
		}else if(i>2 && i%%migFreq !=0){

			genoSelectionTable_afterMigration <- genoSelectionTable
			phenoSelectionTable_afterMigration <- phenoSelectionTable

		}

#####################################################################################################################

## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

 ####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}
	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List,QTL_GenoTable_List,FavorableAlleleFreq,
	  percentHeterozygous_List,FavorableQTLAlleleFreq,percentHeterozygous_QTL_List,Diff_Stats_List_Local,Diff_Stats_List_Global)

	  return(simResults_List)

}

###

## Fn 28:
#######################################################################################
### Select progeny for next gen from a larger set before applying selection


runSimulations20X_IslandModelSelection_GM_Mod_SelProgeny<- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,migrationFrequency,Direction,NAM_LinkMap){



### Assign Variables #########################################################################

  	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues


	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny


	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  nMarkers <- 4289
	  NAM_LinkMap_New <- NAM_LinkMap

############################################################

	  nFamilies <- 20
	  migFreq <- migrationFrequency

      direction <- Direction

	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }

      nCrosses_inFamily <- nSel_inFamily

	  nProgeny <- nProgeny/nSel_inFamily





#### /nSel_inFamily

###############################################################################################
## Get predicted geno and pheno values for cycle1
	  nCyc <-1

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	  newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

########  Split F5 Progeny into list based on families

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }

######## split cycle1_Geno data in to list according to families ########################
######## Get percentHeterozygous & favorableAllele frequency ############################
      if(selectionOnGeno==TRUE){ 
			MarkerEffects <- unlist(PredictionModel_Geno[1])
		} else if(selectionOnGeno==FALSE) {
			MarkerEffects <- unlist(PredictionModel_Pheno[1])
		}
		
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:20){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[nFamily]],1,cycle1_nProgeny)
		}



###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1

		}

## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)

     # percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)


	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)


	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }else if (selectionOnSimulated==TRUE){

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  selectedGenoSimValues <- genoSimValues_List[[nFamily]][GenoTableIndices_topN]

	  attainedGenoValues_List[[nFamily]][1] <- max(selectedGenoSimValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


#########################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- (PhenoVal_NX_N_3c_List[[nFamily]][,1])
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
	}

	GSTable_afterMigration <- migrate2X(genoSelectionTable,1,1,direction)
	PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,1,direction)

#############################################


	genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
	phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
	emigrantGroups <- (GSTable_afterMigration[[2]])
	immigrantGroups<- GSTable_afterMigration[[3]]

	rm(nextGenGenoTable)
	rm(newNextGenGenoTable)

###########################################################################
# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################


	for( i in 2:nCycles){

	    # i<-2
	    nCyc <- i
		print(i)

      	selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		# nSel_inFamily <- no_selected/20

		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }


### Exchange selected geno data among pairs


	  selectedGenoData_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)


### Generate F5 RIL progeny for selected geno data
    if(nSel_inFamily==2){
			nCrosses_inFamily <- 1
	}else{nCrosses_inFamily <- nSel_inFamily}


	Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))),nFamilies)

	Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nCrosses_inFamily*nProgeny))),nFamilies)

	Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nCrosses_inFamily))

	for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nCrosses_inFamily,nMarkers,2,1))

		Cycle_Progeny_F2 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))


###################################################################################


		for( j in 1:(nCrosses_inFamily)){

		  if(j ==1 && nCrosses_inFamily ==1){
			  Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=(nCrosses_inFamily%/%2))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j>(nCrosses_inFamily%/%2))&&(j<=(nCrosses_inFamily-1))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][2,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][(j-(nCrosses_inFamily%/%2)+1),,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
		  }

		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  # for(m in 1:nProgeny){

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1

		  # }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){

			for(nSel in 1:nCrosses_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nCrosses_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}

# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

####################################################################################################################################################################################################
### This section works with full data table

  		# nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,(nSel_inFamily*cycle1_nProgeny),nFamilies,no_selected)


		nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nFamilies,nCrosses_inFamily*nProgeny)
		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)


		for(nFamily in 1:20){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)


		}


		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)


### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	    initIndex <- 1
	    finalIndex <- nProgeny*nCrosses_inFamily


	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()
		genoValues_List.fam <- list()

	for(nFamily in 1:20){


		  genoValues.family <- genoValues[initIndex:finalIndex]
			phenoValues.family <- phenoValues[initIndex:finalIndex]
			genoValSimValues.family<- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues.family <- phenoValSimValues[initIndex:finalIndex]


			initIndex.Cross <-1
			finalIndex.Cross <- nProgeny

			genoValues.selected <- c()
			phenoValues.selected <- c()
			genoValSimValues.selected <- c()
			phenoValSimValues.selected <- c()



				genoValues.sorted <- sort(genoValues.family,decreasing=TRUE,index.return=TRUE)
				genoValues.selFamily <- genoValues.sorted[[1]][1:nProgeny]

				phenoValues.sorted <- sort(phenoValues.family,decreasing=TRUE,index.return=TRUE)
				phenoValues.selFamily <- phenoValues.sorted[[1]][1:nProgeny]

				genoValSimValues.selCross <- genoValSimValues.family[initIndex.Cross:finalIndex.Cross]
				genoValSimValues.sorted <- sort(genoValSimValues.family,decreasing=TRUE,index.return=TRUE)
				genoValSimValues.selFamily <- genoValSimValues.sorted[[1]][1:nProgeny]

				phenoValSimValues.selCross <- phenoValSimValues.family[initIndex.Cross:finalIndex.Cross]
				phenoValSimValues.sorted <- sort(phenoValSimValues.family,decreasing=TRUE,index.return=TRUE)
				phenoValSimValues.selFamily <- phenoValSimValues.sorted[[1]][1:nProgeny]


				genoValues_List[[nFamily]] <- genoValues.selFamily
				phenoValues_List[[nFamily]] <- phenoValues.selFamily

				genoValSimValues_List[[nFamily]] <- genoValSimValues.selFamily
				phenoValSimValues_List[[nFamily]] <- phenoValSimValues.selFamily

				initIndex <- finalIndex+1
				finalIndex <- finalIndex + (nProgeny*nCrosses_inFamily)
	}




### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)

		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)

		}


		if(i>=2 && i%%migFreq==0){

			genoSelectionTable_afterMigration <- migrateN_v1(genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
			phenoSelectionTable_afterMigration <- migrateN_v1(phenoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
		}else if(i>2 && i%%migFreq !=0){

			genoSelectionTable_afterMigration <- genoSelectionTable
			phenoSelectionTable_afterMigration <- phenoSelectionTable

		}

#####################################################################################################################

		## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

 ####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}

	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List)

      return(simResults_List)

}

###################################
## Fn 29: 

runSimulations20X_IslandModelSelection_GM_Mod_SelProgeny_PG <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,migrationFrequency,Direction,NAM_LinkMap){



### Assign Variables #########################################################################

  	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues


	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny


	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  nMarkers <-4289

	  NAM_LinkMap_New <- NAM_LinkMap

############################################################
      # i<-1
	  # k<-1
	  # noCrosses <- 200
	  # no_Progeny<- 100

      # h2 <- h2[i]
	  # noQTLs <- no_QTL[i]

	  # # no_selected <- numberSelected
	  # no_selected <-200

	  # cycle1_nProgeny<-100
	  # nCrosses <- noCrosses
	  # nProgeny <- no_Progeny
	  # nCycles <- noCycles

	  # selectionOnGeno <- TRUE
	  # selectionOnSimulated <- TRUE
	  # nFamilies <- 20
	  # F5_Progeny_List<- F5_Progeny_ListReps[[k]][[i]]

# # #############################################################################

	  # PredictionModel_Pheno <- model_Sol_Pheno_List[[i]][[k]]
	  # PredictionModel_Geno <- model_Sol_Geno_List[[i]][[k]]

	  # cycle1_Geno_Data <- F5_Progeny_ListReps[[k]][[i]]


      # genoSimValues <-  genotypicValuesSimListReps[[k]][[i]]
	  # phenoSimValues <- phenotypeSimListReps[[k]][[i]]

	  # varE <- errorVarListReps[[k]][[i]]

	  # direction <- 1

################

	  nFamilies <- 20
	  migFreq <- migrationFrequency
	  direction <- Direction

	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }

      nCrosses_inFamily <- nSel_inFamily

	  nProgeny <- nProgeny/nSel_inFamily


	  percentHeterozygous_List<- rep(list(rep(0,20)),nFamilies)
	  percentHeterozygous_QTL_List<- rep(list(rep(0,20)),nFamilies)
	  FavorableAlleleFreq <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableAllele <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableQTLAlleleFreq <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableQTLAllele <- rep(list(rep(list(list()),20)),nFamilies)
	  percentHet <- rep(0,20)

	  QTL_GenoTable_List <- list()

	  Diff_Stats_List_Global <- list()
	  Diff_Stats_List_Local <- list()





#### /nSel_inFamily

###############################################################################################
## Get predicted geno and pheno values for cycle1
	  nCyc <-1

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	  percentHet[nCyc] <- getPercentHeterozygous(nextGenGenoTable)

      newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

########  Split F5 Progeny into list based on families

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }

######## split cycle1_Geno data in to list according to families
######## Get percentHeterozygous & favorableAllele frequency
      if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[1])
	  }else if (selectionOnGeno==FALSE){
		MarkerEffects <- unlist(PredictionModel_Pheno[1])
	  }
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:20){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[nFamily]],1,cycle1_nProgeny)
			QTL_GenoTable_List[[nFamily]] <- as.big.matrix(getQTLTable(nextGenGenoTable_List[[nFamily]],no_QTL))

			Freq_Allele_List <-	getFreq_BasePoln(nextGenGenoTable_List[[nFamily]],MarkerEffects)
			FavorableAlleleFreq[[nFamily]][[nCyc]] <- Freq_Allele_List[[1]]
			FavorableAllele[[nFamily]][[nCyc]]  <- Freq_Allele_List[[2]]
			percentHeterozygous_List[[nFamily]][nCyc] <- getPercentHeterozygous(nextGenGenoTable_List[[nFamily]])


			Freq_QTLAllele_List <-	getFreq_BasePoln(QTL_GenoTable_List[[nFamily]],MarkerEffects)
			FavorableQTLAlleleFreq[[nFamily]][[nCyc]] <- Freq_QTLAllele_List[[1]]
			FavorableQTLAllele[[nFamily]][[nCyc]]  <- Freq_QTLAllele_List[[2]]
			percentHeterozygous_QTL_List[[nFamily]][nCyc] <- getPercentHeterozygous(QTL_GenoTable_List[[nFamily]])
	  }

###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1

		}

### get island PG parameters

		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		Diff_Stats <- diff_stats(nextGenGenoTable_GInd)
		Diff_Stats_List_Local[[nCyc]] <- as.big.matrix(Diff_Stats[[1]])
		Diff_Stats_List_Global[[nCyc]] <- Diff_Stats[[2]]

	    rm(nextGenGenoTable_AlleleFormat)
		rm(nextGenGenoTable_GInd)
		rm(Diff_Stats)
		gc()



## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)
	  # percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)


	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)


	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }else if (selectionOnSimulated==TRUE){

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  selectedGenoSimValues <- genoSimValues_List[[nFamily]][GenoTableIndices_topN]

	  attainedGenoValues_List[[nFamily]][1] <- max(selectedGenoSimValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


#########################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- (PhenoVal_NX_N_3c_List[[nFamily]][,1])
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
	}

	GSTable_afterMigration <- migrate2X(genoSelectionTable,1,1,direction)
	PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,1,direction)

#############################################


	genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
	phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
	emigrantGroups <- (GSTable_afterMigration[[2]])
	immigrantGroups<- GSTable_afterMigration[[3]]

	rm(nextGenGenoTable)
	rm(newNextGenGenoTable)

###########################################################################
# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################


	for( i in 2:nCycles){

	    #i<-2
	    nCyc <- i
		print(i)

      	selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		# nSel_inFamily <- no_selected/20

		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }


### Exchange selected geno data among pairs


	  selectedGenoData_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)


### Generate F5 RIL progeny for selected geno data
    if(nSel_inFamily==2){
			nCrosses_inFamily <- 1
	}else { nCrosses_inFamily <- nSel_inFamily}


	Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))),nFamilies)

	Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nCrosses_inFamily*nProgeny))),nFamilies)

	Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nCrosses_inFamily))

	for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nCrosses_inFamily,nMarkers,2,1))

		Cycle_Progeny_F2 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))


###################################################################################


		for( j in 1:(nCrosses_inFamily)){

		  if(j ==1 && nCrosses_inFamily ==1){
			  Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=(nCrosses_inFamily%/%2))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j>(nCrosses_inFamily%/%2))&&(j<=(nCrosses_inFamily-1))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][2,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][(j-(nCrosses_inFamily%/%2)+1),,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
		  }

		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  # for(m in 1:nProgeny){ 		  # }

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){

			for(nSel in 1:nCrosses_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nCrosses_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}

# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

####################################################################################################################################################################################################
### This section works with full data table



		nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nFamilies,nCrosses_inFamily*nProgeny)
		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)


		for(nFamily in 1:20){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)
			QTL_GenoTable_List[[nFamily]] <- getQTLTable(nextGenGenoTable_List[[nFamily]],no_QTL)

			FavorableAlleleFreq[[nFamily]][[nCyc]] <- getFreq(nextGenGenoTable_List[[nFamily]],FavorableAllele[[nFamily]][[1]])
			percentHeterozygous_List[[nFamily]][nCyc] <- getPercentHeterozygous(nextGenGenoTable_List[[nFamily]])

			FavorableQTLAlleleFreq[[nFamily]][[nCyc]] <- getFreq(QTL_GenoTable_List[[nFamily]],FavorableQTLAllele[[nFamily]][[1]])
			percentHeterozygous_QTL_List[[nFamily]][nCyc] <- getPercentHeterozygous(QTL_GenoTable_List[[nFamily]])

		}


		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)


### get island PG parameters


		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		Diff_Stats <- diff_stats(nextGenGenoTable_GInd)
		Diff_Stats_List_Local[[nCyc]] <- as.big.matrix(Diff_Stats[[1]])
		Diff_Stats_List_Global[[nCyc]] <- Diff_Stats[[2]]

	    rm(nextGenGenoTable_AlleleFormat)
		rm(nextGenGenoTable_GInd)
		rm(Diff_Stats)
		gc()


### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	    initIndex <- 1
	    finalIndex <- nProgeny*nCrosses_inFamily


	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()
		genoValues_List.fam <- list()

	for(nFamily in 1:20){


		    genoValues.family <- genoValues[initIndex:finalIndex]
			phenoValues.family <- phenoValues[initIndex:finalIndex]
			genoValSimValues.family<- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues.family <- phenoValSimValues[initIndex:finalIndex]


			initIndex.Cross <-1
			finalIndex.Cross <- nProgeny

			genoValues.selected <- c()
			phenoValues.selected <- c()
			genoValSimValues.selected <- c()
			phenoValSimValues.selected <- c()

				genoValues.sorted <- sort(genoValues.family,decreasing=TRUE,index.return=TRUE)
				genoValues.selFamily <- genoValues.sorted[[1]][1:nProgeny]

				phenoValues.sorted <- sort(phenoValues.family,decreasing=TRUE,index.return=TRUE)
				phenoValues.selFamily <- phenoValues.sorted[[1]][1:nProgeny]

				genoValSimValues.selCross <- genoValSimValues.family[initIndex.Cross:finalIndex.Cross]
				genoValSimValues.sorted <- sort(genoValSimValues.family,decreasing=TRUE,index.return=TRUE)
				genoValSimValues.selFamily <- genoValSimValues.sorted[[1]][1:nProgeny]

				phenoValSimValues.selCross <- phenoValSimValues.family[initIndex.Cross:finalIndex.Cross]
				phenoValSimValues.sorted <- sort(phenoValSimValues.family,decreasing=TRUE,index.return=TRUE)
				phenoValSimValues.selFamily <- phenoValSimValues.sorted[[1]][1:nProgeny]


				genoValues_List[[nFamily]] <- genoValues.selFamily
				phenoValues_List[[nFamily]] <- phenoValues.selFamily

				genoValSimValues_List[[nFamily]] <- genoValSimValues.selFamily
				phenoValSimValues_List[[nFamily]] <- phenoValSimValues.selFamily

				initIndex <- finalIndex+1
				finalIndex <- finalIndex + (nProgeny*nCrosses_inFamily)
	}



### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)

		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)

		}


		if(i>=2 && i%%migFreq==0){

			genoSelectionTable_afterMigration <- migrateN_v1(genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
			phenoSelectionTable_afterMigration <- migrateN_v1(phenoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
		}else if(i>2 && i%%migFreq !=0){

			genoSelectionTable_afterMigration <- genoSelectionTable
			phenoSelectionTable_afterMigration <- phenoSelectionTable

		}

#####################################################################################################################
## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

 ####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}

	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List,QTL_GenoTable_List,FavorableAlleleFreq,percentHeterozygous_List,FavorableQTLAlleleFreq,percentHeterozygous_QTL_List,Diff_Stats_List_Local,Diff_Stats_List_Global)

      return(simResults_List)

}


#######################################################################################
## Fn 30:
## Discrete &Island model with options for weighting (Jannink/DynamicWeighting) ("JN"/"DW"),updating (update frequency)

runSimulations20X_IslandSelection_GM_BayesB_Wght_SelProgeny <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,migrationFrequency,Direction,NAM_LinkMap,ModelRetrain,ModelType,UpdateFrequency,i,k,weighted,wghtMethod){

### Assign Variables #########################################################################

  	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  
## Criterion value  
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues

## Prediction Models

	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

## Geno data 
	  
	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny

## Geno & Pheno values

	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  nMarkers <- 4289
	  NAM_LinkMap_New <- NAM_LinkMap

	  Weighted <- weighted
	  WghtMethod <- wghtMethod
############################################################

	  nFamilies <- 20
	  migFreq <- migrationFrequency
       
      direction <- Direction

	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }

      nCrosses_inFamily <- nSel_inFamily

	  # nProgeny <- nProgeny/nSel_inFamily
	  nProgeny <- nSel_inFamily

      modelType <- ModelType
      modelRetrain <- ModelRetrain
      updateFrequency <- UpdateFrequency

 
	  condition <- i
	  Rep <- k  
	  

#### /nSel_inFamily

###############################################################################################
## Get predicted geno and pheno values for cycle1
####### Predict geno and pheno values ###################################################################

	nCyc <-1

	cycle1GenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable <- generateWeightedGenoTable(cycle1GenoTable,no_QTL)

    nIndividuals <- nrow(cycle1GenoTable) 
	
    if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[3])
	}else if(selectionOnGeno==FALSE){
	    MarkerEffects <- unlist(PredictionModel_Pheno[3])
	}

    Freq_List <- getFreq_BasePoln(cycle1GenoTable,MarkerEffects)
	Freq <- Freq_List[[1]]
	FavAllele <- Freq_List[[2]]


	if(Weighted==TRUE){
		genoValues  <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Geno,Freq)
		phenoValues <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Pheno,Freq)
	}else if(Weighted==FALSE){
		genoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Geno)
		phenoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Pheno)

	}

	 
########  Split F5 Progeny into list based on families

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }

######## split cycle1_Geno data in to list according to families ########################
######## Get percentHeterozygous & favorableAllele frequency ############################
      if(selectionOnGeno==TRUE){ 
			MarkerEffects <- unlist(PredictionModel_Geno[3])
	  } else if(selectionOnGeno==FALSE) {
			MarkerEffects <- unlist(PredictionModel_Pheno[3])
	  }
		
	
		
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:20){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[nFamily]],1,cycle1_nProgeny)
		}



###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1

		}

## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)

     # percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)


	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)


	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }else if (selectionOnSimulated==TRUE){

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  selectedGenoSimValues <- genoSimValues_List[[nFamily]][GenoTableIndices_topN]

	  attainedGenoValues_List[[nFamily]][1] <- max(selectedGenoSimValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


#########################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- (PhenoVal_NX_N_3c_List[[nFamily]][,1])
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
	}

	
		GSTable_afterMigration <- migrate2X(genoSelectionTable,1,1,direction)
		PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,1,direction)
    
#############################################


	genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
	phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
	emigrantGroups <- (GSTable_afterMigration[[2]])
	immigrantGroups<- GSTable_afterMigration[[3]]

	rm(nextGenGenoTable)
	rm(newNextGenGenoTable)

###########################################################################
# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################


	for( i in 2:nCycles){

	    # i<-2
	    nCyc <- i
		print(i)

      	selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		# nSel_inFamily <- no_selected/20

		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }


### Exchange selected geno data among pairs


	  selectedGenoData_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)


### Generate F5 RIL progeny for selected geno data
    if(nSel_inFamily==2){
			nCrosses_inFamily <- 1
	}else{nCrosses_inFamily <- nSel_inFamily}


	Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))),nFamilies)

	Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nCrosses_inFamily*nProgeny))),nFamilies)

	Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nCrosses_inFamily))

	for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nCrosses_inFamily,nMarkers,2,1))

		Cycle_Progeny_F2 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))


###################################################################################


		for( j in 1:(nCrosses_inFamily)){

		  if(j ==1 && nCrosses_inFamily ==1){
			  Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=(nCrosses_inFamily%/%2))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j>(nCrosses_inFamily%/%2))&&(j<=(nCrosses_inFamily-1))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][2,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][(j-(nCrosses_inFamily%/%2)+1),,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
		  }

		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  # for(m in 1:nProgeny){

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1

		  # }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){

			for(nSel in 1:nCrosses_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nCrosses_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}

# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

####################################################################################################################################################################################################
### This section works with full data table

  		# nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,(nSel_inFamily*cycle1_nProgeny),nFamilies,no_selected)


		nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nFamilies,nCrosses_inFamily*nProgeny)
		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)


		for(nFamily in 1:20){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)


		}


		# genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		# phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		# newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
		
		
		
		
########################

    # nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,no_selected,nProgeny)

	Freq <- getFreq(nextGenGenoTable,FavAllele)

### Simulate genotypic and phenotypic values

    genoSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
    phenoSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

    newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)


############################################################################################### 
### Prediction Model Update 	
### Extract training and test data for model build

	IndividualNames <-rep(0,(nCrosses*nProgeny))

    for(nInd in 1:(nCrosses*nProgeny)){

        IndividualNames[nInd]<- paste("Ind",nInd,sep="")

    }

##################################################################################################
	names(phenoSimValues)<- IndividualNames
    Mean_Fixed<- rep(1,nIndividuals)

    phenotypeSimTable<- cbind(phenoSimValues,Mean_Fixed)

#### Table for Simulated Genotypic Values ####################################


    names(genoSimValues)<-IndividualNames
    Mean_Fixed<- rep(1,nIndividuals)

    genotypicValuesSimTable <- cbind(genoSimValues,Mean_Fixed)

#####################

	if(modelRetrain==TRUE){

		indices<-c(1:nIndividuals)

		trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))
		testIndices<- indices[which(!indices %in% trainIndices)]

	#### Training Set table for Genotype and Simulated Phenotype ###################

		# trainGenoTable <- genoTable_Mod[trainIndices,]
		trainGenoNewTable <- newNextGenGenoTable[trainIndices,]
		trainSimPhenoTable <- phenotypeSimTable[trainIndices,]
		trainSimGenoTable <-  genotypicValuesSimTable[trainIndices,]

	#### Validation Set table for Genotype and Simulated Phenotype ###################

		# testGenoTable<- genoTable_Mod[testIndices,]
		testGenoNewTable <- newNextGenGenoTable[testIndices,]
		testSimPhenoTable <- phenotypeSimTable[testIndices,]
		testSimGenoTable <-  genotypicValuesSimTable[testIndices,]

############### Changed RRBLUP to Bayes

	   if(sd(trainSimGenoTable[,1])!=0 && sd(trainSimPhenoTable[,1])!=0 && i%%updateFrequency==0){
	  
	    if(modelType=="BayesB"){ 
	   		PredictionModel_Geno <- (buildBayesBModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable))
			PredictionModel_Pheno <-(buildBayesBModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))
		}else if(modelType=="BL"){ 
		
			PredictionModel_Geno <- (buildBLModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable))
			PredictionModel_Pheno <-(buildBLModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))
		}
	  }

	}

	gc()

#################get predicted geno and pheno values with JM and unweighted GP models
### Get Predicted geno & pheno values for nextGenGenoData	
		
		if(Weighted==TRUE){
			genoValues  <- PredictBayesPhenoValues_Wght(newNextGenGenoTable,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newNextGenGenoTable,PredictionModel_Pheno,Freq)
	    }else if(Weighted==FALSE){
			genoValues <- PredictBayesPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
			phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

	    }
		
### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	    initIndex <- 1
	    finalIndex <- nProgeny*nCrosses_inFamily


	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoSimValues_List <- list()
	    phenoSimValues_List <- list()
		genoValues_List.fam <- list()

	    for(nFamily in 1:20){


		    genoValues.family <- genoValues[initIndex:finalIndex]
			phenoValues.family <- phenoValues[initIndex:finalIndex]
			genoSimValues.family<- genoSimValues[initIndex:finalIndex]
			phenoSimValues.family <- phenoSimValues[initIndex:finalIndex]


			initIndex.Cross <-1
			finalIndex.Cross <- nProgeny

			genoValues.selected <- c()
			phenoValues.selected <- c()
			genoSimValues.selected <- c()
			phenoSimValues.selected <- c()



				genoValues.sorted <- sort(genoValues.family,decreasing=TRUE,index.return=TRUE)
				genoValues.selFamily <- genoValues.sorted[[1]][1:nProgeny]

				phenoValues.sorted <- sort(phenoValues.family,decreasing=TRUE,index.return=TRUE)
				phenoValues.selFamily <- phenoValues.sorted[[1]][1:nProgeny]

				genoSimValues.selCross <- genoSimValues.family[initIndex.Cross:finalIndex.Cross]
				genoSimValues.sorted <- sort(genoSimValues.family,decreasing=TRUE,index.return=TRUE)
				genoSimValues.selFamily <- genoSimValues.sorted[[1]][1:nProgeny]

				phenoSimValues.selCross <- phenoSimValues.family[initIndex.Cross:finalIndex.Cross]
				phenoSimValues.sorted <- sort(phenoSimValues.family,decreasing=TRUE,index.return=TRUE)
				phenoSimValues.selFamily <- phenoSimValues.sorted[[1]][1:nProgeny]


				genoValues_List[[nFamily]] <- genoValues.selFamily
				phenoValues_List[[nFamily]] <- phenoValues.selFamily

				genoSimValues_List[[nFamily]] <- genoSimValues.selFamily
				phenoSimValues_List[[nFamily]] <- phenoSimValues.selFamily

				initIndex <- finalIndex+1
				finalIndex <- finalIndex + (nProgeny*nCrosses_inFamily)
	}




### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)

		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)

		}


		if(i>=2 && i%%migFreq==0){

			genoSelectionTable_afterMigration <- migrateN_v1(genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
			phenoSelectionTable_afterMigration <- migrateN_v1(phenoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
			
			
			
		}else if(i>2 && i%%migFreq !=0){

			genoSelectionTable_afterMigration <- genoSelectionTable
			phenoSelectionTable_afterMigration <- phenoSelectionTable

		}

#####################################################################################################################

		## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

 ####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}

	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List)

      return(simResults_List)

}

###################################
## Fn 31: old version

 runSimulations20X_IslandSelection_GM_BayesB_Wght <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,migrationFrequency,Direction,NAM_LinkMap,ModelRetrain,ModelType,UpdateFrequency,i,k,weighted,wghtMethod,BD){
	 
	  library(adegenet)
	  library(mmod)
	  
### Assign Variables #########################################################################

  	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues


	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny


	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  nMarkers <- 4289

	  NAM_LinkMap_New <- NAM_LinkMap
	  
	  Weighted <- weighted
	  WghtMethod <- wghtMethod
	  
	  BreedMethod <- BD
	
################

	  nFamilies <- 20
	  migFreq <- migrationFrequency
	  direction <- Direction

	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }

      nCrosses_inFamily <- nSel_inFamily

	  nProgeny <- nProgeny/nSel_inFamily


	  modelType <- ModelType
      modelRetrain <- ModelRetrain
      updateFrequency <- UpdateFrequency

 
	  condition <- i
	  Rep <- k  


#### /nSel_inFamily

###############################################################################################
## Get predicted geno and pheno values for cycle1
	  nCyc <-1

	 cycle1GenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	 newcycle1GenoTable <- generateWeightedGenoTable(cycle1GenoTable,no_QTL)
	 
	if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[3])
	}else if(selectionOnGeno==FALSE){
	    MarkerEffects <- unlist(PredictionModel_Pheno[3])
	}

     
     Freq_List <- getFreq_BasePoln(cycle1GenoTable,MarkerEffects)
	 Freq <- Freq_List[[1]]
	 FavAllele <- Freq_List[[2]]


	 if(Weighted==TRUE){
		genoValues  <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Geno,Freq)
		phenoValues <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Pheno,Freq)
	 }else if(Weighted==FALSE){
		genoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Geno)
		phenoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Pheno)
	 }

    
########  Split F5 Progeny into list based on families

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }

######## split cycle1_Geno data in to list according to families
######## Get percentHeterozygous & favorableAllele frequency
      if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[1])
	  }else if(selectionOnGeno==FALSE){ 
	    MarkerEffects <- unlist(PredictionModel_Pheno[1])
	  }
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:20){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[
			nFamily]],1,cycle1_nProgeny)

	  }
###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1

		}

## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)

	  #percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)


	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)


	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }else if (selectionOnSimulated==TRUE){

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  selectedGenoSimValues <- genoSimValues_List[[nFamily]][GenoTableIndices_topN]

	  attainedGenoValues_List[[nFamily]][1] <- max(selectedGenoSimValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


#########################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- (PhenoVal_NX_N_3c_List[[nFamily]][,1])
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
	}

	GSTable_afterMigration <- migrate2X(genoSelectionTable,1,1,direction)
	PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,1,direction)



	genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
	phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
	emigrantGroups <- (GSTable_afterMigration[[2]])
	immigrantGroups<- GSTable_afterMigration[[3]]


###########################################################################
# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################

    # nCycles <- 10

	for( i in 2:nCycles){

	    # i<-2
	    nCyc <- i
		print(i)

        selectedGenoData_List <- list()
		selectedGenoData_List <- list()


		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }


### Exchange selected geno data among pairs


	  selectedGenoData_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)


### Generate F5 RIL progeny for selected geno data
    if(nSel_inFamily==2){
			nCrosses_inFamily <- 1
	}else{nCrosses_inFamily <- nSel_inFamily}


	Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))),nFamilies)

	Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nCrosses_inFamily*nProgeny))),nFamilies)

	Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nCrosses_inFamily))

	for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nCrosses_inFamily,nMarkers,2,1))

		Cycle_Progeny_F2 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))


###################################################################################


		for( j in 1:(nCrosses_inFamily)){

		  if(j ==1 && nCrosses_inFamily ==1){
			  Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=(nCrosses_inFamily%/%2))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j>(nCrosses_inFamily%/%2))&&(j<=(nCrosses_inFamily-1))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][2,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][(j-(nCrosses_inFamily%/%2)+1),,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
		  }

		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  # for(m in 1:nProgeny){

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1

		  # }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){

			for(nSel in 1:nCrosses_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nCrosses_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}

# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

####################################################################################################################################################################################################
### This section works with full data table

  		# nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,(nSel_inFamily*cycle1_nProgeny),nFamilies,no_selected)


		nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nFamilies,nCrosses_inFamily*nProgeny)

		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

		for(nFamily in 1:20){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)

		}



		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)
		
		

		
############################################################################################### 
### Prediction Model Update 	
### Extract training and test data for model build

	IndividualNames <-rep(0,(nCrosses*nProgeny))

    for(nInd in 1:(nCrosses*nProgeny)){

        IndividualNames[nInd]<- paste("Ind",nInd,sep="")

    }

##################################################################################################
	names(phenoSimValues)<- IndividualNames
    Mean_Fixed<- rep(1,nIndividuals)

    phenotypeSimTable<- cbind(phenoSimValues,Mean_Fixed)

#### Table for Simulated Genotypic Values ####################################


    names(genoSimValues)<-IndividualNames
    Mean_Fixed<- rep(1,nIndividuals)

    genotypicValuesSimTable <- cbind(genoSimValues,Mean_Fixed)

#####################

	if(modelRetrain==TRUE){

		indices<-c(1:nIndividuals)

		trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))
		testIndices<- indices[which(!indices %in% trainIndices)]

	#### Training Set table for Genotype and Simulated Phenotype ###################

		# trainGenoTable <- genoTable_Mod[trainIndices,]
		trainGenoNewTable <- newNextGenGenoTable[trainIndices,]
		trainSimPhenoTable <- phenotypeSimTable[trainIndices,]
		trainSimGenoTable <-  genotypicValuesSimTable[trainIndices,]

	#### Validation Set table for Genotype and Simulated Phenotype ###################

		# testGenoTable<- genoTable_Mod[testIndices,]
		testGenoNewTable <- newNextGenGenoTable[testIndices,]
		testSimPhenoTable <- phenotypeSimTable[testIndices,]
		testSimGenoTable <-  genotypicValuesSimTable[testIndices,]

############### Changed RRBLUP to Bayes

	   if(sd(trainSimGenoTable[,1])!=0 && sd(trainSimPhenoTable[,1])!=0 && i%%updateFrequency==0){
	  
	    if(modelType=="BayesB"){ 
	   		PredictionModel_Geno <- (buildBayesBModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable))
			PredictionModel_Pheno <-(buildBayesBModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))
		}else if(modelType=="BL"){ 
		
			PredictionModel_Geno <- (buildBLModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable))
			PredictionModel_Pheno <-(buildBLModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))
		}
	  }

	}

	gc()

#################get predicted geno and pheno values with JM and unweighted GP models
### Get Predicted geno & pheno values for nextGenGenoData	
		
		if(Weighted==TRUE){
			genoValues  <- PredictBayesPhenoValues_Wght(newNextGenGenoTable,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newNextGenGenoTable,PredictionModel_Pheno,Freq)
	    }else if(Weighted==FALSE){
			genoValues <- PredictBayesPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
			phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

	    }
		

### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	    initIndex <- 1
	    finalIndex <- nProgeny*nCrosses_inFamily


	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()

	    for(nFamily in 1:20){


		    genoValues.Family <- genoValues[initIndex:finalIndex]
			phenoValues.Family <- phenoValues[initIndex:finalIndex]
			genoValSimValues.Family<- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues.Family <- phenoValSimValues[initIndex:finalIndex]
			genoValues_List[[nFamily]] <- genoValues.Family
			phenoValues_List[[nFamily]] <- phenoValues.Family

			genoValSimValues_List[[nFamily]] <- genoValSimValues.Family
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues.Family



			initIndex <- finalIndex+1
			finalIndex <- finalIndex + (nProgeny*nCrosses_inFamily)

		}



### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)

		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)

		}


		if(i>=2 && i%%migFreq==0){

			genoSelectionTable_afterMigration <- migrateN_v1(genoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
			phenoSelectionTable_afterMigration <- migrateN_v1(phenoSelectionTable,1,emigrantGroups,immigrantGroups,direction)
		}else if(i>2 && i%%migFreq !=0){

			genoSelectionTable_afterMigration <- genoSelectionTable
			phenoSelectionTable_afterMigration <- phenoSelectionTable

		}

#####################################################################################################################
## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}
		
	

	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List)

	  return(simResults_List)

}


#####################################################################################################################################

#####################
## Fn42:



runSimulations20X_IslandSelection_BD_BayesB_Wght <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateType,UpdateFrequency,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoTable,trainSimGenoValTable,trainSimPhenoValTable){


	options(warn=-1)
	
    library(adegenet)
	library(mmod)
	library(bigmemory)
	library(factoextra)
	library(LTSGenoIM)

### Assign Variables #########################################################################

  	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues


	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny


	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  nMarkers <- 4289

	  NAM_LinkMap_New <- NAM_LinkMap
	  
	  Weighted <- weighted
	  WghtMethod <- wghtMethod
	  
	  BD <- BreedDesign
	    
	  Policy <- policy	
      migPolicyFreq <- policyChangeFrequency
	  
############################################################
 
	  migFreq <- migrationFrequency
	  direction <- Direction

	  nFamilies <- 20
	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }

      nCrosses_inFamily <- nSel_inFamily

	  nProgeny <- nProgeny/nSel_inFamily


	  modelType <- ModelType
      modelRetrain <- ModelRetrain
	  retrainFrequency <- RetrainFrequency
	  
	  modelUpdate <- ModelUpdate
      updateFrequency <- UpdateFrequency
      updateType <- UpdateType
 
	  condition <- i
	  Rep <- k  

      AlleleConversionTable_Combined <-  alleleConversionTable_Combined
	  
	  Migration_List <- list()
	  PCA_List <- list()
	  j<-1

	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
	  nCyc <-1
	  trainGenoNewTable_List[[nCyc]] <-  as.big.matrix(trainGenoNewTable)
      trainSimGenoValTable_List[[nCyc]]<- as.big.matrix(trainSimGenoValTable)
	  trainSimPhenoTable_List[[nCyc]] <- as.big.matrix(trainSimPhenoTable)
		
#### /nSel_inFamily
###############################################################################################
## Get predicted geno and pheno values for cycle1 for combined population set

	cycle1GenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable <- generateWeightedGenoTable(cycle1GenoTable,no_QTL)
	nIndividuals <- nrow(cycle1GenoTable)
	 
	if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[3])
	}else if(selectionOnGeno==FALSE){
	    MarkerEffects <- unlist(PredictionModel_Pheno[3])
	}

     
     Freq_List <- getFreq_BasePoln(cycle1GenoTable,MarkerEffects)
	 Freq <- Freq_List[[1]]
	 FavAllele <- Freq_List[[2]]


	 if(Weighted==TRUE){
		genoValues  <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Geno,Freq)
		phenoValues <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Pheno,Freq)
	 }else if(Weighted==FALSE){
		genoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Geno)
		phenoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Pheno)
	 }

    
########  Split F5 Progeny into list based on families

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:nFamilies){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }

######## split cycle1_Geno data in to list according to families
######## Get percentHeterozygous & favorableAllele frequency
     
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:nFamilies){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[
			nFamily]],1,cycle1_nProgeny)

	  }


## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:nFamilies){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	 

	  #percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)


	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)


	  selectedGenoIndividualIndices2_List <- list()
	
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:nFamilies){


        if(selectionOnSimulated == FALSE){
		 
          if(selectionOnGeno==TRUE){
			  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
			  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
			  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]
			  
			  sortedPhenoValues <- phenoValues_List[[nFamily]][GenoTableIndices_topN]
			  PhenoTableIndices_topN <-  GenoTableIndices_topN
			  topNPhenoValues <- phenoValues_List[[nFamily]][GenoTableIndices_topN]
			
			}else if(selectionOnGeno==FALSE){
		  	  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
			  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
			  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]
			  
			  sortedGenoValues <- genoValues_List[[nFamily]][PhenoTableIndices_topN]
			  GenoTableIndices_topN <- PhenoTableIndices_topN 
			  topNGenoValues <- genoValues_List[[nFamily]][PhenoTableIndices_topN]
	        }
		}
		if (selectionOnSimulated==TRUE){

			if(selectionOnGeno==TRUE){
				sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
				topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
				GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]
				
				sortedPhenoValues <- phenoSimValues_List[[nFamily]][GenoTableIndices_topN]
				PhenoTableIndices_topN <-  GenoTableIndices_topN
				topNPhenoValues <- phenoSimValues_List[[nFamily]][GenoTableIndices_topN]
				
			
			}else if(selectionOnGeno==FALSE){

				sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
				topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
				PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]
				
				GenoTableIndices_topN <- PhenoTableIndices_topN 
				sortedGenoValues <- genoSimValues_List[[nFamily]][PhenoTableIndices_topN]
			    topNGenoValues <- genoSimValues_List[[nFamily]][PhenoTableIndices_topN]
			}
	    }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  selectedGenoSimValues <- genoSimValues_List[[nFamily]][GenoTableIndices_topN]
	  attainedGenoValues_List[[nFamily]][1] <- max(selectedGenoSimValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN

    
#########################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues
	  
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
	}
	
	
	migrationGroups <- migrationPolicy(cycle1GenoTable,genoSelectionTable,Policy,nFamilies)

	
	emigrantGroups <- migrationGroups[[1]]
	immigrantGroups <- migrationGroups[[2]]
	
###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(pop in 1:20){

			final <- pop*100
			populations[init:final]<- rep(pop,100)

			init <- final+1

		}

### get island PG parameters

		cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable,AlleleConversionTable_Combined)
 		cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
		tol <- 1e-07
		cycle1Gen_df <- as.matrix(cycle1GenoTable)
		cycle1Gen_df <- tcrossprod(cycle1Gen_df,cycle1Gen_df) 
		eig1 <- eigen(cycle1Gen_df,symmetric=TRUE) 
		eig <- eig1$values
		rank_X <- sum((eig/eig[1]) > tol)
		NF_Vec <- as.numeric(which(cumsum(eig/sum(eig)) >= 0.99)) 
		NF <- NF_Vec[1]
			
		if(NF == dim(cycle1Gen_df)[1]) { 
			
				NF = rank_X 
		}else if (NF==1){
				NF <-rank_X
		}else if( ! exists("NF")) {
				NF <-rank_X 
		}else if(rank_X== 1){
				NF <- 2
		}
		
			      
		scaledGInd <- scaleGen(cycle1GenoTable_GInd,scale=FALSE)
        
		tryCatch(pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF), error=function(e) paste("Error dudiPCA Cycle No",i,"-",e))
				
		   		
		filename <- paste("Sclass_Plot","-",condition,"-",Rep,"-",nCyc,".png",sep="")
				
		png(filename,width=1024,height=768,pointsize=20)
		tryCatch(s.class(pca.GInd$li, fac=pop(cycle1GenoTable_GInd),col=funky(20),addaxes=TRUE),error=function(e) print(paste("Error Sclass Cycle No",i,"-",e)))
				
								
		tryCatch(title(paste("Cycle - ",i,sep="")),error=function(e) print(paste("Error Title Cycle No",i,"-",e)))
		dev.off()
		tryCatch(system(paste("gzip -f ",filename,sep="")),error=function(e) print(paste("Error gzip Cycle No",i,"-",e)))
				
		tryCatch(PCA_Li_Table <- as.big.matrix(as.matrix(pca.GInd$li)),error=function(e) print(paste("Error BigMatrix Cycle No",i,e)))
		tryCatch(PCA_C1_Table <- as.big.matrix(as.matrix(pca.GInd$c1)),error=function(e) print(paste("Error BigMatrix Cycle No",i,e)))
		tryCatch(PCA_List[[j]] <- list(PCA_Li_Table,PCA_C1_Table),error=function(e)print(paste("Error PCA_List Cycle No",i,e)))
		
		# pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF)
	    # filename <- paste("Sclass_Plot","-",condition,"-",Rep,"-",nCyc,".png",sep="")
		# png(filename,width=1024,height=768,pointsize=20)
		# s.class(pca.GInd$li, fac=pop(cycle1GenoTable_GInd),col=funky(20),addaxes=TRUE)
		# title(paste("Cycle - ",i,sep=""))
		# dev.off()
		# PCA_Li_Table <- as.big.matrix(as.matrix(pca.GInd$li))
		# PCA_C1_Table <- as.big.matrix(as.matrix(pca.GInd$c1)) 
		# PCA_List[[j]] <- list(PCA_Li_Table,PCA_C1_Table)
		# system(paste("gzip ",filename,sep=""))
		
		rm(cycle1GenoTable_GInd)
		rm(scaledGInd)
		rm(pca.GInd) 
		
	    gc()

###########################################################################
# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################
   
    for( nCyc in 2:nCycles){
	   
		print(nCyc)

        selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		for(nFamily in 1:nFamilies){

			if(nCyc==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}
				

			}

			if(nCyc>2){

				# if(selectionOnGeno == TRUE){

					# selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				# }else if(selectionOnGeno == FALSE){
					
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
			 

		    }
	    }


	
### Exchange selected geno data among pairs
    if(migFreq !=0){
      if(nCyc>=2 && nCyc%%migFreq==0){

		output_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction,Policy)
	    
		selectedGenoData_List_afterExchange <- (output_List_afterExchange)[[1]]
		
		familyInfo <- output_List_afterExchange[[2]]
		
		Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   }else if(nCyc>=2 && nCyc%%migFreq !=0){
	  
		selectedGenoData_List_afterExchange <- selectedGenoData_List 
		
	   }
    }else if(migFreq ==0){
		selectedGenoData_List_afterExchange <- selectedGenoData_List 
	}
### Generate F5 RIL progeny for selected geno data
 
	
	F5RILs <- getF5RILs_BD(BD,selectedGenoData_List_afterExchange,nFamilies,nSel_inFamily,cycle1_nProgeny,nMarkers,NAM_LinkMap_New)
	
	
	Cycle_Progeny_F5 <- F5RILs[[1]]
	Cycle_Progeny_F5_List <- F5RILs[[2]]
	
	

 # dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

####################################################################################################################################################################################################
### This section works with full data table

 
		nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nFamilies,nCrosses_inFamily*nProgeny)
		newNextGenGenoTable <-  generateWeightedGenoTable(nextGenGenoTable,no_QTL)
		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

		for(nFamily in 1:nFamilies){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)

		}

		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)
		

		
		
	if(nCyc%%5 == 0 && sd(genoValSimValues !=0)){	

### Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}

### get island PG parameters

		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
      
	    		
		gc()
	   
				
		
			tol <- 1e-07
			nextGen_df <- as.matrix(nextGenGenoTable)
			nextGen_df <- tcrossprod(nextGen_df,nextGen_df) 
			eig1 <- eigen(nextGen_df,symmetric=TRUE) 
			eig <- eig1$values
			rank_X <- sum((eig/eig[1]) > tol)
			NF_Vec <- as.numeric(which(cumsum(eig/sum(eig)) >= 0.99)) 
			NF <- NF_Vec[1]
			
			if(NF == dim(nextGen_df)[1]){ 
				NF = rank_X 
			}else if (NF==1){
				NF <-rank_X
			}else if( ! exists("NF")) {
				NF <-rank_X 
			}else if(rank_X== 1){
				NF <- 2
			}
		
	       
	        if(rank_X > 0 && length(NF)>0){
		    
				j<- j+1
				
				scaledGInd <- scaleGen(nextGenGenoTable_GInd,scale=FALSE)
				pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF)
		        
				tryCatch(pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF), error=function(e) paste("Error dudiPCA Cycle No",nCyc,"-",e))
				
		   		
				filename <- paste("Sclass_Plot","-",condition,"-",Rep,"-",nCyc,".png",sep="")
				
				png(filename,width=1024,height=768,pointsize=20)
				tryCatch(s.class(pca.GInd$li, fac=pop(nextGenGenoTable_GInd),col=funky(20),addaxes=TRUE),error=function(e) print(paste("Error Sclass Cycle No",nCyc,"-",e)))
								
				tryCatch(title(paste("Cycle - ",nCyc,sep="")),error=function(e) print(paste("Error Title Cycle No",nCyc,"-",e)))
				dev.off()
				tryCatch(system(paste("gzip -f ",filename,sep="")),error=function(e) print(paste("Error gzip Cycle No",nCyc,"-",e)))
				
				tryCatch(PCA_Li_Table <- as.big.matrix(as.matrix(pca.GInd$li)),error=function(e) print(paste("Error BigMatrix Cycle No",nCyc,e)))
				tryCatch(PCA_C1_Table <- as.big.matrix(as.matrix(pca.GInd$c1)),error=function(e) print(paste("Error BigMatrix Cycle No",nCyc,e)))
				
				tryCatch(PCA_List[[j]] <- list(PCA_Li_Table,PCA_C1_Table),error=function(e)print(paste("Error PCA_List Cycle No",nCyc,e)))
			}
			rm(nextGenGenoTable_GInd)
			rm(scaledGInd)
			rm(pca.GInd)
		}
			
	
############################################################################################### 
### Prediction Model Update 	
### Extract training and test data for model build

	IndividualNames <-rep(0,(nCrosses*nProgeny))

    for(nInd in 1:(nCrosses*nProgeny)){

        IndividualNames[nInd]<- paste("Ind",nInd,sep="")

    }

##################################################################################################
	names(phenoValSimValues)<- IndividualNames
    Mean_Fixed<- rep(1,nIndividuals)

    phenotypicValuesSimTable<- cbind(phenoValSimValues,Mean_Fixed)

#### Table for Simulated Genotypic Values ####################################


    names(genoValSimValues)<-IndividualNames
    Mean_Fixed<- rep(1,nIndividuals)

    genotypicValuesSimTable <- cbind(genoValSimValues,Mean_Fixed)

#####################

	if(modelRetrain==TRUE){

		indices<-c(1:nIndividuals)

		trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))
		testIndices<- indices[which(!indices %in% trainIndices)]

	#### Training Set table for Genotype and Simulated Phenotype ###################

		# trainGenoTable <- genoTable_Mod[trainIndices,]
		trainGenoNewTable <- newNextGenGenoTable[trainIndices,]
		trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]
		trainSimGenoTable <-  genotypicValuesSimTable[trainIndices,]

	#### Validation Set table for Genotype and Simulated Phenotype ###################

		# testGenoTable<- genoTable_Mod[testIndices,]
		testGenoNewTable <- newNextGenGenoTable[testIndices,]
		testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]
		testSimGenoTable <-  genotypicValuesSimTable[testIndices,]

############### Changed RRBLUP to Bayes

	   if(sd(trainSimGenoTable[,1])!=0 && sd(trainSimPhenoTable[,1])!=0 && nCyc%%updateFrequency==0){
	  
	    if(modelType=="BayesB"){ 
	   		PredictionModel_Geno <- (buildBayesBModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable,h2))
			PredictionModel_Pheno <-(buildBayesBModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
		}else if(modelType=="BL"){ 
		
			PredictionModel_Geno <- (buildBLModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable,h2))
			PredictionModel_Pheno <-(buildBLModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
		}
	  }
	gc()
	  
	}
	
	if(modelUpdate ==TRUE && updateType == "EqPartition"){	
	
        indices<-c(1:nIndividuals)

        trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices<- indices[which(!indices %in% trainIndices)]

#### Training Set table for Genotype and Simulated Phenotype ###################

        trainGenoNewTable <- newNextGenGenoTable[trainIndices,]

        trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]

        trainSimGenoValTable <-  genotypicValuesSimTable[trainIndices,]

#### Validation Set table for Genotype and Simulated Phenotype ###################

        testGenoNewTable <- newNextGenGenoTable[testIndices,]

        testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]

        testSimGenoValTable <-  genotypicValuesSimTable[testIndices,]

		
## Training Set List with equal partitioning among cycles 		
		
		trainGenoNewTable_List[[nCyc]] <- as.big.matrix(trainGenoNewTable)
		trainSimGenoValTable_List[[nCyc]] <- as.big.matrix(trainSimGenoValTable)
		trainSimPhenoTable_List[[nCyc]] <- as.big.matrix(trainSimPhenoTable)
      
	    nTrainingSets <- length(trainGenoNewTable_List)
		nTrainIndividuals <- (0.8*nIndividuals) %/% nCyc
        nTrainIndividualsPlusRemainder <- nTrainIndividuals + ((0.8*nIndividuals) %% nCyc)
		
##### Build prediction model 

        trainGenoNewTableComb <- c() 
		trainSimGenoValTableComb <- c() 
		trainSimPhenoValTableComb <- c()
		 
		
        # trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))
 
        for(nTrainSets in 1:nTrainingSets){
		   
		    if(nTrainSets < nTrainingSets){
			   trainingIndices <- sample(c(1:(0.8*nIndividuals)),nTrainIndividuals)
		    }else if(nTrainSets == nTrainingSets){
			   trainingIndices <- sample(c(1:(0.8*nIndividuals)),nTrainIndividualsPlusRemainder)
		    }   
			   
			   trainGenoTableInCycle <- trainGenoNewTable_List[[nTrainSets]][trainingIndices,]
			   trainSimGenoValTableInCycle <- trainSimGenoValTable_List[[nTrainSets]][trainingIndices,]
			   trainSimPhenoTableInCycle <- trainSimPhenoTable_List[[nTrainSets]][trainingIndices,]
			   
			   trainGenoNewTableComb <- rbind((trainGenoNewTableComb),bigmemory::as.matrix(trainGenoTableInCycle))
			   trainSimGenoValTableComb <- rbind((trainSimGenoValTableComb),bigmemory::as.matrix(trainSimGenoValTableInCycle))
               trainSimPhenoValTableComb <- rbind((trainSimPhenoValTableComb),bigmemory::as.matrix(trainSimPhenoTableInCycle))
		}
		
        print(dim(trainGenoNewTableComb))

        if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){
            if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0 && nCyc%%updateFrequency==0){ 

				if(modelType=="BayesB"){
					PredictionModel_Geno <- (buildBayesBModel((trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					PredictionModel_Pheno <-(buildBayesBModel((trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
				}else if(modelType=="BL"){
			
					PredictionModel_Geno <- (buildBLModel((trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					PredictionModel_Pheno <-(buildBLModel((trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
			    }
		    }
		}
		 
		    # trainGenoNewTablePreCycle <-  trainGenoNewTableComb
            # trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)
            # trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)


            rm(trainGenoNewTableComb)
            rm(trainSimGenoValTableComb)
            rm(trainSimPhenoValTableComb)
            gc()

        }
		 
	if(modelUpdate ==TRUE && updateType == "FullSet"){

	  indices<-c(1:nIndividuals)

        trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices<- indices[which(!indices %in% trainIndices)]

#### Training Set table for Genotype and Simulated Phenotype ###################

        trainGenoNewTable <- newNextGenGenoTable[trainIndices,]

        trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]

        trainSimGenoValTable <-  genotypicValuesSimTable[trainIndices,]

#### Validation Set table for Genotype and Simulated Phenotype ###################

        testGenoNewTable <- newNextGenGenoTable[testIndices,]

        testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]

        testSimGenoValTable <-  genotypicValuesSimTable[testIndices,]


##### Build RRBLUP prediction model every cycle


        trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))

        trainSimGenoValTableComb <- rbind((trainSimGenoValTablePreCycle),trainSimGenoValTable)

        trainSimPhenoValTableComb <- rbind((trainSimPhenoTablePreCycle),trainSimPhenoTable)

        print(dim(trainGenoNewTableComb))

        if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){
            if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0 && nCyc%%updateFrequency==0){ 

				if(modelType=="BayesB"){
					PredictionModel_Geno <- (buildBayesBModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					PredictionModel_Pheno <-(buildBayesBModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
				}else if(modelType=="BL"){
			
					PredictionModel_Geno <- (buildBLModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					PredictionModel_Pheno <-(buildBLModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
			    }
		    }
		 
		    trainGenoNewTablePreCycle <-  trainGenoNewTableComb
            trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)
            trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)


            rm(trainGenoNewTableComb)
            rm(trainSimGenoValTableComb)
            rm(trainSimPhenoValTableComb)
            gc()

        }
		 
		 
	}
	
		
#################get predicted geno and pheno values with JM and unweighted GP models
### Get Predicted geno & pheno values for nextGenGenoData	
		
		if(Weighted==TRUE && WghtMethod=="JM"){
			genoValues  <- PredictBayesPhenoValues_Wght(newNextGenGenoTable,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newNextGenGenoTable,PredictionModel_Pheno,Freq)
	    }else if(Weighted==TRUE && WghtMethod =="DW"){
		    genoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newNextGenGenoTable,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)
		    phenoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newNextGenGenoTable,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

		    genoValues <- genoValues_PredList[[1]]
		    phenoValues <- phenoValues_PredList[[1]]

		    alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]
				
		}else if(Weighted==FALSE){
			genoValues <- PredictBayesPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
			phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

	    }
		

### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	   
	Criterion_List <- getCriterionValueList(BD,genoValues,phenoValues,genoValSimValues,phenoValSimValues,nCrosses_inFamily,nProgeny,selectionOnGeno,selectionOnSimulated,nFamilies)

	genoValues_List <- Criterion_List[[1]]
	phenoValues_List <- Criterion_List[[2]]
	genoValSimValues_List <- Criterion_List[[3]]
	phenoValSimValues_List <- Criterion_List[[4]] 
	 

	 # for(i in 1: length(genoValSimValues_List)){ 
	 
		# print(summary(genoValSimValues_List[[i]]))
	 # }
### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
		    if(selectionOnGeno ==TRUE){
		
				genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
		    
			}else if(selectionOnGeno ==FALSE){ 

				genoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
            }			
			
		}
		if (selectionOnSimulated==FALSE){
			
			if(selectionOnGeno ==TRUE){
			
				genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			}else if(selectionOnGeno ==FALSE){
			
				 genoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			}
		}

		if(nCyc>=2 && nCyc%%migPolicyFreq==0){
			migrationGroups <- migrationPolicy(nextGenGenoTable,genoSelectionTable,Policy,nFamilies)
			emigrantGroups <- migrationGroups[[1]]
			immigrantGroups<- migrationGroups[[2]]
		}
		
#####################################################################################################################
## Assign output variables for cycle n

		for(nFamily in 1:nFamilies){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,nCyc] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,nCyc] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,nCyc]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,nCyc] <- phenoValues_List[[nFamily]]

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
			genoSimSelectedValues <- genoValSimValues_List[[nFamily]][selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues_List[[nFamily]][selectedGenoIndividualIndices]
				
			selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
			attainedGenoValues_List[[nFamily]][nCyc] <- max(genoSimSelectedValues)
			
			if(selectionOnSimulated ==TRUE){
		        if(selectionOnGeno ==TRUE){
		
				GenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- as.vector(genoSelectionTable[[nFamily]][,1])
			    PhenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- phenoValSimValues_List[[nFamily]][selectedGenoIndividualIndices]
				
			   } else if(selectionOnGeno ==FALSE){ 

				GenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- genoValSimValues_List[[nFamily]][selectedGenoIndividualIndices]
				PhenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- as.vector(genoSelectionTable[[nFamily]][,1])
			   }			
			}
			if (selectionOnSimulated==FALSE){
			
				if(selectionOnGeno ==TRUE){
				
				GenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- as.vector(genoSelectionTable[[nFamily]][,1])
			        PhenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- phenoValues_List[[nFamily]][selectedGenoIndividualIndices]
								
				}else if(selectionOnGeno ==FALSE){
				
					GenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- genoValues_List[[nFamily]][selectedGenoIndividualIndices]
					PhenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- as.vector(genoSelectionTable[[nFamily]][,1])
				}
		         }
		}
		
  }	
	
#########

	 simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List,Migration_List,PCA_List)

	  return(simResults_List)

 }




###################### 
## Fn43: 

runSimulations20X_DiscreteSelection_BD_BayesB_Wght <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap,ModelRetrain,ModelType,UpdateFrequency,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined){

    library(bigmemory)
	library(adegenet)
	library(mmod)
### Assign Variables #########################################################################


  	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues


	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny


	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  nMarkers <- 4289

	  NAM_LinkMap_New <- NAM_LinkMap
	  
	  Weighted <- weighted
	  WghtMethod <- wghtMethod
	  
	  BD <-  BreedDesign
	  
	  AlleleConversionTable_Combined <-  alleleConversionTable_Combined

################

	  nFamilies <- 20
	 
	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }

      nCrosses_inFamily <- nSel_inFamily

	  nProgeny <- nProgeny/nSel_inFamily


	  modelType <- ModelType
      modelRetrain <- ModelRetrain
      updateFrequency <- UpdateFrequency

 
	  condition <- i
	  Rep <- k  


#### /nSel_inFamily

###############################################################################################
## Get predicted geno and pheno values for cycle1
	 nCyc <-1

	 cycle1GenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	 newCycle1GenoTable <- generateWeightedGenoTable(cycle1GenoTable,no_QTL)
	 nIndividuals <- nrow(cycle1GenoTable)
	 
	if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[3])
	}else if(selectionOnGeno==FALSE){
	    MarkerEffects <- unlist(PredictionModel_Pheno[3])
	}

     
     Freq_List <- getFreq_BasePoln(cycle1GenoTable,MarkerEffects)
	 Freq <- Freq_List[[1]]
	 FavAllele <- Freq_List[[2]]


	 if(Weighted==TRUE){
		genoValues  <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Geno,Freq)
		phenoValues <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Pheno,Freq)
	 }else if(Weighted==FALSE){
		genoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Geno)
		phenoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Pheno)
	 }

    
########  Split F5 Progeny into list based on families

	  F5_Progeny_List_Family<- cycle1_Geno_Data

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }

######## split cycle1_Geno data in to list according to families
######## Get percentHeterozygous & favorableAllele frequency
      
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:20){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[
			nFamily]],1,cycle1_nProgeny)

	  }
###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1

		}

## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)

	  #percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)


	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)


	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){

		if(selectionOnSimulated == FALSE){
		 
          if(selectionOnGeno==TRUE){
			  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
			  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
			  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]
			  sortedPhenoValues <- phenoValues_List[[nFamily]][GenoTableIndices_topN]
			  PhenoTableIndices_topN <-  GenoTableIndices_topN
			}else if(selectionOnGeno==FALSE){
		  	  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
			  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
			  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]
			  sortedGenoValues <- genoValues_List[[nFamily]][PhenoTableIndices_topN]
			  GenoTableIndices_topN <- PhenoTableIndices_topN 
	        }
		}
		if (selectionOnSimulated==TRUE){

			if(selectionOnGeno==TRUE){
				sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
				topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
				GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]
				sortedPhenoValues <- phenoSimValues_List[[nFamily]][GenoTableIndices_topN]
				PhenoTableIndices_topN <-  GenoTableIndices_topN
			
			}else if(selectionOnGeno==FALSE){

				sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
				topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
				PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]
				GenoTableIndices_topN <- PhenoTableIndices_topN 
				sortedGenoValues <- genoSimValues_List[[nFamily]][PhenoTableIndices_topN]
			
			}
	    }
	
	  # if(selectionOnSimulated == FALSE){
		  # sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  # topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  # GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  # sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  # topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  # PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   # }else if (selectionOnSimulated==TRUE){

		  # sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  # topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  # GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  # sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  # topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  # PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   # }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  selectedGenoSimValues <- genoSimValues_List[[nFamily]][GenoTableIndices_topN]

	  attainedGenoValues_List[[nFamily]][1] <- max(selectedGenoSimValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


#########################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- (PhenoVal_NX_N_3c_List[[nFamily]][,1])
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
	}
	

###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1

		}

### get island PG parameters

		cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable,AlleleConversionTable_Combined)
 		cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
      
		
		
		scaledGInd <- scaleGen(cycle1GenoTable_GInd,scale=FALSE)
        pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=100)
	  
	    filename <- paste("Sclass_Plot","-",condition,"-",Rep,"-",nCyc,".png",sep="")
		
		png(filename,width=1024,height=768,pointsize=20)
		s.class(pca.GInd$li, fac=pop(cycle1GenoTable_GInd),col=funky(20),addaxes=TRUE)
		dev.off()
		system(paste("gzip ",filename,sep=""))
		
		rm(cycle1GenoTable_GInd)
		rm(scaledGInd)
		rm(pca.GInd) 
		
	    gc()
	
###########################################################################

############### Cycles 2*29 ###############################################

    # nCycles <- 10

	for( i in 2:nCycles){

	    # i <- 2
	   
	    nCyc <- i
		print(i)

        selectedGenoData_List <- list()
		selectedGenoData_List <- list()


		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }



### Generate F5 RIL progeny for selected geno data
 

    F5RILs <- getF5RILs_BD(BD,selectedGenoData_List,nFamilies,nSel_inFamily,nProgeny,nMarkers,NAM_LinkMap)

	Cycle_Progeny_F5 <- F5RILs[[1]]
	Cycle_Progeny_F5_List <- F5RILs[[2]]
	
	

 # dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

####################################################################################################################################################################################################
### This section works with full data table

  		# nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,(nSel_inFamily*cycle1_nProgeny),nFamilies,no_selected)


		nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nFamilies,nCrosses_inFamily*nProgeny)
		newNextGenGenoTable <-  generateWeightedGenoTable(nextGenGenoTable,no_QTL)
		
		print(dim(nextGenGenoTable))

		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

		for(nFamily in 1:20){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)

		}
		
		Freq <- getFreq(nextGenGenoTable,FavAllele) 
		genoSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)
		
	
    if(i %% 5 ==0){		
###  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}

### get island PG parameters

		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
      
			
		gc()
	
		
		scaledGInd <- scaleGen(nextGenGenoTable_GInd,scale=FALSE)
        pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=100)
	  
	    filename <- paste("Sclass_Plot","-",condition,"-",Rep,"-",nCyc,".png",sep="")
		
		png(filename,width=1024,height=768,pointsize=20)
		s.class(pca.GInd$li, fac=pop(nextGenGenoTable_GInd),col=funky(20),addaxes=TRUE)
		dev.off()
		system(paste("gzip ",filename,sep=""))
		
		rm(nextGenGenoTable_GInd)
		rm(scaledGInd)
		rm(pca.GInd)
		
	}

	
############################################################################################### 
### Prediction Model Update 	
### Extract training and test data for model build

	IndividualNames <-rep(0,(nCrosses*nProgeny))

    for(nInd in 1:(nCrosses*nProgeny)){

        IndividualNames[nInd]<- paste("Ind",nInd,sep="")

    }

##################################################################################################
	names(phenoSimValues)<- IndividualNames
    Mean_Fixed<- rep(1,nIndividuals)

    phenotypeSimTable<- cbind(phenoSimValues,Mean_Fixed)

#### Table for Simulated Genotypic Values ####################################


    names(genoSimValues)<-IndividualNames
    Mean_Fixed<- rep(1,nIndividuals)

    genotypicValuesSimTable <- cbind(genoSimValues,Mean_Fixed)

#####################

	if(modelRetrain==TRUE){

		indices<-c(1:nIndividuals)

		trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))
		testIndices<- indices[which(!indices %in% trainIndices)]

	#### Training Set table for Genotype and Simulated Phenotype ###################

		# trainGenoTable <- genoTable_Mod[trainIndices,]
		trainGenoNewTable <- newNextGenGenoTable[trainIndices,]
		trainSimPhenoTable <- phenotypeSimTable[trainIndices,]
		trainSimGenoTable <-  genotypicValuesSimTable[trainIndices,]

	#### Validation Set table for Genotype and Simulated Phenotype ###################

		# testGenoTable<- genoTable_Mod[testIndices,]
		testGenoNewTable <- newNextGenGenoTable[testIndices,]
		testSimPhenoTable <- phenotypeSimTable[testIndices,]
		testSimGenoTable <-  genotypicValuesSimTable[testIndices,]

############### Changed RRBLUP to Bayes

	   if(sd(trainSimGenoTable[,1])!=0 && sd(trainSimPhenoTable[,1])!=0 && i%%updateFrequency==0){
	  
	    if(modelType=="BayesB"){ 
	   		PredictionModel_Geno <- (buildBayesBModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable,h2))
			PredictionModel_Pheno <-(buildBayesBModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
		}else if(modelType=="BL"){ 
		
			PredictionModel_Geno <- (buildBLModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable))
			PredictionModel_Pheno <-(buildBLModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))
		}
	  }

	}

	gc()

#################get predicted geno and pheno values with JM and unweighted GP models
### Get Predicted geno & pheno values for nextGenGenoData	
		
		if(Weighted==TRUE){
			genoValues  <- PredictBayesPhenoValues_Wght(newNextGenGenoTable,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newNextGenGenoTable,PredictionModel_Pheno,Freq)
	    }else if(Weighted==FALSE){
			genoValues <- PredictBayesPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
			phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

	    }
		

### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	   
Criterion_List <-getCriterionValueList(BD,genoValues,phenoValues,genoSimValues,phenoSimValues,nCrosses_inFamily,nProgeny,selectionOnGeno,selectionOnSimulated,nFamilies)

 genoValues_List <- Criterion_List[[1]]
 phenoValues_List <- Criterion_List[[2]]
 genoValSimValues_List <- Criterion_List[[3]]
 phenoValSimValues_List <- Criterion_List[[4]] 
 

### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)

		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)

		}


		
#####################################################################################################################
## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}
		
	}
	
	
	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List)

	  return(simResults_List)

}


############################ 
### Fn 44: 


runSimulations20X_IslandSelection_BD_Wght <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,SclassPlot){


	options(warn=-1)
	
### Assign Variables #########################################################################

  	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues


	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny


	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  nMarkers <- 4289

	  NAM_LinkMap_New <- NAM_LinkMap
	  
	  Weighted <- weighted
	  WghtMethod <- wghtMethod
	  
	  BD <- BreedDesign
	    
	  Policy <- policy	
      migPolicyFreq <- policyChangeFrequency
	  
############################################################
 
	  migFreq <- migrationFrequency
	  migrationSize <- MigrationSize
	  direction <- Direction

	  nFamilies <- 20
	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }

      nCrosses_inFamily <- nSel_inFamily

	  nProgeny <- nProgeny/nSel_inFamily


	  modelType <- ModelType
      modelRetrain <- ModelRetrain
	  retrainFrequency <- RetrainFrequency
	  
	  modelUpdate <- ModelUpdate
      updateFrequency <- UpdateFrequency
      updateType <- UpdateType
 
	  condition <- i
	  Rep <- k 
       
      sclassPlot <- SclassPlot 

      AlleleConversionTable_Combined <-  alleleConversionTable_Combined
	 
### Initialize variables
	 
	  Migration_List <- list()
	  PCA_List <- list()
	  j<-1

	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
	  nCyc <-1
	  trainGenoNewTable_List[[nCyc]] <-  as.big.matrix(trainGenoNewTable)
      trainSimGenoValTable_List[[nCyc]] <- as.big.matrix(trainSimGenoValTable)
	  trainSimPhenoTable_List[[nCyc]] <- as.big.matrix(trainSimPhenoValTable)
	  
	  trainGenoNewTablePreCycle <- as.big.matrix(trainGenoNewTable)
      trainSimPhenoTablePreCycle <- trainSimPhenoValTable
      trainSimGenoValTablePreCycle <- trainSimGenoValTable
		
#### /nSel_inFamily
###############################################################################################
## Get predicted geno and pheno values for cycle1 for combined population set

	cycle1GenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable <- generateWeightedGenoTable(cycle1GenoTable,no_QTL)
	nIndividuals <- nrow(cycle1GenoTable)
	
    if(modelType=="BayesB" || modelType =="BL"){	
		if(selectionOnGeno==TRUE){
			MarkerEffects <- unlist(PredictionModel_Geno[3])
		}else if(selectionOnGeno==FALSE){
			MarkerEffects <- unlist(PredictionModel_Pheno[3])
		}
    }	
    if(modelType=="RRBLUP" || modelType =="RRBLUP_REML"){	
		if(selectionOnGeno==TRUE){
			MarkerEffects <- unlist(PredictionModel_Geno[[1]])
		}else if(selectionOnGeno==FALSE){
			MarkerEffects <- unlist(PredictionModel_Pheno[[1]])
		}
    }
     
     Freq_List <- getFreq_BasePoln(cycle1GenoTable,MarkerEffects)
	 Freq <- Freq_List[[1]]
	 FavAllele <- Freq_List[[2]]


     if(modelType=="BayesB" || modelType=="BL"){

         if(Weighted==TRUE){
                genoValues  <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Geno,Freq)
                phenoValues <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Pheno,Freq)
         }else if(Weighted==FALSE){
                genoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Geno)
                phenoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Pheno)
         }
     }

     if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newCycle1GenoTable,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newCycle1GenoTable,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newCycle1GenoTable,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newCycle1GenoTable,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

            genoValues  <- PredictRRBLUPPhenoValues(newCycle1GenoTable,PredictionModel_Geno)
            phenoValues <- PredictRRBLUPPhenoValues(newCycle1GenoTable,PredictionModel_Pheno)

        }
}


########  Split F5 Progeny into list based on families

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:nFamilies){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }

######## split cycle1_Geno data in to list according to families
######## Get percentHeterozygous & favorableAllele frequency
     
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:nFamilies){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[
			nFamily]],1,cycle1_nProgeny)

	  }


## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()


	  for(nFamily in 1:nFamilies){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }


### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	 

	  #percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)


	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)


	  selectedGenoIndividualIndices2_List <- list()
	
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()

######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:nFamilies){


        if(selectionOnSimulated == FALSE){
		 
          if(selectionOnGeno==TRUE){
			  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
			  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
			  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]
			  
			  sortedPhenoValues <- phenoValues_List[[nFamily]][GenoTableIndices_topN]
			  PhenoTableIndices_topN <-  GenoTableIndices_topN
			  topNPhenoValues <- phenoValues_List[[nFamily]][GenoTableIndices_topN]
			
			}else if(selectionOnGeno==FALSE){
		  	  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
			  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
			  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]
			  
			  sortedGenoValues <- genoValues_List[[nFamily]][PhenoTableIndices_topN]
			  GenoTableIndices_topN <- PhenoTableIndices_topN 
			  topNGenoValues <- genoValues_List[[nFamily]][PhenoTableIndices_topN]
	        }
		}
		if (selectionOnSimulated==TRUE){

			if(selectionOnGeno==TRUE){
				sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
				topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
				GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]
				
				sortedPhenoValues <- phenoSimValues_List[[nFamily]][GenoTableIndices_topN]
				PhenoTableIndices_topN <-  GenoTableIndices_topN
				topNPhenoValues <- phenoSimValues_List[[nFamily]][GenoTableIndices_topN]
				
			
			}else if(selectionOnGeno==FALSE){

				sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
				topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
				PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]
				
				GenoTableIndices_topN <- PhenoTableIndices_topN 
				sortedGenoValues <- genoSimValues_List[[nFamily]][PhenoTableIndices_topN]
			    topNGenoValues <- genoSimValues_List[[nFamily]][PhenoTableIndices_topN]
			}
	    }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  selectedGenoSimValues <- genoSimValues_List[[nFamily]][GenoTableIndices_topN]
	  attainedGenoValues_List[[nFamily]][1] <- max(selectedGenoSimValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN

    
#########################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues
	  
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
	}
	
	
	migrationGroups <- migrationPolicy(cycle1GenoTable,genoSelectionTable,Policy,nFamilies,no_QTL)

	
	emigrantGroups <- migrationGroups[[1]]
	immigrantGroups <- migrationGroups[[2]]
	
###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(pop in 1:20){

			final <- pop*100
			populations[init:final]<- rep(pop,100)

			init <- final+1

		}


if(sclassPlot==TRUE){
### get island PG parameters

		cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable,AlleleConversionTable_Combined)
 		cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
		tol <- 1e-07
		cycle1Gen_df <- as.matrix(cycle1GenoTable)
		cycle1Gen_df <- tcrossprod(cycle1Gen_df,cycle1Gen_df) 
		eig1 <- eigen(cycle1Gen_df,symmetric=TRUE) 
		eig <- eig1$values
		rank_X <- sum((eig/eig[1]) > tol)
		NF_Vec <- as.numeric(which(cumsum(eig/sum(eig)) >= 0.99)) 
		NF <- NF_Vec[1]
			
		if(NF == dim(cycle1Gen_df)[1]) { 
			
				NF = rank_X 
		}else if (NF==1){
				NF <-rank_X
		}else if( ! exists("NF")) {
				NF <-rank_X 
		}else if(rank_X== 1){
				NF <- 2
		}
		
			      
		scaledGInd <- scaleGen(cycle1GenoTable_GInd,scale=FALSE)
        
		tryCatch(pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF), error=function(e) paste("Error dudiPCA Cycle No",i,"-",e))
				
		   		
		filename <- paste("Sclass_Plot","-",condition,"-",Rep,"-",nCyc,".png",sep="")
				
		png(filename,width=1024,height=768,pointsize=20)
		tryCatch(s.class(pca.GInd$li, fac=pop(cycle1GenoTable_GInd),col=funky(20),addaxes=TRUE),error=function(e) print(paste("Error Sclass Cycle No",i,"-",e)))
				
								
		tryCatch(title(paste("Cycle - ",i,sep="")),error=function(e) print(paste("Error Title Cycle No",i,"-",e)))
		dev.off()
		tryCatch(system(paste("gzip -f ",filename,sep="")),error=function(e) print(paste("Error gzip Cycle No",i,"-",e)))
				
		tryCatch(PCA_Li_Table <- as.big.matrix(as.matrix(pca.GInd$li)),error=function(e) print(paste("Error BigMatrix Cycle No",i,e)))
		tryCatch(PCA_C1_Table <- as.big.matrix(as.matrix(pca.GInd$c1)),error=function(e) print(paste("Error BigMatrix Cycle No",i,e)))
		tryCatch(PCA_List[[j]] <- list(PCA_Li_Table,PCA_C1_Table),error=function(e)print(paste("Error PCA_List Cycle No",i,e)))
		
		# pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF)
	    # filename <- paste("Sclass_Plot","-",condition,"-",Rep,"-",nCyc,".png",sep="")
		# png(filename,width=1024,height=768,pointsize=20)
		# s.class(pca.GInd$li, fac=pop(cycle1GenoTable_GInd),col=funky(20),addaxes=TRUE)
		# title(paste("Cycle - ",i,sep=""))
		# dev.off()
		# PCA_Li_Table <- as.big.matrix(as.matrix(pca.GInd$li))
		# PCA_C1_Table <- as.big.matrix(as.matrix(pca.GInd$c1)) 
		# PCA_List[[j]] <- list(PCA_Li_Table,PCA_C1_Table)
		# system(paste("gzip ",filename,sep=""))
		
		rm(cycle1GenoTable_GInd)
		rm(scaledGInd)
		rm(pca.GInd) 
		
	    gc()
}

                rm(cycle1GenoTable) 
                rm(newCycle1GenoTable)

###########################################################################
# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################
   
    for( nCyc in 2:nCycles){
	   
		print(nCyc)

        selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		for(nFamily in 1:nFamilies){

			if(nCyc==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}
				

			}

			if(nCyc>2){

									
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
			 

		         }
	    }


	
### Exchange selected geno data among pairs
    if(migFreq !=0){
      if(nCyc>=2 && nCyc%%migFreq==0){

		output_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
		selectedGenoData_List_afterExchange <- (output_List_afterExchange)[[1]]
		
		familyInfo <- output_List_afterExchange[[2]]
		
		Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   }else if(nCyc>=2 && nCyc%%migFreq !=0){
	  
		selectedGenoData_List_afterExchange <- selectedGenoData_List 
		
	   }
    }else if(migFreq ==0){
		selectedGenoData_List_afterExchange <- selectedGenoData_List 
	}
### Generate F5 RIL progeny for selected geno data
 
	
	F5RILs <- getF5RILs_BD(BD,selectedGenoData_List_afterExchange,nFamilies,nSel_inFamily,cycle1_nProgeny,nMarkers,NAM_LinkMap_New)
	
	
	Cycle_Progeny_F5 <- F5RILs[[1]]
	Cycle_Progeny_F5_List <- F5RILs[[2]]
	
	

 # dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

####################################################################################################################################################################################################
### This section works with full data table

 
		nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nFamilies,nCrosses_inFamily*nProgeny)
		newNextGenGenoTable <-  generateWeightedGenoTable(nextGenGenoTable,no_QTL)
		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

		for(nFamily in 1:nFamilies){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)

		}

		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)
		
if(sclassPlot==TRUE){		
		
	if(nCyc%%5 == 0 && sd(genoValSimValues !=0)){	

### Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}

### get island PG parameters

		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
      
	    		
		gc()
	   
				
		
			tol <- 1e-07
			nextGen_df <- as.matrix(nextGenGenoTable)
			nextGen_df <- tcrossprod(nextGen_df,nextGen_df) 
			eig1 <- eigen(nextGen_df,symmetric=TRUE) 
			eig <- eig1$values
			rank_X <- sum((eig/eig[1]) > tol)
			NF_Vec <- as.numeric(which(cumsum(eig/sum(eig)) >= 0.99)) 
			NF <- NF_Vec[1]
			
			if(NF == dim(nextGen_df)[1]){ 
				NF = rank_X 
			}else if (NF==1){
				NF <-rank_X
			}else if( ! exists("NF")) {
				NF <-rank_X 
			}else if(rank_X== 1){
				NF <- 2
			}
		
	       
	        if(rank_X > 0 && length(NF)>0){
		    
				j<- j+1
				
				scaledGInd <- scaleGen(nextGenGenoTable_GInd,scale=FALSE)
				pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF)
		        
				tryCatch(pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF), error=function(e) paste("Error dudiPCA Cycle No",nCyc,"-",e))
				
		   		
				filename <- paste("Sclass_Plot","-",condition,"-",Rep,"-",nCyc,".png",sep="")
				
				png(filename,width=1024,height=768,pointsize=20)
				tryCatch(s.class(pca.GInd$li, fac=pop(nextGenGenoTable_GInd),col=funky(20),addaxes=TRUE),error=function(e) print(paste("Error Sclass Cycle No",nCyc,"-",e)))
								
				tryCatch(title(paste("Cycle - ",nCyc,sep="")),error=function(e) print(paste("Error Title Cycle No",nCyc,"-",e)))
				dev.off()
				tryCatch(system(paste("gzip -f ",filename,sep="")),error=function(e) print(paste("Error gzip Cycle No",nCyc,"-",e)))
				
				tryCatch(PCA_Li_Table <- as.big.matrix(as.matrix(pca.GInd$li)),error=function(e) print(paste("Error BigMatrix Cycle No",nCyc,e)))
				tryCatch(PCA_C1_Table <- as.big.matrix(as.matrix(pca.GInd$c1)),error=function(e) print(paste("Error BigMatrix Cycle No",nCyc,e)))
				
				tryCatch(PCA_List[[j]] <- list(PCA_Li_Table,PCA_C1_Table),error=function(e)print(paste("Error PCA_List Cycle No",nCyc,e)))
			}
			rm(nextGenGenoTable_GInd)
			rm(scaledGInd)
			rm(pca.GInd)
		}
			
}	
############################################################################################### 
### Prediction Model Update 	
### Extract training and test data for model build


if(modelUpdate == TRUE && nCyc%%updateFrequency==0){

	IndividualNames <-rep(0,(nCrosses*nProgeny))

    for(nInd in 1:(nCrosses*nProgeny)){

        IndividualNames[nInd]<- paste("Ind",nInd,sep="")

    }

##################################################################################################
	names(phenoValSimValues)<- IndividualNames
    Mean_Fixed<- rep(1,nIndividuals)

    phenotypicValuesSimTable<- cbind(phenoValSimValues,Mean_Fixed)

#### Table for Simulated Genotypic Values ####################################


    names(genoValSimValues)<-IndividualNames
    Mean_Fixed<- rep(1,nIndividuals)

    genotypicValuesSimTable <- cbind(genoValSimValues,Mean_Fixed)

#####################

	if(modelRetrain==TRUE){

		indices<-c(1:nIndividuals)

		trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))
		testIndices<- indices[which(!indices %in% trainIndices)]

	#### Training Set table for Genotype and Simulated Phenotype ###################

		# trainGenoTable <- genoTable_Mod[trainIndices,]
		trainGenoNewTable <- newNextGenGenoTable[trainIndices,]
		trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]
		trainSimGenoTable <-  genotypicValuesSimTable[trainIndices,]

	#### Validation Set table for Genotype and Simulated Phenotype ###################

		# testGenoTable<- genoTable_Mod[testIndices,]
		testGenoNewTable <- newNextGenGenoTable[testIndices,]
		testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]
		testSimGenoTable <-  genotypicValuesSimTable[testIndices,]

############### Changed RRBLUP to Bayes

	   if(sd(trainSimGenoTable[,1])!=0 && sd(trainSimPhenoTable[,1])!=0 && nCyc%%updateFrequency==0){
	  
	    if(modelType=="BayesB"){ 
	   		PredictionModel_Geno <- (buildBayesBModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable,h2))
			PredictionModel_Pheno <-(buildBayesBModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
		}else if(modelType=="BL"){ 
		
			PredictionModel_Geno <- (buildBLModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable,h2))
			PredictionModel_Pheno <-(buildBLModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
		}else if(modelType=="RRBLUP"){

                        PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableiComb),trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))

        }else if(modelType=="RRBLUP_REML"){


                        PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        PredictionModel_Pheno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

        }
		
	  }
	  gc()
	  
	}
	
	if(modelUpdate ==TRUE && updateType == "EqPartition"){	
	
	
	PredictionModel_Geno_PreCycle <- PredictionModel_Geno
        PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno
        
	indices<-c(1:nIndividuals)

        trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices<- indices[which(!indices %in% trainIndices)]

#### Training Set table for Genotype and Simulated Phenotype ###################

        trainGenoNewTable <- newNextGenGenoTable[trainIndices,]

        trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]

        trainSimGenoValTable <-  genotypicValuesSimTable[trainIndices,]

#### Validation Set table for Genotype and Simulated Phenotype ###################

        testGenoNewTable <- newNextGenGenoTable[testIndices,]

        testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]

        testSimGenoValTable <-  genotypicValuesSimTable[testIndices,]

		
## Training Set List with equal partitioning among cycles 		
		
	trainGenoNewTable_List[[nCyc]] <- as.big.matrix(trainGenoNewTable)
	trainSimGenoValTable_List[[nCyc]] <- as.big.matrix(trainSimGenoValTable)
	trainSimPhenoTable_List[[nCyc]] <- as.big.matrix(trainSimPhenoTable)
      
	nTrainingSets <- length(trainGenoNewTable_List)
	nTrainIndividuals <- (0.8*nIndividuals) %/% nCyc
        nTrainIndividualsPlusRemainder <- nTrainIndividuals + ((0.8*nIndividuals) %% nCyc)
		
##### Build prediction model 

        trainGenoNewTableComb <- c() 
		trainSimGenoValTableComb <- c() 
		trainSimPhenoValTableComb <- c()
		 
		
        # trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))
 
        for(nTrainSets in 1:nTrainingSets){
		   
		    if(nTrainSets < nTrainingSets){
			   trainingIndices <- sample(c(1:(0.8*nIndividuals)),nTrainIndividuals)
		    }else if(nTrainSets == nTrainingSets){
			   trainingIndices <- sample(c(1:(0.8*nIndividuals)),nTrainIndividualsPlusRemainder)
		    }   
			   
			   trainGenoTableInCycle <- trainGenoNewTable_List[[nTrainSets]][trainingIndices,]
			   trainSimGenoValTableInCycle <- trainSimGenoValTable_List[[nTrainSets]][trainingIndices,]
			   trainSimPhenoTableInCycle <- trainSimPhenoTable_List[[nTrainSets]][trainingIndices,]
			   
			   trainGenoNewTableComb <- rbind((trainGenoNewTableComb),bigmemory::as.matrix(trainGenoTableInCycle))
			   trainSimGenoValTableComb <- rbind((trainSimGenoValTableComb),bigmemory::as.matrix(trainSimGenoValTableInCycle))
                           trainSimPhenoValTableComb <- rbind((trainSimPhenoValTableComb),bigmemory::as.matrix(trainSimPhenoTableInCycle))
		}
		
        print(dim(trainGenoNewTableComb))

        if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){
            if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0 && nCyc%%updateFrequency==0){ 

				if(modelType=="BayesB"){
					PredictionModel_Geno <- (buildBayesBModel((trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					PredictionModel_Pheno <-(buildBayesBModel((trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
				}else if(modelType=="BL"){
			
					PredictionModel_Geno <- (buildBLModel((trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					PredictionModel_Pheno <-(buildBLModel((trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
			    }else if(modelType=="RRBLUP"){

                        PredictionModel_Geno <- (buildRRBLUPModel((trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel((trainGenoNewTableiComb),trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))

                }else if(modelType=="RRBLUP_REML"){



                        PredictionModel_Geno <- buildRRBLUPModel_REML((trainGenoNewTableComb), trainSimGenoValTableComb,(testGenoNewTable),testSimGenoValTable)

                        PredictionModel_Pheno <- buildRRBLUPModel_REML((trainGenoNewTableComb), trainSimPhenoValTableComb,(testGenoNewTable),testSimPhenoTable)

                }
		    }
	}
		 
		    if((length(levels(factor(PredictionModel_Pheno[[1]]))) == 1 && sd(PredictionModel_Pheno[[1]])==0) ||(length(levels(factor(PredictionModel_Geno[[1]]))) ==1 && sd(PredictionModel_Geno[[1]])==0)){

                         PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         PredictionModel_Geno <- PredictionModel_Geno_PreCycle

                     }
		   

            rm(trainGenoNewTableComb)
            rm(trainSimGenoValTableComb)
            rm(trainSimPhenoValTableComb)
            gc()

        }
		 
	if(modelUpdate ==TRUE && updateType == "FullSet"){
	
	   PredictionModel_Geno_PreCycle <- PredictionModel_Geno
       PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno

	    indices<-c(1:nIndividuals)

        trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices<- indices[which(!indices %in% trainIndices)]

#### Training Set table for Genotype and Simulated Phenotype ###################

        trainGenoNewTable <- newNextGenGenoTable[trainIndices,]

        trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]

        trainSimGenoValTable <-  genotypicValuesSimTable[trainIndices,]

#### Validation Set table for Genotype and Simulated Phenotype ###################

        testGenoNewTable <- newNextGenGenoTable[testIndices,]

        testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]

        testSimGenoValTable <-  genotypicValuesSimTable[testIndices,]


##### Build RRBLUP prediction model every cycle


        trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))

        trainSimGenoValTableComb <- rbind((trainSimGenoValTablePreCycle),trainSimGenoValTable)

        trainSimPhenoValTableComb <- rbind((trainSimPhenoTablePreCycle),trainSimPhenoTable)

        print(dim(trainGenoNewTableComb))

        if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){
            if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0 && nCyc%%updateFrequency==0){ 

				if(modelType=="BayesB"){
					PredictionModel_Geno <- (buildBayesBModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					PredictionModel_Pheno <-(buildBayesBModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
				}else if(modelType=="BL"){
			
					PredictionModel_Geno <- (buildBLModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					PredictionModel_Pheno <-(buildBLModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
			    }else if(modelType=="RRBLUP"){

                        PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableiComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable))

                }else if(modelType=="RRBLUP_REML"){



                        PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        PredictionModel_Pheno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

                }
				
		    }
		 
		    trainGenoNewTablePreCycle <-  trainGenoNewTableComb
            trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)
            trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
			
	if((length(levels(factor(PredictionModel_Pheno[[1]]))) == 1 && sd(PredictionModel_Pheno[[1]])==0) ||(length(levels(factor(PredictionModel_Geno[[1]]))) ==1 && sd(PredictionModel_Geno[[1]])==0)){

                         PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         PredictionModel_Geno <- PredictionModel_Geno_PreCycle

       }
			
			

            rm(trainGenoNewTableComb)
            rm(trainSimGenoValTableComb)
            rm(trainSimPhenoValTableComb)
            gc()

        }
}
	

}	



#################get predicted geno and pheno values with JM and unweighted GP models
### Get Predicted geno & pheno values for nextGenGenoData	

if(modelType=="BayesB" || modelType=="BL"){
		
		if(Weighted==TRUE && WghtMethod=="JM"){
			genoValues  <- PredictBayesPhenoValues_Wght(newNextGenGenoTable,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newNextGenGenoTable,PredictionModel_Pheno,Freq)
	        }else if(Weighted==TRUE && WghtMethod =="DW"){
		    genoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newNextGenGenoTable,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)
		    phenoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newNextGenGenoTable,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

		    genoValues <- genoValues_PredList[[1]]
		    phenoValues <- phenoValues_PredList[[1]]

		    alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]
				
		}else if(Weighted==FALSE){
			genoValues <- PredictBayesPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
			phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

	        }	
}

if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

            genoValues  <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
            phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

        }
}


### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	   
	Criterion_List <- getCriterionValueList(BD,genoValues,phenoValues,genoValSimValues,phenoValSimValues,nCrosses_inFamily,nProgeny,selectionOnGeno,selectionOnSimulated,nFamilies)

	genoValues_List <- Criterion_List[[1]]
	phenoValues_List <- Criterion_List[[2]]
	genoValSimValues_List <- Criterion_List[[3]]
	phenoValSimValues_List <- Criterion_List[[4]] 
	 
### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
		    if(selectionOnGeno ==TRUE){
		
				genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
		    
			}else if(selectionOnGeno ==FALSE){ 

				genoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
            }			
			
		}
		if (selectionOnSimulated==FALSE){
			
			if(selectionOnGeno ==TRUE){
			
				genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			}else if(selectionOnGeno ==FALSE){
			
				 genoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			}
		}

		if(nCyc>=2 && nCyc%%migPolicyFreq==0){
			migrationGroups <- migrationPolicy(nextGenGenoTable,genoSelectionTable,Policy,nFamilies,no_QTL)
			emigrantGroups <- migrationGroups[[1]]
			immigrantGroups<- migrationGroups[[2]]
		}
		
#####################################################################################################################
## Assign output variables for cycle n

	for(nFamily in 1:nFamilies){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,nCyc] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,nCyc] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,nCyc]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,nCyc] <- phenoValues_List[[nFamily]]

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
			genoSimSelectedValues <- genoValSimValues_List[[nFamily]][selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues_List[[nFamily]][selectedGenoIndividualIndices]
				
			selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
			attainedGenoValues_List[[nFamily]][nCyc] <- max(genoSimSelectedValues)
			
			if(selectionOnSimulated ==TRUE){
		        if(selectionOnGeno ==TRUE){
		
				GenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- as.vector(genoSelectionTable[[nFamily]][,1])
			    PhenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- phenoValSimValues_List[[nFamily]][selectedGenoIndividualIndices]
				
			   } else if(selectionOnGeno ==FALSE){ 

				GenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- genoValSimValues_List[[nFamily]][selectedGenoIndividualIndices]
				PhenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- as.vector(genoSelectionTable[[nFamily]][,1])
			   }			
			}
			if (selectionOnSimulated==FALSE){
			
				if(selectionOnGeno ==TRUE){
				
				GenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- as.vector(genoSelectionTable[[nFamily]][,1])
			        PhenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- phenoValues_List[[nFamily]][selectedGenoIndividualIndices]
								
				}else if(selectionOnGeno ==FALSE){
				
					GenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- genoValues_List[[nFamily]][selectedGenoIndividualIndices]
					PhenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- as.vector(genoSelectionTable[[nFamily]][,1])
				}
		    }
	}
		

   rm(nextGenGenoTable) 
   rm(newNextGenGenoTable)
   rm(nextGenGenoTable_List)


  }	
	
#########
     if(sclassPlot==TRUE){
	    simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List,Migration_List,PCA_List)
	 }else if(sclassPlot==FALSE){
	    simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List)
	 }

	  return(simResults_List)

 }


### Island Selection GM Frontier



runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3_Complete_nPar1 <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach){

	 
         options(warn=-1)
### Assign Variables #########################################################################
## Selection Parameters 

  	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues

	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno
	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny
	  
	  nFamilies <- 20
	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }

      nCrosses_inFamily <- nSel_inFamily
	  nProgeny <- nProgeny/nSel_inFamily

## QG Parameters	 
	  
	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList
	  varE <- varEList

	  no_QTL <- noQTLs
	  h2 <- H2
	  
	  nMarkers <- 4289
	  NAM_LinkMap_New <- NAM_LinkMap

	  condition <- i
	  Rep <- k  

      AlleleConversionTable_Combined <-  alleleConversionTable_Combined
 	
## Model Update Parameters	
	  modelType <- ModelType
      modelRetrain <- ModelRetrain
	  retrainFrequency <- RetrainFrequency
	  
	  modelUpdate <- ModelUpdate
      updateFrequency <- UpdateFrequency
	  updateType <- UpdateType

	  trainGenoNewTablePreCycle <- as.big.matrix(trainGenoNewTable)
      trainSimPhenoTablePreCycle <- trainSimPhenoValTable
      trainSimGenoValTablePreCycle <- trainSimGenoValTable

	  trainTableWindowSize <- TrainTableWindowSize
	  
## Weighting Parameters	  
	  Weighted <- weighted
	  WghtMethod <- wghtMethod

## BD Parameters	  
	  BD <- BreedDesign

## Migration Parameters	  
	  
	  Policy <- policy	
      migPolicyChangeFrequency <- policyChangeFrequency
	  migFreq <- migrationFrequency
	  direction <- Direction
	  migSize <- MigrationSize
	  migrationSize <- MigrationSize 
	  sClassPlot <- SClassPlot

## GM Parameters 
	  
	  GM_Param_List <- GM_Frontier_Param_List
	  no_cores <- GM_Frontier_Param_List$no_cores
	  
	  WorkspaceName <- WorkSpaceName
	  startCycle <- StartCycle
	  
	  nCores_FE <- noCores_ForEach
	  
	  cycle_nProgeny_Fam <- cycle1_nProgeny
	  
###########################################################################################
### Initialize variables

  if(startCycle ==1){
	
   	  Migration_List <- list()
	  PCA_List <- list()
	  j<-1
 	  selectedParents_Geno_List <- list() 
	  selectedParents_Pheno_List <- list()
	  selectedParents_NumProgeny_Geno_List <- list()
	  selectedParents_NumProgeny_Pheno_List <- list() 
	  
	  familyPairs_List <- list()
	  familyInd_List <- list()

	  combinedTableCount <- 1
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
# ###############################################################################################	  	  	  
# ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	
	#cl <- makeCluster(nCores_FE,outfile=paste("log_IM_GM",nCores_FE,sep=""))
	registerDoParallel(nCores_FE)
	
## Get genotype table from for cycle 1 genotype data
				
	nIndividuals <- nrow(cycle1GenoTable_GM)
	
    cycle1GenoTable_GM <-  generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny_Fam)
    newCycle1GenoTable_GM <- generateWeightedGenoTable(cycle1GenoTable_GM,no_QTL)
    cycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(cycle1GenoTable_GM)

## Get simulated geno and pheno values 	
 
    genoValSimValues <- simulateGenotypicValues(cycle1GenoTable_GM,no_QTL) 
    phenoValSimValues <- simulatePhenotypicValues(cycle1GenoTable_GM,varE,no_QTL) 

### Marker Effects, Frequency of Favorable allele 
	
      if(modelType=="BayesB" || modelType=="BL"){
	 if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[3])
	 }else if(selectionOnGeno==FALSE){
	    MarkerEffects <- unlist(PredictionModel_Pheno[3])
	 }
        }
	
	if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){
	 if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[1])
	 }else if(selectionOnGeno==FALSE){
	    MarkerEffects <- unlist(PredictionModel_Pheno[1])
	 }
	}
	
     Freq_List <- getFreq_BasePoln(cycle1GenoTable_GM,MarkerEffects)
	 Freq <- Freq_List[[1]]
	 FavAllele <- Freq_List[[2]]

### Get predicted geno and pheno values for unweighted and weighted genomic selection 

     if(modelType=="BayesB" || modelType=="BL"){

         if(Weighted==TRUE){
                genoValues  <- PredictBayesPhenoValues_Wght(newCycle1GenoTable_GM,PredictionModel_Geno,Freq)
                phenoValues <- PredictBayesPhenoValues_Wght(newCycle1GenoTable_GM,PredictionModel_Pheno,Freq)
         }else if(Weighted==FALSE){
                genoValues <- PredictBayesPhenoValues(newCycle1GenoTable_GM,PredictionModel_Geno)
                phenoValues <- PredictBayesPhenoValues(newCycle1GenoTable_GM,PredictionModel_Pheno)
         }
     }

    if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newCycle1GenoTable_GM,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newCycle1GenoTable_GM,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newCycle1GenoTable_GM,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newCycle1GenoTable_GM,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

            genoValues  <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM,PredictionModel_Geno)
            phenoValues <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM,PredictionModel_Pheno)

        }
    }

# # #######################################################################################################
# # ### Split cycle1  geno data in to families		

	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	 	  
	  for(nFamily in 1:20){ 
		
		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,] 
		
	  } 

# # ## Split Geno and Pheno values according to families  ########################


	 initIndex <- 1
	 finalIndex <- cycle1_nProgeny 
	 genoValues_List <- rep(list(list()),nFamilies)
     phenoValues_List <- rep(list(list()),nFamilies)
     genoValSimValues_List <- rep(list(list()),nFamilies)
     phenoValSimValues_List <- rep(list(list()),nFamilies)
 
	 cycle1_GenoTable_GM_List <- rep(list(matrix(rep(0,cycle1_nProgeny*nMarkers),nrow=cycle1_nProgeny,ncol=nMarkers)),nFamilies)
	 cycle1_GenoTableMod_GM_List <- rep(list(matrix(rep(0,cycle1_nProgeny*nMarkers),nrow=cycle1_nProgeny,ncol=nMarkers)),nFamilies)
	  
	 
	 for(nFamily in 1:nFamilies){ 
		
		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoValSimValues_List[[nFamily]] <- (genoValSimValues)[initIndex:finalIndex]
		phenoValSimValues_List[[nFamily]] <- (phenoValSimValues)[initIndex:finalIndex]
		
		cycle1_GenoTable_GM_List[[nFamily]] <- cycle1GenoTable_GM[initIndex:finalIndex,]
		cycle1_GenoTableMod_GM_List[[nFamily]] <- cycle1GenoTableMod_GM[initIndex:finalIndex,]
		
		initIndex <- finalIndex+1 
		finalIndex <- finalIndex+cycle1_nProgeny
	 
	 } 
	 
###### Select island pairs based on migration policy 
### Split criterion list 
	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- genoValSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- genoValues_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- phenoValSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- phenoValues_List
		} 
      } 



	families <-  c(1:nFamilies)
	
## Apply selection to cycle 1 population 
### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
			selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
            selectedPhenoIndividualIndices <- lapply(phenoSelectionTable,function(x) as.vector(x[,2]))
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
            selectedPhenoIndividualIndices <- lapply(phenoSelectionTable,function(x) as.vector(x[,2]))
	}
	
## selected Line indices to select lines from cycle genotable
    if(selectionOnGeno == TRUE) {
      selectedLineIndices <- selectedGenoIndividualIndices
    } else if(selectionOnGeno==FALSE){ 
      selectedLineIndices <- selectedPhenoIndividualIndices
    }   
	
	
### Migration Policy to identify immigration and emigration groups
	
	if(Policy != "GeneticClusters"){
	  migrationGroups <- migrationPolicy_GM_V2(cycle1GenoTable_GM,genoSelectionTable,Policy,nFamilies,no_QTL)
	
	  emigrantGroups <- migrationGroups[[1]]
	  immigrantGroups <- migrationGroups[[2]]
	
	}else if (Policy == "GeneticClusters"){ 
	
	  migrationGroups <- migrationPolicy_GM_V2(cycle1GenoTable_GM,genoSelectionTable,Policy,nFamilies,no_QTL)
	
	  emigrantGroups <- migrationGroups[[1]]
	  immigrantGroups <- migrationGroups[[2]]
	  emigrantGroup_Prob <- migrationGroups[[3]]
	
	}
	
	
### Exchange selected geno data among pairs ##############
    Cycle_Progeny_F5_List_Fam <- cycle1_GenoData_List

    if(migFreq !=0){
       
 	    if(nCyc==1 && nCyc%%migFreq==0){
	  
	    if(Policy != "GeneticClusters"){
		 
		 output_List_afterExchange <- exchangeGenoData_GM(Cycle_Progeny_F5_List,genoValSimValues_List,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
		}else if(Policy == "GeneticClusters"){ 
		
		 
		 output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		}
		  
		Cycle_Progeny_F5_List <- (output_List_afterExchange)[[1]]
		
		familyInfo <- output_List_afterExchange[[2]]
		
		Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 
	   
        Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,cycle1_nProgeny))
        nProgIndex <-  dim(Cycle_Progeny_F5_List[[nFamily]])[3]
	    nProg <-1
				
	    for(nFamily in 1:length(families)){ 
					 
		    Cycle_Progeny_F5[nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_List[[nFamily]][,,nProg:nProgIndex] 
	    } 
	
    }
	   
	
### i) Get cycle 1 genotype table and criterion value list for sets of singleislands and pairs of islands	


### Get list of cycle1 geno table for family ind and family pairs ######################################

		
	cycle1GenoTable_GM_List <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
	cycle1GenoTableMod_GM_List <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
       
	
	initIndex <- 1
	finalIndex <- cycle1_nProgeny_Fam
	
	for(nFam in 1:nFamilies){ 
		
		cycle1GenoTable_GM_List[[nFam]] <- cycle1GenoTable_GM[initIndex:finalIndex,]
		cycle1GenoTableMod_GM_List[[nFam]] <- cycle1GenoTableMod_GM[initIndex:finalIndex,]
            
	    initIndex <- finalIndex+1
		finalIndex <- finalIndex + cycle1_nProgeny_Fam
		
	}

######################################################################################   

#### GM method for parent selection in cycle 1 ################################################
# #### GM method for isolated families
    nParameters <- 30
 			 
 	selectedParents_Geno_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles) 
	CriterionValue_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  system.time({

    genoSelected_List <- foreach(k=1:nFamilies) %dopar% (getIM_GM_Frontier_selectedGenoList_SingleFams_V2(cycle1GenoTable_GM_List[[k]],cycle1GenoTableMod_GM_List[[k]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))
  })
 
   for(nFamily in 1:nFamilies){
	  for(nParameter in 1:nParameters){
		selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <- genoSelected_List[[nFamily]][[1]][[nParameter]]
		CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- genoSelected_List[[nFamily]][[2]][[nParameter]]
		selectedParents_Geno <- genoSelected_List[[nFamily]][[1]][[nParameter]]
      }
   }


#######################################################################################

    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
	

	Cycle_GenoData_List <- cycle1_GenoData_List
## Get F5 RILs for selected parents in cycle 1
   
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% (getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(),nFamilies)),nParameters)
    nCrosses_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
		

	for(nFamily in 1:nFamilies){
          for(nParameter in 1:nParameters){
          	Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]] <- F5RILs_GM_Fam[[nFamily]][[nParameter]][[1]]
          	nCrosses_List[[nFamily]][[nParameter]] <- F5RILs_GM_Fam[[nFamily]][[nParameter]][[2]]
 	
         }
    }


    gc()
 
## Initiate Output list for all families 	  
	  
	  GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
	  GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
	  GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
	
	  PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
	  PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
	  PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

      attainedPhenoValues_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
      attainedGenoValues_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
	
	  topNGenoValues_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
      GenoTableIndices_topN_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),,nCycles)

      topNPhenoValues_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
      PhenoTableIndices_topN_List_Cycle <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
	  
	
 	 
##
for(nFamily in 1:nFamilies){
 
   for(nParameter in 1:nParameters){
   
        uniqueParentIDs <-  unique(as.vector(selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues_List[[nFamily]][uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues_List[[nFamily]][uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues_List[[nFamily]][uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues_List[[nFamily]][uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs

		
### GM Output in selected individual ids 
############################

        GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- genoValSimValues_List[[nFamily]]
		
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]]<- genoValues_List[[nFamily]] 
	    attainedGenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- max(genoValSimValues_List[[nFamily]][GenoTableIndices_topN])
		
	    topNGenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- topNGenoValues_List
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- GenoTableIndices_topN_List
	    GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]]
	   
	    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- phenoValSimValues_List[[nFamily]]
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <-phenoValues_List[[nFamily]] 
	    attainedPhenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <-  max(phenoValSimValues_List[[nFamily]][PhenoTableIndices_topN])
		
	    topNPhenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <-  PhenoTableIndices_topN
	    PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]]
		
		
	}
}

###########################################################################
if(sClassPlot== TRUE){

## Sclass plot code and get PG parameters for islands	
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}

# ### get island PG parameters

		cycle1GenoTable_AlleleFormat_GM_FamPair <- getAlleleFormat(cycle1GenoTable_GM_FamPair,AlleleConversionTable_Combined)
 		cycle1GenoTable_GInd_GM_FamPair <- df2genind(cycle1GenoTable_AlleleFormat_GM_FamPair,pop=populations,sep="")
      			
		gc()
		
# #################################################################################	 
	

    		tol <- 1e-07
			cycle1_df <- as.matrix(cycle1GenoTable_GM)
			cycle1_df <- tcrossprod(cycle1_df,cycle1_df) 
			eig1 <- eigen(cycle1_df,symmetric=TRUE) 
			eig <- eig1$values
			rank_X <- sum((eig/eig[1]) > tol)
			NF_Vec <- as.numeric(which(cumsum(eig/sum(eig)) >= 0.99)) 
			NF <- NF_Vec[1]
			
						
			if(NF == dim(cycle1_df)[1]) { 
				NF = rank_X 
			}else if (NF==1){
				NF <-rank_X
			}else if( ! exists("NF")) {
				NF <-rank_X 
			}else if(rank_X== 1){
				NF <- 2
			}
			
			
	       
	        if(rank_X > 0 && length(NF)>0){
		    
				j<- j+1
				
			    scaledGInd <- scaleGen(cycle1GenoTable_GInd,scale=FALSE)
				pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF)
		        
				tryCatch(pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF), error=function(e) paste("Error dudiPCA Cycle No",i,"-",e))
				
		   		
				filename <- paste("Sclass_Plot","-",condition,"-",Rep,"-",nCyc,".png",sep="")
				
				png(filename,width=1024,height=768,pointsize=20)
				tryCatch(s.class(pca.GInd$li, fac=pop(cycle1GenoTable_GInd),col=funky(20),addaxes=TRUE),error=function(e) print(paste("Error Sclass Cycle No",i,"-",e)))
								
				tryCatch(title(paste("Cycle - ",i,sep="")),error=function(e) print(paste("Error Title Cycle No",i,"-",e)))
				dev.off()
				tryCatch(system(paste("gzip -f ",filename,sep="")),error=function(e) print(paste("Error gzip Cycle No",i,"-",e)))
				
				tryCatch(PCA_Li_Table <- as.big.matrix(as.matrix(pca.GInd$li)),error=function(e) print(paste("Error BigMatrix Cycle No",i,e)))
				tryCatch(PCA_C1_Table <- as.big.matrix(as.matrix(pca.GInd$c1)),error=function(e) print(paste("Error BigMatrix Cycle No",i,e)))
				
				tryCatch(PCA_List[[j]] <- list(PCA_Li_Table,PCA_C1_Table),error=function(e)print(paste("Error PCA_List Cycle No",i,e)))
			
			
			}
			rm(cycle1GenoTable_GInd)
			rm(scaledGInd)
			rm(pca.GInd)
  }
#######################################################

	startCycle <- startCycle +1 
	save.image(WorkspaceName) 
	gc() 

}

###########################################################################

  startCycle <- 2
  nextGenGenoTable_GM_List <- rep(list(rep(list(),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(),nParameters)),nFamilies)
  PredictionModel_Geno_List <- rep(list(rep(list(),nParameters)),nFamilies)
  PredictionModel_Pheno_List <- rep(list(rep(list(),nParameters)),nFamilies)
  genoValues_List <- rep(list(rep(list(),nParameters)),nFamilies)
  phenoValues_List <- rep(list(rep(list(),nParameters)),nFamilies)
  genoValSimValues_List <- rep(list(rep(list(),nParameters)),nFamilies)
  phenoValSimValues_List <- rep(list(rep(list(),nParameters)),nFamilies)
 
 nCyc <- 2
  
 for(nCyc in startCycle:nCycles){
   
    print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
    for(nFamily in 1:nFamilies){
    
	  for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <-   Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]]
		nCrosses <-  nCrosses_List[[nFamily]][[nParameter]] 
		
		nextGenGenoTable_GM <- generateMinus101GenoFormat_IM_GM_Frontier_V1(Cycle_Progeny_F5,nCrosses)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)

        nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTable_GM
        nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableMod_GM
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <- newNextGenGenoTable_GM 
		
		
		genoValSimValues_List[[nFamily]][[nParameter]] <- genoValSimValues 
		phenoValSimValues_List[[nFamily]][[nParameter]] <- phenoValSimValues 

	  }
	}
	
	
	PredictionModel_Geno_Parameter_List <- list() 
	PredictionModel_Pheno_Parameter_List <- list()
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	for(nParameter in 1:nParameters){ 
	
	   genoValSimValues_Fam <- c()
	   phenoValSimValues_Fam <- c()
	   nextGenGenoTable_GM_Fam <- c()
	   nextGenGenoTableMod_GM_Fam <- c()
	   newNextGenGenoTable_GM_Fam <- c() 
	   
	    for(nFamily in 1:nFamilies){
	
			nextGenGenoTable_GM_Fam <- rbind(nextGenGenoTable_GM_Fam,nextGenGenoTable_GM_List[[nFamily]][[nParameter]])
		
			nextGenGenoTableMod_GM_Fam <- rbind(nextGenGenoTableMod_GM_Fam,nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]]) 
			newNextGenGenoTable_GM_Fam <- rbind(newNextGenGenoTable_GM_Fam,newNextGenGenoTable_GM_List[[nFamily]][[nParameter]]) 
			genoValSimValues_Fam <- c(genoValSimValues_Fam,genoValSimValues_List[[nFamily]][[nParameter]])
			phenoValSimValues_Fam <- c(phenoValSimValues_Fam,phenoValSimValues_List[[nFamily]][[nParameter]]) 
        }
############################################################################################### 
        if(modelUpdate ==TRUE && nCyc%%updateFrequency==0 && nParameter ==1){

			PredictionModels <- getUpdatedGSModel(newNextGenGenoTable_GM_Fam,genoValSimValues_Fam,phenoValSimValues_Fam,PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,trainGenoNewTablePreCycle,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount) 
 
			PredictionModel_Geno <- PredictionModels[[1]]
			PredictionModel_Pheno <- PredictionModels[[2]]
 
	    }
		
		PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModel_Geno
        PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModel_Pheno
 
#################get predicted geno and pheno values with JM and unweighted GP models
### Get Predicted geno & pheno values for nextGenGenoData	

	    if(modelType=="BayesB" || modelType=="BL"){
		
		   if(Weighted==TRUE && WghtMethod=="JM"){
			genoValues  <- PredictBayesPhenoValues_Wght(newNextGenGenoTable_GM,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq)
	        }else if(Weighted==TRUE && WghtMethod =="DW"){
		    genoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)
		    phenoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

		    genoValues <- genoValues_PredList[[1]]
		    phenoValues <- phenoValues_PredList[[1]]

		    alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]
				
		}else if(Weighted==FALSE){
		
		    genoValues_FamPair <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
			phenoValues_FamPair <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
			
			genoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
			
			genoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
			phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)

	    }	
	    }


	    if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM_Fam,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM_Fam,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

            				
			genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
		
        }
     }

  	
			genoValues_Parameter_List[[nParameter]] <- genoValues
			phenoValues_Parameter_List[[nParameter]] <- phenoValues
}

####
   	 
    for(nParameter in 1:nParameters){ 
	
	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny

	  for(nFamily in 1:nFamilies){
	  
	    genoValues_List[[nFamily]][[nParameter]] <- genoValues_Parameter_List[[nParameter]][initIndex:finalIndex]
		phenoValues_List[[nFamily]][[nParameter]] <- phenoValues_Parameter_List[[nParameter]][initIndex:finalIndex]
		initIndex <- finalIndex+1 
		finalIndex <- finalIndex + cycle1_nProgeny
	 
	  }
	}
	
### Apply selection to cycle 1 population 
### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete_Frontier(genoValSimValues_List,no_selected,nFamilies,nParameters)
			phenoSelectionTable <- ApplySelection_Discrete_Frontier(phenoValSimValues_List,no_selected,nFamilies,nParameters)
			genoSimSelectionTable <- genoSelectionTable
			selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) lapply(x,function(y) as.vector(y[,2]))) 
            selectedPhenoIndividualIndices <- lapply(phenoSelectionTable,function(x) lapply(x,function(y) as.vector(y[,2])))
		


	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete_Frontier(genoValues_List,no_selected,nFamilies,nParameters)
			phenoSelectionTable <- ApplySelection_Discrete_Frontier(phenoValues_List,no_selected,nFamilies,nParameters)
			genoSimSelectionTable <- ApplySelection_Discrete_Frontier(genoValSimValues_List,no_selected,nFamilies,nParameters)
			selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) lapply(x,function(y) as.vector(y[,2])))
            selectedPhenoIndividualIndices <- lapply(phenoSelectionTable,function(x) lapply(x,function(y) as.vector(y[,2])))
	}
	
## selected Line indices to select lines from cycle genotable
    if(selectionOnGeno == TRUE) {
	    selectedLineIndices <- selectedGenoIndividualIndices
    } else if(selectionOnGeno==FALSE){ 
        selectedLineIndices <- selectedPhenoIndividualIndices
    }  


######################################################################################   
## get correct table format for migration and exchange genodata functions

      Cycle_Progeny_F5_List_Fam_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
      nextGenGenoTable_GM_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters) 
      genoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      genoValSimValues_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)

      for(nFamily in 1:nFamilies){
    
	    for(nParameter in 1:nParameters){
	
	        Cycle_Progeny_F5 <-   Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]]
	 
	        Cycle_Progeny_F5_List_Fam_Reverse[[nParameter]][[nFamily]] <- Cycle_Progeny_F5
                genoValSimValues_List_Reverse[[nParameter]][[nFamily]] <- genoValSimValues_List[[nFamily]][[nParameter]]

			         
	    }
	}
 
       
     for(nParameter in 1:nParameters){

	   for(nFamily in 1:nFamilies){

        	nextGenGenoTable_GM_List_Reverse[[nParameter]][[nFamily]] <- nextGenGenoTable_GM_List[[nFamily]][[nParameter]] 

                genoSelectionTable_Parameter[[nParameter]][[nFamily]] <- genoSelectionTable[[nFamily]][[nParameter]]

           }
     }

################################################

        Cycle_Progeny_F5_Fam_ParamList <- rep(list(array(0,c(nFamilies,nMarkers,2,cycle1_nProgeny))),nParameters)

        Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)),nParameters)
        
	nProgInit <- 1
	nProgFinal <-  dim(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[1]])[4]
        nCrosses_inFamily <- nSel_inFamily

        nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
    for(nParameter in 1:nParameters){
	 		
	   for(nFamily in 1:length(families)){ 

              nProg <- nProgInit 
	      nProgIndex <- nProgFinal

              for(nCross in 1:nCrosses_inFamily){
					 
		        Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,nProgInit:nProgFinal] 
                        Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]][,,nProg:nProgIndex]<- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,nProgInit:nProgFinal] 
		    
		        nProg <- nProgIndex+1 
                        nProgIndex <- nProgIndex + nProgFinal


			}
	    } 
	}
	


#################################################################################################################
### Migration Policy to identify immigration and emigration groups ##############################################

	
	if(Policy != "GeneticClusters"){
	       migrationGroups <- migrationPolicy_GM_V2(nextGenGenoTable_GM_List_Reverse[[nParameter]],genoSelectionTable_Parameter[[nParameter]],Policy,nFamilies,no_QTL)
	
	       emigrantGroups <- migrationGroups[[1]]
	       immigrantGroups <- migrationGroups[[2]]
	
	}else if (Policy == "GeneticClusters"){ 
	
	       migrationGroups <- migrationPolicy_GM_V2(nextGenGenoTable_GM_List_Reverse[[nParameter]],genoSelectionTable_Parameter[[nParameter]],Policy,nFamilies,no_QTL)
	
	       emigrantGroups <- migrationGroups[[1]]
	       immigrantGroups <- migrationGroups[[2]]
	       emigrantGroup_Prob <- migrationGroups[[3]]
	
	}
	
### Exchange selected geno data among pairs ##############
     
################################

    if(migFreq !=0){
   
 	    if(nCyc>=2 && nCyc%%migFreq==0){
	  
	         if(Policy != "GeneticClusters"){
		 
			output_List_afterExchange <- foreach(k=1:nParameters) %dopar% (exchangeGenoData_GM(Cycle_Progeny_F5_Fam_CrossList[[k]],genoValSimValues_List_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
		  }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		  }
	  
		  Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
		
		  familyInfo <- output_List_afterExchange[[2]]
		
		  Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 
	}   
	   
      
############################################################################################### 
  
	 if(selectionOnGeno==TRUE){ 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- genoValSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- genoValues_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- phenoValSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- phenoValues_List
		} 
	 }
#### GM method for isolated families

 system.time({

	genoSelected_List <- foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% 
(getIM_GM_Frontier_selectedGenoList_SingleFams_V2(nextGenGenoTable_GM_List[[k]][[j]],nextGenGenoTableMod_GM_List[[k]][[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))

 })


## Reduce Parameter set to 30  

 
  selectedParents_Geno_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
  CriterionValue_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  for(nFamily in 1:nFamilies){
    for(nParameter in 1:nParameters){
   
      selectedParents_Geno_Fam_List[[nCyc]][[nFamily]][[nParameter]] <- genoSelected_List[[nFamily]][[nParameter]][[1]]
      CriterionValue_Fam_List[[nCyc]][[nFamily]][[nParameter]] <-genoSelected_List[[nFamily]][[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nFamily]][[nParameter]][[1]]
      
    }
  }  
 
   
    
    selectedParents_Geno_Par_List <- rep(list(rep(list(list()),nParameters*nParameters)),nFamilies)
    GainValues <- rep(list(c()),nFamilies)
    UsefulnessValues <- rep(list(c()),nFamilies)
    InbreedValues <- rep(list(c()),nFamilies)
   
    for(nFamily in 1:nFamilies){ 


      count <- 1
     
     for (nPar1 in 1: nParameters){ 
 
      for(nPar2 in 1:nParameters){ 
	    GainValues[[nFamily]][count] <- CriterionValue_Fam_List[[nCyc]][[nFamily]][[nPar1]][[nPar2]][1]
	    UsefulnessValues[[nFamily]][count] <- CriterionValue_Fam_List[[nCyc]][[nFamily]][[nPar1]][[nPar2]][2]
	    InbreedValues[[nFamily]][count] <- CriterionValue_Fam_List[[nCyc]][[nFamily]][[nPar1]][[nPar2]][3]
	    selectedParents_Geno_Par_List[[nFamily]][[count]]  <-  selectedParents_Geno_Fam_List[[nCyc]][[nFamily]][[nPar1]][[nPar2]]
	 
	    count <- count+1
	  }
    }
    
      
    Quartile3 <- summary(GainValues[[nFamily]])[5]
    Quartile1 <- summary(GainValues[[nFamily]])[2]
			
    indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
    indicesRange1 <- indicesSortRange1[[2]][1:10]
			
    indices2 <- which((GainValues[[nFamily]] <= Quartile3) & (GainValues[[nFamily]] >= Quartile1))
    indicesRange2 <- sample(indices2,10) 
			
    indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
    indicesRange3 <- indicesSortRange3[[2]][1:10]
			
    indicesRange <- c(indicesRange1,indicesRange2,indicesRange3)
    nParameters <- length(indicesRange)

    			
    for(nParameter in 1:nParameters){
                nParameterIndex <- indicesRange[nParameter]
	        
			   selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <-  selectedParents_Geno_Par_List[[nFamily]][[nParameterIndex]]
			   CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- c(GainValues[[nFamily]][nParameterIndex],UsefulnessValues[[nFamily]][nParameterIndex],InbreedValues[[nFamily]][nParameterIndex])
		
	} 
   }

  save.image(WorkspaceName)
  gc() 	
   	
###########################################################################
## Get F5 RILs for selected parents in cycle 1

        Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
	nProgInit <- 1
	nProgFinal <-  dim(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[1]])[4]
        nCrosses_inFamily <- nSel_inFamily

         nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	   for(nFamily in 1:length(families)){ 

           nProg <- nProgInit 
	       nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					 
		        Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,nProgInit:nProgFinal] 
		    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
	    } 
	}

######################################################################################################## 

    Cycle_GenoData_List <- Cycle_Progeny_F5_Fam_CrossList
		
 
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% (getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]][[j]],selectedParents_Geno_Fam_List[[nCyc]][[k]][[j]][[1]],selectionOnGeno,nMarkers,NAM_LinkMap_New))

	for(nFamily in 1:nFamilies){
        
           for(nParameter in 1:nParameters){

          	Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]] <- F5RILs_GM_Fam[[nFamily]][[nParameter]][[1]]
          	nCrosses_List[[nFamily]][[nParameter]] <- F5RILs_GM_Fam[[nFamily]][[nParameter]][[2]]
 	   }

    }

    gc()

 
########################################################################
### Output variables for isolated families

for(nFamily in 1:nFamilies){
 
   for(nParameter in 1:nParameters){
   
        uniqueParentIDs <- unique(as.vector(selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues_List[[nFamily]][[nParameter]][uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues_List[[nFamily]][[nParameter]]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues_List[[nFamily]][[nParameter]][uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues_List[[nFamily]][[nParameter]][uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs

        topNGenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <-  topNGenoValues
        topNPhenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- topNPhenoValues
		
### GM Output in selected individual ids 
############################

	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <-  genoValSimValues_List[[nFamily]][[nParameter]]
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]]<- genoValues_List[[nFamily]][[nParameter]]
	    attainedGenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- max(genoValSimValues_List[[nFamily]][[nParameter]][GenoTableIndices_topN])
		
	    topNGenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- GenoTableIndices_topN_List
	    GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]]
	   
	    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <-  phenoValSimValues_List[[nFamily]][[nParameter]]
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <-  phenoValues_List[[nFamily]][[nParameter]]
	    attainedPhenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <-  max(phenoValSimValues_List[[nFamily]][[nParameter]][PhenoTableIndices_topN])
		
	    topNPhenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <-  PhenoTableIndices_topN
	    PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nFamily]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nFamily]][[nParameter]]
		
		
	}
 }

###########################################################################

	print(paste(nCyc,"-Fam20-",mean(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]]),sep=""))

		 
## Sclass plot for family clustering 	
####################################################################################################

 if(sClassPlot == TRUE){
		
	if(i%%5 ==0){	

###  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}

### get island PG parameters for family pairs

		nextGenGenoTable_AlleleFormat_FamPair <- getAlleleFormat(nextGenGenoTable_GM_FamPair,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat_FamPair,pop=populations,sep="")
      	
		gc()
		
#################################################################################		
	 				
			tol <- 1e-07
			nextGen_df_FP <- as.matrix(nextGenGenoTable_GM_FamPair)
			nextGen_df_FP <- tcrossprod(nextGen_df_FP,nextGen_df_FP) 
			eig1 <- eigen(nextGen_df_FP,symmetric=TRUE) 
			eig <- eig1$values
			rank_X <- sum((eig/eig[1]) > tol)
			NF_Vec <- as.numeric(which(cumsum(eig/sum(eig)) >= 0.99)) 
			NF <- NF_Vec[1]
			
						
			if(NF == dim(nextGen_df)[1]) { 
				NF = rank_X 
			}else if (NF==1){
				NF <-rank_X
			}else if( ! exists("NF")) {
				NF <-rank_X 
			}else if(rank_X== 1){
				NF <- 2
			}
			
			
	       
	        if(rank_X > 0 && length(NF)>0){
		    
				j<- j+1
				
			    scaledGInd <- scaleGen(nextGenGenoTable_GInd,scale=FALSE)
				pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF)
		        
				tryCatch(pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF), error=function(e) paste("Error dudiPCA Cycle No",i,"-",e))
				
		   		
				filename <- paste("Sclass_Plot_FamPair","-",condition,"-",Rep,"-",nCyc,".png",sep="")
				
				png(filename,width=1024,height=768,pointsize=20)
				tryCatch(s.class(pca.GInd$li, fac=pop(nextGenGenoTable_GInd),col=funky(20),addaxes=TRUE),error=function(e) print(paste("Error Sclass Cycle No",i,"-",e)))
								
				tryCatch(title(paste("Cycle - ",i,sep="")),error=function(e) print(paste("Error Title Cycle No",i,"-",e)))
				dev.off()
				tryCatch(system(paste("gzip -f ",filename,sep="")),error=function(e) print(paste("Error gzip Cycle No",i,"-",e)))
				
				tryCatch(PCA_Li_Table <- as.big.matrix(as.matrix(pca.GInd$li)),error=function(e) print(paste("Error BigMatrix Cycle No",i,e)))
				tryCatch(PCA_C1_Table <- as.big.matrix(as.matrix(pca.GInd$c1)),error=function(e) print(paste("Error BigMatrix Cycle No",i,e)))
				
				tryCatch(PCA_List[[j]] <- list(PCA_Li_Table,PCA_C1_Table),error=function(e)print(paste("Error PCA_List Cycle No",i,e)))
			
			
			}
			rm(nextGenGenoTable_GInd)
			rm(scaledGInd)
			rm(pca.GInd)
	
### Sclass plots for single families
	
	        tol <- 1e-07
			nextGen_df<- as.matrix(nextGenGenoTable_GM_Fam)
			nextGen_df <- tcrossprod(nextGen_df,nextGen_df) 
			eig1 <- eigen(nextGen_df,symmetric=TRUE) 
			eig <- eig1$values
			rank_X <- sum((eig/eig[1]) > tol)
			NF_Vec <- as.numeric(which(cumsum(eig/sum(eig)) >= 0.99)) 
			NF <- NF_Vec[1]
			
						
			if(NF == dim(nextGen_df)[1]) { 
				NF = rank_X 
			}else if (NF==1){
				NF <-rank_X
			}else if( ! exists("NF")) {
				NF <-rank_X 
			}else if(rank_X== 1){
				NF <- 2
			}
			
			
	       
	        if(rank_X > 0 && length(NF)>0){ 		    
				j<- j+1
				
			    scaledGInd <- scaleGen(nextGenGenoTable_GInd,scale=FALSE)
				pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF)
		        
				tryCatch(pca.GInd <- dudi.pca(scaledGInd,center=FALSE,scale=FALSE,scannf=FALSE,nf=NF), error=function(e) paste("Error dudiPCA Cycle No",i,"-",e))
				
		   		
				filename <- paste("Sclass_Plot_Fam","-",condition,"-",Rep,"-",nCyc,".png",sep="")
				
				png(filename,width=1024,height=768,pointsize=20)
				tryCatch(s.class(pca.GInd$li, fac=pop(nextGenGenoTable_GInd),col=funky(20),addaxes=TRUE),error=function(e) print(paste("Error Sclass Cycle No",i,"-",e)))
								
				tryCatch(title(paste("Cycle - ",i,sep="")),error=function(e) print(paste("Error Title Cycle No",i,"-",e)))
				dev.off()
				tryCatch(system(paste("gzip -f ",filename,sep="")),error=function(e) print(paste("Error gzip Cycle No",i,"-",e)))
				
				tryCatch(PCA_Li_Table <- as.big.matrix(as.matrix(pca.GInd$li)),error=function(e) print(paste("Error BigMatrix Cycle No",i,e)))
				tryCatch(PCA_C1_Table <- as.big.matrix(as.matrix(pca.GInd$c1)),error=function(e) print(paste("Error BigMatrix Cycle No",i,e)))
				
				tryCatch(PCA_List[[j]] <- list(PCA_Li_Table,PCA_C1_Table),error=function(e)print(paste("Error PCA_List Cycle No",i,e)))
			
			
			}
			rm(nextGenGenoTable_GInd)
			rm(scaledGInd)
	     	rm(pca.GInd)
	
	}
   }		 
	 
}
 
######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle)

	   return(simResults_List)

}
