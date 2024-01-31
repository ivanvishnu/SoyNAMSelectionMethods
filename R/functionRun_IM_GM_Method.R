### Simulations for IM GM Method 
 

runSimulations20X_IslandSelection_BD_Wght_GM <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Param_List){

	 
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
	 
	  sClassPlot <- SClassPlot

## GM Parameters 
	  
	  GM_Param_List <- GM_Param_List
	  no_cores <- GM_Param_List$no_cores
	  
	  WorkspaceName <- WorkSpaceName
	  startCycle <- StartCycle

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
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <-  generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable_GM <- generateWeightedGenoTable(cycle1GenoTable_GM,no_QTL)
	cycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(cycle1GenoTable_GM)
	
		
	nIndividuals <- nrow(cycle1GenoTable_GM)
	 

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
# # ### Split geno and pheno data tables in to families		
	 
	  cycle1_nProgeny <- 100
	  nProgeny <- 100
	  
	  nSel_inFamily <- 1
	  
	  F5_Progeny_List_Family<- F5_Progeny_List
	  
	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  
	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)
	  
	  for(nFamily in 1:20){ 
		
		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,] 
	  
	  }
	  
# ####################		
	 	  
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	 	  
	  for(nFamily in 1:20){ 
		
		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,] 
		
	  } 
	  
	  
# # ## Split Geno and Pheno values according to families  ########################


	 initIndex <- 1
	 finalIndex <- cycle1_nProgeny 
	 genoValues_List <- list()
	 phenoValues_List <- list()
	 genoSimValues_List <- list()
	 phenoSimValues_List <- list()
	 cycle1_GenoTable_GM_List <- rep(list(matrix(rep(0,cycle1_nProgeny*nMarkers),nrow=cycle1_nProgeny,ncol=nMarkers)),nFamilies)
	 	  
	 cycle1_GenoTableMod_GM_List <- rep(list(matrix(rep(0,cycle1_nProgeny*nMarkers),nrow=cycle1_nProgeny,ncol=nMarkers)),nFamilies)
	  
	 
	 for(nFamily in 1:nFamilies){ 
		
		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]
		
		cycle1_GenoTable_GM_List[[nFamily]] <-cycle1GenoTable_GM[initIndex:finalIndex,]
		cycle1_GenoTableMod_GM_List[[nFamily]] <-cycle1GenoTableMod_GM[initIndex:finalIndex,]
		
		initIndex <- finalIndex+1 
		finalIndex <- finalIndex+cycle1_nProgeny
	 
	 } 

	 
	
###### Select island pairs based on migration policy 
### Split criterion list 
	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- genoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- genoValues_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- phenoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- phenoValues_List
		} 
     } 


### i)  Get pairs of families and isolated families based on migration policy 
 
	
	familyPairs <-  migrationPolicy_GM(cycle1GenoTable_GM,SelCriterion_List,Policy,nFamilies)
	nFamilyPairs <- dim(familyPairs)[1]
	nUniqueFam_in_Pairs <-   length(unique(as.vector(familyPairs)))
	uniqueFamilies_in_Pairs <- unique(as.vector(familyPairs))
	
	families <-  setdiff(c(1:nFamilies),unique(as.vector(familyPairs)))
    nFamilyInd <- length(families) 
	
	
	familyPairs_List[[nCyc]] <- familyPairs
	familyInd_List[[nCyc]] <- families
	
	
	
	cycle_nProgeny_FamPair <- cycle1_nProgeny*2
	cycle_nProgeny_Fam <- cycle1_nProgeny
	
	selectedParents_Geno_FamPair_List <- rep(list(rep(list(),nFamilyPairs)),nCycles)
	selected_Parents_Geno_FamPair_Ind_List <-  rep(list(rep(list(),nFamilyPairs)),nCycles)
	selectedParents_NumProgeny_Geno_FamPair_List <- rep(list(rep(list(),nFamilyPairs)),nCycles)
	 
 	selectedParents_Geno_Fam_List <- rep(list(rep(list(),nFamilyInd)),nCycles) 
	selectedParents_NumProgeny_Geno_Fam_List <- rep(list(rep(list(),nFamilyInd)),nCycles)
	

### ii) Split cycle 1 genotype data in to sets for singleislands and pairs of islands	
 
	 
	cycle1_GenoData_FamPair <- array(0,c(nFamilyPairs,nMarkers,2,cycle1_nProgeny_FamPair)) 
	cycle1_GenoData_Fam <- array(0,c(nFamilyInd,nMarkers,2,cycle1_nProgeny_Fam)) 
	#cycle1_GenoData_Fam_in_Pair <- array(0,c(nUniqueFam_in_Pairs,nMarkers,2,cycle1_nProgeny_Fam)) 
	
	cycle1_GenoData_List_FamPair <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny_FamPair))),nFamilyPairs) 
	cycle1_GenoData_List_Fam <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny_Fam))),nFamilyInd) 
	#cycle1_GenoData_List_Fam_in_FamPair <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny_Fam))),nUniqueFam_in_Pairs) 

	for(nFamPairs in 1:nFamilyPairs){

           fam1 <- familyPairs[nFamPairs,1]
           fam2 <- familyPairs[nFamPairs,2]

           cycle1_GenoData_List_FamPair[[nFamPairs]][,,1:cycle1_nProgeny_Fam]<- cycle1_Geno_Data[fam1,,,]
           cycle1_GenoData_List_FamPair[[nFamPairs]][,,(cycle1_nProgeny_Fam+1):cycle1_nProgeny_FamPair] <- cycle1_Geno_Data[fam2,,,]
		  
		   cycle1_GenoData_FamPair[nFamPairs,,,1:cycle1_nProgeny_Fam] <- cycle1_Geno_Data[fam1,,,] 
		   cycle1_GenoData_FamPair[nFamPairs,,,(cycle1_nProgeny_Fam+1):cycle1_nProgeny_FamPair] <- cycle1_Geno_Data[fam2,,,] 
    }

	for(nFamily in 1:nFamilyInd){

             fam <- families[nFamily]

             cycle1_GenoData_List_Fam[[nFamily]]<- cycle1_Geno_Data[fam,,,]
             cycle1_GenoData_Fam[nFamily,,,]<- cycle1_Geno_Data[fam,,,]
	}
	
	# for(nUniqueFamily in 1: nUniqueFam_in_Pairs){
	
	         # fam <- uniqueFamilies_in_Pairs[nUniqueFamily]
	         # cycle1_GenoData_List_Fam_in_FamPair[[nUniqueFamily]] <- cycle1_Geno_Data[fam,,,]
	         # cycle1_GenoData_Fam_in_Pair[nUniqueFamily,,,] <- cycle1_Geno_Data[fam,,,]
	
	# }
	
### iii) Get cycle 1 genotype table and criterion value list for sets of singleislands and pairs of islands	
	
    cycle1GenoTable_GM_FamPair<-  generateMinus101GenoFormat_V1(cycle1_GenoData_FamPair,nFamilyPairs,cycle1_nProgeny_FamPair)
    newCycle1GenoTable_GM_FamPair <- generateWeightedGenoTable(cycle1GenoTable_GM_FamPair,no_QTL)
    cycle1GenoTableMod_GM_FamPair <- generate012GenoFormat_Inverse(cycle1GenoTable_GM_FamPair)
	
    cycle1GenoTable_GM_Fam <-  generateMinus101GenoFormat_V1(cycle1_GenoData_Fam,nFamilyInd,cycle1_nProgeny_Fam)
    newCycle1GenoTable_GM_Fam<- generateWeightedGenoTable(cycle1GenoTable_GM_Fam,no_QTL)
    cycle1GenoTableMod_GM_Fam <- generate012GenoFormat_Inverse(cycle1GenoTable_GM_Fam)

    # cycle1GenoTable_GM_Fam_in_Pair<-  generateMinus101GenoFormat_V1(cycle1_GenoData_Fam_in_Pair,nUniqueFam_in_Pairs,cycle1_nProgeny_Fam) 
    # newCycle1GenoTable_GM_Fam_in_Pair <- generateWeightedGenoTable(cycle1GenoTable_GM_Fam_in_Pair,no_QTL)
	# cycle1GenoTableMod_GM_Fam_in_Pair<- generate012GenoFormat_Inverse(cycle1GenoTable_GM_Fam_in_Pair)

    genoValSimValues_FamPair <- simulateGenotypicValues(cycle1GenoTable_GM_FamPair,no_QTL)
    phenoValSimValues_FamPair <- simulatePhenotypicValues(cycle1GenoTable_GM_FamPair,varE,no_QTL)

    genoValSimValues_Fam <- simulateGenotypicValues(cycle1GenoTable_GM_Fam,no_QTL) 
    phenoValSimValues_Fam <- simulatePhenotypicValues(cycle1GenoTable_GM_Fam,varE,no_QTL) 

# genoValSimValues_Fam_in_Pair <- simulateGenotypicValues(cycle1GenoTable_GM_Fam_in_Pair,no_QTL) 
# phenoValSimValues_Fam_in_pair <- simulatePhenotypicValues(cycle1GenoTable_GM_Fam_in_Pair,varE,no_QTL) 


if(modelType=="BayesB" || modelType=="BL"){
		
		if(Weighted==TRUE && WghtMethod=="JM"){
			genoValues  <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Pheno,Freq)
	        }else if(Weighted==TRUE && WghtMethod =="DW"){
		    genoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newCycle1GenoTable,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)
		    phenoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newCycle1GenoTable,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

		    genoValues <- genoValues_PredList[[1]]
		    phenoValues <- phenoValues_PredList[[1]]

		    alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]
				
		}else if(Weighted==FALSE){
		
		    genoValues_FamPair <- PredictBayesPhenoValues(newCycle1GenoTable_GM_FamPair,PredictionModel_Geno)
			phenoValues_FamPair <- PredictBayesPhenoValues(newCycle1GenoTable_GM_FamPair,PredictionModel_Pheno)
			
			genoValues_Fam <- PredictBayesPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Pheno)
			
			# genoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Geno)
			# phenoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Pheno)

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

            			
			genoValues_FamPair <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM_FamPair,PredictionModel_Geno)
			phenoValues_FamPair <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM_FamPair,PredictionModel_Pheno)
			
			genoValues_Fam <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Pheno)
			
			# genoValues_Fam_in_Pair <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM_Fam_in_Pair,PredictionModel_Geno)
			# phenoValues_Fam_in_Pair <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM_Fam_in_Pair,PredictionModel_Pheno)
		
        }
}


### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

Criterion_List <- getCriterionValueList_GM(genoValues_FamPair,phenoValues_FamPair,genoValSimValues_FamPair,phenoValSimValues_FamPair,cycle_nProgeny_FamPair,cycle_nProgeny_Fam,nFamilyPairs,nFamilyInd)


 genoValues_FamPair_List <- Criterion_List[[1]]
 phenoValues_FamPair_List <- Criterion_List[[2]]
 genoValSimValues_FamPair_List <- Criterion_List[[3]]
 phenoValSimValues_FamPair_List <- Criterion_List[[4]]

 genoValues_Fam_List <- Criterion_List[[5]]
 phenoValues_Fam_List <- Criterion_List[[6]]
 genoValSimValues_Fam_List <- Criterion_List[[7]]
 phenoValSimValues_Fam_List <- Criterion_List[[8]]



### Get list of cycle1 geno table for family ind and family pairs ######################################

	
	cycle1GenoTable_GM_List_FamPair <- rep(list(matrix(0,nrow=cycle1_nProgeny_FamPair,ncol=nMarkers)),nFamilyPairs)
	cycle1GenoTableMod_GM_List_FamPair <- rep(list(matrix(0,nrow=cycle1_nProgeny_FamPair,ncol=nMarkers)),nFamilyPairs)
	
	initIndex <- 1
	finalIndex <- cycle1_nProgeny_FamPair
	
	for(nFamPairs in 1:nFamilyPairs){ 
	
		cycle1GenoTable_GM_List_FamPair[[nFamPairs]] <- cycle1GenoTable_GM_FamPair[initIndex:finalIndex,]
		cycle1GenoTableMod_GM_List_FamPair[[nFamPairs]] <- cycle1GenoTableMod_GM_FamPair[initIndex:finalIndex,]
	
	    initIndex <- finalIndex+1
		finalIndex <- finalIndex + cycle1_nProgeny_FamPair
	} 
	
	
	cycle1GenoTable_GM_List_Fam <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilyInd)
	cycle1GenoTableMod_GM_List_Fam <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilyInd)
	
	initIndex <- 1
	finalIndex <- cycle1_nProgeny_Fam
	
	for(nFam in 1:nFamilyInd){ 
		
		cycle1GenoTable_GM_List_Fam[[nFam]] <- cycle1GenoTable_GM_Fam[initIndex:finalIndex,]
		cycle1GenoTableMod_GM_List_Fam[[nFam]] <- cycle1GenoTableMod_GM_Fam[initIndex:finalIndex,]
	
	    initIndex <- finalIndex+1
		finalIndex <- finalIndex + cycle1_nProgeny_Fam
		
	}

    # cycle1GenoTable_GM_List_Fam_in_Pair <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nUniqueFam_in_Pairs)
	# cycle1GenoTableMod_GM_List_Fam_in_Pair <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nUniqueFam_in_Pairs)
	
	# initIndex <- 1
	# finalIndex <- cycle1_nProgeny_Fam
	
	# for(nFam in 1:nUniqueFam_in_Pairs){ 
		
		# cycle1GenoTable_GM_List_Fam_in_Pair[[nFam]] <- cycle1GenoTable_GM_Fam_in_Pair[initIndex:finalIndex,]
		# cycle1GenoTableMod_GM_List_Fam_in_Pair[[nFam]] <- cycle1GenoTableMod_GM_Fam_in_Pair[initIndex:finalIndex,]
	
	    # initIndex <- finalIndex+1
		# finalIndex <- finalIndex + cycle1_nProgeny_Fam
		
	# }
	
    selectedGenoIndividualIndices2 <- rep(0,no_selected) 
	attainedGenoValues <- rep(0,nCycles)
	attainedPhenoValues <- rep(0,nCycles)
	 
	   
######################################################################################   
## Output list for all families 
	
	 GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	 GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	 GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 

	 PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles)  
	 PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	 PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	
     attainedPhenoValues_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
     attainedGenoValues_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 

## Output list for family pairs 	

 	  	  
	  selectedGenoIndividualIndices_List_FamPair <-  rep(list(list()),nFamilyPairs)
	  selectedPhenoIndividualIndices_List_FamPair <- rep(list(list()),nFamilyPairs)

## Output list for single fam 	  
	 	
		 	 	    
      selectedGenoIndividualIndices_List_Fam <- rep(list(list()),nFamilyInd)
      selectedPhenoIndividualIndices_List_Fam <- rep(list(list()),nFamilyInd)

	
#### GM method for parent selection in cycle 1 ################################################
#### GM method for parent selection for pairs of families ################################################
  system.time({
		
	genoSelected_List_Pairs <- getGM_selectedGenoList_FamilyPairs(cycle1GenoTable_GM_List_FamPair,cycle1GenoTableMod_GM_List_FamPair,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,familyPairs,GM_Param_List)
	
	selectedParents_Geno_FamPair_List[[nCyc]] <- genoSelected_List_Pairs[[1]]
	selectedParents_Geno <- genoSelected_List_Pairs[[1]]
	selectedParents_NumProgeny_Geno_FamPair_List[[nCyc]] <- genoSelected_List_Pairs[[2]]
 	selected_Parents_Geno_FamPair_Ind_List[[nCyc]] <- genoSelected_List_Pairs[[3]]
   })
	
# save.image(WorkspaceName)

    gc()
	
## Assign family pair values to variables 

		
	selectedGenoIndividualIndices_List_FamPair <-  rep(list(list()),nFamilyPairs)
	selectedPhenoIndividualIndices_List_FamPair <- rep(list(list()),nFamilyPairs)
	  
	topNGenoValues_List_FamPair <- rep(list(list()),nFamilyPairs)
	GenoTableIndices_topN_List_FamPair <- rep(list(list()),nFamilyPairs)
	  
	topNPhenoValues_List_FamPair <- rep(list(list()),nFamilyPairs)
	PhenoTableIndices_topN_List_FamPair <-rep(list(list()),nFamilyPairs)
	 
	
	for(nFamPairs in 1:nFamilyPairs){
	
	    initIndex <- 1
	    finalIndex <- cycle1_nProgeny
	
	    for(nPairs in 1:2){
	   
	    nfamily <- familyPairs[nFamPairs,nPairs]
	
	 	fam1 <- familyPairs[nFamPairs,1]
		fam2 <- familyPairs[nFamPairs,2]
		
		genoValSimValues_FamPair_List[[nFamPairs]] <- c(genoSimValues_List[[fam1]],genoSimValues_List[[fam2]]) 
		phenoValSimValues_FamPair_List[[nFamPairs]] <- c(phenoSimValues_List[[fam1]],phenoSimValues_List[[fam2]])
	    genoValues_FamPair_List[[nFamPairs]] <- c(genoValues_List[[fam1]],genoValues_List[[fam2]])
		phenoValues_FamPair_List[[nFamPairs]] <- c(phenoValues_List[[fam1]],phenoValues_List[[fam2]])
	  
	         
        genoIndices  <- selectedParents_Geno_FamPair_List[[nCyc]][[nFamPairs]]
		
		topNGenoSimValues <- genoValSimValues_FamPair_List[[nFamPairs]][genoIndices]
		topNPhenoSimValues <- phenoValSimValues_FamPair_List[[nFamPairs]][genoIndices]
		
		GenoTableIndices_topN <- genoIndices 
		PhenoTableIndices_topN <- genoIndices 
		
		topNGenoValues <- genoValues_FamPair_List[[nFamPairs]][genoIndices]
		topNPhenoValues <- phenoValues_FamPair_List[[nFamPairs]][genoIndices]
		
		
### GM Output in selected individual ids and Num Progeny from selected parents
### Genotypic values for nextGen Geno data  
 	   	  	  
	    topNGenoValues_List_FamPair[[nFamPairs]] <- topNGenoValues
	    GenoTableIndices_topN_List_FamPair[[nFamPairs]] <- GenoTableIndices_topN
	    selectedGenoIndividualIndices_List_FamPair[[nFamPairs]] <- GenoTableIndices_topN 
	 		
#### Phenotypic values for nextGen Geno data 
  	   
	    topNPhenoValues_List_FamPair[[nFamPairs]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_FamPair[[nFamPairs]] <- PhenoTableIndices_topN
		selectedPhenoIndividualIndices_List_FamPair[[nFamPairs]] <- PhenoTableIndices_topN
		
 ####
    	
		
		GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <- genoValSimValues_FamPair_List[[nFamPairs]][initIndex:finalIndex]
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <-  genoValues_FamPair_List[[nFamPairs]][initIndex:finalIndex]
		
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfamily]] <- topNGenoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nfamily]] <- max(genoValSimValues_FamPair_List[[nFamPairs]][GenoTableIndices_topN])
	
        PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <-  phenoValSimValues_FamPair_List[[nFamPairs]][initIndex:finalIndex]
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <- phenoValues_FamPair_List[[nFamPairs]][initIndex:finalIndex]
		
		PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfamily]] <- topNPhenoValues
		attainedPhenoValues_List_Cycle[[nCyc]][[nfamily]] <- max(phenoValSimValues_FamPair_List[[nFamPairs]][PhenoTableIndices_topN])
	
		initIndex <- finalIndex+1
		finalIndex <- finalIndex + cycle1_nProgeny
	
	}	
	
  }
 
						
# #################################################################################
# #### GM method for isolated families
 
selectedParents_Geno_Fam_List <- list()
 selectedParents_NumProgeny_Geno_Fam_List <- list()
 selectedParents_Geno_Fam <- list()

 system.time({
 
 genoSelected_List_Ind <- getGM_selectedGenoList_SingleFams(cycle1GenoTable_GM_List_Fam,cycle1GenoTableMod_GM_List_Fam,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,families,GM_Param_List)

 selectedParents_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[1]]
 selectedParents_NumProgeny_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[2]]
 selectedParents_Geno_Fam <- genoSelected_List_Ind[[1]]
 
})

 # save.image(WorkspaceName) 	
 #load(WorkspaceName) 
 gc()
 
### Output variables for isolated families
	
	  selectedGenoIndividualIndices_List_Fam <-  rep(list(list()),nFamilyInd)
	  selectedPhenoIndividualIndices_List_Fam <- rep(list(list()),nFamilyInd)
	  
	  topNGenoValues_List_Fam <- rep(list(list()),nFamilyInd)
	  GenoTableIndices_topN_List_Fam <- rep(list(list()),nFamilyInd)
	  
	  topNPhenoValues_List_Fam <- rep(list(list()),nFamilyInd)
	  PhenoTableIndices_topN_List_Fam <-rep(list(list()),nFamilyInd)


	for(nfam in 1:nFamilyInd){

			nfamily <- families[nfam]
			
			topNGenoSimValues <- genoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoSimValues <- phenoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
			topNGenoValues <- genoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoValues <- phenoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
### GM Output in selected individual ids and Num Progeny from selected parents
		
   	  	  
			topNGenoValues_List_Fam[[nfam]] <- topNGenoValues
			GenoTableIndices_topN_List_Fam[[nfam]] <- GenoTableIndices_topN
			selectedGenoIndividualIndices_List_Fam[[nfam]] <- GenoTableIndices_topN 
		 
#### Phenotypic values for nextGen Geno data
  
	   		topNPhenoValues_List_Fam[[nfam]] <- topNPhenoValues
			PhenoTableIndices_topN_List_Fam[[nfam]] <- PhenoTableIndices_topN
			selectedPhenoIndividualIndices_List_Fam[[nfam]] <- PhenoTableIndices_topN
			
					
#### Output data for cycle- nCyc			
		 
			
	       GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <- genoValSimValues_Fam_List[[nfam]]
	       GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <-  genoValues_Fam_List[[nfam]]
		   GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfamily]] <- topNGenoValues
	       attainedGenoValues_List_Cycle[[nCyc]][[nfamily]] <-max(genoValSimValues_Fam_List[[nfam]][GenoTableIndices_topN])

		   PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <-  phenoValSimValues_Fam_List[[nfam]]
	       PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <- phenoValues_Fam_List[[nfam]]
		   PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfamily]] <- topNPhenoValues
		   attainedPhenoValues_List_Cycle[[nCyc]][[nfamily]] <- max(phenoValSimValues_Fam_List[[nfam]][PhenoTableIndices_topN])
 		 
	}


###################################################################################	

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
## Get F5 RILs for selected parents in cycle 1
 
 
    Cycle_GenoData_List <- cycle1_GenoData_List

    Cycle_GenoData_List_FamPair <- cycle1_GenoData_List_FamPair
    Cycle_GenoData_List_Fam <- cycle1_GenoData_List_Fam


    F5RILs_GM_FamPair <- getF5RILs_GM_FamPair(Cycle_GenoData_List_FamPair,selectedParents_Geno_FamPair_List[[nCyc]],selected_Parents_Geno_FamPair_Ind_List[[nCyc]],selectedParents_Geno_Fam_List[[nCyc]],selectedParents_NumProgeny_Geno_FamPair_List[[nCyc]],selectedParents_NumProgeny_Geno_Fam_List[[nCyc]],nFamilies,nMarkers,NAM_LinkMap,familyPairs,families,cycle1_nProgeny)

    F5RILs_GM_Fam <- getF5RILs_GM_Fam(Cycle_GenoData_List_Fam,selectedParents_Geno_FamPair_List[[nCyc]],selected_Parents_Geno_FamPair_Ind_List[[nCyc]],selectedParents_Geno_Fam_List[[nCyc]],selectedParents_NumProgeny_Geno_FamPair_List[[nCyc]],selectedParents_NumProgeny_Geno_Fam_List[[nCyc]],nFamilies,nMarkers,NAM_LinkMap,familyPairs,families,cycle1_nProgeny)
	
	
    Cycle_Progeny_F5_FamPair <- F5RILs_GM_FamPair[[1]]
    Cycle_Progeny_F5_List_FamPair <- F5RILs_GM_FamPair[[2]]

    
	Cycle_Progeny_F5_Fam <- F5RILs_GM_Fam[[1]]
    Cycle_Progeny_F5_List_Fam <- F5RILs_GM_Fam[[2]]


	save.image(WorkspaceName) 
	gc() 

	startCycle <- startCycle +1 
}
###########################################################################


 for(nCyc in startCycle:nCycles){

        print(nCyc)

### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
     		
		nextGenGenoTable_GM_FamPair <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5_FamPair,nFamilyPairs,cycle1_nProgeny_FamPair) 
		nextGenGenoTable_GM_Fam <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5_Fam,nFamilyInd,cycle1_nProgeny_Fam)
		
		nextGenGenoTableMod_GM_FamPair <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_FamPair)
		nextGenGenoTableMod_GM_Fam <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_Fam)
		
		
		newNextGenGenoTable_GM_FamPair <- generateWeightedGenoTable(nextGenGenoTable_GM_FamPair,no_QTL)	
        newNextGenGenoTable_GM_Fam <- generateWeightedGenoTable(nextGenGenoTable_GM_Fam,no_QTL)
		
			
		nextGenGenoTable_GM_List_FamPair <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_FamPair,ncol=nMarkers))),nFamilyPairs)
        nextGenGenoTable_GM_List_Fam <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilyInd)
		
		nextGenGenoTableMod_GM_List_FamPair <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_FamPair,ncol=nMarkers))),nFamilyPairs)
        nextGenGenoTableMod_GM_List_Fam <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilyInd)
		
		for(nFamPairs in 1:nFamilyPairs){

				nextGenGenoTable_GM_List_FamPair[[nFamPairs]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List_FamPair[[nFamPairs]],1,cycle1_nProgeny_FamPair)
                
				nextGenGenoTableMod_GM_List_FamPair[[nFamPairs]] <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_List_FamPair[[nFamPairs]])
		
		}
		for(nfam in 1:nFamilyInd){ 
		
				nextGenGenoTable_GM_List_Fam[[nfam]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List_Fam[[nfam]],1,cycle1_nProgeny_Fam)
				nextGenGenoTableMod_GM_List_Fam[[nfam]] <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_List_Fam[[nfam]])
		
		}
			

		genoValSimValues_FamPair <- simulateGenotypicValues(nextGenGenoTable_GM_FamPair,no_QTL)
		phenoValSimValues_FamPair <- simulatePhenotypicValues(nextGenGenoTable_GM_FamPair,varE,no_QTL)
			
		genoValSimValues_Fam <- simulateGenotypicValues(nextGenGenoTable_GM_Fam,no_QTL)
		phenoValSimValues_Fam <- simulatePhenotypicValues(nextGenGenoTable_GM_Fam,varE,no_QTL)
		
  			
#########################
## Sclass plot for family clustering 	
#########################	
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
	
############################################################################################### 
###### Get updated or retrained prediction model for geno and pheno values

if(modelUpdate== TRUE && nCyc%%updateFrequency==0){

    nFamilyPair_List <- length(unique(as.vector(familyPairs)))
	nFamilyInd_List <- length(families)
	nFamilies_List <- nFamilyPair_List + nFamilyInd_List
	nextGenGeno_nProgeny_FamPair <- cycle1_nProgeny_Fam
	
	nextGenGenoTable_GM <- rbind(nextGenGenoTable_GM_FamPair,nextGenGenoTable_GM_Fam)
	nextGenGenoTableMod_GM <- rbind(nextGenGenoTableMod_GM_FamPair,nextGenGenoTableMod_GM_Fam)
	newNextGenGenoTable_GM <- rbind(newNextGenGenoTable_GM_FamPair,newNextGenGenoTable_GM_Fam)
	
	genoValSimValues <- c(genoValSimValues_FamPair,genoValSimValues_Fam)
	phenoValSimValues <- c(phenoValSimValues_FamPair,phenoValSimValues_Fam)

    nIndividuals <- length(genoValSimValues)
	IndividualNames <-rep(0,nIndividuals)

    for(nInd in 1:(nIndividuals)){

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
		trainGenoNewTable <- newNextGenGenoTable_GM[trainIndices,]
		trainSimPhenoTable <- phenotypeSimTable[trainIndices,]
		trainSimGenoTable <-  genotypicValuesSimTable[trainIndices,]

	#### Validation Set table for Genotype and Simulated Phenotype ###################

		# testGenoTable<- genoTable_Mod[testIndices,]
		testGenoNewTable <- newNextGenGenoTable_GM[testIndices,]
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
	
    if(modelUpdate ==TRUE && updateType == "FullSet"){
	
	   PredictionModel_Geno_PreCycle <- PredictionModel_Geno
       PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno

	    indices<-c(1:nIndividuals)

        trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices<- indices[which(!indices %in% trainIndices)]

#### Training Set table for Genotype and Simulated Phenotype ###################

        trainGenoNewTable <- newNextGenGenoTable_GM[trainIndices,]

        trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]

        trainSimGenoValTable <-  genotypicValuesSimTable[trainIndices,]

#### Validation Set table for Genotype and Simulated Phenotype ###################

        testGenoNewTable <- newNextGenGenoTable_GM[testIndices,]

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

                        PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable))

                }else if(modelType=="RRBLUP_REML"){



                        PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        PredictionModel_Pheno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

                }
				
		    }
		 
		    trainGenoNewTablePreCycle <-  trainGenoNewTableComb
            trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)
            trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
			
	if((length(levels(factor(PredictionModel_Pheno[[1]]))) == 1 && sd()==0) ||(length(levels(factor(PredictionModel_Geno[[1]]))) ==1 && sd()==0)){

                         PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         PredictionModel_Geno <- PredictionModel_Geno_PreCycle

       }
			
			

            rm(trainGenoNewTableComb)
            rm(trainSimGenoValTableComb)
            rm(trainSimPhenoValTableComb)
            gc()

        }
}	

    if(modelUpdate ==TRUE && updateType=="TrainingCycleWindow"){

		
			
		PredictionModel_Geno_PreCycle <- PredictionModel_Geno
		PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno
        
	    indices<-c(1:nIndividuals)

        trainIndices <-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices <- indices[which(!indices %in% trainIndices)]
		
		nTrainIndices <- length(trainIndices)

#### Training Set table for Genotype and Simulated Phenotype ###################

        trainGenoNewTable <- newNextGenGenoTable_GM[trainIndices,]
        trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]
        trainSimGenoValTable <- genotypicValuesSimTable[trainIndices,]

#### Validation Set table for Genotype and Simulated Phenotype ###################
        testGenoNewTable <- newNextGenGenoTable_GM[testIndices,]
        testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]
        testSimGenoValTable <-  genotypicValuesSimTable[testIndices,]

##### Build RRBLUP prediction model every cycle


        trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))
        trainSimGenoValTableComb <- rbind((trainSimGenoValTablePreCycle),trainSimGenoValTable)
        trainSimPhenoValTableComb <- rbind((trainSimPhenoTablePreCycle),trainSimPhenoTable)
        print(dim(trainGenoNewTableComb))
		Avg_GenoSim_Prev <- mean(trainSimGenoValTableComb[,1])
		Avg_GenoSim_Current <- mean(trainSimGenoValTable[,1])

		change_GenoSim <- Avg_GenoSim_Prev - Avg_GenoSim_Current
		combinedTableCount <- combinedTableCount+1
         print(change_GenoSim)

        if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){

           if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0 && change_GenoSim !=0 ){ 

           if(modelType=="RRBLUP"){

                        PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable))

           }else if(modelType=="RRBLUP_REML"){


                        PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        PredictionModel_Pheno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

           }

	     } 

           print(sd(PredictionModel_Pheno[[1]])) 
           print(sd(PredictionModel_Geno[[1]])) 
           print(length(levels(factor(PredictionModel_Pheno[[1]]))))
           print(length(levels(factor(PredictionModel_Geno[[1]]))))	 

           genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
           phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
           print(sd(phenoValues))
           print(sd(genoValues))

 
           if(length(levels(factor(PredictionModel_Pheno[[1]]))) ==1 || length(levels(factor(PredictionModel_Geno[[1]]))) ==1 || (is.na(sd(phenoValues)) || sd(phenoValues)==0) || (is.na(sd(genoValues)) || sd(genoValues)==0) || change_GenoSim==0){

                         PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         PredictionModel_Geno <- PredictionModel_Geno_PreCycle

           }
	 }
		
		if(combinedTableCount < trainTableWindowSize){
			
			    trainGenoNewTablePreCycle <-  trainGenoNewTableComb

                trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)

                trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
					
		}else if(combinedTableCount >= trainTableWindowSize){
			
			    trainGenoNewTablePreCycle <-  trainGenoNewTableComb[-c(1:nTrainIndices),]

                trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb[-c(1:nTrainIndices),])

                trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb[-c(1:nTrainIndices),])
				
				combinedTableCount <- (trainTableWindowSize-1)
					
		}
		
	
                rm(trainGenoNewTableComb)
                rm(trainSimGenoValTableComb)
                rm(trainSimPhenoValTableComb)
                gc()
       }
}
################################################
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
		
		    genoValues_FamPair <- PredictBayesPhenoValues(newNextGenGenoTable_GM_FamPair,PredictionModel_Geno)
			phenoValues_FamPair <- PredictBayesPhenoValues(newNextGenGenoTable_GM_FamPair,PredictionModel_Pheno)
			
			genoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
			# genoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
			# phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)

	    }	
}

if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

            			
			genoValues_FamPair <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_FamPair,PredictionModel_Geno)
			phenoValues_FamPair <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_FamPair,PredictionModel_Pheno)
			
			genoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
		
        }
}
	

### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

Criterion_List <- getCriterionValueList_GM (genoValues_FamPair,phenoValues_FamPair,genoValSimValues_FamPair,phenoValSimValues_FamPair,cycle_nProgeny_FamPair,cycle_nProgeny_Fam,nFamilyPairs,nFamilyInd)
	

 genoValues_FamPair_List <- Criterion_List[[1]]
 phenoValues_FamPair_List <- Criterion_List[[2]]
 genoValSimValues_FamPair_List <- Criterion_List[[3]]
 phenoValSimValues_FamPair_List <- Criterion_List[[4]] 
  
 genoValues_Fam_List <- Criterion_List[[5]]
 phenoValues_Fam_List <- Criterion_List[[6]]
 genoValSimValues_Fam_List <- Criterion_List[[7]]
 phenoValSimValues_Fam_List <- Criterion_List[[8]] 
 

### Migration policy change (nextGenGenoTable_GM_FamPair has to be changed for GC policy) 
 
 if( (migPolicyChangeFrequency !=0) && (nCyc %% migPolicyChangeFrequency ==0)){

 	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_FamPair_List <- genoValSimValues_FamPair_List
			SelCriterion_Fam_List <- genoValSimValues_Fam_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_FamPair_List <- genoValues_FamPair_List
			SelCriterion_Fam_List <- genoValues_Fam_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_FamPair_List <- phenoValSimValues_FamPair_List
			SelCriterion_Fam_List <- phenoValSimValues_Fam_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_FamPair_List <- phenoValues_FamPair_List
			SelCriterion_Fam_List <- phenoValues_Fam_List
		} 
     } 
 
    SelCriterion_List <- append(SelCriterion_FamPair_List,SelCriterion_Fam_List)
		
	nFamilyPair_List <- length(unique(as.vector(familyPairs)))
	nFamilyInd_List <- length(families)
	
	nFamilies_List <- nFamilyPair_List + nFamilyInd_List
	
	nextGenGeno_nProgeny_FamPair <- cycle1_nProgeny_Fam
	
	
 	# genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
	# phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
	
	# newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)

    nextGenGenoTable_GM <- rbind(nextGenGenoTable_GM_FamPair,nextGenGenoTable_GM_Fam)
	nextGenGenoTableMod_GM <- rbind(nextGenGenoTableMod_GM_FamPair,nextGenGenoTableMod_GM_Fam)
	newNextGenGenoTable_GM <- rbind(newNextGenGenoTable_GM_FamPair,newNextGenGenoTable_GM_Fam)
	
	genoValSimValues <- c(genoValSimValues_FamPair,genoValSimValues_Fam)
	phenoValSimValues <- c(phenoValSimValues_FamPair,phenoValSimValues_Fam)

    nIndividuals <- length(genoValSimValues)
	IndividualNames <-rep(0,nIndividuals)

    for(nInd in 1:(nIndividuals)){

        IndividualNames[nInd]<- paste("Ind",nInd,sep="")

    }
	
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
		
		   			
			genoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
			phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
			
		}	
    }

    if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

     						
			genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
			phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
				
        }
    }

 
	 initIndex <- 1
	 finalIndex <- cycle1_nProgeny 
	 genoValues_List <- list()
	 phenoValues_List <- list()
	 genoSimValues_List <- list()
	 phenoSimValues_List <- list()
	
	 
	 for(nFamily in 1:nFamilies_List){ 
		
		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoValSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoValSimValues)[initIndex:finalIndex]
						
		initIndex <- finalIndex+1 
		finalIndex <- finalIndex+cycle1_nProgeny
	 
	 } 
 	
###### Select island pairs based on migration policy 
### SPlit criterion list 
	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- genoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- genoValues_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- phenoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- phenoValues_List
		} 
     } 

#### Migration Policy change from the previous combined set ###################### 	

	familyPairs <-  migrationPolicy_GM(nextGenGenoTable_GM,SelCriterion_List,Policy,nFamilies_List)
	nFamilyPairs <- dim(familyPairs)[1]
 
    families <-  setdiff(c(1:nFamilies_List),unique(as.vector(familyPairs)))
	nFamilyInd <- length(families) 
	
	familyPairs_List[[nCyc]] <- familyPairs
	familyInd_List[[nCyc]] <- families
 
   
    nextGenGenoTable_GM_List_FamPair_MP_Change<- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_FamPair,ncol=nMarkers))),nFamilyPairs)
	nextGenGenoTableMod_GM_List_FamPair_MP_Change <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_FamPair,ncol=nMarkers))),nFamilyPairs)
    
	nextGenGenoTable_GM_List_Fam_MP_Change <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilyInd)
    nextGenGenoTableMod_GM_List_Fam_MP_Change <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilyInd)
	
	
	for(nFamPairs in 1:nFamilyPairs){
	  
	   FamPairGenoTable <- c()
	   FamPairGenoTableMod <- c()
	
	    for(nPairs in 1:2){
	   
	      nfam <- familyPairs[nFamPairs,nPairs]
		
		  FamPairGenoTable <- rbind(FamPairGenoTable,nextGenGenoTable_GM_List_FamPair[[nfam]])
		  FamPairGenoTableMod <- rbind(FamPairGenoTableMod,nextGenGenoTableMod_GM_List_FamPair[[nfam]])
		  
		}  
	 
	    nextGenGenoTable_GM_List_FamPair_MP_Change[[nFamPairs]] <-  FamPairGenoTable
		nextGenGenoTableMod_GM_List_FamPair_MP_Change[[nFamPairs]] <- FamPairGenoTableMod
	
	}
	
	
	for(nFamInd in 1:nFamilyInd){
	
	   fam <- families[nFamInd]
	   nextGenGenoTable_GM_List_Fam_MP_Change[[nFamInd]] <- nextGenGenoTable_GM_List_Fam[[fam]]
	   nextGenGenoTableMod_GM_List_Fam_MP_Change[[nFamInd]] <- nextGenGenoTableMod_GM_List_Fam[[fam]]
	
	}
				
	nextGenGenoTable_GM_List_FamPair <- nextGenGenoTable_GM_List_FamPair_MP_Change
	nextGenGenoTableMod_GM_List_FamPair <- nextGenGenoTableMod_GM_List_FamPair_MP_Change
	
	nextGenGenoTable_GM_List_Fam <- nextGenGenoTable_GM_List_Fam_MP_Change
	nextGenGenoTableMod_GM_List_Fam <- nextGenGenoTableMod_GM_List_Fam_MP_Change
	
###########			
	nextGenGenoTable_GM_FamPair <- c()
	nextGenGenoTable_GM_Fam	<- c()	
		
	for(nFamPairs in 1:nFamilyPairs){
		  nextGenGenoTable_GM_FamPair <- rbind(nextGenGenoTable_GM_FamPair,nextGenGenoTable_GM_List_FamPair[[nFamPairs]])
        	
	}
	
	for(nfam in 1:nFamilyInd){ 
		  nextGenGenoTable_GM_Fam <-	rbind(nextGenGenoTable_GM_Fam,nextGenGenoTable_GM_List_Fam[[nfam]]) 
	}	
	
	nextGenGenoTableMod_GM_FamPair <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_FamPair)
	nextGenGenoTableMod_GM_Fam <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_Fam)
		
	newNextGenGenoTable_GM_FamPair <- generateWeightedGenoTable(nextGenGenoTable_GM_FamPair,no_QTL)	
    newNextGenGenoTable_GM_Fam <- generateWeightedGenoTable(nextGenGenoTable_GM_Fam,no_QTL)
     
	genoValSimValues_FamPair <- simulateGenotypicValues(nextGenGenoTable_GM_FamPair,no_QTL)
	phenoValSimValues_FamPair <- simulatePhenotypicValues(nextGenGenoTable_GM_FamPair,varE,no_QTL)
			
	genoValSimValues_Fam <- simulateGenotypicValues(nextGenGenoTable_GM_Fam,no_QTL)
	phenoValSimValues_Fam <- simulatePhenotypicValues(nextGenGenoTable_GM_Fam,varE,no_QTL)
		

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
		
		    genoValues_FamPair <- PredictBayesPhenoValues(newNextGenGenoTable_GM_FamPair,PredictionModel_Geno)
			phenoValues_FamPair <- PredictBayesPhenoValues(newNextGenGenoTable_GM_FamPair,PredictionModel_Pheno)
			
			genoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
			

	    }	
    }

    if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

            			
			genoValues_FamPair <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_FamPair,PredictionModel_Geno)
			phenoValues_FamPair <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_FamPair,PredictionModel_Pheno)
			
			genoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
		
        }
    }

  
Criterion_List <- getCriterionValueList_GM (genoValues_FamPair,phenoValues_FamPair,genoValSimValues_FamPair,phenoValSimValues_FamPair,cycle_nProgeny_FamPair,cycle_nProgeny_Fam,nFamilyPairs,nFamilyInd)
	

 genoValues_FamPair_List <- Criterion_List[[1]]
 phenoValues_FamPair_List <- Criterion_List[[2]]
 genoValSimValues_FamPair_List <- Criterion_List[[3]]
 phenoValSimValues_FamPair_List <- Criterion_List[[4]] 
  
 genoValues_Fam_List <- Criterion_List[[5]]
 phenoValues_Fam_List <- Criterion_List[[6]]
 genoValSimValues_Fam_List <- Criterion_List[[7]]
 phenoValSimValues_Fam_List <- Criterion_List[[8]] 		
		
} 



### GM method for nextGen Geno data set family pairs


 genoSelected_List_Pairs <- getGM_selectedGenoList_FamilyPairs(nextGenGenoTable_GM_List_FamPair,nextGenGenoTableMod_GM_List_FamPair,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,familyPairs,GM_Param_List)
	
	
  selectedParents_Geno_FamPair_List[[nCyc]] <- genoSelected_List_Pairs[[1]]
  selectedParents_Geno <- genoSelected_List_Pairs[[1]]
  selectedParents_NumProgeny_Geno_FamPair_List[[nCyc]] <- genoSelected_List_Pairs[[2]]
  selected_Parents_Geno_FamPair_Ind_List[[nCyc]]<- genoSelected_List_Pairs[[3]]
 
   
 save.image(WorkspaceName)
 gc()
 
  
    selectedGenoIndividualIndices_List_FamPair <-  rep(list(list()),nFamilyPairs)
    selectedPhenoIndividualIndices_List_FamPair <- rep(list(list()),nFamilyPairs)
	  
    topNGenoValues_List_FamPair <- rep(list(list()),nFamilyPairs)
	GenoTableIndices_topN_List_FamPair <- rep(list(list()),nFamilyPairs)
	  
	topNPhenoValues_List_FamPair <- rep(list(list()),nFamilyPairs)
	PhenoTableIndices_topN_List_FamPair <-rep(list(list()),nFamilyPairs)
 
  
### Assign values to ouput variables 
### Output variables for family pairs

  
    familyPairs <- familyPairs_List[[nCyc]]
    families <- familyInd_List[[nCyc]]

	for(nFamPairs in 1:nFamilyPairs){
	
	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	
	  for(nPairs in 1:2){
	   
	    nfamily <- familyPairs[nFamPairs,nPairs]
		      
        genoIndices  <- selectedParents_Geno_FamPair_List[[nCyc]][[nFamPairs]]
		
		topNGenoSimValues <- genoValSimValues_FamPair_List[[nFamPairs]][genoIndices]
		topNPhenoSimValues <- phenoValSimValues_FamPair_List[[nFamPairs]][genoIndices]
		
		GenoTableIndices_topN <- genoIndices 
		PhenoTableIndices_topN <- genoIndices 
		
		topNGenoValues <- genoValues_FamPair_List[[nFamPairs]][genoIndices]
		topNPhenoValues <- phenoValues_FamPair_List[[nFamPairs]][genoIndices]
		
### GM Output in selected individual ids and Num Progeny from selected parents
### Genotypic values for nextGen Geno data  
 	   	  	  
	    topNGenoValues_List_FamPair[[nFamPairs]] <- topNGenoValues
	   	GenoTableIndices_topN_List_FamPair[[nFamPairs]] <- GenoTableIndices_topN
	    selectedGenoIndividualIndices_List_FamPair[[nFamPairs]] <- GenoTableIndices_topN 
		
	
#### Phenotypic values for nextGen Geno data 
  	   
	   	topNPhenoValues_List_FamPair[[nFamPairs]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_FamPair[[nFamPairs]] <- PhenoTableIndices_topN
	    selectedPhenoIndividualIndices_List_FamPair[[nFamPairs]] <- PhenoTableIndices_topN
				

#### Output data for cycle- nCyc					
	 
	   		 
		GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <- genoValSimValues_FamPair_List[[nFamPairs]][initIndex:finalIndex]
	     GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <-  genoValues_FamPair_List[[nFamPairs]][initIndex:finalIndex]
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfamily]] <- topNGenoValues
	     attainedGenoValues_List_Cycle[[nCyc]][[nfamily]] <- max(genoValSimValues_FamPair_List[[nFamPairs]][GenoTableIndices_topN]) 
	
            PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <-  phenoValSimValues_FamPair_List[[nFamPairs]][initIndex:finalIndex]
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <- phenoValues_FamPair_List[[nFamPairs]][initIndex:finalIndex]
		PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfamily]] <- topNPhenoValues
		attainedPhenoValues_List_Cycle[[nCyc]][[nfamily]] <-  max(phenoValSimValues_FamPair_List[[nFamPairs]][PhenoTableIndices_topN])
		
		initIndex <- finalIndex+1
		finalIndex <- finalIndex + cycle1_nProgeny
	
 
	   }
	}


#### GM method for isolated families


 genoSelected_List_Ind <- getGM_selectedGenoList_SingleFams(nextGenGenoTable_GM_List_Fam,nextGenGenoTableMod_GM_List_Fam,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,families,GM_Param_List)

 
  selectedParents_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[1]]
  selectedParents_NumProgeny_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[2]]
  selectedParents_Geno_Fam <- genoSelected_List_Ind[[1]]
 
  save.image(WorkspaceName)
  gc() 	
  
########################################################################
### Output variables for isolated families
	
  selectedGenoIndividualIndices_List_Fam <-  rep(list(list()),nFamilyInd)
  selectedPhenoIndividualIndices_List_Fam <- rep(list(list()),nFamilyInd)
	  
  topNGenoValues_List_Fam <- rep(list(list()),nFamilyInd)
  GenoTableIndices_topN_List_Fam <- rep(list(list()),nFamilyInd)
	  
  topNPhenoValues_List_Fam <- rep(list(list()),nFamilyInd)
  PhenoTableIndices_topN_List_Fam <-rep(list(list()),nFamilyInd)
  
  
	for(nfam in 1:nFamilyInd){

			nfamily <- families[nfam]
			
			topNGenoSimValues <- genoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoSimValues <- phenoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
			topNGenoValues <- genoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoValues <- phenoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
### GM Output in selected individual ids and Num Progeny from selected parents
#### Genotypic values for nextGen Geno data		
   	  	  
			topNGenoValues_List_Fam[[nfam]] <- topNGenoValues
			GenoTableIndices_topN_List_Fam[[nfam]] <- GenoTableIndices_topN
			selectedGenoIndividualIndices_List_Fam[[nfam]] <- GenoTableIndices_topN 
#### Phenotypic values for nextGen Geno data
			
			topNPhenoValues_List_Fam[[nfam]] <- topNPhenoValues
			PhenoTableIndices_topN_List_Fam[[nfam]] <- PhenoTableIndices_topN
			selectedPhenoIndividualIndices_List_Fam[[nfam]] <- PhenoTableIndices_topN
			

#### Output data for cycle- nCyc			
		 
			
	       GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <- genoValSimValues_Fam_List[[nfam]]
	       GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <-  genoValues_Fam_List[[nfam]]
		   GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfamily]] <- topNGenoValues
	       attainedGenoValues_List_Cycle[[nCyc]][[nfamily]] <-max(genoValSimValues_Fam_List[[nfam]][GenoTableIndices_topN])
		   
		   
		   PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <-  phenoValSimValues_Fam_List[[nfam]]
	       PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfamily]] <- phenoValues_Fam_List[[nfam]]
		   PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfamily]] <- topNPhenoValues
		   attainedPhenoValues_List_Cycle[[nCyc]][[nfamily]] <- max(phenoValSimValues_Fam_List[[nfam]][PhenoTableIndices_topN])
 
				 
	}

	print(paste(nCyc,"-Fam20-",mean(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]]),sep=""))
#########################

 	Cycle_GenoData_List_FamPair <- Cycle_Progeny_F5_List_FamPair
	Cycle_GenoData_List_Fam <- Cycle_Progeny_F5_List_Fam
	
	F5RILs_GM_FamPair <- getF5RILs_GM_FamPair(Cycle_GenoData_List_FamPair,selectedParents_Geno_FamPair_List[[nCyc]],selected_Parents_Geno_FamPair_Ind_List[[nCyc]],selectedParents_Geno_Fam_List[[nCyc]],selectedParents_NumProgeny_Geno_FamPair_List[[nCyc]],selectedParents_NumProgeny_Geno_Fam_List[[nCyc]],nFamilies,nMarkers,NAM_LinkMap,familyPairs,families,cycle1_nProgeny)	

    F5RILs_GM_Fam <- getF5RILs_GM_Fam(Cycle_GenoData_List_Fam,selectedParents_Geno_FamPair_List[[nCyc]],selected_Parents_Geno_FamPair_Ind_List[[nCyc]],selectedParents_Geno_Fam_List[[nCyc]],selectedParents_NumProgeny_Geno_FamPair_List[[nCyc]],selectedParents_NumProgeny_Geno_Fam_List[[nCyc]],nFamilies,nMarkers,NAM_LinkMap,familyPairs,families,cycle1_nProgeny)	

	Cycle_Progeny_F5_FamPair <- F5RILs_GM_FamPair[[1]]
	Cycle_Progeny_F5_List_FamPair <- F5RILs_GM_FamPair[[2]]
	Cycle_Progeny_F5_Fam <- F5RILs_GM_Fam[[1]]
	Cycle_Progeny_F5_List_Fam <- F5RILs_GM_Fam[[2]]
	
	gc()
	
	
	 
}


######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,attainedPhenoValues_List_Cycle,familyPairs_List,familyInd_List)

	   return(simResults_List)

}

#### 


runSimulations20X_IslandSelection_BD_Wght_GM_V2 <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Param_List){

	 
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
	  
	  GM_Param_List <- GM_Param_List
	  no_cores <- GM_Param_List$no_cores
	  
	  WorkspaceName <- WorkSpaceName
	  startCycle <- StartCycle

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
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <-  generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable_GM <- generateWeightedGenoTable(cycle1GenoTable_GM,no_QTL)
	cycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(cycle1GenoTable_GM)
	
		
	nIndividuals <- nrow(cycle1GenoTable_GM)
	 

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
# # ### Split geno and pheno data tables in to families		
	 
	  cycle1_nProgeny <- 100
	  nProgeny <- 100
	  
	  nSel_inFamily <- 1
	  
	  F5_Progeny_List_Family<- F5_Progeny_List
	  
	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  
	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)
	  
	  for(nFamily in 1:20){ 
		
		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,] 
	  
	  }
	  
# ####################		
	 	  
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	 	  
	  for(nFamily in 1:20){ 
		
		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,] 
		
	  } 
	  
	  
# # ## Split Geno and Pheno values according to families  ########################


	 initIndex <- 1
	 finalIndex <- cycle1_nProgeny 
	 genoValues_List <- list()
	 phenoValues_List <- list()
	 genoSimValues_List <- list()
	 phenoSimValues_List <- list()
	 cycle1_GenoTable_GM_List <- rep(list(matrix(rep(0,cycle1_nProgeny*nMarkers),nrow=cycle1_nProgeny,ncol=nMarkers)),nFamilies)
	 	  
	 cycle1_GenoTableMod_GM_List <- rep(list(matrix(rep(0,cycle1_nProgeny*nMarkers),nrow=cycle1_nProgeny,ncol=nMarkers)),nFamilies)
	  
	 
	 for(nFamily in 1:nFamilies){ 
		
		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]
		
		cycle1_GenoTable_GM_List[[nFamily]] <- cycle1GenoTable_GM[initIndex:finalIndex,]
		cycle1_GenoTableMod_GM_List[[nFamily]] <-cycle1GenoTableMod_GM[initIndex:finalIndex,]
		
		initIndex <- finalIndex+1 
		finalIndex <- finalIndex+cycle1_nProgeny
	 
	 } 

	 
	
###### Select island pairs based on migration policy 
### Split criterion list 
	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- genoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- genoValues_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- phenoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- phenoValues_List
		} 
         } 


### i)  Get pairs of families and isolated families based on migration policy 
 
	families <-  c(1:nFamilies)
    
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
	
		
	cycle_nProgeny_Fam <- cycle1_nProgeny
			 
 	selectedParents_Geno_Fam_List <- rep(list(rep(list(),nFamilies)),nCycles) 
	selectedParents_NumProgeny_Geno_Fam_List <- rep(list(rep(list(),nFamilies)),nCycles)
	
### i) Get cycle 1 genotype table and criterion value list for sets of singleislands and pairs of islands	
		
    cycle1GenoTable_GM_Fam <-  generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny_Fam)
    newCycle1GenoTable_GM_Fam <- generateWeightedGenoTable(cycle1GenoTable_GM_Fam,no_QTL)
    cycle1GenoTableMod_GM_Fam <- generate012GenoFormat_Inverse(cycle1GenoTable_GM_Fam)

   
   
    genoValSimValues_Fam <- simulateGenotypicValues(cycle1GenoTable_GM_Fam,no_QTL) 
    phenoValSimValues_Fam <- simulatePhenotypicValues(cycle1GenoTable_GM_Fam,varE,no_QTL) 



if(modelType=="BayesB" || modelType=="BL"){
		
		if(Weighted==TRUE && WghtMethod=="JM"){
			genoValues  <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Pheno,Freq)
	        }else if(Weighted==TRUE && WghtMethod =="DW"){
		    genoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newCycle1GenoTable,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)
		    phenoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newCycle1GenoTable,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

		    genoValues <- genoValues_PredList[[1]]
		    phenoValues <- phenoValues_PredList[[1]]

		    alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]
				
		}else if(Weighted==FALSE){
		
		    genoValues_FamPair <- PredictBayesPhenoValues(newCycle1GenoTable_GM_FamPair,PredictionModel_Geno)
			phenoValues_FamPair <- PredictBayesPhenoValues(newCycle1GenoTable_GM_FamPair,PredictionModel_Pheno)
			
			genoValues_Fam <- PredictBayesPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Pheno)
			
			# genoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Geno)
			# phenoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Pheno)

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

            			
			genoValues_Fam <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Pheno)
			
		
        }
}


### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

Criterion_List <- getCriterionValueList_GM_V2(genoValues_Fam,phenoValues_Fam,genoValSimValues_Fam,phenoValSimValues_Fam,cycle_nProgeny_Fam,nFamilies)

 genoValues_Fam_List <- Criterion_List[[1]]
 phenoValues_Fam_List <- Criterion_List[[2]]
 genoValSimValues_Fam_List <- Criterion_List[[3]]
 phenoValSimValues_Fam_List <- Criterion_List[[4]]



### Get list of cycle1 geno table for family ind and family pairs ######################################

		
	cycle1GenoTable_GM_List_Fam <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
	cycle1GenoTableMod_GM_List_Fam <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
	
	initIndex <- 1
	finalIndex <- cycle1_nProgeny_Fam
	
	for(nFam in 1:nFamilies){ 
		
		cycle1GenoTable_GM_List_Fam[[nFam]] <- cycle1GenoTable_GM_Fam[initIndex:finalIndex,]
		cycle1GenoTableMod_GM_List_Fam[[nFam]] <- cycle1GenoTableMod_GM_Fam[initIndex:finalIndex,]
                           
	        initIndex <- finalIndex+1
		finalIndex <- finalIndex + cycle1_nProgeny_Fam
		
	}

 	
        selectedGenoIndividualIndices2 <- rep(0,no_selected) 
	attainedGenoValues <- rep(0,nCycles)
	attainedPhenoValues <- rep(0,nCycles)
	 
	   
######################################################################################   
## Output list for all families 
	
	 GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	 GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	 GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 

	 PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles)  
	 PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	 PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	
     attainedPhenoValues_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
     attainedGenoValues_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 


## Output list for single fam 	  
	 	
		 	 	    
      selectedGenoIndividualIndices_List_Fam <- rep(list(list()),nFamilies)
      selectedPhenoIndividualIndices_List_Fam <- rep(list(list()),nFamilies)

	
#### GM method for parent selection in cycle 1 ################################################
# #### GM method for isolated families
 
 selectedParents_Geno_Fam_List <- list()
 selectedParents_NumProgeny_Geno_Fam_List <- list()
 selectedParents_Geno_Fam <- list()

 system.time({
 
 genoSelected_List_Ind <- getGM_selectedGenoList_SingleFams_V2(cycle1GenoTable_GM_List_Fam,cycle1GenoTableMod_GM_List_Fam,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,families,GM_Param_List)

 selectedParents_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[1]]
 selectedParents_NumProgeny_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[2]]
 selectedParents_Geno_Fam <- genoSelected_List_Ind[[1]]
 
})

 # save.image(WorkspaceName) 	
 #load(WorkspaceName) 
 gc()
 
### Output variables for isolated families
	
	  selectedGenoIndividualIndices_List_Fam <-  rep(list(list()),nFamilies)
	  selectedPhenoIndividualIndices_List_Fam <- rep(list(list()),nFamilies)
	  
	  topNGenoValues_List_Fam <- rep(list(list()),nFamilies)
	  GenoTableIndices_topN_List_Fam <- rep(list(list()),nFamilies)
	  
	  topNPhenoValues_List_Fam <- rep(list(list()),nFamilies)
	  PhenoTableIndices_topN_List_Fam <-rep(list(list()),nFamilies)


	for(nfam in 1:nFamilies){

			#nfamily <- families[nfam]
			
			topNGenoSimValues <- genoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoSimValues <- phenoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
			topNGenoValues <- genoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoValues <- phenoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
### GM Output in selected individual ids and Num Progeny from selected parents
		
   	  	  
			topNGenoValues_List_Fam[[nfam]] <- topNGenoValues
			GenoTableIndices_topN_List_Fam[[nfam]] <- GenoTableIndices_topN
			selectedGenoIndividualIndices_List_Fam[[nfam]] <- GenoTableIndices_topN 
		 
#### Phenotypic values for nextGen Geno data
  
	   		topNPhenoValues_List_Fam[[nfam]] <- topNPhenoValues
			PhenoTableIndices_topN_List_Fam[[nfam]] <- PhenoTableIndices_topN
			selectedPhenoIndividualIndices_List_Fam[[nfam]] <- PhenoTableIndices_topN
			
					
#### Output data for cycle- nCyc			
		 
			
	                GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <- genoValSimValues_Fam_List[[nfam]]
	                GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <-  genoValues_Fam_List[[nfam]]
		        GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfam]] <- topNGenoValues
	                attainedGenoValues_List_Cycle[[nCyc]][[nfam]] <-max(genoValSimValues_Fam_List[[nfam]][GenoTableIndices_topN])

		        PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <-  phenoValSimValues_Fam_List[[nfam]]
	                PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <- phenoValues_Fam_List[[nfam]]
		        PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfam]] <- topNPhenoValues
		        attainedPhenoValues_List_Cycle[[nCyc]][[nfam]] <- max(phenoValSimValues_Fam_List[[nfam]][PhenoTableIndices_topN])
 		 
	}


###################################################################################	

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
## Get F5 RILs for selected parents in cycle 1
 
 
    Cycle_GenoData_List <- cycle1_GenoData_List
      
    F5RILs_GM_Fam <- getF5RILs_GM_Fam_V2(Cycle_GenoData_List,selectedParents_Geno_Fam_List[[nCyc]],selectedParents_NumProgeny_Geno_Fam_List[[nCyc]],nFamilies,nMarkers,NAM_LinkMap,families,cycle1_nProgeny)
	 
	Cycle_Progeny_F5_Fam <- F5RILs_GM_Fam[[1]]
        Cycle_Progeny_F5_List_Fam <- F5RILs_GM_Fam[[2]]


	save.image(WorkspaceName) 
	gc() 

	startCycle <- startCycle +1 
}
###########################################################################

 startCycle <- 2
 for(nCyc in startCycle:nCycles){

        
        print(nCyc)
		
 #exchangeGenoData_GM <- function(Cycle_Progeny_F5_List_Fam,genoValSimValues_Fam_List,MigrationSize,EmigrantGroups,ImmigrantGroups,Direction,Policy){
			
### Exchange selected geno data among pairs
    if(migFreq !=0){
          if(nCyc>=2 && nCyc%%migFreq==0){
	  
	    if(Policy != "GeneticClusters"){
		 
		 output_List_afterExchange <- exchangeGenoData_GM(Cycle_Progeny_F5_List_Fam,genoValSimValues_Fam_List,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
		}else if(Policy == "GeneticClusters"){ 
		
		 
		 output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		}
		  
		Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
		
		familyInfo <- output_List_afterExchange[[2]]
		
		Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 
	   
	   
           Cycle_Progeny_F5_Fam <- array(0,c(nFamilies,nMarkers,2,cycle1_nProgeny))
           nProgIndex <-  dim(Cycle_Progeny_F5_List_Fam[[nFamily]])[3]
	   nProg <-1
				
	   for(nFamily in 1:length(families)){ 
					 
		Cycle_Progeny_F5_Fam[nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][,,nProg:nProgIndex] 

	   } 
	
      }
	   
	      
		
##    Cycle_Progeny_F5_Fam <-   Cycle_Progeny_F5_List_Fam 
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
     		
	
	    nextGenGenoTable_GM_Fam <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5_Fam,nFamilies,cycle1_nProgeny_Fam)
	    nextGenGenoTableMod_GM_Fam <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_Fam)
	    newNextGenGenoTable_GM_Fam <- generateWeightedGenoTable(nextGenGenoTable_GM_Fam,no_QTL)
		
	    nextGenGenoTable_GM_List_Fam <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilies)
        nextGenGenoTableMod_GM_List_Fam <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilies)
		
		
		for(nfam in 1:nFamilies){ 
		
				nextGenGenoTable_GM_List_Fam[[nfam]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List_Fam[[nfam]],1,cycle1_nProgeny_Fam)
				nextGenGenoTableMod_GM_List_Fam[[nfam]] <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_List_Fam[[nfam]])
		
		}
		
		genoValSimValues_Fam <- simulateGenotypicValues(nextGenGenoTable_GM_Fam,no_QTL)
		phenoValSimValues_Fam <- simulatePhenotypicValues(nextGenGenoTable_GM_Fam,varE,no_QTL)
		
  			
#########################
## Sclass plot for family clustering 	
#########################	
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
	
############################################################################################### 
###### Get updated or retrained prediction model for geno and pheno values

if(modelUpdate== TRUE && nCyc%%updateFrequency==0){

	nFamilies_List <- length(families) 
	nextGenGeno_nProgeny_FamPair <- cycle1_nProgeny_Fam
	
	nextGenGenoTable_GM <- nextGenGenoTable_GM_Fam
	nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_Fam
	newNextGenGenoTable_GM <- newNextGenGenoTable_GM_Fam
	
	genoValSimValues <- genoValSimValues_Fam
	phenoValSimValues <- phenoValSimValues_Fam

        nIndividuals <- length(genoValSimValues)
	IndividualNames <-rep(0,nIndividuals)

        for(nInd in 1:(nIndividuals)){

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
		trainGenoNewTable <- newNextGenGenoTable_GM[trainIndices,]
		trainSimPhenoTable <- phenotypeSimTable[trainIndices,]
		trainSimGenoTable <-  genotypicValuesSimTable[trainIndices,]

	#### Validation Set table for Genotype and Simulated Phenotype ###################

		# testGenoTable<- genoTable_Mod[testIndices,]
		testGenoNewTable <- newNextGenGenoTable_GM[testIndices,]
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
	
    if(modelUpdate ==TRUE && updateType == "FullSet"){
	
	   PredictionModel_Geno_PreCycle <- PredictionModel_Geno
       PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno

	    indices<-c(1:nIndividuals)

        trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices<- indices[which(!indices %in% trainIndices)]

#### Training Set table for Genotype and Simulated Phenotype ###################

        trainGenoNewTable <- newNextGenGenoTable_GM[trainIndices,]

        trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]

        trainSimGenoValTable <-  genotypicValuesSimTable[trainIndices,]

#### Validation Set table for Genotype and Simulated Phenotype ###################

        testGenoNewTable <- newNextGenGenoTable_GM[testIndices,]

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

                        PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable))

                }else if(modelType=="RRBLUP_REML"){



                        PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        PredictionModel_Pheno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

                }
				
		    }
		 
		    trainGenoNewTablePreCycle <-  trainGenoNewTableComb
            trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)
            trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
			
	if((length(levels(factor(PredictionModel_Pheno[[1]]))) == 1 && sd()==0) ||(length(levels(factor(PredictionModel_Geno[[1]]))) ==1 && sd()==0)){

                         PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         PredictionModel_Geno <- PredictionModel_Geno_PreCycle

       }
			
			

            rm(trainGenoNewTableComb)
            rm(trainSimGenoValTableComb)
            rm(trainSimPhenoValTableComb)
            gc()

        }
}	

    if(modelUpdate ==TRUE && updateType=="TrainingCycleWindow"){

		
			
		PredictionModel_Geno_PreCycle <- PredictionModel_Geno
		PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno
        
	    indices<-c(1:nIndividuals)

        trainIndices <-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices <- indices[which(!indices %in% trainIndices)]
		
		nTrainIndices <- length(trainIndices)

#### Training Set table for Genotype and Simulated Phenotype ###################

        trainGenoNewTable <- newNextGenGenoTable_GM[trainIndices,]
        trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]
        trainSimGenoValTable <- genotypicValuesSimTable[trainIndices,]

#### Validation Set table for Genotype and Simulated Phenotype ###################
        testGenoNewTable <- newNextGenGenoTable_GM[testIndices,]
        testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]
        testSimGenoValTable <-  genotypicValuesSimTable[testIndices,]

##### Build RRBLUP prediction model every cycle


        trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))
        trainSimGenoValTableComb <- rbind((trainSimGenoValTablePreCycle),trainSimGenoValTable)
        trainSimPhenoValTableComb <- rbind((trainSimPhenoTablePreCycle),trainSimPhenoTable)
        print(dim(trainGenoNewTableComb))
		Avg_GenoSim_Prev <- mean(trainSimGenoValTableComb[,1])
		Avg_GenoSim_Current <- mean(trainSimGenoValTable[,1])

		change_GenoSim <- Avg_GenoSim_Prev - Avg_GenoSim_Current
		combinedTableCount <- combinedTableCount+1
         print(change_GenoSim)

        if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){

           if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0 && change_GenoSim !=0 ){ 

           if(modelType=="RRBLUP"){

                        PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable))

           }else if(modelType=="RRBLUP_REML"){


                        PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        PredictionModel_Pheno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

           }

	     } 

           print(sd(PredictionModel_Pheno[[1]])) 
           print(sd(PredictionModel_Geno[[1]])) 
           print(length(levels(factor(PredictionModel_Pheno[[1]]))))
           print(length(levels(factor(PredictionModel_Geno[[1]]))))	 

           genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
           phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
           print(sd(phenoValues))
           print(sd(genoValues))

 
           if(length(levels(factor(PredictionModel_Pheno[[1]]))) ==1 || length(levels(factor(PredictionModel_Geno[[1]]))) ==1 || (is.na(sd(phenoValues)) || sd(phenoValues)==0) || (is.na(sd(genoValues)) || sd(genoValues)==0) || change_GenoSim==0){

                         PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         PredictionModel_Geno <- PredictionModel_Geno_PreCycle

           }
	 }
		
		if(combinedTableCount < trainTableWindowSize){
			
			    trainGenoNewTablePreCycle <-  trainGenoNewTableComb

                trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)

                trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
					
		}else if(combinedTableCount >= trainTableWindowSize){
			
			    trainGenoNewTablePreCycle <-  trainGenoNewTableComb[-c(1:nTrainIndices),]

                trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb[-c(1:nTrainIndices),])

                trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb[-c(1:nTrainIndices),])
				
				combinedTableCount <- (trainTableWindowSize-1)
					
		}
		
	
                rm(trainGenoNewTableComb)
                rm(trainSimGenoValTableComb)
                rm(trainSimPhenoValTableComb)
                gc()
       }
}
################################################
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
				
			genoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
	
	    }	
}

if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

        
			genoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
		
        }
}
	

### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

Criterion_List <- getCriterionValueList_GM_V2(genoValues_Fam,phenoValues_Fam,genoValSimValues_Fam,phenoValSimValues_Fam,cycle_nProgeny_Fam,nFamilies)
	
  
 genoValues_Fam_List <- Criterion_List[[1]]
 phenoValues_Fam_List <- Criterion_List[[2]]
 genoValSimValues_Fam_List <- Criterion_List[[3]]
 phenoValSimValues_Fam_List <- Criterion_List[[4]] 
 

### Migration policy change (nextGenGenoTable_GM_FamPair has to be changed for GC policy) 
 
 if( (migPolicyChangeFrequency !=0) && (nCyc %% migPolicyChangeFrequency ==0)){

 	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_Fam_List <- genoValSimValues_Fam_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_Fam_List <- genoValues_Fam_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_Fam_List <- phenoValSimValues_Fam_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_Fam_List <- phenoValues_Fam_List
		} 
     } 
 
    SelCriterion_List <- SelCriterion_Fam_List
	
	nFamilies_List <- length(families)
	
	nextGenGeno_nProgeny_FamPair <- cycle1_nProgeny
	
	
 	# genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
	# phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
	# newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)

        nextGenGenoTable_GM <- nextGenGenoTable_GM_Fam
	nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_Fam
	newNextGenGenoTable_GM <- newNextGenGenoTable_GM_Fam
	
	genoValSimValues <- genoValSimValues_Fam
	phenoValSimValues <- phenoValSimValues_Fam

        nIndividuals <- length(genoValSimValues)
	IndividualNames <-rep(0,nIndividuals)

        for(nInd in 1:(nIndividuals)){

          IndividualNames[nInd]<- paste("Ind",nInd,sep="")

        }
	
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
		
		   			
			genoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
			phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
			
		}	
    }

    if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

     						
			genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
			phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
				
        }
    }

 
	 initIndex <- 1
	 finalIndex <- cycle1_nProgeny 
	 genoValues_List <- list()
	 phenoValues_List <- list()
	 genoSimValues_List <- list()
	 phenoSimValues_List <- list()
	
	 
	 for(nFamily in 1:nFamilies_List){ 
		
		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoValSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoValSimValues)[initIndex:finalIndex]
						
		initIndex <- finalIndex+1 
		finalIndex <- finalIndex+cycle1_nProgeny
	 
	 } 
 	
###### Select island pairs based on migration policy 
### SPlit criterion list 
	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- genoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- genoValues_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- phenoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- phenoValues_List
		} 
     } 

#### Migration Policy change from the previous combined set ###################### 	

	#familyPairs <-  migrationPolicy_GM(nextGenGenoTable_GM,SelCriterion_List,Policy,nFamilies_List)
	#nFamilyPairs <- dim(familyPairs)[1]

        if(Policy != "GeneticClusters"){
	  migrationGroups <- migrationPolicy_GM_V2(nextGenGenoTable_GM,genoSelectionTable,Policy,nFamilies,no_QTL)
	
	  emigrantGroups <- migrationGroups[[1]]
	  immigrantGroups <- migrationGroups[[2]]
	
	}else if (Policy == "GeneticClusters"){ 
	
	  migrationGroups <- migrationPolicy_GM_V2(nextGenGenoTable_GM,genoSelectionTable,Policy,nFamilies,no_QTL)
	
	  emigrantGroups <- migrationGroups[[1]]
	  immigrantGroups <- migrationGroups[[2]]
	  emigrantGroup_Prob <- migrationGroups[[3]]
	
	}


        families <-  c(1:nFamilies_List)
	# nFamilyInd <- length(families) 
	
	    
	nextGenGenoTable_GM_List_Fam_MP_Change <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilies)
       nextGenGenoTableMod_GM_List_Fam_MP_Change <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilies)
	
			
	
	for(nFamInd in 1:nFamilies){
	
	   nextGenGenoTable_GM_List_Fam_MP_Change[[nFamInd]] <- nextGenGenoTable_GM_List_Fam[[nFamInd]]
	   nextGenGenoTableMod_GM_List_Fam_MP_Change[[nFamInd]] <- nextGenGenoTableMod_GM_List_Fam[[nFamInd]]
	
	}
				
		
	nextGenGenoTable_GM_List_Fam <- nextGenGenoTable_GM_List_Fam_MP_Change
	nextGenGenoTableMod_GM_List_Fam <- nextGenGenoTableMod_GM_List_Fam_MP_Change
	
###########			
	
	nextGenGenoTable_GM_Fam	<- c()	
		
		
	for(nfam in 1:nFamilies){ 
		  nextGenGenoTable_GM_Fam <-	rbind(nextGenGenoTable_GM_Fam,nextGenGenoTable_GM_List_Fam[[nfam]]) 
	}	
	
	nextGenGenoTableMod_GM_Fam <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_Fam)
	newNextGenGenoTable_GM_Fam <- generateWeightedGenoTable(nextGenGenoTable_GM_Fam,no_QTL)
     
	genoValSimValues_Fam <- simulateGenotypicValues(nextGenGenoTable_GM_Fam,no_QTL)
	phenoValSimValues_Fam <- simulatePhenotypicValues(nextGenGenoTable_GM_Fam,varE,no_QTL)
		

#################get predicted geno and pheno values with JM and unweighted GP models
### Get Predicted geno & pheno values for nextGenGenoData	

    if(modelType=="BayesB" || modelType=="BL"){
		
		if(Weighted==TRUE && WghtMethod=="JM"){
			genoValues  <- PredictBayesPhenoValues_Wght(newNextGenGenoTable_GM_Fam,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno,Freq)
	        }else if(Weighted==TRUE && WghtMethod =="DW"){
		    genoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newNextGenGenoTable_GM_Fam,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)
		    phenoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

		    genoValues <- genoValues_PredList[[1]]
		    phenoValues <- phenoValues_PredList[[1]]

		    alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]
				
		}else if(Weighted==FALSE){
		
		    			
			genoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
			

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

            			
			genoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
		}
    }

 

Criterion_List <- getCriterionValueList_GM_V2(genoValues_Fam,phenoValues_Fam,genoValSimValues_Fam,phenoValSimValues_Fam,cycle_nProgeny_Fam,nFamilies)
	
 
 genoValues_Fam_List <- Criterion_List[[1]]
 phenoValues_Fam_List <- Criterion_List[[2]]
 genoValSimValues_Fam_List <- Criterion_List[[3]]
 phenoValSimValues_Fam_List <- Criterion_List[[4]] 		
		
} 



#### GM method for isolated families


 genoSelected_List_Ind <- getGM_selectedGenoList_SingleFams_V2(nextGenGenoTable_GM_List_Fam,nextGenGenoTableMod_GM_List_Fam,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,families,GM_Param_List)

 
  selectedParents_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[1]]
  selectedParents_NumProgeny_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[2]]
  selectedParents_Geno_Fam <- genoSelected_List_Ind[[1]]
 
  save.image(WorkspaceName)
  gc() 	
  
########################################################################
### Output variables for isolated families
	
  selectedGenoIndividualIndices_List_Fam <-  rep(list(list()),nFamilies)
  selectedPhenoIndividualIndices_List_Fam <- rep(list(list()),nFamilies)
	  
  topNGenoValues_List_Fam <- rep(list(list()),nFamilies)
  GenoTableIndices_topN_List_Fam <- rep(list(list()),nFamilies)
	  
  topNPhenoValues_List_Fam <- rep(list(list()),nFamilies)
  PhenoTableIndices_topN_List_Fam <-rep(list(list()),nFamilies)
  
  
	for(nfam in 1:nFamilies){

			#nfamily <- families[nfam]
			
			topNGenoSimValues <- genoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoSimValues <- phenoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
			topNGenoValues <- genoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoValues <- phenoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
### GM Output in selected individual ids and Num Progeny from selected parents
#### Genotypic values for nextGen Geno data		
   	  	  
			topNGenoValues_List_Fam[[nfam]] <- topNGenoValues
			GenoTableIndices_topN_List_Fam[[nfam]] <- GenoTableIndices_topN
			selectedGenoIndividualIndices_List_Fam[[nfam]] <- GenoTableIndices_topN 
#### Phenotypic values for nextGen Geno data
			
			topNPhenoValues_List_Fam[[nfam]] <- topNPhenoValues
			PhenoTableIndices_topN_List_Fam[[nfam]] <- PhenoTableIndices_topN
			selectedPhenoIndividualIndices_List_Fam[[nfam]] <- PhenoTableIndices_topN
			

#### Output data for cycle- nCyc			
		 
			
	       GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <- genoValSimValues_Fam_List[[nfam]]
	       GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <-  genoValues_Fam_List[[nfam]]
 	       GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfam]] <- topNGenoValues
	       attainedGenoValues_List_Cycle[[nCyc]][[nfam]] <-max(genoValSimValues_Fam_List[[nfam]][GenoTableIndices_topN])
		   
		   
	       PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <-  phenoValSimValues_Fam_List[[nfam]]
	       PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <- phenoValues_Fam_List[[nfam]]
	       PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfam]] <- topNPhenoValues
	       #attainedPhenoValues_List_Cycle[[nCyc]][[nfam]] <- max(phenoValSimValues_Fam_List[[nfam]][PhenoTableIndices_topN])
 
				 
	}

	print(paste(nCyc,"-Fam20-",mean(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]]),sep=""))
#########################

       Cycle_GenoData_List <- Cycle_Progeny_F5_List_Fam
		
     F5RILs_GM_Fam <- getF5RILs_GM_Fam_V2(Cycle_GenoData_List,selectedParents_Geno_Fam_List[[nCyc]],selectedParents_NumProgeny_Geno_Fam_List[[nCyc]],nFamilies,nMarkers,NAM_LinkMap,families,cycle1_nProgeny)
	
		
	
	Cycle_Progeny_F5_Fam <- F5RILs_GM_Fam[[1]]
	Cycle_Progeny_F5_List_Fam <- F5RILs_GM_Fam[[2]]
	
	gc()
	
	
	 
}


######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,attainedPhenoValues_List_Cycle)

	   return(simResults_List)

}

####################### 



runSimulations20X_IslandSelection_BD_Wght_GM_V3 <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Param_List){

	 
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
	  
	  GM_Param_List <- GM_Param_List
	  no_cores <- GM_Param_List$no_cores
	  
	  WorkspaceName <- WorkSpaceName
	  startCycle <- StartCycle

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
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <-  generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable_GM <- generateWeightedGenoTable(cycle1GenoTable_GM,no_QTL)
	cycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(cycle1GenoTable_GM)
	
		
	nIndividuals <- nrow(cycle1GenoTable_GM)
	 

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
# # ### Split geno and pheno data tables in to families		
	 
	  cycle1_nProgeny <- 100
	  
	  
	  	  
	  F5_Progeny_List_Family<- F5_Progeny_List
	  
	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  
	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  
	  for(nFamily in 1:20){ 
		
		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,] 
	  
	  }
	  
# ####################		
	 	  
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	 	  
	  for(nFamily in 1:20){ 
		
		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,] 
		
	  } 
	  
	  
# # ## Split Geno and Pheno values according to families  ########################


	 initIndex <- 1
	 finalIndex <- cycle1_nProgeny 
	 genoValues_List <- list()
	 phenoValues_List <- list()
	 genoSimValues_List <- list()
	 phenoSimValues_List <- list()
	 cycle1_GenoTable_GM_List <- rep(list(matrix(rep(0,cycle1_nProgeny*nMarkers),nrow=cycle1_nProgeny,ncol=nMarkers)),nFamilies)
	 	  
	 cycle1_GenoTableMod_GM_List <- rep(list(matrix(rep(0,cycle1_nProgeny*nMarkers),nrow=cycle1_nProgeny,ncol=nMarkers)),nFamilies)
	  
	 
	 for(nFamily in 1:nFamilies){ 
		
		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]
		
		cycle1_GenoTable_GM_List[[nFamily]] <- cycle1GenoTable_GM[initIndex:finalIndex,]
		cycle1_GenoTableMod_GM_List[[nFamily]] <-cycle1GenoTableMod_GM[initIndex:finalIndex,]
		
		initIndex <- finalIndex+1 
		finalIndex <- finalIndex+cycle1_nProgeny
	 
	 } 

	 
	
###### Select island pairs based on migration policy 
### Split criterion list 
	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- genoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- genoValues_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- phenoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- phenoValues_List
		} 
         } 


### i)  Get pairs of families and isolated families based on migration policy 
 
	families <-  c(1:nFamilies)
    
### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
			selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
            selectedPhenoIndividualIndices <- lapply(phenoSelectionTable,function(x) as.vector(x[,2]))
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
            selectedPhenoIndividualIndices <- lapply(phenoSelectionTable,function(x) as.vector(x[,2]))
	}
	
## selected Line indices to select lines from cycle genotable
    if(selectionOnGeno == TRUE) {
      selectedLineIndices <- selectedGenoIndividualIndices
    } else if(selectionOnGeno==FALSE){ 
      selectedLineIndices <- selectedPhenoIndividualIndices
    }   
	
	
### Migration Policy
	
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
	
###
	
	cycle_nProgeny_Fam <- cycle1_nProgeny
			 
 	selectedParents_Geno_Fam_List <- rep(list(rep(list(),nFamilies)),nCycles) 
	selectedParents_NumProgeny_Geno_Fam_List <- rep(list(rep(list(),nFamilies)),nCycles)
	
### i) Get cycle 1 genotype table and criterion value list for sets of singleislands and pairs of islands	
		
    cycle1GenoTable_GM_Fam <-  generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny_Fam)
    newCycle1GenoTable_GM_Fam <- generateWeightedGenoTable(cycle1GenoTable_GM_Fam,no_QTL)
    cycle1GenoTableMod_GM_Fam <- generate012GenoFormat_Inverse(cycle1GenoTable_GM_Fam)

   
   
    genoValSimValues_Fam <- simulateGenotypicValues(cycle1GenoTable_GM_Fam,no_QTL) 
    phenoValSimValues_Fam <- simulatePhenotypicValues(cycle1GenoTable_GM_Fam,varE,no_QTL) 



if(modelType=="BayesB" || modelType=="BL"){
		
		if(Weighted==TRUE && WghtMethod=="JM"){
			genoValues  <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newCycle1GenoTable,PredictionModel_Pheno,Freq)
	        }else if(Weighted==TRUE && WghtMethod =="DW"){
		    genoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newCycle1GenoTable,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)
		    phenoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newCycle1GenoTable,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

		    genoValues <- genoValues_PredList[[1]]
		    phenoValues <- phenoValues_PredList[[1]]

		    alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]
				
		}else if(Weighted==FALSE){
		
		    genoValues_FamPair <- PredictBayesPhenoValues(newCycle1GenoTable_GM_FamPair,PredictionModel_Geno)
			phenoValues_FamPair <- PredictBayesPhenoValues(newCycle1GenoTable_GM_FamPair,PredictionModel_Pheno)
			
			genoValues_Fam <- PredictBayesPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Pheno)
			
			# genoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Geno)
			# phenoValues <- PredictBayesPhenoValues(newCycle1GenoTable,PredictionModel_Pheno)

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

            			
			genoValues_Fam <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictRRBLUPPhenoValues(newCycle1GenoTable_GM_Fam,PredictionModel_Pheno)
			
		
        }
}


### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

Criterion_List <- getCriterionValueList_GM_V2(genoValues_Fam,phenoValues_Fam,genoValSimValues_Fam,phenoValSimValues_Fam,cycle_nProgeny_Fam,nFamilies)

 genoValues_Fam_List <- Criterion_List[[1]]
 phenoValues_Fam_List <- Criterion_List[[2]]
 genoValSimValues_Fam_List <- Criterion_List[[3]]
 phenoValSimValues_Fam_List <- Criterion_List[[4]]


### Get list of cycle1 geno table for family ind and family pairs ######################################

		
	cycle1GenoTable_GM_List_Fam <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
	cycle1GenoTableMod_GM_List_Fam <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
        selectedCycle1GenoTable_GM_List_Fam <- rep(list(matrix(0,nrow=nSel_inFamily,ncol=nMarkers)),nFamilies)
        selectedCycle1GenoTableMod_GM_List_Fam <- rep(list(matrix(0,nrow=nSel_inFamily,ncol=nMarkers)),nFamilies) 

	
	initIndex <- 1
	finalIndex <- cycle1_nProgeny_Fam
	
	for(nFam in 1:nFamilies){ 
		
		cycle1GenoTable_GM_List_Fam[[nFam]] <- cycle1GenoTable_GM_Fam[initIndex:finalIndex,]
		cycle1GenoTableMod_GM_List_Fam[[nFam]] <- cycle1GenoTableMod_GM_Fam[initIndex:finalIndex,]
                selectedCycle1GenoTable_GM_List_Fam[[nFam]] <- cycle1GenoTable_GM_List_Fam[[nFam]][selectedLineIndices[[nFam]],] 
                selectedCycle1GenoTableMod_GM_List_Fam[[nFam]] <- cycle1GenoTableMod_GM_List_Fam[[nFam]][selectedLineIndices[[nFam]],]
	        initIndex <- finalIndex+1
		finalIndex <- finalIndex + cycle1_nProgeny_Fam
		
	}

    
######################################################################################   
## Output list for all families 
	 attainedGenoValues <- rep(0,nCycles)
	 attainedPhenoValues <- rep(0,nCycles)

	 GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	 GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	 GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 

	 PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles)  
	 PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	 PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
	
     attainedPhenoValues_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 
     attainedGenoValues_List_Cycle <- rep(list(rep(list(list()),nFamilies)),nCycles) 


## Output list for single fam 	  
	 	
		 	 	    
     selectedGenoIndividualIndices_List_Fam <- rep(list(list()),nFamilies)
     selectedPhenoIndividualIndices_List_Fam <- rep(list(list()),nFamilies)

	
#### GM method for parent selection in cycle 1 ################################################
# #### GM method for isolated families
 
 selectedParents_Geno_Fam_List <- list()
 selectedParents_NumProgeny_Geno_Fam_List <- list()
 selectedParents_Geno_Fam <- list()

 system.time({
 
 genoSelected_List_Ind <- getGM_selectedGenoList_SingleFams_V2(selectedCycle1GenoTable_GM_List_Fam,selectedCycle1GenoTableMod_GM_List_Fam,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,families,GM_Param_List)

 selectedParents_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[1]]
 selectedParents_NumProgeny_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[2]]
 selectedParents_Geno_Fam <- genoSelected_List_Ind[[1]]
 
})

 # save.image(WorkspaceName) 	
 #load(WorkspaceName) 
 gc()
 
### Output variables for isolated families
	
	  selectedGenoIndividualIndices_List_Fam <-  rep(list(list()),nFamilies)
	  selectedPhenoIndividualIndices_List_Fam <- rep(list(list()),nFamilies)
	  
	  topNGenoValues_List_Fam <- rep(list(list()),nFamilies)
	  GenoTableIndices_topN_List_Fam <- rep(list(list()),nFamilies)
	  
	  topNPhenoValues_List_Fam <- rep(list(list()),nFamilies)
	  PhenoTableIndices_topN_List_Fam <-rep(list(list()),nFamilies)


	for(nfam in 1:nFamilies){

			#nfamily <- families[nfam]
			
			topNGenoSimValues <- genoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoSimValues <- phenoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
			topNGenoValues <- genoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoValues <- phenoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
### GM Output in selected individual ids and Num Progeny from selected parents
		
   	  	  
			topNGenoValues_List_Fam[[nfam]] <- topNGenoValues
			GenoTableIndices_topN_List_Fam[[nfam]] <- GenoTableIndices_topN
			selectedGenoIndividualIndices_List_Fam[[nfam]] <- GenoTableIndices_topN 
		 
#### Phenotypic values for nextGen Geno data
  
	   		topNPhenoValues_List_Fam[[nfam]] <- topNPhenoValues
			PhenoTableIndices_topN_List_Fam[[nfam]] <- PhenoTableIndices_topN
			selectedPhenoIndividualIndices_List_Fam[[nfam]] <- PhenoTableIndices_topN
			
					
#### Output data for cycle- nCyc			
		 
			
	        GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <- genoValSimValues_Fam_List[[nfam]]
	        GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <-  genoValues_Fam_List[[nfam]]
		    GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfam]] <- topNGenoValues
	        attainedGenoValues_List_Cycle[[nCyc]][[nfam]] <-max(genoValSimValues_Fam_List[[nfam]][GenoTableIndices_topN])

		    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <-  phenoValSimValues_Fam_List[[nfam]]
	        PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <- phenoValues_Fam_List[[nfam]]
		    PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfam]] <- topNPhenoValues
		    attainedPhenoValues_List_Cycle[[nCyc]][[nfam]] <- max(phenoValSimValues_Fam_List[[nfam]][PhenoTableIndices_topN])
 		 
	}


###################################################################################	

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
## Get F5 RILs for selected parents in cycle 1
 
    selectedCycleGenoData_List <- list()
    Cycle_GenoData_List <- cycle1_GenoData_List
		
    for(nFamily in 1:nFamilies){ 
	
	  	selectedCycleGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedLineIndices[[nFamily]],Cycle_GenoData_List[[nFamily]],nSel_inFamily)
	 
    } 
	
	
    F5RILs_GM_Fam <- getF5RILs_GM_Fam_V2(selectedCycleGenoData_List,selectedParents_Geno_Fam_List[[nCyc]],selectedParents_NumProgeny_Geno_Fam_List[[nCyc]],nFamilies,nMarkers,NAM_LinkMap,families,cycle1_nProgeny)
	 
    Cycle_Progeny_F5_Fam <- F5RILs_GM_Fam[[1]]
    Cycle_Progeny_F5_List_Fam <- F5RILs_GM_Fam[[2]]


	save.image(WorkspaceName) 
	gc() 

	startCycle <- startCycle +1 
}
###########################################################################

 startCycle <- 2
 for(nCyc in startCycle:nCycles){

    print(nCyc)
### Exchange selected geno data among pairs
    if(migFreq !=0){
        if(nCyc>=2 && nCyc%%migFreq==0){
	  
	    if(Policy != "GeneticClusters"){
		 
		 output_List_afterExchange <- exchangeGenoData_GM(Cycle_Progeny_F5_List_Fam,genoValSimValues_Fam_List,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
		}else if(Policy == "GeneticClusters"){ 
		
		 
		 output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		}
		  
		Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
		
		familyInfo <- output_List_afterExchange[[2]]
		
		Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 
	   
	   
        Cycle_Progeny_F5_Fam <- array(0,c(nFamilies,nMarkers,2,cycle1_nProgeny))
        nProgIndex <-  dim(Cycle_Progeny_F5_List_Fam[[nFamily]])[3]
	    nProg <-1
				
	    for(nFamily in 1:length(families)){ 
					 
		    Cycle_Progeny_F5_Fam[nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][,,nProg:nProgIndex] 
	    } 
	
      }
	   
	
##    Cycle_Progeny_F5_Fam <-   Cycle_Progeny_F5_List_Fam 

### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
     		
	
	    nextGenGenoTable_GM_Fam <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5_Fam,nFamilies,cycle1_nProgeny_Fam)
	    nextGenGenoTableMod_GM_Fam <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_Fam)
	    newNextGenGenoTable_GM_Fam <- generateWeightedGenoTable(nextGenGenoTable_GM_Fam,no_QTL)
		
	    nextGenGenoTable_GM_List_Fam <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilies)
           nextGenGenoTableMod_GM_List_Fam <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilies)
		
		
		for(nfam in 1:nFamilies){ 
		
				nextGenGenoTable_GM_List_Fam[[nfam]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List_Fam[[nfam]],1,cycle1_nProgeny_Fam)
				nextGenGenoTableMod_GM_List_Fam[[nfam]] <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_List_Fam[[nfam]])
		
		}
		
		genoValSimValues_Fam <- simulateGenotypicValues(nextGenGenoTable_GM_Fam,no_QTL)
		phenoValSimValues_Fam <- simulatePhenotypicValues(nextGenGenoTable_GM_Fam,varE,no_QTL)
		
#########################
## Sclass plot for family clustering 	
#########################	
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
	
############################################################################################### 
###### Get updated or retrained prediction model for geno and pheno values

if(modelUpdate== TRUE && nCyc%%updateFrequency==0){

	nFamilies_List <- length(families) 
	nextGenGeno_nProgeny_FamPair <- cycle1_nProgeny_Fam
	
	nextGenGenoTable_GM <- nextGenGenoTable_GM_Fam
	nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_Fam
	newNextGenGenoTable_GM <- newNextGenGenoTable_GM_Fam
	
	genoValSimValues <- genoValSimValues_Fam
	phenoValSimValues <- phenoValSimValues_Fam

        nIndividuals <- length(genoValSimValues)
	IndividualNames <-rep(0,nIndividuals)

        for(nInd in 1:(nIndividuals)){

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
		trainGenoNewTable <- newNextGenGenoTable_GM[trainIndices,]
		trainSimPhenoTable <- phenotypeSimTable[trainIndices,]
		trainSimGenoTable <-  genotypicValuesSimTable[trainIndices,]

	#### Validation Set table for Genotype and Simulated Phenotype ###################

		# testGenoTable<- genoTable_Mod[testIndices,]
		testGenoNewTable <- newNextGenGenoTable_GM[testIndices,]
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
	
    if(modelUpdate ==TRUE && updateType == "FullSet"){
	
	   PredictionModel_Geno_PreCycle <- PredictionModel_Geno
       PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno

	    indices<-c(1:nIndividuals)

        trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices<- indices[which(!indices %in% trainIndices)]

#### Training Set table for Genotype and Simulated Phenotype ###################

        trainGenoNewTable <- newNextGenGenoTable_GM[trainIndices,]

        trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]

        trainSimGenoValTable <-  genotypicValuesSimTable[trainIndices,]

#### Validation Set table for Genotype and Simulated Phenotype ###################

        testGenoNewTable <- newNextGenGenoTable_GM[testIndices,]

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

                        PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable))

                }else if(modelType=="RRBLUP_REML"){



                        PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        PredictionModel_Pheno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

                }
				
		    }
		 
		    trainGenoNewTablePreCycle <-  trainGenoNewTableComb
            trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)
            trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
			
	if((length(levels(factor(PredictionModel_Pheno[[1]]))) == 1 && sd()==0) ||(length(levels(factor(PredictionModel_Geno[[1]]))) ==1 && sd()==0)){

                         PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         PredictionModel_Geno <- PredictionModel_Geno_PreCycle

       }
			
			

            rm(trainGenoNewTableComb)
            rm(trainSimGenoValTableComb)
            rm(trainSimPhenoValTableComb)
            gc()

        }
}	

    if(modelUpdate ==TRUE && updateType=="TrainingCycleWindow"){

		
			
	PredictionModel_Geno_PreCycle <- PredictionModel_Geno
	PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno
        
	indices<-c(1:nIndividuals)

        trainIndices <-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices <- indices[which(!indices %in% trainIndices)]
		
		nTrainIndices <- length(trainIndices)

#### Training Set table for Genotype and Simulated Phenotype ###################

        trainGenoNewTable <- newNextGenGenoTable_GM[trainIndices,]
        trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]
        trainSimGenoValTable <- genotypicValuesSimTable[trainIndices,]

#### Validation Set table for Genotype and Simulated Phenotype ###################
        testGenoNewTable <- newNextGenGenoTable_GM[testIndices,]
        testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]
        testSimGenoValTable <-  genotypicValuesSimTable[testIndices,]

##### Build RRBLUP prediction model every cycle


        trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))
        trainSimGenoValTableComb <- rbind((trainSimGenoValTablePreCycle),trainSimGenoValTable)
        trainSimPhenoValTableComb <- rbind((trainSimPhenoTablePreCycle),trainSimPhenoTable)
        print(dim(trainGenoNewTableComb))
		Avg_GenoSim_Prev <- mean(trainSimGenoValTableComb[,1])
		Avg_GenoSim_Current <- mean(trainSimGenoValTable[,1])

		change_GenoSim <- Avg_GenoSim_Prev - Avg_GenoSim_Current
		combinedTableCount <- combinedTableCount+1
         print(change_GenoSim)

        if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){

           if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0 && change_GenoSim !=0 ){ 

           if(modelType=="RRBLUP"){

                        PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable))

           }else if(modelType=="RRBLUP_REML"){


                        PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        PredictionModel_Pheno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

           }

	     } 

           print(sd(PredictionModel_Pheno[[1]])) 
           print(sd(PredictionModel_Geno[[1]])) 
           print(length(levels(factor(PredictionModel_Pheno[[1]]))))
           print(length(levels(factor(PredictionModel_Geno[[1]]))))	 

           genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
           phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
           print(sd(phenoValues))
           print(sd(genoValues))

 
           if(length(levels(factor(PredictionModel_Pheno[[1]]))) ==1 || length(levels(factor(PredictionModel_Geno[[1]]))) ==1 || (is.na(sd(phenoValues)) || sd(phenoValues)==0) || (is.na(sd(genoValues)) || sd(genoValues)==0) || change_GenoSim==0){

                         PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         PredictionModel_Geno <- PredictionModel_Geno_PreCycle

           }
	 }
		
		if(combinedTableCount < trainTableWindowSize){
			
			    trainGenoNewTablePreCycle <-  trainGenoNewTableComb

                trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)

                trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
					
		}else if(combinedTableCount >= trainTableWindowSize){
			
			    trainGenoNewTablePreCycle <-  trainGenoNewTableComb[-c(1:nTrainIndices),]

                trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb[-c(1:nTrainIndices),])

                trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb[-c(1:nTrainIndices),])
				
				combinedTableCount <- (trainTableWindowSize-1)
					
		}
		
	
                rm(trainGenoNewTableComb)
                rm(trainSimGenoValTableComb)
                rm(trainSimPhenoValTableComb)
                gc()
       }
}
################################################
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
				
			genoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
	
	    }	
}

if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

        
			genoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
		
        }
}
	

### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

Criterion_List <- getCriterionValueList_GM_V2(genoValues_Fam,phenoValues_Fam,genoValSimValues_Fam,phenoValSimValues_Fam,cycle_nProgeny_Fam,nFamilies)
	
  
 genoValues_Fam_List <- Criterion_List[[1]]
 phenoValues_Fam_List <- Criterion_List[[2]]
 genoValSimValues_Fam_List <- Criterion_List[[3]]
 phenoValSimValues_Fam_List <- Criterion_List[[4]] 
 
### Migration policy change (nextGenGenoTable_GM_FamPair has to be changed for GC policy) 
 
 if( (migPolicyChangeFrequency !=0) && (nCyc %% migPolicyChangeFrequency ==0)){

 	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_Fam_List <- genoValSimValues_Fam_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_Fam_List <- genoValues_Fam_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_Fam_List <- phenoValSimValues_Fam_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_Fam_List <- phenoValues_Fam_List
		} 
         } 
 
        SelCriterion_List <- SelCriterion_Fam_List
	
	nFamilies_List <- length(families)
	
	nextGenGeno_nProgeny_FamPair <- cycle1_nProgeny
	
	
 # genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
# phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
# newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)

        nextGenGenoTable_GM <- nextGenGenoTable_GM_Fam
	nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_Fam
	newNextGenGenoTable_GM <- newNextGenGenoTable_GM_Fam
	
	genoValSimValues <- genoValSimValues_Fam
	phenoValSimValues <- phenoValSimValues_Fam

        nIndividuals <- length(genoValSimValues)
	IndividualNames <-rep(0,nIndividuals)

        for(nInd in 1:(nIndividuals)){

          IndividualNames[nInd]<- paste("Ind",nInd,sep="")

        }
	
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
		
		   			
			genoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
			phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
			
		}	
    }

    if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){

        if(Weighted==TRUE && WghtMethod=="JM"){

                        genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Geno,Freq)

                        phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq)

        }else if (Weighted==TRUE && WghtMethod =="DW"){

                        cycleNumber <- nCyc

                        genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable_GM,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

                        genoValues <- genoValues_PredList[[1]]
                        phenoValues <- phenoValues_PredList[[1]]
                        alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]

        }else if(Weighted==FALSE){

     						
			genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
			phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
				
        }
    }

 
	 initIndex <- 1
	 finalIndex <- cycle1_nProgeny 
	 genoValues_List <- list()
	 phenoValues_List <- list()
	 genoSimValues_List <- list()
	 phenoSimValues_List <- list()
	
	 
	 for(nFamily in 1:nFamilies_List){ 
		
		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoValSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoValSimValues)[initIndex:finalIndex]
						
		initIndex <- finalIndex+1 
		finalIndex <- finalIndex+cycle1_nProgeny
	 
	 } 
 	
###### Select island pairs based on migration policy 
### SPlit criterion list 
	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- genoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- genoValues_List
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- phenoSimValues_List
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- phenoValues_List
		} 
         } 
	 
	 
### Genoselection table & Island migration for cycle1

	if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
			genoSimSelectionTable <- genoSelectionTable
			selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
            selectedPhenoIndividualIndices <- lapply(phenoSelectionTable,function(x) as.vector(x[,2]))
	}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			genoSimSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
			selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
            selectedPhenoIndividualIndices <- lapply(phenoSelectionTable,function(x) as.vector(x[,2]))
	}
	
## selected Line indices to select lines from cycle genotable
    if(selectionOnGeno == TRUE) {
      selectedLineIndices <- selectedGenoIndividualIndices
    } else if(selectionOnGeno==FALSE){ 
      selectedLineIndices <- selectedPhenoIndividualIndices
    }   

#### Migration Policy change from the previous combined set ###################### 	

	#familyPairs <-  migrationPolicy_GM(nextGenGenoTable_GM,SelCriterion_List,Policy,nFamilies_List)
	#nFamilyPairs <- dim(familyPairs)[1]

    if(Policy != "GeneticClusters"){
	  migrationGroups <- migrationPolicy_GM_V2(nextGenGenoTable_GM,genoSelectionTable,Policy,nFamilies,no_QTL)
	
	  emigrantGroups <- migrationGroups[[1]]
	  immigrantGroups <- migrationGroups[[2]]
	
	}else if (Policy == "GeneticClusters"){ 
	
	  migrationGroups <- migrationPolicy_GM_V2(nextGenGenoTable_GM,genoSelectionTable,Policy,nFamilies,no_QTL)
	
	  emigrantGroups <- migrationGroups[[1]]
	  immigrantGroups <- migrationGroups[[2]]
	  emigrantGroup_Prob <- migrationGroups[[3]]
	
	}


    families <-  c(1:nFamilies_List)
	
		
	    
	 nextGenGenoTable_GM_List_Fam_MP_Change <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilies)
         nextGenGenoTableMod_GM_List_Fam_MP_Change <- rep(list(matrix(rep(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers))),nFamilies)
	 selectedNextGenGenoTable_GM_List_Fam_MP_Change <- rep(list(matrix(rep(0,nrow=nSel_inFamily,ncol=nMarkers))),nFamilies)
	 selectedNextGenGenoTableMod_GM_List_Fam_MP_Change <- rep(list(matrix(rep(0,nrow=nSel_inFamily,ncol=nMarkers))),nFamilies)
	
	for(nFamInd in 1:nFamilies){
	
	   nextGenGenoTable_GM_List_Fam_MP_Change[[nFamInd]] <- nextGenGenoTable_GM_List_Fam[[nFamInd]]
	   nextGenGenoTableMod_GM_List_Fam_MP_Change[[nFamInd]] <- nextGenGenoTableMod_GM_List_Fam[[nFamInd]]
	   selectedNextGenGenoTable_GM_List_Fam_MP_Change[[nFamInd]] <- nextGenGenoTable_GM_List_Fam_MP_Change[[nFamInd]][selectedLineIndices[[nFamily]],]
	   selectedNextGenGenoTableMod_GM_List_Fam_MP_Change[[nFamInd]] <-  nextGenGenoTableMod_GM_List_Fam_MP_Change[[nFamInd]][selectedLineIndices[[nFamily]],]
	}
	
		
	nextGenGenoTable_GM_List_Fam <- nextGenGenoTable_GM_List_Fam_MP_Change
	nextGenGenoTableMod_GM_List_Fam <- nextGenGenoTableMod_GM_List_Fam_MP_Change
	selectedNextGenGenoTable_GM_List_Fam <- selectedNextGenGenoTable_GM_List_Fam_MP_Change
	selectedNextGenGenoTableMod_GM_List_Fam <- selectedNextGenGenoTableMod_GM_List_Fam_MP_Change

###########			
	
	nextGenGenoTable_GM_Fam	<- c()	
		
		
	for(nfam in 1:nFamilies){ 
		  nextGenGenoTable_GM_Fam <- rbind(nextGenGenoTable_GM_Fam,nextGenGenoTable_GM_List_Fam[[nfam]]) 
	}	
	
	nextGenGenoTableMod_GM_Fam <- generate012GenoFormat_Inverse(nextGenGenoTable_GM_Fam)
	newNextGenGenoTable_GM_Fam <- generateWeightedGenoTable(nextGenGenoTable_GM_Fam,no_QTL)
     
	genoValSimValues_Fam <- simulateGenotypicValues(nextGenGenoTable_GM_Fam,no_QTL)
	phenoValSimValues_Fam <- simulatePhenotypicValues(nextGenGenoTable_GM_Fam,varE,no_QTL)
		

#################get predicted geno and pheno values with JM and unweighted GP models
### Get Predicted geno & pheno values for nextGenGenoData	

    if(modelType=="BayesB" || modelType=="BL"){
		
		if(Weighted==TRUE && WghtMethod=="JM"){
			genoValues  <- PredictBayesPhenoValues_Wght(newNextGenGenoTable_GM_Fam,PredictionModel_Geno,Freq)
			phenoValues <- PredictBayesPhenoValues_Wght(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno,Freq)
	        }else if(Weighted==TRUE && WghtMethod =="DW"){
		    genoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newNextGenGenoTable_GM_Fam,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)
		    phenoValues_PredList <- PredictBayesPhenoValues_WGS_DWModel(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

		    genoValues <- genoValues_PredList[[1]]
		    phenoValues <- phenoValues_PredList[[1]]

		    alphaPar_Vector[nCyc] <- genoValues_PredList[[2]]
				
		}else if(Weighted==FALSE){
		
		    			
			genoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
			

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

            			
			genoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
		}
    }

 

Criterion_List <- getCriterionValueList_GM_V2(genoValues_Fam,phenoValues_Fam,genoValSimValues_Fam,phenoValSimValues_Fam,cycle_nProgeny_Fam,nFamilies)
	
 
 genoValues_Fam_List <- Criterion_List[[1]]
 phenoValues_Fam_List <- Criterion_List[[2]]
 genoValSimValues_Fam_List <- Criterion_List[[3]]
 phenoValSimValues_Fam_List <- Criterion_List[[4]] 		
		
} 


#### GM method for isolated families


 genoSelected_List_Ind <- getGM_selectedGenoList_SingleFams_V2(selectedNextGenGenoTable_GM_List_Fam,selectedNextGenGenoTableMod_GM_List_Fam,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,families,GM_Param_List)

 
  selectedParents_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[1]]
  selectedParents_NumProgeny_Geno_Fam_List[[nCyc]] <- genoSelected_List_Ind[[2]]
  selectedParents_Geno_Fam <- genoSelected_List_Ind[[1]]
 
  save.image(WorkspaceName)
  gc() 	
  
########################################################################
### Output variables for isolated families
	
  selectedGenoIndividualIndices_List_Fam <-  rep(list(list()),nFamilies)
  selectedPhenoIndividualIndices_List_Fam <- rep(list(list()),nFamilies)
	  
  topNGenoValues_List_Fam <- rep(list(list()),nFamilies)
  GenoTableIndices_topN_List_Fam <- rep(list(list()),nFamilies)
	  
  topNPhenoValues_List_Fam <- rep(list(list()),nFamilies)
  PhenoTableIndices_topN_List_Fam <-rep(list(list()),nFamilies)
  
  
	for(nfam in 1:nFamilies){

			#nfamily <- families[nfam]
			
			topNGenoSimValues <- genoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoSimValues <- phenoValSimValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
			topNGenoValues <- genoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			GenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]] 
			   
			topNPhenoValues <- phenoValues_Fam_List[[nfam]][selectedParents_Geno_Fam[[nfam]]]
			PhenoTableIndices_topN <- selectedParents_Geno_Fam[[nfam]]
		
### GM Output in selected individual ids and Num Progeny from selected parents
#### Genotypic values for nextGen Geno data		
   	  	  
			topNGenoValues_List_Fam[[nfam]] <- topNGenoValues
			GenoTableIndices_topN_List_Fam[[nfam]] <- GenoTableIndices_topN
			selectedGenoIndividualIndices_List_Fam[[nfam]] <- GenoTableIndices_topN 
#### Phenotypic values for nextGen Geno data
			
			topNPhenoValues_List_Fam[[nfam]] <- topNPhenoValues
			PhenoTableIndices_topN_List_Fam[[nfam]] <- PhenoTableIndices_topN
			selectedPhenoIndividualIndices_List_Fam[[nfam]] <- PhenoTableIndices_topN
			

#### Output data for cycle- nCyc			
		 
			
	       GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <- genoValSimValues_Fam_List[[nfam]]
	       GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <-  genoValues_Fam_List[[nfam]]
 	       GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfam]] <- topNGenoValues
	       attainedGenoValues_List_Cycle[[nCyc]][[nfam]] <-max(genoValSimValues_Fam_List[[nfam]][GenoTableIndices_topN])
		   
		   
	       PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <-  phenoValSimValues_Fam_List[[nfam]]
	       PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nfam]] <- phenoValues_Fam_List[[nfam]]
	       PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nfam]] <- topNPhenoValues
	       #attainedPhenoValues_List_Cycle[[nCyc]][[nfam]] <- max(phenoValSimValues_Fam_List[[nfam]][PhenoTableIndices_topN])
 
				 
	}

	print(paste(nCyc,"-Fam20-",mean(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]]),sep=""))
#########################

       Cycle_GenoData_List <- Cycle_Progeny_F5_List_Fam
       selectedCycleGenoData_List <- list()
   
	 for(nFamily in 1:nFamilies){ 
	
	   	selectedCycleGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedLineIndices[[nFamily]],Cycle_GenoData_List[[nFamily]],nSel_inFamily)
		
          } 
	 
			
      F5RILs_GM_Fam <- getF5RILs_GM_Fam_V2(selectedCycleGenoData_List,selectedParents_Geno_Fam_List[[nCyc]],selectedParents_NumProgeny_Geno_Fam_List[[nCyc]],nFamilies,nMarkers,NAM_LinkMap,families,cycle1_nProgeny)
	
		
	
	Cycle_Progeny_F5_Fam <- F5RILs_GM_Fam[[1]]
	Cycle_Progeny_F5_List_Fam <- F5RILs_GM_Fam[[2]]
	
	gc()
	
	
	 
}


######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle)

	   return(simResults_List)

}












