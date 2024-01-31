
#####################################################################

runSimulations20X_BD_Wght_GM <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,WorkSpaceName,StartCycle,GM_Param_List){

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
       
## BD Parameter	  
	  BD <- BreedDesign
	  

##  Model Update Parameters
	  
      modelRetrain <- ModelRetrain
	  retrainFrequency <-  RetrainFrequency
	  modelUpdate <- ModelUpdate
      updateFrequency <- UpdateFrequency
	  modelType <- ModelType
      updateType <- UpdateType
	  
	   
	  trainGenoNewTablePreCycle <- as.big.matrix(trainGenoNewTable)
      trainSimPhenoTablePreCycle <- trainSimPhenoValTable
      trainSimGenoValTablePreCycle <- trainSimGenoValTable

	  trainTableWindowSize <- TrainTableWindowSize
	  
	  
## Weighting Parameters	  
	  Weighted <- weighted
	  WghtMethod <- wghtMethod

	  
## GM Parameters 
	  
	  GM_Param_List <- GM_Param_List
	  no_cores <- GM_Param_List$no_cores
	  startCycle <- StartCycle
	 
#################################################################################################
#################################################################################################
### Initialize variables	
	
   	  j<-1
 	  selectedParents_Geno_List <- list() 
	  selectedParents_Pheno_List <- list()
	  selectedParents_NumProgeny_Geno_List <- list()
	  selectedParents_NumProgeny_Pheno_List <- list() 

# ###############################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
	  combinedTableCount <- 1
# ###############################################################################################	  	  	  
# ### Get predicted geno and pheno values for cycle1 ###################################################
    
	  nCyc <-1
	  print(nCyc)
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <-  generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable_GM <- generateWeightedGenoTable(cycle1GenoTable_GM,no_QTL)
	cycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(cycle1GenoTable_GM)
	  
	   
	genoValSimValues <- simulateGenotypicValues(cycle1GenoTable_GM,no_QTL) 
	phenoValSimValues <- simulatePhenotypicValues(cycle1GenoTable_GM,varE,no_QTL) 


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

#########################################################################################################

# ##### Select island pairs based on migration policy 
### Split criterion list 
	  
	 if(selectionOnGeno==TRUE){ 
	 
	    if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- genoValSimValues
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <- genoValues
		} 
	 }else if(selectionOnGeno==FALSE){ 
		if(selectionOnSimulated==TRUE){ 
			SelCriterion_List <- phenoValSimValues
		}else if (selectionOnSimulated==FALSE){
			SelCriterion_List <-  phenoValues
		} 
     } 

 
	  selectedGenoIndividualIndices2 <- rep(0,no_selected) 
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)
	 
################################################################################### 
## Output list 	  
	  
	  GenoVal_NX_N_3c_List_Cycle <- rep(list(list()),nCycles) 
	  GenoVal_NX_2k_3c_List_Cycle <- rep(list(list()),nCycles)
	  GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(list()),nCycles)
	  
	 
	  PhenoVal_NX_N_3c_List_Cycle <- rep(list(list()),nCycles) 
	  PhenoVal_NX_2k_3c_List_Cycle <- rep(list(list()),nCycles)
	  PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(list()),nCycles)
	  
	  
      attainedPhenoValues_List_Cycle <- rep(list(list(),nCycles))
      attainedGenoValues_List_Cycle <- rep(list(list(),nCycles))
	
	  topNGenoValues_List_Cycle <- rep(list(list(),nCycles))
      GenoTableIndices_topN_List_Cycle <- rep(list(list(),nCycles))

      topNPhenoValues_List_Cycle <- rep(list(list(),nCycles))
      PhenoTableIndices_topN_List_Cycle <- rep(list(list(),nCycles))
	 
	
# #################################################################################
#### GM method for parent selection in cycle 1 #####################################
 selectedParents_Geno_List <- list()
 selectedParents_NumProgeny_Geno_List <- list()
 selectedParents_Geno <- list()

 system.time({
 
 genoSelected_List <- getGM_FS_selectedGenoList(cycle1GenoTable_GM,cycle1GenoTableMod_GM,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List)

 selectedParents_Geno_List[[nCyc]] <- genoSelected_List[[1]]
 selectedParents_NumProgeny_Geno_List[[nCyc]] <- genoSelected_List[[2]]
 selectedParents_Geno <- genoSelected_List[[1]]
 
})

 save.image(WorkspaceName) 	
 gc()
 
     
 
	  topNGenoSimValues <- genoValSimValues[selectedParents_Geno]
	  GenoTableIndices_topN <- selectedParents_Geno 
			   
	  topNPhenoSimValues <- phenoValSimValues[selectedParents_Geno]
	  PhenoTableIndices_topN <- selectedParents_Geno
		
	  topNGenoValues <- genoValues[selectedParents_Geno]
	  GenoTableIndices_topN <- selectedParents_Geno 
			   
	  topNPhenoValues <- phenoValues[selectedParents_Geno]
	  PhenoTableIndices_topN <- selectedParents_Geno
		
### GM Output in selected individual ids and Num Progeny from selected parents

############################

	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]] <- genoValues
	    attainedGenoValues_List_Cycle[[nCyc]] <- max(genoValSimValues[GenoTableIndices_topN])
		
	    topNGenoValues_List_Cycle[[nCyc]] <- topNGenoValues_List
	    GenoTableIndices_topN_List_Cycle[[nCyc]] <- GenoTableIndices_topN_List
	   
	    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]] <-  phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
	    topNPhenoValues_List_Cycle[[nCyc]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]] <-  PhenoTableIndices_topN
	    
			  

###################################################################################	
## Get F5 RILs for selected parents in cycle 1
 
   
    nCrosses <- dim(cycle1_Geno_Data)[1]
 
    Cycle_GenoData  <- rep(list(array(0,c(1,nMarkers,2,cycle1_nProgeny))),nCrosses)

    for(nCrs in 1: nCrosses){

	Cycle_GenoData[[nCrs]][1,,,] <- cycle1_Geno_Data[nCrs,,,]
    }

   
    F5RILs_GM <- getF5RILs_FS_GM(Cycle_GenoData,selectedParents_Geno_List[[nCyc]],selectedParents_NumProgeny_Geno_List[[nCyc]],selectionOnGeno,nMarkers,NAM_LinkMap)

    Cycle_Progeny_F5 <- F5RILs_GM[[1]]
    nCrosses <- F5RILs_GM[[2]]
    nProgenyVector  <- F5RILs_GM[[3]]

    startCycle <- startCycle + 1 	 
	save.image(WorkspaceName) 
	gc() 
	
###########################################################################
nCrosses <- length(nProgenyVector)
startCycle <- 2
for(nCyc in startCycle:nCycles){
        #nCyc <- 2

        print(nCyc)

### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
    		 
		nextGenGenoTable_GM <- generateMinus101GenoFormat_GM_V1(Cycle_Progeny_F5,nCrosses,nProgenyVector)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		
        newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
	
############################################################################################### 
#####################

if(nCyc%%updateFrequency==0){

    nIndividuals <- 2000
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

        trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        testIndices<- indices[which(!indices %in% trainIndices)]
		
		nTrainIndices <- length(trainIndices)

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

 
#	   if((sd(phenoValues)==0) ||(sd(genoValues)==0) || change_GenoSim==0 ){
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

		

### Split Geno and Pheno values according to families  ########################
#### GM method for parent selection in cycle nCyc #####################################
 

 genoSelected_List <- getGM_FS_selectedGenoList(nextGenGenoTable_GM,nextGenGenoTableMod_GM,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List)

 
  selectedParents_Geno_List[[nCyc]] <- genoSelected_List[[1]]
  selectedParents_NumProgeny_Geno_List[[nCyc]] <- genoSelected_List[[2]]
  selectedParents_Geno <- genoSelected_List[[1]]
 
  save.image(WorkspaceName)
  gc() 	
  
#####################################################################################################################

	   topNGenoSimValues <- genoValSimValues[selectedParents_Geno]
	   GenoTableIndices_topN <- selectedParents_Geno 
			   
	   topNPhenoSimValues <- phenoValSimValues[selectedParents_Geno]
	   PhenoTableIndices_topN <- selectedParents_Geno
		
	   topNGenoValues <- genoValues[selectedParents_Geno]
	   GenoTableIndices_topN <- selectedParents_Geno 
			   
	   topNPhenoValues <- phenoValues[selectedParents_Geno]
	   PhenoTableIndices_topN <- selectedParents_Geno
		
### GM Output in selected individual ids and Num Progeny from selected parents
		   		
	
	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]] <-   genoValues
	    attainedGenoValues_List_Cycle[[nCyc]] <-  max(genoValSimValues[GenoTableIndices_topN])
		
		topNGenoValues_List_Cycle[[nCyc]] <-  topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]] <- GenoTableIndices_topN
		
		PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]] <- phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
		topNPhenoValues_List_Cycle[[nCyc]] <-  topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]] <- PhenoTableIndices_topN
		
			
             print(paste(nCyc,"-",mean(genoValSimValues),sep=""))	
#########################
	
	Cycle_GenoData <- Cycle_Progeny_F5
     	
    F5RILs_GM <- getF5RILs_FS_GM(Cycle_GenoData,selectedParents_Geno_List[[nCyc]],selectedParents_NumProgeny_Geno_List[[nCyc]],selectionOnGeno,nMarkers,NAM_LinkMap)	

	
	Cycle_Progeny_F5 <- F5RILs_GM[[1]]
	nCrosses <- F5RILs_GM[[2]]
	nProgenyVector  <- F5RILs_GM[[3]]
	
	if(nCyc %% 1==0){
	
		save.image(WorkspaceName) 
		gc()
	}
 
}

######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,attainedPhenoValues_List_Cycle)

	   return(simResults_List)

} 



 getGM_FS_selectedGenoList <- function(nextGenGenoTable_GM_List,nextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List){

	nextGenGenoTable_GM <- nextGenGenoTable_GM_List 
	nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_List 
	
	selectedParents_Geno_List <- list()
	selectedParents_Pheno_List <- list()
	  
	selectedParents_NumProgeny_Geno_List <- list()
	selectedParents_NumProgeny_Pheno_List <- list() 
	
### GM_Parameters
	Minparents= GM_Param_List$minparents
	Impinbreedstepsize= GM_Param_List$impinbreedstepsize 
	Impvar= GM_Param_List$impvar
	Mc.cores= GM_Param_List$no_cores
	No_cores = GM_Param_List$no_cores
	Mutprob= GM_Param_List$mutprob
	NpopGA= GM_Param_List$npopGA
	NitGA= GM_Param_List$nitGA
	Impforinbreed= GM_Param_List$impforinbreed
	Plotiters= GM_Param_List$plotiters
	Nelite= GM_Param_List$nelite
	Method= GM_Param_List$method
	PlotMates= GM_Param_List$plotMates
	
#######
 
	M <- nextGenGenoTable_GM
	K <- GenomicMatingV2::Amat.pieces(M,pieces=15,mc.cores=No_cores) 
	    
	Markers <- nextGenGenoTableMod_GM
	
	if(modelType=="BayesB" || modelType=="BL"){
		   markerEffects_Geno <- as.vector(unlist(PredictionModel_Geno[[3]]))
		   markerEffects_Pheno <- as.vector(unlist(PredictionModel_Pheno[[3]]))
    } 
	if(modelType =="RRBLUP" || modelType =="RRBLUP_REML"){
		   markerEffects_Geno <- as.vector(unlist(PredictionModel_Geno[[1]]))
		   markerEffects_Pheno <- as.vector(unlist(PredictionModel_Pheno[[1]]))
	} 
		
					
						
## selection on genotypic values 
		
		if(selectionOnGeno ==TRUE){
			baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutions(Markers,Markers2=NULL,K,markerEffects_Geno,minparents=Minparents,impinbreedstepsize=Impinbreedstepsize,impvar=Impvar,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,impforinbreed=Impforinbreed,plotiters=Plotiters,nelite=Nelite,method=Method,plotMates=PlotMates)
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Geno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
			
			
			tableselected_Geno <-table(factor(unlist(c(baseGASolns_GenoModel[[1]])), levels=1:N))
			tableselected2_Geno<-table(factor(paste(baseGASolns_GenoModel[[1]][,1],baseGASolns_GenoModel[[1]][,2], sep="_")))
			tableselected2_Geno<- sort(tableselected2_Geno[tableselected2_Geno>0], decreasing=T)
	 	 
			selectedParents_Geno <- as.numeric(unlist(strsplit(names(tableselected2_Geno),"_")))
		
			selectedParents_NumProgeny_Geno <- c() 
		
			for(nGeno in 1:length(tableselected2_Geno)){
				selectedParents_NumProgeny_Geno <- c(selectedParents_NumProgeny_Geno,rep(tableselected2_Geno[nGeno],2))
			} 

			selectedParents_Geno_List <- selectedParents_Geno
		
			selectedParents_NumProgeny_Geno_List <- as.vector(selectedParents_NumProgeny_Geno)
			
		}else if(selectionOnGeno==FALSE){
		
## selection on phenotypic values 
			baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutions(Markers,Markers2=NULL,K,markerEffects_Pheno,minparents=Minparents,impinbreedstepsize=Impinbreedstepsize,impvar=Impvar,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,impforinbreed=Impforinbreed,plotiters=Plotiters,nelite=Nelite,method=Method,plotMates=PlotMates) 
			
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Pheno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
			
## gets the list of pairs 
			tableselected_Geno <- table(factor(unlist(c(baseGASolns_GenoModel[[1]])), levels=1:N))
		    tableselected2_Geno<- table(factor(paste(baseGASolns_GenoModel[[1]][,1],baseGASolns_GenoModel[[1]][,2], sep="_")))
		    tableselected2_Geno<- sort(tableselected2_Geno[tableselected2_Geno>0], decreasing=T)
	 
			selectedParents_Geno <- as.numeric(unlist(strsplit(names(tableselected2_Geno),"_")))
			selectedParents_NumProgeny_Geno <- c() 
		
## gets the list of num_progeny per pair
			for(nGeno in 1:length(tableselected2_Geno)){
				selectedParents_NumProgeny_Geno <- c(selectedParents_NumProgeny_Geno,rep(tableselected2_Geno[nGeno],2))
			} 
		    selectedParents_Geno_List <- selectedParents_Geno
			selectedParents_NumProgeny_Geno_List <- as.vector(selectedParents_NumProgeny_Geno)
			
			
		}
		
	
		return(list(selectedParents_Geno_List,selectedParents_NumProgeny_Geno_List))
	
	
}


 
getF5RILs_FS_GM <- function(GenoData_List,SelectedParents_Geno_List,SelectedParents_NumProgeny_Geno_List,SelectionOnGeno,nMarkers,NAM_LinkMap){

	    
                BD <- "GM"
		        NAM_LinkMap_New <- NAM_LinkMap
                GenoData_List <- GenoData_List
                selectedParents_Geno_List <- SelectedParents_Geno_List
		        selectedParents_NumProgeny_Geno_List <- SelectedParents_NumProgeny_Geno_List
                selectionOnGeno <- SelectionOnGeno

##### Generate F5 RIL Progeny from selected parent set for all families 
#####	Initialize lists for progeny data for all families every cylce 
		
		nProgeny_list <- list()
		nSel_inFamily_list <- list()
		
		selectedGenoData_List <- extractGenoData_GM(selectedParents_Geno_List,GenoData_List)


                if(selectionOnGeno==TRUE){

				nProgeny_list <- (selectedParents_NumProgeny_Geno_List)

				nSel_inFamily_list <-  (length(selectedParents_NumProgeny_Geno_List))/2
				
		}else if(selectionOnGeno==FALSE){

				nProgeny_list <- (selectedParents_NumProgeny_Geno_List)

				nSel_inFamily_list <-  (length(selectedParents_NumProgeny_Geno_List))/2
		}	

					
		nSel_inFamily <- nSel_inFamily_list
		 
		nParents <- nSel_inFamily*2
		 
		nProgenyVec <- nProgeny_list
		nProgenyVecPerPair <- nProgenyVec[seq(1,length(nProgenyVec),by=2)]
		
		nPairs <- 1
			
### Loop througn selected parent pairs to generate progeny for each family
			
		Cycle_Progeny_F5_temp1_List <- rep(list(array(0,c(1,nMarkers,2,nProgeny))),nSel_inFamily)	

                              while(nPairs < (nParents)) {
			
				nProgeny <- nProgenyVec[nPairs]
													 
				Cycle_Progeny_F1 <- array(0,c(1,nMarkers,2,1)) 
				
				Cycle_Progeny_F2 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F3 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F4 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F5_temp1 <- array(0,c(1,nMarkers,2,nProgeny)) 

			
	###################################################################################
		 
		### Get Parent pair geno data
				Parent1<-  selectedGenoData_List[nPairs,,]
				Parent2<-  selectedGenoData_List[nPairs+1,,]
			  
				j<-1
		
		###  F1	progeny
		
				Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			  
		
		###  F2 progeny
				Parent1<- Cycle_Progeny_F1[j,,,1]  
				Parent2<- Cycle_Progeny_F1[j,,,1] 
			
				progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New) 
			 
				Cycle_Progeny_F2[j,,,]  <- progeny1 
			  
		###  F3 Progeny	  
			  for(m in 1:nProgeny){ 
				
				Parent1<- Cycle_Progeny_F2[j,,,m]  
				Parent2<- Cycle_Progeny_F2[j,,,m] 
				progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New) 
				Cycle_Progeny_F3[j,,,m] <- progeny1 
			  }
			  
		### F4 Progeny   
			  for(m in 1:nProgeny){ 
				
				Parent1<- Cycle_Progeny_F3[j,,,m]  
				Parent2<- Cycle_Progeny_F3[j,,,m] 
				progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New) 
				Cycle_Progeny_F4[j,,,m] <- progeny1 
			  }
			  
		### F5 Progeny	  
			  
			  for(m in 1:nProgeny){ 
				
				Parent1<- Cycle_Progeny_F4[j,,,m]  
				Parent2<- Cycle_Progeny_F4[j,,,m] 
				progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New) 
				Cycle_Progeny_F5_temp1[j,,,m] <- progeny1 
			  }
			  
			   Cycle_Progeny_F5_temp1_List[[(nPairs+1)/2]] <- Cycle_Progeny_F5_temp1
			
				nPairs <- nPairs+2
			 
			}
			
			   Cycle_Progeny_F5_List_Temp <- Cycle_Progeny_F5_temp1_List

	          
	   
	   
### Get Cycle_Progeny F5 data in correct format 
      
   	  
	                   Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Temp
       	 
	 		

	    return(list(Cycle_Progeny_F5,nSel_inFamily,nProgenyVecPerPair)) 
  
  }

##################################################################################################

generateMinus101GenoFormat_GM_V1 <- function(Cycle_Progeny_Data,number_of_Crosses,n_ProgenyVector){

  Cycle_Progeny <- Cycle_Progeny_Data
  nCrosses <- number_of_Crosses
  n_ProgenyVec <- n_ProgenyVector
 
  nMarkers <- 4289
 
  Cycle_Progeny_table <-c()
  
  for(i in 1:nCrosses){
  
    nProgeny <- n_ProgenyVec[i]
  
    cycle_family_table<- array(0,c(1,nMarkers,nProgeny))
   
    for(k in 1:nProgeny){


		cycle_family_table[1,,k] <- (Cycle_Progeny[[i]][1,,1,k] + Cycle_Progeny[[i]][1,,2,k])

	}

    Cycle_Progeny_table <-cbind(Cycle_Progeny_table,cycle_family_table[1,,])

  }

  ##IndividualNames <-rep(0,(nCrosses*nProgeny))

  #for(i in 1:(nCrosses*nProgeny)){
  #  IndividualNames[i]<- paste("Ind",i,sep="")
  # }

  # IndividualNames <- c("Animal_ID",IndividualNames)
  #### RowNames
  markerNames <-rep(0,nMarkers)

  for(i in 1:(nMarkers)){
    markerNames[i]<- paste("m",i,sep="")
  }

  #colnames(Cycle_Progeny_table)<- IndividualNames
  row.names(Cycle_Progeny_table)<- markerNames

  genoTable<- t(Cycle_Progeny_table)

  ########## Translate 0-1-2 code to -1-0-1 code in genotable

  genoTable_Mod <- apply(genoTable,2,function(x)(gsub("\\b0\\b","minus1",as.matrix(x))))
  genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b1\\b","0",as.matrix(x))))
  genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b2\\b","1",as.matrix(x))))
  genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\bminus1\\b","-1",as.matrix(x))))
  genoTable_Mod <- apply(genoTable_Mod,2,as.numeric)


  return(genoTable_Mod)
}




extractGenoData_GM <- function(GenoTableIndices1,ProgenyGenoData){
  
        GenoTableIndices <- GenoTableIndices1
	Cycle_GenoData <- ProgenyGenoData
	
        nCrosses <- length(Cycle_GenoData)

## combine progeny data from all crosses

        nPoplnSize <- 2000

        comb_GenoData <- array(0,c(1,nMarkers,2,nPoplnSize))
  
        initIndex <- 1
        finalIndex <- dim(Cycle_GenoData[[1]])[4]
 
       for(nCrs in 1:nCrosses){
   
          comb_GenoData[1,,,initIndex:finalIndex] <- Cycle_GenoData[[nCrs]][1,,,] 

          initIndex <- finalIndex +1
          
        if(nCrs < nCrosses){
           finalIndex <- initIndex+(dim(Cycle_GenoData[[(nCrs+1)]])[4])-1
          }
       }

## select lines from combined progeny geno data 
	nSelectedIndividuals <- length(GenoTableIndices) 
  
	selected_GenoData <- array(0,c(nSelectedIndividuals,nMarkers,2))
	
					
	for(j in 1:nSelectedIndividuals){
		
		individualIndex <- GenoTableIndices[j]
			
		selected_GenoData[j,,] <- comb_GenoData[1,,,individualIndex] 
		
	}

## return seleced geno data in chr format 
		
	# dimnames(selected_GenoData)[[1]] <- as.vector(GenoTableIndices)
		
	  
	return(selected_GenoData) 
	
}
        
################################################################################################3
