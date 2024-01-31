## Run simulations for GM parameter set

runSimulations20X_BD_Wght_GM_Frontier <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,WorkSpaceName,StartCycle,GM_Frontier_Param_List){

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
	  
	  GM_Param_List <- GM_Frontier_Param_List
	  no_cores <- GM_Param_List$no_cores
	  startCycle <- StartCycle
	  WorkspaceName <- WorkSpaceName
	 
#################################################################################################
#################################################################################################
### Initialize variables	
	
   	  j<-1
 	 
################################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
	  combinedTableCount <- 1
	  
################################################################################################	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
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
	 

# #################################################################################
    nCrosses <- dim(cycle1_Geno_Data)[1]
    nProgeny <- dim(cycle1_Geno_Data)[4]
	
    Cycle_GenoData  <- (array(0,c(1,nMarkers,2,cycle1_nProgeny*nCrosses)))
    initIndex <- 1
    finalIndex <- nProgeny
   
    for(nCrs in 1:nCrosses){

	
	 Cycle_GenoData[1,,,initIndex:finalIndex] <- cycle1_Geno_Data[nCrs,,,] 
	 initIndex <- initIndex + nProgeny
	 finalIndex <- finalIndex + nProgeny
	}

#### GM method for parent selection in cycle 1 #####################################
 nParameters <- 30 
 selectedParents_Geno_List <- rep(list(rep(list(list()),nParameters)),nCycles)
 CriterionValue_List <- rep(list(rep(list(list()),nParameters)),nCycles)

 system.time({
 
 genoSelected_List <- getGM_Frontier_FS_selectedGenoList(cycle1GenoTable_GM,cycle1GenoTableMod_GM,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List)

 selectedParents_Geno_List[[nCyc]] <- genoSelected_List[[1]]
 CriterionValue_List[[nCyc]] <- genoSelected_List[[2]]
 selectedParents_Geno <- genoSelected_List[[1]]
 nParameters <- genoSelected_List[[3]]
 
})

 save.image(WorkspaceName) 	
 gc()
 
################################################################################### 
## Initiate Output list 	  
	  
	  GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	 
	  PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  
      attainedPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      attainedGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	
	  topNGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      GenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)

      topNPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      PhenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  Cycle_Progeny_F5_List <- list()
      nCrosses_List <- list()
	  
	  
	 
##
 
   for(nParameter in 1:nParameters){
   
        uniqueParentIDs <-  unique(as.vector(selectedParents_Geno[[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
### GM Output in selected individual ids 

############################

	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]]<- genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <- max(genoValSimValues[GenoTableIndices_topN])
		
	    topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
	   
	    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
	    topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <-  PhenoTableIndices_topN
		PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		
	    
		 print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
###################################################################################	
## Get F5 RILs for selected parents in cycle 1
 
      
        F5RILs_GM <- getF5RILs_FS_GM_Frontier(Cycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New)

        Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
    
	}

    startCycle <- startCycle + 1 	 
	save.image(WorkspaceName) 
	gc() 
	
###########################################################################

    startCycle <- 2
    nCores <- 10
    cl <- makeCluster(nCores,outfile='log.txt')
    registerDoParallel(cl)

    nextGenGenoTable_GM_List <- list()
    nextGenGenoTableMod_GM_List <- list()
    newNextGenGenoTable_GM_List <- list()
	PredictionModel_Geno_List <- list()
	PredictionModel_Pheno_List <- list()
	genoValues_List <- list()
    phenoValues_List <- list()
	genoValSimValues_List <- list()
	phenoValSimValues_List <- list()


for(nCyc in startCycle:nCycles){
        #nCyc <- 2

    print(nCyc)

### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
   for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <-  Cycle_Progeny_F5_List[[nParameter]]
		nCrosses <-  nCrosses_List[[nParameter]] 
		
		nextGenGenoTable_GM <- generateMinus101GenoFormat_GM_Frontier_V1(Cycle_Progeny_F5,nCrosses)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
	
############################################################################################### 

    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0 && nParameter ==1){

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



                        PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

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

        nextGenGenoTable_GM_List[[nParameter]] <- nextGenGenoTable_GM
        nextGenGenoTableMod_GM_List[[nParameter]] <- nextGenGenoTableMod_GM
		newNextGenGenoTable_GM_List[[nParameter]] <- newNextGenGenoTable_GM 
		PredictionModel_Geno_List[[nParameter]] <- PredictionModel_Geno
		PredictionModel_Pheno_List[[nParameter]] <- PredictionModel_Pheno
		genoValues_List[[nParameter]] <- genoValues
        phenoValues_List[[nParameter]] <- phenoValues
		genoValSimValues_List[[nParameter]] <- genoValSimValues 
		phenoValSimValues_List[[nParameter]] <- phenoValSimValues 
		
 
 }

### Split Geno and Pheno values according to families  ########################
#### GM method for parent selection in cycle nCyc #####################################
 
 
 
 genoSelected_List <- foreach(k=1:nParameters) %dopar% (getGM_Frontier_FS_selectedGenoList(nextGenGenoTable_GM_List[[k]],nextGenGenoTableMod_GM_List[[k]],PredictionModel_Geno_List[[k]],PredictionModel_Pheno_List[[k]],modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List))
  
## Reduce Parameter set to 30  
 
    for(nParameter in 1:nParameters){
   
      selectedParents_Geno_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[1]]
      CriterionValue_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nParameter]][[1]]
      nParameters <- genoSelected_List[[nParameter]][[3]]
    
    }  
 
    count <- 1
	GainValues <- c()
	UsefulnessValues <- c()
	InbreedValues <- c()
	selectedParents_Geno_Par_List <- list()
 
    for (nPar1 in 1: nParameters){ 
 
      for(nPar2 in 1:nParameters){ 
	    GainValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][1]
	    UsefulnessValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][2]
	    InbreedValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][3]
	    selectedParents_Geno_Par_List[[count]]  <-  selectedParents_Geno_List[[nCyc]][[nPar1]][[nPar2]]
	 
	    count <- count+1
	  }
    }
 
      
    Quartile3 <- summary(GainValues)[5]
	Quartile1 <- summary(GainValues)[2]
			
	indices1 <- which(GainValues >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1 <- indicesSortRange1[[2]][1:10]
			
	indices2 <- which((GainValues <= Quartile3) & (GainValues >= Quartile1))
	indicesRange2 <- sample(indices2,10) 
			
    indices3 <- which(GainValues <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3 <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1,indicesRange2,indicesRange3)
	nParameters <- length(indicesRange)
			
	for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	   
			   selectedParents_Geno_List[[nCyc]][[nParameter]] <-  selectedParents_Geno_Par_List[[nParameterIndex]]
			   CriterionValue_List[[nCyc]][[nParameter]] <- c(GainValues[nParameterIndex],UsefulnessValues[nParameterIndex],InbreedValues[nParameterIndex])
		
	} 
   
    save.image(WorkspaceName)
    gc() 	
  
#####################################################################################################################
    for(nParameter in 1:nParameters){
	
	    genoValues <- genoValues_List[[nParameter]] 
        phenoValues <- phenoValues_List[[nParameter]] 
		genoValSimValues <- genoValSimValues_List[[nParameter]] 
		phenoValSimValues <- phenoValSimValues_List[[nParameter]] 
		
	       
    	uniqueParentIDs <-  unique(as.vector(selectedParents_Geno_List[[nCyc]][[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs

   		
### GM Output in selected individual ids and Num Progeny from selected parents
		   		
	
	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-   genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(genoValSimValues[GenoTableIndices_topN])
		
		topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
		topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- PhenoTableIndices_topN
		PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
#########################
	
	Cycle_Progeny_F5 <- Cycle_Progeny_F5_List[[nParameter]]
	
	nCrosses <- length(Cycle_Progeny_F5)
	nProgeny_per_Cross <- 10 
	Cycle_GenoData  <- (array(0,c(1,nMarkers,2,nCrosses*nProgeny_per_Cross)))
	
	nCrsInit <- 1
	nCrsFinal <- nProgeny_per_Cross
   
    for(nCrs in 1:nCrosses){

	 Cycle_GenoData[1,,,nCrsInit:nCrsFinal] <- Cycle_Progeny_F5[[nCrs]]
	 
	 nCrsInit <- nCrsFinal+1
	 nCrsFinal <-  nCrsFinal+ nProgeny_per_Cross 
	}
	
     	
    F5RILs_GM <- getF5RILs_FS_GM_Frontier(Cycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New)
	
	Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
    nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
	
	}
	

		save.image(WorkspaceName) 
		gc()
	 
}

######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	   return(simResults_List)

}

	
## try getF5RILs_FS_GM with equal number of progeny and select 
	
getF5RILs_FS_GM_Frontier <- function(GenoData_List,SelectedParents_Geno_List,SelectionOnGeno,nMarkers,NAM_LinkMap,nProgeny_per_Cross){

	    
        BD <- "GM"
		NAM_LinkMap_New <- NAM_LinkMap
        GenoData_List <- GenoData_List
        selectedParents_Geno_List <- SelectedParents_Geno_List
		selectionOnGeno <- SelectionOnGeno

##### Generate F5 RIL Progeny from selected parent set for all families 
#####	Initialize lists for progeny data for all families every cylce 
				
		nParents <- dim(selectedParents_Geno_List)[1]
		nUniqueParents <-  length(unique(as.vector(selectedParents_Geno_List)))		
		nProgeny <- nProgeny_per_Cross
		nPairs <- 1
			
### Loop througn selected parent pairs to generate progeny for each family
			
		Cycle_Progeny_F5_temp1_List <- rep(list(array(0,c(1,nMarkers,2,nProgeny))),nParents)	

        while(nPairs <=(nParents)) {
															 
				Cycle_Progeny_F1 <- array(0,c(1,nMarkers,2,1)) 
				
				Cycle_Progeny_F2 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F3 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F4 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F5_temp1 <- array(0,c(1,nMarkers,2,nProgeny)) 

			
###################################################################################
		 
### Get Parent pair geno data
		        parent1Index <- selectedParents_Geno_List[nPairs,1]
				parent2Index <- selectedParents_Geno_List[nPairs,2]
				Parent1<-  GenoData_List[1,,,parent1Index]
				Parent2<-  GenoData_List[1,,,parent2Index]
			  
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
			  
			   Cycle_Progeny_F5_temp1_List[[nPairs]] <- Cycle_Progeny_F5_temp1
			
				nPairs <- nPairs+1
			 
			}
			
			   Cycle_Progeny_F5_List_Temp <- Cycle_Progeny_F5_temp1_List

	          
	   
	   
### Get Cycle_Progeny F5 data in correct format 
      
   	  
	            Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Temp
       	 
	 		

	    return(list(Cycle_Progeny_F5,nParents)) 
  
  }


	
getF5RILs_FS_GM_Frontier_Complete <- function(GenoData_List,SelectedParents_Geno_List,SelectionOnGeno,nMarkers,NAM_LinkMap,nProgeny_per_Cross){

	    
        BD <- "GM"
		NAM_LinkMap_New <- NAM_LinkMap
        GenoData_List <- GenoData_List
        selectedParents_Geno_List <- SelectedParents_Geno_List
		selectionOnGeno <- SelectionOnGeno

##### Generate F5 RIL Progeny from selected parent set for all families 
#####	Initialize lists for progeny data for all families every cylce 
				
		nParents <- dim(selectedParents_Geno_List)[1]
		nUniqueParents <-  length(unique(as.vector(selectedParents_Geno_List)))		
		nProgeny <- nProgeny_per_Cross
		nPairs <- 1
			
### Loop througn selected parent pairs to generate progeny for each family
			
		Cycle_Progeny_F5_temp1_List <- rep(list(array(0,c(1,nMarkers,2,nProgeny))),nParents)	

        while(nPairs <=(nParents)) {
															 
				Cycle_Progeny_F1 <- array(0,c(1,nMarkers,2,1)) 
				
				Cycle_Progeny_F2 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F3 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F4 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F5_temp1 <- array(0,c(1,nMarkers,2,nProgeny)) 

			
###################################################################################
		 
### Get Parent pair geno data
		        parent1Index <- selectedParents_Geno_List[nPairs,1]
				parent2Index <- selectedParents_Geno_List[nPairs,2]
				Parent1<-  GenoData_List[1,,,parent1Index]
				Parent2<-  GenoData_List[1,,,parent2Index]
			  
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
			  
			   Cycle_Progeny_F5_temp1_List[[nPairs]] <- Cycle_Progeny_F5_temp1
			
				nPairs <- nPairs+1
			 
			}
			
			   Cycle_Progeny_F5_List_Temp <- Cycle_Progeny_F5_temp1_List

	          
	   
	   
### Get Cycle_Progeny F5 data in correct format 
      
   	  
	            Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Temp
       	 
	 		

	    return(list(Cycle_Progeny_F5,nParents)) 
  
  }





##################################################################################################


 getGM_Frontier_FS_selectedGenoList <- function(nextGenGenoTable_GM_List,nextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List){

	nextGenGenoTable_GM <- nextGenGenoTable_GM_List 
	nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_List 
	
	selectedParents_Geno_List <- list()
	selectedCriterion_List <- list()
	
    GM_Param_List <- GM_Frontier_Param_List
	
### GM_Parameters
	
	Mc.cores= GM_Param_List$no_cores
	No_cores = GM_Param_List$no_cores
	Mutprob= GM_Param_List$mutprob
	NpopGA= GM_Param_List$npopGA
	NitGA= GM_Param_List$nitGA
	Plotiters= GM_Param_List$plotiters
	Nelite= GM_Param_List$nelite
	Method= GM_Param_List$method
	PlotMates= GM_Param_List$plotMates
    NMates= GM_Param_List$nMates
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
			
			baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutionsFrontier(Markers,Markers2=NULL,K,markerEffects_Geno,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,plotiters=Plotiters,method=Method,nmates=NMates)
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Geno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
			
			
			Quartile3 <- summary(baseGASolns_GenoModel[[1]][,1])[5]
			Quartile1 <- summary(baseGASolns_GenoModel[[1]][,1])[2]
			
			
			indices1 <- which(baseGASolns_GenoModel[[1]][,1]>= Quartile3)
            indicesSortRange1 <- sort.int(baseGASolns_GenoModel[[1]][indices1,1],decreasing=TRUE,index.return=TRUE)
			indicesRange1 <- indicesSortRange1[[2]][1:10]
			
			indices2 <- which((baseGASolns_GenoModel[[1]][,1]<= Quartile3) & (baseGASolns_GenoModel[[1]][,1]>= Quartile1))
			indicesRange2 <- sample(indices2,10) 
			
            indices3 <- which(baseGASolns_GenoModel[[1]][,1]<= Quartile1)
            indicesSortRange3 <- sort.int(baseGASolns_GenoModel[[1]][indices3,1],decreasing=FALSE,index.return=TRUE) 
			indicesRange3 <- indicesSortRange3[[2]][1:10]
			
			indicesRange <- c(indicesRange1,indicesRange2,indicesRange3)
			nParameters <- length(indicesRange)
			for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	        	       
			   selectedParents_Geno_List[[nParameter]] <-  baseGASolns_GenoModel[[2]][[nParameterIndex]]
			   selectedCriterion_List[[nParameter]] <- baseGASolns_GenoModel[[1]][nParameterIndex,]
		
			} 
			
			
		}else if(selectionOnGeno==FALSE){
		
## selection on phenotypic values 
			baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutionsFrontier(Markers,Markers2=NULL,K,markerEffects_Pheno,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,plotiters=Plotiters,method=Method,nmates=NMates) 
			
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Pheno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
			
			
			Quartile3 <- summary(baseGASolns_GenoModel[[1]][,1])[5]
			Quartile1 <- summary(baseGASolns_GenoModel[[1]][,1])[2]
			
			
			indices1 <- which(baseGASolns_GenoModel[[1]][,1]>= Quartile3)
            indicesSortRange1 <- sort.int(baseGASolns_GenoModel[[1]][indices1,1],decreasing=TRUE,index.return=TRUE)
			indicesRange1 <- indicesSortRange1[[2]][1:10]
			
			indices2 <- which((baseGASolns_GenoModel[[1]][,1]<= Quartile3) & (baseGASolns_GenoModel[[1]][,1]>= Quartile1))
			indicesRange2 <- sample(indices2,10) 
			
            indices3 <- which(baseGASolns_GenoModel[[1]][,1]<= Quartile1)
            indicesSortRange3 <- sort.int(baseGASolns_GenoModel[[1]][indices3,1],decreasing=FALSE,index.return=TRUE) 
			indicesRange3 <- indicesSortRange3[[2]][1:10]
			
			indicesRange <- c(indicesRange1,indicesRange2,indicesRange3)
			nParameters <- length(indicesRange)
			for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	        	       
			   selectedParents_Geno_List[[nParameter]] <-  baseGASolns_GenoModel[[2]][[nParameterIndex]]
			   selectedCriterion_List[[nParameter]] <- baseGASolns_GenoModel[[1]][nParameterIndex,]
		
			} 
			
		}
		
	
	return(list(selectedParents_Geno_List,selectedCriterion_List,nParameters))
	
	
} 






 getGM_Frontier_FS_selectedGenoList_Cyc1 <- function(nextGenGenoTable_GM_List,nextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List){

	nextGenGenoTable_GM <- nextGenGenoTable_GM_List 
	nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_List 
	
	selectedParents_Geno_List <- list()
	selectedCriterion_List <- list()
	
    GM_Param_List <- GM_Frontier_Param_List
	
### GM_Parameters
	
	Mc.cores= GM_Param_List$no_cores
	No_cores = GM_Param_List$no_cores
	Mutprob= GM_Param_List$mutprob
	NpopGA= GM_Param_List$npopGA
	NitGA= GM_Param_List$nitGA
	Plotiters= GM_Param_List$plotiters
	Nelite= GM_Param_List$nelite
	Method= GM_Param_List$method
	PlotMates= GM_Param_List$plotMates
    NMates= GM_Param_List$nMates
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
			
			baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutionsFrontier(Markers,Markers2=NULL,K,markerEffects_Geno,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,plotiters=Plotiters,method=Method,nmates=NMates)
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Geno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
			
			
			Quartile3 <- summary(baseGASolns_GenoModel[[1]][,1])[5]
			Quartile1 <- summary(baseGASolns_GenoModel[[1]][,1])[2]
			
			
			indices1 <- which(baseGASolns_GenoModel[[1]][,1]>= Quartile3)
            indicesSortRange1 <- sort.int(baseGASolns_GenoModel[[1]][indices1,1],decreasing=TRUE,index.return=TRUE)
			indicesRange1 <- indicesSortRange1[[2]][1:10]
			
			indices2 <- which((baseGASolns_GenoModel[[1]][,1]<= Quartile3) & (baseGASolns_GenoModel[[1]][,1]>= Quartile1))
			indicesRange2 <- sample(indices2,10) 
			
            indices3 <- which(baseGASolns_GenoModel[[1]][,1]<= Quartile1)
            indicesSortRange3 <- sort.int(baseGASolns_GenoModel[[1]][indices3,1],decreasing=FALSE,index.return=TRUE) 
			indicesRange3 <- indicesSortRange3[[2]][1:10]
			
			indicesRange <- c(indicesRange1,indicesRange2,indicesRange3)
			nParameters <- length(indicesRange)
			for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	        	       
			   selectedParents_Geno_List[[nParameter]] <-  baseGASolns_GenoModel[[2]][[nParameterIndex]]
			   selectedCriterion_List[[nParameter]] <- baseGASolns_GenoModel[[1]][nParameterIndex,]
		
			} 
			
			
		}else if(selectionOnGeno==FALSE){
		
## selection on phenotypic values 
			baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutionsFrontier(Markers,Markers2=NULL,K,markerEffects_Pheno,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,plotiters=Plotiters,method=Method,nmates=NMates) 
			
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Pheno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
			
			
			Quartile3 <- summary(baseGASolns_GenoModel[[1]][,1])[5]
			Quartile1 <- summary(baseGASolns_GenoModel[[1]][,1])[2]
			
			
			indices1 <- which(baseGASolns_GenoModel[[1]][,1]>= Quartile3)
            indicesSortRange1 <- sort.int(baseGASolns_GenoModel[[1]][indices1,1],decreasing=TRUE,index.return=TRUE)
			indicesRange1 <- indicesSortRange1[[2]][1:10]
			
			indices2 <- which((baseGASolns_GenoModel[[1]][,1]<= Quartile3) & (baseGASolns_GenoModel[[1]][,1]>= Quartile1))
			indicesRange2 <- sample(indices2,10) 
			
            indices3 <- which(baseGASolns_GenoModel[[1]][,1]<= Quartile1)
            indicesSortRange3 <- sort.int(baseGASolns_GenoModel[[1]][indices3,1],decreasing=FALSE,index.return=TRUE) 
			indicesRange3 <- indicesSortRange3[[2]][1:10]
			
			indicesRange <- c(indicesRange1,indicesRange2,indicesRange3)
			nParameters <- length(indicesRange)
			for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	        	       
			   selectedParents_Geno_List[[nParameter]] <-  baseGASolns_GenoModel[[2]][[nParameterIndex]]
			   selectedCriterion_List[[nParameter]] <- baseGASolns_GenoModel[[1]][nParameterIndex,]
		
			} 
			
		}
		
	
	return(list(selectedParents_Geno_List,selectedCriterion_List,nParameters,baseGASolns_GenoModel))
	
	
} 

########### 


########### 

 getGM_Frontier_FS_selectedGenoList_Cyc_2_N<- function(nextGenGenoTable_GM_List,nextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List){

	nextGenGenoTable_GM <- nextGenGenoTable_GM_List 
	nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_List 
	
	selectedParents_Geno_List <- list()
	selectedCriterion_List <- list()
	
    GM_Param_List <- GM_Frontier_Param_List
	
### GM_Parameters
	
	Mc.cores= GM_Param_List$no_cores
	No_cores = GM_Param_List$no_cores
	Mutprob= GM_Param_List$mutprob
	NpopGA= GM_Param_List$npopGA
	NitGA= GM_Param_List$nitGA
	Plotiters= GM_Param_List$plotiters
	Nelite= GM_Param_List$nelite
	Method= GM_Param_List$method
	PlotMates= GM_Param_List$plotMates
    NMates= GM_Param_List$nMates
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
			
			baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutionsFrontier(Markers,Markers2=NULL,K,markerEffects_Geno,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,plotiters=Plotiters,method=Method,nmates=NMates)
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Geno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
				
			selectedParents_Geno_List <-  baseGASolns_GenoModel[[2]]
			selectedCriterion_List <- baseGASolns_GenoModel[[1]]
			nParameters <- length(baseGASolns_GenoModel[[2]])
				
		}else if(selectionOnGeno==FALSE){
		
## selection on phenotypic values 
			baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutionsFrontier(Markers,Markers2=NULL,K,markerEffects_Pheno,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,plotiters=Plotiters,method=Method,nmates=NMates) 
					
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Pheno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
			
		
			selectedParents_Geno_List <-  baseGASolns_GenoModel[[2]]
			selectedCriterion_List <- baseGASolns_GenoModel[[1]] 
			nParameters <- length(baseGASolns_GenoModel[[2]])
		
			
		}
		
	
	return(list(selectedParents_Geno_List,selectedCriterion_List,nParameters,baseGASolns_GenoModel))
	
	
} 




### Generate Minus101 format genotype table for frontier progeny set

generateMinus101GenoFormat_GM_Frontier_V1<- function(Cycle_Progeny_Data,number_of_Crosses,nProgeny_per_Cross){

  Cycle_Progeny <- Cycle_Progeny_Data
  nCrosses <- number_of_Crosses
  
  nMarkers <- 4289
  nProgeny <- nProgeny_per_Cross 
 
  Cycle_Progeny_table <-c()
  
  for(i in 1:nCrosses){
   cycle_family_table<- array(0,c(1,nMarkers,nProgeny))
  
    for(k in 1:nProgeny){
   
	 cycle_family_table[1,,k] <- (Cycle_Progeny[[i]][1,,1,k] + Cycle_Progeny[[i]][1,,2,k])

	}
   
    Cycle_Progeny_table <-cbind(Cycle_Progeny_table,cycle_family_table[1,,])

  }
 IndividualNames <-rep(0,(nCrosses*nProgeny))

  for(i in 1:(nCrosses*nProgeny)){
    IndividualNames[i]<- paste("Ind",i,sep="")
  }

  # IndividualNames <- c("Animal_ID",IndividualNames)
  #### RowNames

  markerNames <-rep(0,nMarkers)


  for(i in 1:(nMarkers)){
    markerNames[i]<- paste("m",i,sep="")
  }
  
  colnames(Cycle_Progeny_table)<- IndividualNames
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



################################################################################################## 

generateMinus101GenoFormat_GM_Frontier_V1_Complete<- function(Cycle_Progeny_Data,number_of_Crosses){

  Cycle_Progeny <- Cycle_Progeny_Data
  nCrosses <- number_of_Crosses
  
  nMarkers <- 4289
  nProgeny <- 1 
 
  Cycle_Progeny_table <-c()
  k <- 1
  
  for(i in 1:nCrosses){
    cycle_family_table<- array(0,c(1,nMarkers,nProgeny))
  
    #for(k in 1:nProgeny){
	   
	 cycle_family_table[1,,k] <- (Cycle_Progeny[[i]][1,,1,k] + Cycle_Progeny[[i]][1,,2,k])

	#}
   
    Cycle_Progeny_table <-cbind(Cycle_Progeny_table,cycle_family_table[1,,])

  }
 IndividualNames <-rep(0,(nCrosses*nProgeny))

  for(i in 1:(nCrosses*nProgeny)){
    IndividualNames[i]<- paste("Ind",i,sep="")
  }

  # IndividualNames <- c("Animal_ID",IndividualNames)
  #### RowNames

  markerNames <-rep(0,nMarkers)


  for(i in 1:(nMarkers)){
    markerNames[i]<- paste("m",i,sep="")
  }
  
  colnames(Cycle_Progeny_table)<- IndividualNames
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


####

runSimulations20X_BD_Wght_GM_Frontier_V2 <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,WorkSpaceName,StartCycle,GM_Frontier_Param_List,nCore_ForEach){

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
	  
	  GM_Param_List <- GM_Frontier_Param_List
	  no_cores <- GM_Param_List$no_cores
	  startCycle <- StartCycle
	  WorkspaceName <- WorkSpaceName
	  
	  nCores <- nCore_ForEach
	 
#################################################################################################
#################################################################################################
### Initialize variables	
	
   	  j<-1
 	 
################################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
	  combinedTableCount <- 1
	  
################################################################################################	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
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
	 

# #################################################################################
    nCrosses <- dim(cycle1_Geno_Data)[1]
    nProgeny <- dim(cycle1_Geno_Data)[4]
    no_selected <- 200

	
    Cycle_GenoData  <- (array(0,c(1,nMarkers,2,cycle1_nProgeny*nCrosses)))
    selectedCycle_GenoData  <- (array(0,c(1,nMarkers,2,no_selected))) 
   
    initIndex <- 1
    finalIndex <- nProgeny
   
    for(nCrs in 1:nCrosses){

	
	 Cycle_GenoData[1,,,initIndex:finalIndex] <- cycle1_Geno_Data[nCrs,,,] 
	 initIndex <- initIndex + nProgeny
	 finalIndex <- finalIndex + nProgeny
     }
	
	
	
############  
### ApplySelection to select top nSel lines 
	
	if(selectionOnSimulated ==TRUE){

		genoSelectionTable <- ApplySelection(genoValSimValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)


		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedGenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedGenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]

		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedPhenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}

		
	}
    if(selectionOnSimulated==FALSE){
		genoSelectionTable <- ApplySelection(genoValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

    	print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
		print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
		print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedGenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedGenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]
		
		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedPhenoIndividualIndices,]
		
		    selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}
		
	} 
	
## 

  

#### GM method for parent selection in cycle 1 #####################################
 nParameters <- 30 
 selectedParents_Geno_List <- rep(list(rep(list(list()),nParameters)),nCycles)
 CriterionValue_List <- rep(list(rep(list(list()),nParameters)),nCycles) 

 system.time({
 
 genoSelected_List <- getGM_Frontier_FS_selectedGenoList(selectedCycle1GenoTable_GM,selectedCycle1GenoTableMod_GM,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List)

 selectedParents_Geno_List[[nCyc]] <- genoSelected_List[[1]]
 CriterionValue_List[[nCyc]] <- genoSelected_List[[2]]
 selectedParents_Geno <- genoSelected_List[[1]]
 nParameters <- genoSelected_List[[3]]
 
})

# save.image(WorkspaceName) 	
 gc()
 
################################################################################### 
## Initiate Output list 	  
	  
	  GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	 
	  PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  
      attainedPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      attainedGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	
	  topNGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      GenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)

      topNPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      PhenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  Cycle_Progeny_F5_List <- list()
      nCrosses_List <- list()
	  
	  
	 
##
 
   for(nParameter in 1:nParameters){
   
        uniqueParentIDs <-  unique(as.vector(selectedParents_Geno[[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
### GM Output in selected individual ids 

############################

	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]]<- genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <- max(genoValSimValues[GenoTableIndices_topN])
		
	    topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN_List
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
	   
	    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
	    topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <-  PhenoTableIndices_topN
	    PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		
	    
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
###################################################################################	
## Get F5 RILs for selected parents in cycle 1
  
      
        F5RILs_GM <- getF5RILs_FS_GM_Frontier(selectedCycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New)

        Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
    
	}

	simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)


    startCycle <- startCycle + 1 

    save(simResults_List_Cycle,file=WorkspaceName)
		
	#save.image(WorkspaceName) 
	gc() 
	
###########################################################################

    startCycle <- 2
    #nCores <- 10
    cl <- makeCluster(nCores,outfile='log.txt')
    registerDoParallel(cl)

    nextGenGenoTable_GM_List <- list()
    nextGenGenoTableMod_GM_List <- list()
    newNextGenGenoTable_GM_List <- list()
	PredictionModel_Geno_List <- list()
	PredictionModel_Pheno_List <- list()
	genoValues_List <- list()
    phenoValues_List <- list()
	genoValSimValues_List <- list()
	phenoValSimValues_List <- list()
 
    selectedNextGenGenoTable_GM_List <- list() 
    selectedNextGenGenoTableMod_GM_List <- list() 
#	selectedNewNextGenGenoTable_GM_List <- list() 

for(nCyc in startCycle:nCycles){
        #nCyc <- 2

    print(nCyc)

### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
   for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <-  Cycle_Progeny_F5_List[[nParameter]]
		nCrosses <-  nCrosses_List[[nParameter]] 
		
		nextGenGenoTable_GM <- generateMinus101GenoFormat_GM_Frontier_V1(Cycle_Progeny_F5,nCrosses)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
	
############################################################################################### 

    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0 && nParameter ==1){

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

                        PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

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

        nextGenGenoTable_GM_List[[nParameter]] <- nextGenGenoTable_GM
        nextGenGenoTableMod_GM_List[[nParameter]] <- nextGenGenoTableMod_GM
		newNextGenGenoTable_GM_List[[nParameter]] <- newNextGenGenoTable_GM 
		PredictionModel_Geno_List[[nParameter]] <- PredictionModel_Geno
		PredictionModel_Pheno_List[[nParameter]] <- PredictionModel_Pheno
		genoValues_List[[nParameter]] <- genoValues
        phenoValues_List[[nParameter]] <- phenoValues
		genoValSimValues_List[[nParameter]] <- genoValSimValues 
		phenoValSimValues_List[[nParameter]] <- phenoValSimValues 

		
### ApplySelection to nextGenGenoTable_GM to select top 200 lines 

    no_selected <- GM_Param_List$nMates
	
	if(selectionOnSimulated ==TRUE){

		genoSelectionTable <- ApplySelection(genoValSimValues ,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)


		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedGenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedGenoIndividualIndices,]

		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		}
	}
    if(selectionOnSimulated==FALSE){

		genoSelectionTable <- ApplySelection(genoValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

    	print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
		print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
		print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		
		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		}
		
	}
	 
	    selectedNextGenGenoTable_GM_List[[nParameter]] <- selectedNextGenGenoTable_GM
        selectedNextGenGenoTableMod_GM_List[[nParameter]] <- selectedNextGenGenoTableMod_GM
 
 }

### Split Geno and Pheno values according to families  ########################
#### GM method for parent selection in cycle nCyc #####################################
  
 
 genoSelected_List <- foreach(k=1:nParameters) %dopar% (getGM_Frontier_FS_selectedGenoList(selectedNextGenGenoTable_GM_List[[k]],selectedNextGenGenoTableMod_GM_List[[k]],PredictionModel_Geno_List[[k]],PredictionModel_Pheno_List[[k]],modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List))
  
## Reduce Parameter set to 30  
 
    for(nParameter in 1:nParameters){
   
      selectedParents_Geno_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[1]]
      CriterionValue_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nParameter]][[1]]
      nParameters <- genoSelected_List[[nParameter]][[3]]
    
    }  
 
    count <- 1
	GainValues <- c()
	UsefulnessValues <- c()
	InbreedValues <- c()
	selectedParents_Geno_Par_List <- list()
 
    for (nPar1 in 1: nParameters){ 
 
      for(nPar2 in 1:nParameters){ 
	    GainValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][1]
	    UsefulnessValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][2]
	    InbreedValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][3]
	    selectedParents_Geno_Par_List[[count]]  <-  selectedParents_Geno_List[[nCyc]][[nPar1]][[nPar2]]
	 
	    count <- count+1
	  }
    }
 
      
    Quartile3 <- summary(GainValues)[5]
	Quartile1 <- summary(GainValues)[2]
			
	indices1 <- which(GainValues >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1 <- indicesSortRange1[[2]][1:10]
			
	indices2 <- which((GainValues <= Quartile3) & (GainValues >= Quartile1))
	indicesRange2 <- sample(indices2,10) 
			
    indices3 <- which(GainValues <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3 <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1,indicesRange2,indicesRange3)
	nParameters <- length(indicesRange)
			
	for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	   
			   selectedParents_Geno_List[[nCyc]][[nParameter]] <-  selectedParents_Geno_Par_List[[nParameterIndex]]
			   CriterionValue_List[[nCyc]][[nParameter]] <- c(GainValues[nParameterIndex],UsefulnessValues[nParameterIndex],InbreedValues[nParameterIndex])
		
	} 
   
    gc() 	
  
#####################################################################################################################
    for(nParameter in 1:nParameters){
	
	    genoValues <- genoValues_List[[nParameter]] 
        phenoValues <- phenoValues_List[[nParameter]] 
        genoValSimValues <- genoValSimValues_List[[nParameter]] 
	    phenoValSimValues <- phenoValSimValues_List[[nParameter]] 
	       
    	uniqueParentIDs <-  unique(as.vector(selectedParents_Geno_List[[nCyc]][[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs

### GM Output in selected individual ids and Num Progeny from selected parents
		   		
	
	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-   genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(genoValSimValues[GenoTableIndices_topN])
		
		topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
		topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- PhenoTableIndices_topN
		PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
#########################
	
	Cycle_Progeny_F5 <- Cycle_Progeny_F5_List[[nParameter]]
	
	nCrosses <- length(Cycle_Progeny_F5)
	nProgeny_per_Cross <- 10
        no_selected <- 200 
	Cycle_GenoData  <- (array(0,c(1,nMarkers,2,nCrosses*nProgeny_per_Cross)))
        selectedCycle_GenoData  <- (array(0,c(1,nMarkers,2,no_selected)))	
	nCrsInit <- 1
	nCrsFinal <- nProgeny_per_Cross
   
    for(nCrs in 1:nCrosses){

	  Cycle_GenoData[1,,,nCrsInit:nCrsFinal] <- Cycle_Progeny_F5[[nCrs]]
	 
	  nCrsInit <- nCrsFinal+1
	  nCrsFinal <-  nCrsFinal+ nProgeny_per_Cross 
	
	}
	
	if(selectionOnGeno == TRUE){
	   selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]
	}else if(selectionOnGeno == FALSE){ 
	   selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
	}
     	
        F5RILs_GM <- getF5RILs_FS_GM_Frontier(selectedCycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New)
	
	    Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
	
	}
	
	simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

  
    save(simResults_List_Cycle,file=WorkspaceName)
		
	
	gc()
	 
}

######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	   return(simResults_List)

}


# testSelection <- ApplySelection(genoValues,20)

 

runSimulations20X_BD_Wght_GM_Frontier_V3 <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach){

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
	  
	  GM_Param_List <- GM_Frontier_Param_List
	  no_cores <- GM_Param_List$no_cores
	  startCycle <- StartCycle
	  WorkspaceName <- WorkSpaceName
	  nCores_ForEach <- noCores_ForEach 
	  
	  save.image(WorkspaceName)
#################################################################################################
#################################################################################################
### Initialize variables	
	
   	  j<-1
	  CombinedTableCount <- 1
################################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
	
	  
################################################################################################	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
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
	 

# #################################################################################
    nCrosses <- dim(cycle1_Geno_Data)[1]
    nProgeny <- dim(cycle1_Geno_Data)[4]
    no_selected <- 200

	
    Cycle_GenoData  <- (array(0,c(1,nMarkers,2,cycle1_nProgeny*nCrosses)))
    selectedCycle_GenoData  <- (array(0,c(1,nMarkers,2,no_selected))) 
   
    initIndex <- 1
    finalIndex <- nProgeny
   
    for(nCrs in 1:nCrosses){

	
	 Cycle_GenoData[1,,,initIndex:finalIndex] <- cycle1_Geno_Data[nCrs,,,] 
	 initIndex <- initIndex + nProgeny
	 finalIndex <- finalIndex + nProgeny
     }
	
	
	
############  
### ApplySelection to select top 200 lines 
	
	if(selectionOnSimulated ==TRUE){

		genoSelectionTable <- ApplySelection(genoValSimValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)


		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedGenoIndividualIndices,]	
                        selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedGenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]

		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedPhenoIndividualIndices,]	
                        selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedPhenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}

		
	}
    if(selectionOnSimulated==FALSE){
		genoSelectionTable <- ApplySelection(genoValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

    	print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
		print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
		print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
		        genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedGenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedGenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]
		
		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedPhenoIndividualIndices,]
		
		    selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}
		
	} 
	
## 

#### GM method for parent selection in cycle 1 #####################################
 nParameters <- 30 
 selectedParents_Geno_List <- rep(list(rep(list(list()),nParameters)),nCycles)
 CriterionValue_List <- rep(list(rep(list(list()),nParameters)),nCycles) 

 system.time({
 
 genoSelected_List <- getGM_Frontier_FS_selectedGenoList_Cyc1(selectedCycle1GenoTable_GM,selectedCycle1GenoTableMod_GM,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List)

 selectedParents_Geno_List[[nCyc]] <- genoSelected_List[[1]]
 CriterionValue_List[[nCyc]] <- genoSelected_List[[2]]
 selectedParents_Geno <- genoSelected_List[[1]]
 nParameters <- genoSelected_List[[3]]
 
})

save(genoSelected_List,file=WorkspaceName)

 #save.image(WorkspaceName) 	
 gc()
 
################################################################################### 
## Initiate Output list 	  
	  
	  GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	 
	  PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  
      attainedPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      attainedGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	
	  topNGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      GenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)

      topNPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      PhenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  Cycle_Progeny_F5_List <- list()
      nCrosses_List <- list()
	  
	  
	 
##
 
   for(nParameter in 1:nParameters){
   
        uniqueParentIDs <-  unique(as.vector(selectedParents_Geno[[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
### GM Output in selected individual ids 

############################

	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]]<- genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <- max(genoValSimValues[GenoTableIndices_topN])
		
	    topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN_List
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
	   
	    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
	    topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <-  PhenoTableIndices_topN
	    PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		
	    
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
###################################################################################	
## Get F5 RILs for selected parents in cycle 1
  
      
        F5RILs_GM <- getF5RILs_FS_GM_Frontier(selectedCycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New)

        Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
    
	}

	simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)
	
	save(simResults_List_Cycle,file=WorkspaceName)
	
    startCycle <- startCycle + 1 	 
	 
	gc() 
	
###########################################################################

    startCycle <- 2
    nCores <- nCores_ForEach 
    cl <- makeCluster(nCores,outfile='log.txt')
    registerDoParallel(cl)

    nextGenGenoTable_GM_List <- list()
    nextGenGenoTableMod_GM_List <- list()
    newNextGenGenoTable_GM_List <- list()
	PredictionModel_Geno_List <- list()
	PredictionModel_Pheno_List <- list()
	genoValues_List <- list()
    phenoValues_List <- list()
	genoValSimValues_List <- list()
	phenoValSimValues_List <- list()
 
    selectedNextGenGenoTable_GM_List <- list() 
    selectedNextGenGenoTableMod_GM_List <- list() 
#	selectedNewNextGenGenoTable_GM_List <- list() 

for(nCyc in startCycle:nCycles){
        #nCyc <- 2

    print(nCyc)

### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
   for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <-  Cycle_Progeny_F5_List[[nParameter]]
		nCrosses <-  nCrosses_List[[nParameter]] 
		
		nextGenGenoTable_GM <- generateMinus101GenoFormat_GM_Frontier_V1(Cycle_Progeny_F5,nCrosses)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
	
############################################################################################### 
		if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){

			PredictionModels <- getUpdatedGSModel(newNextGenGenoTable_GM,genoValSimValues,phenoValSimValues,PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,trainTableWindowSize,trainGenoNewTablePreCycle,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount,nCyc) 
 
			PredictionModel_Geno <- PredictionModels[[1]]
			PredictionModel_Pheno <- PredictionModels[[2]]
			CombinedTableCount <- PredictionModels[[3]]
			trainGenoNewTablePreCycle <- PredictionModels[[4]]
			trainSimPhenoTablePreCycle <- PredictionModels[[5]]
			trainSimGenoValTablePreCycle <- PredictionModels[[6]]
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

        nextGenGenoTable_GM_List[[nParameter]] <- nextGenGenoTable_GM
        nextGenGenoTableMod_GM_List[[nParameter]] <- nextGenGenoTableMod_GM
		newNextGenGenoTable_GM_List[[nParameter]] <- newNextGenGenoTable_GM 
		PredictionModel_Geno_List[[nParameter]] <- PredictionModel_Geno
		PredictionModel_Pheno_List[[nParameter]] <- PredictionModel_Pheno
		genoValues_List[[nParameter]] <- genoValues
        phenoValues_List[[nParameter]] <- phenoValues
		genoValSimValues_List[[nParameter]] <- genoValSimValues 
		phenoValSimValues_List[[nParameter]] <- phenoValSimValues 

		
### ApplySelection to nextGenGenoTable_GM to select top 200 lines 

    no_selected <- GM_Param_List$nMates
	
	if(selectionOnSimulated ==TRUE){

		genoSelectionTable <- ApplySelection(genoValSimValues ,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)


		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedGenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedGenoIndividualIndices,]

		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		}
	}
    if(selectionOnSimulated==FALSE){

		genoSelectionTable <- ApplySelection(genoValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

    	print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
		print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
		print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		
		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		}
		
	}
	 
	    selectedNextGenGenoTable_GM_List[[nParameter]] <- selectedNextGenGenoTable_GM
        selectedNextGenGenoTableMod_GM_List[[nParameter]] <- selectedNextGenGenoTableMod_GM
 
 }

### Split Geno and Pheno values according to families  ########################
#### GM method for parent selection in cycle nCyc #####################################
 
       
	genoSelected_List <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (getGM_Frontier_FS_selectedGenoList_Cyc_2_N(selectedNextGenGenoTable_GM_List[[k]],selectedNextGenGenoTableMod_GM_List[[k]],PredictionModel_Geno_List[[k]],PredictionModel_Pheno_List[[k]],modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List))

    save(genoSelected_List,file=WorkspaceName)
  
## Reduce Parameter set to 30  
   
    for(nParameter in 1:nParameters){
    
      selectedParents_Geno_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[1]]
      CriterionValue_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nParameter]][[1]]
          
    }  
 
    count <- 1
	GainValues <- c()
	UsefulnessValues <- c()
	InbreedValues <- c()
	selectedParents_Geno_Par_List <- list()
	
    for (nPar1 in 1: nParameters){ 
	
	   nParameters2 <- length(selectedParents_Geno_List[[nCyc]][[nPar1]])
 
      for(nPar2 in 1:nParameters2){ 
	    GainValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][1]
	    UsefulnessValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][2]
	    InbreedValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][3]
	    selectedParents_Geno_Par_List[[count]]  <-  selectedParents_Geno_List[[nCyc]][[nPar1]][[nPar2]]
	 
	    count <- count+1
	  }
    }
 
	
	### 
	
    Quartile3 <- summary(GainValues)[5]
	Quartile1 <- summary(GainValues)[2]
	MedianG <- 	summary(GainValues)[3]	
	
    IntG1 <- Quartile3 - MedianG
	IntG2 <- MedianG - Quartile1

    
	
	indices1 <- which(GainValues >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:5]
			
	indices2A <- which((GainValues <= MedianG) & (GainValues >= Quartile1))
	indices2B <- which((GainValues >= MedianG) & (GainValues <= Quartile3))
	indices2CG <- (which(GainValues <= MedianG+1 & GainValues >= MedianG-1))[1]
	indices2DG <- which((GainValues <= MedianG+(0.1*IntG1)) & (GainValues >= MedianG- (0.1*IntG2)))
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2DG,5),sample(indices2B,5))

	indices3 <- which(GainValues <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
###########
				
	for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	   
			   selectedParents_Geno_List[[nCyc]][[nParameter]] <-  selectedParents_Geno_Par_List[[nParameterIndex]]
			   CriterionValue_List[[nCyc]][[nParameter]] <- c(GainValues[nParameterIndex],UsefulnessValues[nParameterIndex],InbreedValues[nParameterIndex])
	} 
   
  #  save.image(WorkspaceName)
    gc() 	
  
#####################################################################################################################
    for(nParameter in 1:nParameters){
	
	    genoValues <- genoValues_List[[nParameter]] 
        phenoValues <- phenoValues_List[[nParameter]] 
        genoValSimValues <- genoValSimValues_List[[nParameter]] 
	    phenoValSimValues <- phenoValSimValues_List[[nParameter]] 
	       
    	uniqueParentIDs <-  unique(as.vector(selectedParents_Geno_List[[nCyc]][[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs

### GM Output in selected individual ids and Num Progeny from selected parents
		
	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-   genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(genoValSimValues[GenoTableIndices_topN])
		
		topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
		topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- PhenoTableIndices_topN
		PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 

#########################
	
	Cycle_Progeny_F5 <- Cycle_Progeny_F5_List[[nParameter]]
	
	nCrosses <- length(Cycle_Progeny_F5)
	nProgeny_per_Cross <- 10
    no_selected <- 200 
	Cycle_GenoData  <- (array(0,c(1,nMarkers,2,nCrosses*nProgeny_per_Cross)))
    selectedCycle_GenoData  <- (array(0,c(1,nMarkers,2,no_selected)))	
	nCrsInit <- 1
	nCrsFinal <- nProgeny_per_Cross
   
    for(nCrs in 1:nCrosses){

	  Cycle_GenoData[1,,,nCrsInit:nCrsFinal] <- Cycle_Progeny_F5[[nCrs]]
	 
	  nCrsInit <- nCrsFinal+1
	  nCrsFinal <-  nCrsFinal+ nProgeny_per_Cross 
	
	}
	
	if(selectionOnGeno == TRUE){
	   selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]
	}else if(selectionOnGeno == FALSE){ 
	   selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
	}
     	
        F5RILs_GM <- getF5RILs_FS_GM_Frontier(selectedCycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New)
	
	    Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
	
	}
	

	#	save.image(WorkspaceName) 
		gc() 
		
    simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)
	
	save(simResults_List_Cycle,file=WorkSpaceName)
	 
}

######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	   return(simResults_List)

}

##### 




runSimulations20X_BD_Wght_GM_Frontier_V4 <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach){

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
	  
	  GM_Param_List <- GM_Frontier_Param_List
	  no_cores <- GM_Param_List$no_cores
	  startCycle <- StartCycle
	  WorkspaceName <- WorkSpaceName
	  nCores_ForEach <- noCores_ForEach 
	  
	  save.image(WorkspaceName)
#################################################################################################
#################################################################################################
### Initialize variables	
	
   	  j<-1
 	  CombinedTableCount <- 1
################################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
	  
	  
################################################################################################	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
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
	 

# #################################################################################
    nCrosses <- dim(cycle1_Geno_Data)[1]
    nProgeny <- dim(cycle1_Geno_Data)[4]
    no_selected <- 200

	
    Cycle_GenoData  <- (array(0,c(1,nMarkers,2,cycle1_nProgeny*nCrosses)))
    selectedCycle_GenoData  <- (array(0,c(1,nMarkers,2,no_selected))) 
   
    initIndex <- 1
    finalIndex <- nProgeny
   
    for(nCrs in 1:nCrosses){

	
	 Cycle_GenoData[1,,,initIndex:finalIndex] <- cycle1_Geno_Data[nCrs,,,] 
	 initIndex <- initIndex + nProgeny
	 finalIndex <- finalIndex + nProgeny
     }
	
	
	
############  
### ApplySelection to select top 200 lines 
	
	if(selectionOnSimulated ==TRUE){

		genoSelectionTable <- ApplySelection(genoValSimValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)


		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedGenoIndividualIndices,]	
                        selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedGenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]

		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedPhenoIndividualIndices,]	
                        selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedPhenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}

		
	}
    if(selectionOnSimulated==FALSE){
		genoSelectionTable <- ApplySelection(genoValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

    	print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
		print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
		print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
		        genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedGenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedGenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]
		
		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedPhenoIndividualIndices,]
		
		    selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}
		
	} 
	
## 

  

#### GM method for parent selection in cycle 1 #####################################
 nParameters <- 30 
 selectedParents_Geno_List <- rep(list(rep(list(list()),nParameters)),nCycles)
 CriterionValue_List <- rep(list(rep(list(list()),nParameters)),nCycles) 

 system.time({
 
 genoSelected_List <- getGM_Frontier_FS_selectedGenoList_Cyc1(selectedCycle1GenoTable_GM,selectedCycle1GenoTableMod_GM,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List)

 selectedParents_Geno_List[[nCyc]] <- genoSelected_List[[1]]
 CriterionValue_List[[nCyc]] <- genoSelected_List[[2]]
 selectedParents_Geno <- genoSelected_List[[1]]
 nParameters <- genoSelected_List[[3]]
 
})

save(genoSelected_List,file=WorkspaceName)

 #save.image(WorkspaceName) 	
 gc()
 
################################################################################### 
## Initiate Output list 	  
	  
	  GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	 
	  PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  
      attainedPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      attainedGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	
	  topNGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      GenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)

      topNPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      PhenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  Cycle_Progeny_F5_List <- list()
      nCrosses_List <- list()
	  
	  
	 
##
 
   for(nParameter in 1:nParameters){
   
        uniqueParentIDs <-  unique(as.vector(selectedParents_Geno[[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
### GM Output in selected individual ids 

############################

	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]]<- genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <- max(genoValSimValues[GenoTableIndices_topN])
		
	    topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN_List
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
	   
	    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
	    topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <-  PhenoTableIndices_topN
	    PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		
	    
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
###################################################################################	
## Get F5 RILs for selected parents in cycle 1
  
      
        F5RILs_GM <- getF5RILs_FS_GM_Frontier(selectedCycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New)

        Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
    
	}

	simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)
	
	save(simResults_List_Cycle,file=WorkspaceName)
	
    startCycle <- startCycle + 1 	 
	 
	gc() 
	
###########################################################################

    startCycle <- 2
    nCores <- nCores_ForEach 
    cl <- makeCluster(nCores,outfile='log.txt')
    registerDoParallel(cl)

    nextGenGenoTable_GM_List <- list()
    nextGenGenoTableMod_GM_List <- list()
    newNextGenGenoTable_GM_List <- list()
	PredictionModel_Geno_List <- list()
	PredictionModel_Pheno_List <- list()
	genoValues_List <- list()
    phenoValues_List <- list()
	genoValSimValues_List <- list()
	phenoValSimValues_List <- list()
 
    selectedNextGenGenoTable_GM_List <- list() 
    selectedNextGenGenoTableMod_GM_List <- list() 
#	selectedNewNextGenGenoTable_GM_List <- list() 

for(nCyc in startCycle:nCycles){
        #nCyc <- 2

    print(nCyc)

### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
   for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <-  Cycle_Progeny_F5_List[[nParameter]]
		nCrosses <-  nCrosses_List[[nParameter]] 
		
		nextGenGenoTable_GM <- generateMinus101GenoFormat_GM_Frontier_V1(Cycle_Progeny_F5,nCrosses)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
	
###############################################################################################
		if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){

			PredictionModels <- getUpdatedGSModel(newNextGenGenoTable_GM,genoValSimValues,phenoValSimValues,PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,trainTableWindowSize,bigmemory::as.matrix(trainGenoNewTablePreCycle),trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount,nCyc)
 
			PredictionModel_Geno <- PredictionModels[[1]]
			PredictionModel_Pheno <- PredictionModels[[2]]
			CombinedTableCount <- PredictionModels[[3]]
			trainGenoNewTablePreCycle <- PredictionModels[[4]]
			trainSimPhenoTablePreCycle <- PredictionModels[[5]]
			trainSimGenoValTablePreCycle <- PredictionModels[[6]]
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

        nextGenGenoTable_GM_List[[nParameter]] <- nextGenGenoTable_GM
        nextGenGenoTableMod_GM_List[[nParameter]] <- nextGenGenoTableMod_GM
		newNextGenGenoTable_GM_List[[nParameter]] <- newNextGenGenoTable_GM 
		PredictionModel_Geno_List[[nParameter]] <- PredictionModel_Geno
		PredictionModel_Pheno_List[[nParameter]] <- PredictionModel_Pheno
		genoValues_List[[nParameter]] <- genoValues
        phenoValues_List[[nParameter]] <- phenoValues
		genoValSimValues_List[[nParameter]] <- genoValSimValues 
		phenoValSimValues_List[[nParameter]] <- phenoValSimValues 

		
### ApplySelection to nextGenGenoTable_GM to select top 200 lines 

    no_selected <- GM_Param_List$nMates
	
	if(selectionOnSimulated ==TRUE){

		genoSelectionTable <- ApplySelection(genoValSimValues ,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)


		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedGenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedGenoIndividualIndices,]

		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		}
	}
    if(selectionOnSimulated==FALSE){

		genoSelectionTable <- ApplySelection(genoValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

    	print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
		print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
		print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		
		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		}
		
	}
	 
	    selectedNextGenGenoTable_GM_List[[nParameter]] <- selectedNextGenGenoTable_GM
        selectedNextGenGenoTableMod_GM_List[[nParameter]] <- selectedNextGenGenoTableMod_GM
 
 }

### Split Geno and Pheno values according to families  ########################
#### GM method for parent selection in cycle nCyc #####################################
 
  
	genoSelected_List <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (getGM_Frontier_FS_selectedGenoList_Cyc_2_N(selectedNextGenGenoTable_GM_List[[k]],selectedNextGenGenoTableMod_GM_List[[k]],PredictionModel_Geno_List[[k]],PredictionModel_Pheno_List[[k]],modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List))
 
 
    save(genoSelected_List,file=WorkspaceName)
  
## Reduce Parameter set to 30  
   
    for(nParameter in 1:nParameters){
    
      selectedParents_Geno_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[1]]
      CriterionValue_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nParameter]][[1]]
          
    } 

    
    count <- 1
	GainValues <- c()
	UsefulnessValues <- c()
	InbreedValues <- c()
	selectedParents_Geno_Par_List <- list()
	
    for (nPar1 in 1: nParameters){ 
	
	   nParameters2 <- length(selectedParents_Geno_List[[nCyc]][[nPar1]])
 
      for(nPar2 in 1:nParameters2){ 
	    GainValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][nPar2,1]
	    UsefulnessValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][nPar2,2]
	    InbreedValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][nPar2,3]
	    selectedParents_Geno_Par_List[[count]]  <-  selectedParents_Geno_List[[nCyc]][[nPar1]][[nPar2]]
	 
	    count <- count+1
	  }
    }
 
##### 
   
    Quartile3 <- summary(GainValues)[5]
	Quartile1 <- summary(GainValues)[2]
	MedianG <- 	summary(GainValues)[3]	
	
	
	indices1 <- which(GainValues >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues <= MedianG) & (GainValues >= Quartile1))
	indices2B <- which((GainValues >= MedianG) & (GainValues <= Quartile3))
	indices2CG <- (which(GainValues <= MedianG+1 & GainValues >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,5),sample(indices2B,5))

	indices3 <- which(GainValues <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange_G <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
###########

    Quartile3 <- summary(InbreedValues)[5]
	Quartile1 <- summary(InbreedValues)[2]
	MedianI <- summary(InbreedValues)[3]
			
	indices1 <- which(InbreedValues >= Quartile3)
    indicesSortRange1 <- sort.int(InbreedValues[indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1I <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((InbreedValues <= MedianI) & (InbreedValues >= Quartile1))
	indices2B <- which((InbreedValues >= MedianI) & (InbreedValues <= Quartile3))
	indices2CI <- (which(InbreedValues <= MedianI+1 & InbreedValues >= MedianI-1))[1]
	indicesRange2I <- c(sample(indices2A,5),sample(indices2B,5))
			
    indices3 <- which(InbreedValues <= Quartile1)
    indicesSortRange3 <- sort.int(InbreedValues[indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3I <- indicesSortRange3[[2]][1:10]
			
	indicesRange_I <- c(indicesRange1I,indicesRange2I,indicesRange3I)
	
	
	indicesRange_GI <- unique(c(indices2CG,indices2CI,indicesRange1I,indicesRange3G,sample(indicesRange2G,5),sample(indicesRange2I,5)))
	if(length(indicesRange_GI) >=30){
	  indicesRange <- indicesRange_GI[1:30]
	}else if(length(indicesRange_GI) <30){ 
	 
	  indicesRange <- c(indices2CG,indices2CI,indicesRange1I,indicesRange3G,sample(c(indicesRange2G,indicesRange2I),8))
	
	}

##########	
	
	
	nParameters <- length(indicesRange)
			
	for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	   
			   selectedParents_Geno_List[[nCyc]][[nParameter]] <-  selectedParents_Geno_Par_List[[nParameterIndex]]
			   CriterionValue_List[[nCyc]][[nParameter]] <- c(GainValues[nParameterIndex],UsefulnessValues[nParameterIndex],InbreedValues[nParameterIndex])
		
	} 
   
  #  save.image(WorkspaceName)
    gc() 	
  
#####################################################################################################################
    for(nParameter in 1:nParameters){
	
	    genoValues <- genoValues_List[[nParameter]] 
        phenoValues <- phenoValues_List[[nParameter]] 
        genoValSimValues <- genoValSimValues_List[[nParameter]] 
	    phenoValSimValues <- phenoValSimValues_List[[nParameter]] 
	       
    	uniqueParentIDs <-  unique(as.vector(selectedParents_Geno_List[[nCyc]][[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs

   		
### GM Output in selected individual ids and Num Progeny from selected parents
		   		
	
	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-   genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(genoValSimValues[GenoTableIndices_topN])
		
		topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
		topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- PhenoTableIndices_topN
		PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
#########################
	
	Cycle_Progeny_F5 <- Cycle_Progeny_F5_List[[nParameter]]
	
	nCrosses <- length(Cycle_Progeny_F5)
	nProgeny_per_Cross <- 10
        no_selected <- 200 
	Cycle_GenoData  <- (array(0,c(1,nMarkers,2,nCrosses*nProgeny_per_Cross)))
        selectedCycle_GenoData  <- (array(0,c(1,nMarkers,2,no_selected)))	
	nCrsInit <- 1
	nCrsFinal <- nProgeny_per_Cross
   
    for(nCrs in 1:nCrosses){

	  Cycle_GenoData[1,,,nCrsInit:nCrsFinal] <- Cycle_Progeny_F5[[nCrs]]
	 
	  nCrsInit <- nCrsFinal+1
	  nCrsFinal <-  nCrsFinal+ nProgeny_per_Cross 
	
	}
	
	if(selectionOnGeno == TRUE){
	   selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]
	}else if(selectionOnGeno == FALSE){ 
	   selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
	}
     	
        F5RILs_GM <- getF5RILs_FS_GM_Frontier(selectedCycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New)
	
	    Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
	
	}
	

	#	save.image(WorkspaceName) 
		gc() 
		
    simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)
	
	save(simResults_List_Cycle,file=WorkSpaceName)
	 
}

######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	   return(simResults_List)

}


   

 getUpdatedGSModel <- function(NewNextGenGenoTable_GM,genoValSimValues,phenoValSimValues,PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,TrainGenoNewTablePreCycle,TrainSimPhenoTablePreCycle,TrainSimGenoValTablePreCycle,CombinedTableCount,nCyc){
       
	  	newNextGenGenoTable_GM <- NewNextGenGenoTable_GM
		combinedTableCount <- CombinedTableCount
		trainGenoNewTablePreCycle <- TrainGenoNewTablePreCycle
		trainSimPhenoTablePreCycle <- TrainSimPhenoTablePreCycle
		trainSimGenoValTablePreCycle <- TrainSimGenoValTablePreCycle
						
		trainTableWindowSize <- TrainTableWindowSize
		nIndividuals <- length(genoValSimValues)
		IndividualNames <-rep(0,nIndividuals)
		IndividualNames<- paste("Ind",c(1:(nIndividuals)),sep="")

		

##################################################################################################
		names(phenoValSimValues)<- IndividualNames
		Mean_Fixed<- rep(1,nIndividuals)

		phenotypicValuesSimTable <- cbind(phenoValSimValues,Mean_Fixed)

#### Table for Simulated Genotypic Values ####################################

		names(genoValSimValues)<-IndividualNames
		Mean_Fixed<- rep(1,nIndividuals)

		genotypicValuesSimTable <- cbind(genoValSimValues,Mean_Fixed)

#############################################################################################
	
    if(updateType == "FullSet"){
	
	   PredictionModel_Geno_PreCycle <- PredictionModel_Geno
       PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno

	    indices<-c(1:nIndividuals)
		initSeed <- 25+nCyc
        set.seed(initSeed)
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

		TrainGenoNewTableComb <- as.matrix(rbind(as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))

       # trainSimGenoValTableComb <- rbind((trainSimGenoValTablePreCycle),trainSimGenoValTable)

       # trainSimPhenoValTableComb <- rbind((trainSimPhenoTablePreCycle),trainSimPhenoTable)

        print(dim(TrainGenoNewTableComb))

        TrainSimGenoValTableComb <- as.matrix(rbind(as.matrix(trainSimGenoValTablePreCycle),as.matrix(trainSimGenoValTable)))
        TrainSimPhenoValTableComb <- as.matrix(rbind(as.matrix(trainSimPhenoTablePreCycle),as.matrix(trainSimPhenoTable)))
      
        trainGenoNewTableComb <- apply(TrainGenoNewTableComb,2,function(x) as.numeric(as.character(x))) 
	    trainSimGenoValTableComb <- apply(TrainSimGenoValTableComb,2,function(x) as.numeric(as.character(x)))
        trainSimPhenoValTableComb <- apply(TrainSimPhenoValTableComb,2,function(x) as.numeric(as.character(x))) 


        if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){
            if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0){

				if(modelType=="BayesB"){
					PredictionModel_Geno <- (buildBayesBModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					PredictionModel_Pheno <-(buildBayesBModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
				}else if(modelType=="BL"){
			
					PredictionModel_Geno <- (buildBLModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					PredictionModel_Pheno <-(buildBLModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
			    }else if(modelType=="RRBLUP"){

                        PredictionModel_Geno <- (buildRRBLUPModel((trainGenoNewTableComb),trainSimGenoValTableComb,(testGenoNewTable),testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable))

                }else if(modelType=="RRBLUP_REML"){

                        PredictionModel_Geno <- buildRRBLUPModel_REML((trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        PredictionModel_Pheno <- buildRRBLUPModel_REML((trainGenoNewTableComb), trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

                }
				
		    }
		 
		    trainGenoNewTablePreCycle <-  trainGenoNewTableComb
            trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)
            trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
			
		   print(sd(PredictionModel_Pheno[[1]])) 
           print(sd(PredictionModel_Geno[[1]])) 
           print(length(levels(factor(PredictionModel_Pheno[[1]]))))
           print(length(levels(factor(PredictionModel_Geno[[1]]))))	 

           genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Geno)
           phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable_GM,PredictionModel_Pheno)
           print(sd(phenoValues))
           print(sd(genoValues))
			
	
	if(length(levels(factor(PredictionModel_Pheno[[1]]))) ==1 || length(levels(factor(PredictionModel_Geno[[1]]))) ==1 || (is.na(sd(phenoValues)) || sd(phenoValues)==0) || (is.na(sd(genoValues)) || sd(genoValues)==0)){

                PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                PredictionModel_Geno <- PredictionModel_Geno_PreCycle

    }
			
	        rm(trainGenoNewTableComb)
            rm(trainSimGenoValTableComb)
            rm(trainSimPhenoValTableComb)
            gc()

    }
}

    if(updateType=="TrainingCycleWindow"){

			
	PredictionModel_Geno_PreCycle <- PredictionModel_Geno
	PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno
        
        indices<-c(1:nIndividuals)
        indices<-c(1:nIndividuals)
		initSeed <- 25+nCyc
        set.seed(initSeed)
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
         # bigmemory::as.matrix

        trainGenoNewTableComb <- (rbind(as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))
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

                        PredictionModel_Geno <- (buildRRBLUPModel((trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable))

                        PredictionModel_Pheno <-(buildRRBLUPModel((trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable))

           }else if(modelType=="RRBLUP_REML"){
				
						PredictionModel_Geno <- tryCatch(buildRRBLUPModel_REML((trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable),error=function(e) print(paste("Error build RRREML",nCyc,nrep,"-",e)))

                        PredictionModel_Pheno <- tryCatch(buildRRBLUPModel_REML((trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable),error=function(e) print(paste("Error build RRREML",nCyc,nrep,"-",e)))

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
	
	nLines_trainGeno <- nrow(trainGenoNewTablePreCycle)
	rownames(trainGenoNewTablePreCycle) <- paste("Ind",c(1:nLines_trainGeno),sep="")
    rownames(trainSimPhenoTablePreCycle) <- paste("Ind",c(1:nLines_trainGeno),sep="")
	rownames(trainSimGenoValTablePreCycle) <- paste("Ind",c(1:nLines_trainGeno),sep="")
	
	return(list(PredictionModel_Geno,PredictionModel_Pheno,combinedTableCount,trainGenoNewTablePreCycle,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle))

}



 runSimulations20X_BD_Wght_GM_Frontier_V3A <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach){

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
	  nLines <- nFamilies * cycle1_nProgeny 
      nProgeny_per_Cross <- nLines/nCrosses		
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
	  
	   
	  #trainGenoNewTablePreCycle <- as.big.matrix(trainGenoNewTable)
	 
 	  trainGeno_backingFileName <- paste("trainTable_Cyc1",condition,"_",Rep,".bin",sep="") 
      trainGeno_descriptorFileName <- paste("trainTable_Cyc1",condition,"_",Rep,".desc",sep="")
		 			
	  trainGenoNewTablePreCycle <- deepcopy(as.big.matrix(trainGenoNewTable),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
  
	  
      trainSimPhenoTablePreCycle <- trainSimPhenoValTable
      trainSimGenoValTablePreCycle <- trainSimGenoValTable

	  trainTableWindowSize <- TrainTableWindowSize
	  
	  
## Weighting Parameters	  
	  Weighted <- weighted
	  WghtMethod <- wghtMethod

	  
## GM Parameters 
	  
	  GM_Param_List <- GM_Frontier_Param_List
	  no_cores <- GM_Param_List$no_cores
	  startCycle <- StartCycle
	  WorkspaceName <- WorkSpaceName
	  nCores_ForEach <- noCores_ForEach 
	  
	  #save.image(WorkspaceName)
#################################################################################################
#################################################################################################
### Initialize variables	
	
   	  j<-1
 	  CombinedTableCount <- 1
################################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
################################################################################################	  	  
#### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable_GM <- generateWeightedGenoTable(cycle1GenoTable_GM,no_QTL)
	cycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(cycle1GenoTable_GM)
	  
	genoValSimValues <- simulateGenotypicValues(cycle1GenoTable_GM,no_QTL) 
	phenoValSimValues <- simulatePhenotypicValues(cycle1GenoTable_GM,varE,no_QTL) 

	nIndividuals <- nrow(cycle1GenoTable_GM)
	
	
## Write Genotable DF 
    populations <- rep(0,nFamilies*100)
	init <- 1
	for(nPop in 1:nFamilies){

		final <- nPop*100
		populations[init:final]<- rep(nPop,100)

		init <- final+1

	}
	
	cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 
	
	
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

###### Select island pairs based on migration policy 
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
	 

# #################################################################################
   
    no_selected <- numberSelected
	nCrosses <- noCrosses
	nProgeny <- no_Progeny
	cycle1_nProgeny <- nProgeny
	nFamilies <- 20
	nLines <- nFamilies*cycle1_nProgeny 
    nProgeny_per_Cross <- nLines/nCrosses	
   	
    Cycle_GenoData  <- (array(0,c(1,nMarkers,2,cycle1_nProgeny*nFamilies)))
    selectedCycle_GenoData  <- (array(0,c(1,nMarkers,2,no_selected))) 
   
    initIndex <- 1
    finalIndex <- nProgeny
   
    for(nCrs in 1:nFamilies){

	
	 Cycle_GenoData[1,,,initIndex:finalIndex] <- cycle1_Geno_Data[nCrs,,,] 
	 initIndex <- initIndex + nProgeny
	 finalIndex <- finalIndex + nProgeny
    }

############  
### ApplySelection to select top nSel lines 
	
	if(selectionOnSimulated ==TRUE){

		genoSelectionTable <- ApplySelection(genoValSimValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)


		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedGenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedGenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]

		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedPhenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}
	}
    if(selectionOnSimulated==FALSE){
		genoSelectionTable <- ApplySelection(genoValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

    	print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
		print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
		print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
		        genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedGenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedGenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]
		
		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedPhenoIndividualIndices,]
		
		    selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}
		
	} 

#### GM method for parent selection in cycle 1 #####################################
 nParameters <- 30 
 selectedParents_Geno_List <- rep(list(rep(list(list()),nParameters)),nCycles)
 CriterionValue_List <- rep(list(rep(list(list()),nParameters)),nCycles) 

 
 genoSelected_List <- getGM_Frontier_FS_selectedGenoList_Cyc1(selectedCycle1GenoTable_GM,selectedCycle1GenoTableMod_GM,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List)

 selectedParents_Geno_List[[nCyc]] <- genoSelected_List[[1]]
 CriterionValue_List[[nCyc]] <- genoSelected_List[[2]]
 selectedParents_Geno <- genoSelected_List[[1]]
 nParameters <- genoSelected_List[[3]]
 


save(genoSelected_List,file=WorkspaceName)

 #save.image(WorkspaceName) 	
 gc()
 
################################################################################### 
## Initiate Output list 	  
	  
	  GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	 
	  PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  
      attainedPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      attainedGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	
	  topNGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      GenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)

      topNPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      PhenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  Cycle_Progeny_F5_List <- list()
      nCrosses_List <- list()
	  
	  
	 
##
 
   for(nParameter in 1:nParameters){
   
        uniqueParentIDs <-  unique(as.vector(selectedParents_Geno[[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
### GM Output in selected individual ids 

############################

	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]]<- genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <- max(genoValSimValues[GenoTableIndices_topN])
		
	    topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN_List
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
	   
	    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
	    topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <-  PhenoTableIndices_topN
	    PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		
	    
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
###################################################################################	
## Get F5 RILs for selected parents in cycle 1
        
         
        F5RILs_GM <- getF5RILs_FS_GM_Frontier(selectedCycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New,nProgeny_per_Cross)

        Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
		
	   
	 	     
	}

	simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)
	
	save(simResults_List_Cycle,file=WorkspaceName)
	
    startCycle <- startCycle + 1 	 
	 
	gc() 
	
###########################################################################

    startCycle <- 2
    nCores <- nCores_ForEach 
    #registerDoParallel(nCores)
    nextGenGenoTable_GM_List <- list()
    nextGenGenoTableMod_GM_List <- list()
    newNextGenGenoTable_GM_List <- list()
	
	genoValues_List <- list()
    phenoValues_List <- list()
	genoValSimValues_List <- list()
	phenoValSimValues_List <- list()
	
	PredictionModel_Geno_List <- list()
	PredictionModel_Pheno_List <- list()
		   
	trainGenoNewTablePreCycle_List <- list()
	trainGenoNewTablePreCycle_Mat_List <- list()
	
	trainSimPhenoTablePreCycle_List <- list()
	trainSimGenoValTablePreCycle_List <- list()
	
	trainGeno_backingFileName_List <- list()
    trainGeno_descriptorFileName_List <- list()	
 
    selectedNextGenGenoTable_GM_List <- list() 
    selectedNextGenGenoTableMod_GM_List <- list() 

for(nCyc in startCycle:nCycles){
      
    print(nCyc)

### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
    for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <-  Cycle_Progeny_F5_List[[nParameter]]
		nCrosses <-  nCrosses_List[[nParameter]] 
		
		nextGenGenoTable_GM <- generateMinus101GenoFormat_GM_Frontier_V1(Cycle_Progeny_F5,nCrosses,nProgeny_per_Cross)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
		
		nextGenGenoTable_GM_List[[nParameter]] <- nextGenGenoTable_GM
        nextGenGenoTableMod_GM_List[[nParameter]] <- nextGenGenoTableMod_GM
		newNextGenGenoTable_GM_List[[nParameter]] <- newNextGenGenoTable_GM 
		genoValSimValues_List[[nParameter]] <- genoValSimValues 
		phenoValSimValues_List[[nParameter]] <- phenoValSimValues 
		
## Write genotable df

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1
		}
		
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
		rm(nextGenGenoTable_AlleleFormat)
		
	}
	
############################################################################################### 
	if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
		
		  if(nCyc ==2){
		   
		    trainGenoNewTablePreCycle <- attach.big.matrix(trainGeno_descriptorFileName)
   		    trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)

			PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_List[[j]],genoValSimValues_List[[j]],phenoValSimValues_List[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,trainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount,nCyc))
		  }else if(nCyc >2){
		   
			for(nParameter in 1:nParameters){
			  
			  trainGeno_FileName <- trainGeno_backingFileName_List[[nParameter]] 
			  trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- read.table(trainGeno_FileName,sep="\t")
			}
		    PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_List[[j]],genoValSimValues_List[[j]],phenoValSimValues_List[[j]],PredictionModel_Geno_List[[j]],PredictionModel_Pheno_List[[j]],modelType,updateType,trainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount,nCyc))
 
		  }
		  gc()
		 		   		
		for(nParameter in 1:nParameters){
   
			PredictionModel_Geno_List[[nParameter]] <- PredictionModels[[nParameter]][[1]]
			PredictionModel_Pheno_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
			CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		   
		    PredictionModel_Geno <- PredictionModel_Geno_List[[nParameter]] 
		    PredictionModel_Pheno <- PredictionModel_Pheno_List[[nParameter]]  
			
								
			trainGeno_FileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".txt",sep="") 
			trainGeno_backingFileName_List[[nParameter]] <- trainGeno_FileName
			
			trainGenoNewTablePreCycle_List <- PredictionModels[[nParameter]][[4]]
			trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
			trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
			
		   	write.table(trainGenoNewTablePreCycle_List,trainGeno_FileName,sep="\t") 
			rm(trainGenoNewTablePreCycle_List)
			
		    nextGenGenoTable_GM <- nextGenGenoTable_GM_List[[nParameter]] 
            nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_List[[nParameter]] 
		    newNextGenGenoTable_GM <- newNextGenGenoTable_GM_List[[nParameter]]
			
			Freq <- getFreq(nextGenGenoTable_GM,FavAllele)
		  
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

				}else if(Weighted==TRUE && WghtMethod =="DW"){

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

			genoValues_List[[nParameter]] <- genoValues
			phenoValues_List[[nParameter]] <- phenoValues
		
			genoValSimValues <- genoValSimValues_List[[nParameter]]  
			phenoValSimValues <- phenoValSimValues_List[[nParameter]]   
		
				
### ApplySelection to nextGenGenoTable_GM to select top 200 lines 

  
		no_selected <- GM_Param_List$nMates
	
		if(selectionOnSimulated ==TRUE){

			genoSelectionTable <- ApplySelection(genoValSimValues,no_selected)
			phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)

		  if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedGenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedGenoIndividualIndices,]

		  }else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		  }
		}
        if(selectionOnSimulated==FALSE){

			genoSelectionTable <- ApplySelection(genoValues,no_selected)
			phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

			print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
			print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
			print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

			if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		
		    }else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		  }
		
	    }
	 
	    selectedNextGenGenoTable_GM_List[[nParameter]] <- selectedNextGenGenoTable_GM
        selectedNextGenGenoTableMod_GM_List[[nParameter]] <- selectedNextGenGenoTableMod_GM
 
    }

}


    if(modelUpdate ==FALSE || selectionOnSimulated==TRUE){
	
	    for(nParameter in 1:nParameters){
		
		    PredictionModel_Geno_List[[nParameter]] <-  PredictionModel_Geno 
		    PredictionModel_Pheno_List[[nParameter]] <-  PredictionModel_Pheno  
   			
		    nextGenGenoTable_GM <- nextGenGenoTable_GM_List[[nParameter]] 
            nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_List[[nParameter]] 
		    newNextGenGenoTable_GM <- newNextGenGenoTable_GM_List[[nParameter]]
			
			Freq <- getFreq(nextGenGenoTable_GM,FavAllele)
  
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

				}else if(Weighted==TRUE && WghtMethod =="DW"){

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

			genoValues_List[[nParameter]] <- genoValues
			phenoValues_List[[nParameter]] <- phenoValues
		
			genoValSimValues <- genoValSimValues_List[[nParameter]]  
			phenoValSimValues <- phenoValSimValues_List[[nParameter]]   
		
### ApplySelection to nextGenGenoTable_GM to select top 200 lines 

			no_selected <- GM_Param_List$nMates
		
			if(selectionOnSimulated ==TRUE){

				genoSelectionTable <- ApplySelection(genoValSimValues,no_selected)
				phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)

			  if(selectionOnGeno==TRUE){

				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
				genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
				selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedGenoIndividualIndices,]	
				selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedGenoIndividualIndices,]

			  }else if(selectionOnGeno==FALSE){

				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
				genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
				
				selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
				selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
			  }
			}
			if(selectionOnSimulated==FALSE){

				genoSelectionTable <- ApplySelection(genoValues,no_selected)
				phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

				print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
				print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
				print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

				if(selectionOnGeno==TRUE){

				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
				genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
				phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
				genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
				
				selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
				selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
			
				}else if(selectionOnGeno==FALSE){

				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
				genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
				phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
				genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
				
				selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
				selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
			  }
		
			}
	 
			selectedNextGenGenoTable_GM_List[[nParameter]] <- selectedNextGenGenoTable_GM
			selectedNextGenGenoTableMod_GM_List[[nParameter]] <- selectedNextGenGenoTableMod_GM
		}
	}


### Split Geno and Pheno values according to families  ########################
#### GM method for parent selection in cycle nCyc #####################################
       
	genoSelected_List <- foreach(k=1:nParameters) %dopar% (LTSGenoIMV2::getGM_Frontier_FS_selectedGenoList_Cyc_2_N(selectedNextGenGenoTable_GM_List[[k]],selectedNextGenGenoTableMod_GM_List[[k]],PredictionModel_Geno_List[[k]],PredictionModel_Pheno_List[[k]],modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List))

 
    save(genoSelected_List,simResults_List_Cycle,file=WorkSpaceName)
	rm(simResults_List_Cycle)
	gc()
## Reduce Parameter set to 30  
   
    for(nParameter in 1:nParameters){
    
      selectedParents_Geno_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[1]]
      CriterionValue_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nParameter]][[1]]
          
    }  
 
    count <- 1
	GainValues <- c()
	UsefulnessValues <- c()
	InbreedValues <- c()
	selectedParents_Geno_Par_List <- list()
	
    for (nPar1 in 1: nParameters){ 
	
	  nParameters2 <- length(selectedParents_Geno_List[[nCyc]][[nPar1]])
 
      for(nPar2 in 1:nParameters2){ 
	    GainValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][1]
	    UsefulnessValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][2]
	    InbreedValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][3]
	    selectedParents_Geno_Par_List[[count]]  <-  selectedParents_Geno_List[[nCyc]][[nPar1]][[nPar2]]
	 
	    count <- count+1
	  }
    }
 
	
	### 
	
    Quartile3 <- summary(GainValues)[5]
	Quartile1 <- summary(GainValues)[2]
	MedianG <- 	summary(GainValues)[3]	
		
	indices1 <- which(GainValues >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues <= MedianG) & (GainValues >= Quartile1))
	indices2B <- which((GainValues >= MedianG) & (GainValues <= Quartile3))
	indices2CG <- (which(GainValues <= MedianG+1 & GainValues >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
###########
				
	for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	   
			   selectedParents_Geno_List[[nCyc]][[nParameter]] <-  selectedParents_Geno_Par_List[[nParameterIndex]]
			   CriterionValue_List[[nCyc]][[nParameter]] <- c(GainValues[nParameterIndex],UsefulnessValues[nParameterIndex],InbreedValues[nParameterIndex])
	} 
   
  #  save.image(WorkspaceName)
    gc() 	
  
#####################################################################################################################
    for(nParameter in 1:nParameters){
	
	    genoValues <- genoValues_List[[nParameter]] 
        phenoValues <- phenoValues_List[[nParameter]] 
        genoValSimValues <- genoValSimValues_List[[nParameter]] 
	    phenoValSimValues <- phenoValSimValues_List[[nParameter]] 
	       
    	uniqueParentIDs <-  unique(as.vector(selectedParents_Geno_List[[nCyc]][[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs

### GM Output in selected individual ids and Num Progeny from selected parents
		
	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-   genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(genoValSimValues[GenoTableIndices_topN])
		
		topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
		topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- PhenoTableIndices_topN
		PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 

#########################
	
		
		no_selected <- numberSelected
		nCrosses <- noCrosses
		nProgeny <- no_Progeny
		cycle1_nProgeny <- nProgeny
		nFamilies <- 20
		nLines <- nFamilies*cycle1_nProgeny 
		nProgeny_per_Cross <- nLines/nCrosses			
	
				
		Cycle_GenoData  <- (array(0,c(1,nMarkers,2,nCrosses*nProgeny_per_Cross)))
		selectedCycle_GenoData  <- (array(0,c(1,nMarkers,2,no_selected)))	
		nCrsInit <- 1
		nCrsFinal <- nProgeny_per_Cross
		Cycle_Progeny_F5 <- Cycle_Progeny_F5_List[[nParameter]]
   
		for(nCrs in 1:nCrosses){

		  Cycle_GenoData[1,,,nCrsInit:nCrsFinal] <- Cycle_Progeny_F5[[nCrs]]
		 
		  nCrsInit <- nCrsFinal+1
		  nCrsFinal <-  nCrsFinal+ nProgeny_per_Cross 
		
		}
		
		if(selectionOnGeno == TRUE){
		   selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]
		}else if(selectionOnGeno == FALSE){ 
		   selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}
			
			F5RILs_GM <- getF5RILs_FS_GM_Frontier(selectedCycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New,nProgeny_per_Cross)
		
			Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
			nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
		
	}
	

	#	save.image(WorkspaceName) 
		gc() 
		
    simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)
	
}

######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	   return(simResults_List)

}

###### 


 runSimulations20X_BD_Wght_GM_Frontier_V3A_Split <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,WorkSpaceName,StartCycle,FinishCycle,GM_Frontier_Param_List,noCores_ForEach){

	  options(warn=-1) 	
### Assign Variables #########################################################################
### Selection Parameters 

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
	  nLines <- nFamilies * cycle1_nProgeny 
      nProgeny_per_Cross <- nLines/nCrosses		
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
	 
 	  trainGeno_backingFileName <- paste("trainTable_Cyc1",condition,"_",Rep,".bin",sep="") 
      trainGeno_descriptorFileName <- paste("trainTable_Cyc1",condition,"_",Rep,".desc",sep="")
		 			
	  #trainGenoNewTablePreCycle <- deepcopy(as.big.matrix(trainGenoNewTable),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
  
      	  
      trainSimPhenoTablePreCycle <- trainSimPhenoValTable
      trainSimGenoValTablePreCycle <- trainSimGenoValTable

	  trainTableWindowSize <- TrainTableWindowSize
	  
	  
## Weighting Parameters	  
	  Weighted <- weighted
	  WghtMethod <- wghtMethod

	  
## GM Parameters 
	  
	  GM_Param_List <- GM_Frontier_Param_List
	  no_cores <- GM_Param_List$no_cores
	  startCycle <- StartCycle
	  finishCycle <- FinishCycle
	  WorkspaceName <- WorkSpaceName
	  nCores_ForEach <- noCores_ForEach 
	  
	  #save.image(WorkspaceName)
#################################################################################################
#################################################################################################
### Initialize variables	
	
   	  j<-1
 	  CombinedTableCount <- 1
################################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
################################################################################################	  	  
#### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable_GM <- generateWeightedGenoTable(cycle1GenoTable_GM,no_QTL)
	cycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(cycle1GenoTable_GM)
	  
	genoValSimValues <- simulateGenotypicValues(cycle1GenoTable_GM,no_QTL) 
	phenoValSimValues <- simulatePhenotypicValues(cycle1GenoTable_GM,varE,no_QTL) 

	nIndividuals <- nrow(cycle1GenoTable_GM)
	
	
## Write Genotable DF 
    populations <- rep(0,nFamilies*100)
	init <- 1
	for(nPop in 1:nFamilies){

		final <- nPop*100
		populations[init:final]<- rep(nPop,100)

		init <- final+1

	}
	
	cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 
	
	
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

###### Select island pairs based on migration policy 
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
	 

# #################################################################################
   
    no_selected <- numberSelected
	nCrosses <- noCrosses
	nProgeny <- no_Progeny
	cycle1_nProgeny <- nProgeny
	nFamilies <- 20
	nLines <- nFamilies*cycle1_nProgeny 
    nProgeny_per_Cross <- nLines/nCrosses	
   	
    Cycle_GenoData  <- (array(0,c(1,nMarkers,2,cycle1_nProgeny*nFamilies)))
    selectedCycle_GenoData  <- (array(0,c(1,nMarkers,2,no_selected))) 
   
    initIndex <- 1
    finalIndex <- nProgeny
   
    for(nCrs in 1:nFamilies){

	
	 Cycle_GenoData[1,,,initIndex:finalIndex] <- cycle1_Geno_Data[nCrs,,,] 
	 initIndex <- initIndex + nProgeny
	 finalIndex <- finalIndex + nProgeny
    }

############  
### ApplySelection to select top nSel lines 
	
	if(selectionOnSimulated ==TRUE){

		genoSelectionTable <- ApplySelection(genoValSimValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)


		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedGenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedGenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]

		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedPhenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}
	}
    if(selectionOnSimulated==FALSE){
		genoSelectionTable <- ApplySelection(genoValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

    	print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
		print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
		print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
		        genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedGenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedGenoIndividualIndices,]
			selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]
		
		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
			
			selectedCycle1GenoTable_GM <- cycle1GenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedCycle1GenoTableMod_GM <-cycle1GenoTableMod_GM[selectedPhenoIndividualIndices,]
		
		    selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}
		
	} 

#### GM method for parent selection in cycle 1 #####################################
 nParameters <- 30 
 selectedParents_Geno_List <- rep(list(rep(list(list()),nParameters)),nCycles)
 CriterionValue_List <- rep(list(rep(list(list()),nParameters)),nCycles) 

 
 genoSelected_List <- getGM_Frontier_FS_selectedGenoList_Cyc1(selectedCycle1GenoTable_GM,selectedCycle1GenoTableMod_GM,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List)

 selectedParents_Geno_List[[nCyc]] <- genoSelected_List[[1]]
 CriterionValue_List[[nCyc]] <- genoSelected_List[[2]]
 selectedParents_Geno <- genoSelected_List[[1]]
 nParameters <- genoSelected_List[[3]]
 


save(genoSelected_List,file=WorkspaceName)

 #save.image(WorkspaceName) 	
 gc()
 
################################################################################### 
## Initiate Output list 	  
	  
	  GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	 
	  PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  
      attainedPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      attainedGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	
	  topNGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      GenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)

      topNPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      PhenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  Cycle_Progeny_F5_List <- list()
      nCrosses_List <- list()
	  
	  
	 
##
 
   for(nParameter in 1:nParameters){
   
        uniqueParentIDs <-  unique(as.vector(selectedParents_Geno[[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
### GM Output in selected individual ids 

############################

	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]]<- genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <- max(genoValSimValues[GenoTableIndices_topN])
		
	    topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN_List
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
	   
	    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
	    topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <-  PhenoTableIndices_topN
	    PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		
	    
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
###################################################################################	
## Get F5 RILs for selected parents in cycle 1
        
         
        F5RILs_GM <- getF5RILs_FS_GM_Frontier(selectedCycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New,nProgeny_per_Cross)

        Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
		
	   
	 	     
	}

	simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)
	
	save(simResults_List_Cycle,file=WorkspaceName)
	
    startCycle <- startCycle + 1 	 
	 
	gc() 
	
###########################################################################

    startCycle <- 2
    nCores <- nCores_ForEach 
    #registerDoParallel(nCores)
    nextGenGenoTable_GM_List <- list()
    nextGenGenoTableMod_GM_List <- list()
    newNextGenGenoTable_GM_List <- list()
	
	genoValues_List <- list()
    phenoValues_List <- list()
	genoValSimValues_List <- list()
	phenoValSimValues_List <- list()
	
	PredictionModel_Geno_List <- list()
	PredictionModel_Pheno_List <- list()
		   
	trainGenoNewTablePreCycle_List <- list()
	trainGenoNewTablePreCycle_Mat_List <- list()
	
	trainSimPhenoTablePreCycle_List <- list()
	trainSimGenoValTablePreCycle_List <- list()
	
	trainGeno_backingFileName_List <- list()
    trainGeno_descriptorFileName_List <- list()	
 
    selectedNextGenGenoTable_GM_List <- list() 
    selectedNextGenGenoTableMod_GM_List <- list() 

for(nCyc in startCycle:finishCycle){
      
    print(nCyc)

### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
    for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <- Cycle_Progeny_F5_List[[nParameter]]
		nCrosses <-  nCrosses_List[[nParameter]] 
						
		
		nextGenGenoTable_GM <- generateMinus101GenoFormat_GM_Frontier_V1(Cycle_Progeny_F5,nCrosses,nProgeny_per_Cross)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
		
		nextGenGenoTable_GM_List[[nParameter]] <- nextGenGenoTable_GM
        nextGenGenoTableMod_GM_List[[nParameter]] <- nextGenGenoTableMod_GM
		newNextGenGenoTable_GM_List[[nParameter]] <- newNextGenGenoTable_GM 
		genoValSimValues_List[[nParameter]] <- genoValSimValues 
		phenoValSimValues_List[[nParameter]] <- phenoValSimValues 
		
## Write genotable df

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1
		}
		
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
		rm(nextGenGenoTable_AlleleFormat)
		
	}
	
		
############################################################################################### 
	if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
		
		  if(nCyc ==2){
		   
		    # trainGenoNewTablePreCycle <- attach.big.matrix(trainGeno_descriptorFileName)
   		    
			trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)

			PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_List[[j]],genoValSimValues_List[[j]],phenoValSimValues_List[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,trainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount,nCyc))
		  }else if(nCyc >2){
		   
			for(nParameter in 1:nParameters){
			  
			  trainGeno_FileName <- trainGeno_backingFileName_List[[nParameter]] 
			  trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- read.table(trainGeno_FileName,sep="\t")
			}
		    PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_List[[j]],genoValSimValues_List[[j]],phenoValSimValues_List[[j]],PredictionModel_Geno_List[[j]],PredictionModel_Pheno_List[[j]],modelType,updateType,trainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount,nCyc))
 
		  }
		  gc()
		  
				 		   		
		for(nParameter in 1:nParameters){
   
			PredictionModel_Geno_List[[nParameter]] <- PredictionModels[[nParameter]][[1]]
			PredictionModel_Pheno_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
			CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		   
		    PredictionModel_Geno <- PredictionModel_Geno_List[[nParameter]] 
		    PredictionModel_Pheno <- PredictionModel_Pheno_List[[nParameter]]  
			
								
			trainGeno_FileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".txt",sep="") 
			trainGeno_backingFileName_List[[nParameter]] <- trainGeno_FileName
			
			trainGenoNewTablePreCycle_List <- PredictionModels[[nParameter]][[4]]
			trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
			trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
			
		   	write.table(trainGenoNewTablePreCycle_List,trainGeno_FileName,sep="\t") 
			rm(trainGenoNewTablePreCycle_List)
			
		    nextGenGenoTable_GM <- nextGenGenoTable_GM_List[[nParameter]] 
            nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_List[[nParameter]] 
		    newNextGenGenoTable_GM <- newNextGenGenoTable_GM_List[[nParameter]]
			
			Freq <- getFreq(nextGenGenoTable_GM,FavAllele)
		  
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

				}else if(Weighted==TRUE && WghtMethod =="DW"){

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

			genoValues_List[[nParameter]] <- genoValues
			phenoValues_List[[nParameter]] <- phenoValues
		
			genoValSimValues <- genoValSimValues_List[[nParameter]]  
			phenoValSimValues <- phenoValSimValues_List[[nParameter]]   
		
				
### ApplySelection to nextGenGenoTable_GM to select top 200 lines 

  
		no_selected <- GM_Param_List$nMates
	
		if(selectionOnSimulated ==TRUE){

			genoSelectionTable <- ApplySelection(genoValSimValues,no_selected)
			phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)

		  if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedGenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedGenoIndividualIndices,]

		  }else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		  }
		}
        if(selectionOnSimulated==FALSE){

			genoSelectionTable <- ApplySelection(genoValues,no_selected)
			phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

			print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
			print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
			print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

			if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		
		    }else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
			
			selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
            selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
		  }
		
	    }
	 
	    selectedNextGenGenoTable_GM_List[[nParameter]] <- selectedNextGenGenoTable_GM
        selectedNextGenGenoTableMod_GM_List[[nParameter]] <- selectedNextGenGenoTableMod_GM
 
    }

}


    if(modelUpdate ==FALSE || selectionOnSimulated==TRUE){
	
	    for(nParameter in 1:nParameters){
		
		    PredictionModel_Geno_List[[nParameter]] <-  PredictionModel_Geno 
		    PredictionModel_Pheno_List[[nParameter]] <-  PredictionModel_Pheno  
   			
		    nextGenGenoTable_GM <- nextGenGenoTable_GM_List[[nParameter]] 
            nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_List[[nParameter]] 
		    newNextGenGenoTable_GM <- newNextGenGenoTable_GM_List[[nParameter]]
			
			Freq <- getFreq(nextGenGenoTable_GM,FavAllele)
  
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

				}else if(Weighted==TRUE && WghtMethod =="DW"){

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

			genoValues_List[[nParameter]] <- genoValues
			phenoValues_List[[nParameter]] <- phenoValues
		
			genoValSimValues <- genoValSimValues_List[[nParameter]]  
			phenoValSimValues <- phenoValSimValues_List[[nParameter]]   
		
### ApplySelection to nextGenGenoTable_GM to select top 200 lines 

			no_selected <- GM_Param_List$nMates
		
			if(selectionOnSimulated ==TRUE){

				genoSelectionTable <- ApplySelection(genoValSimValues,no_selected)
				phenoSelectionTable <- ApplySelection(phenoValSimValues,no_selected)

			  if(selectionOnGeno==TRUE){

				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
				genoSimSelectedValues <- genoValSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoValSimValues[selectedGenoIndividualIndices]
				selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedGenoIndividualIndices,]	
				selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedGenoIndividualIndices,]

			  }else if(selectionOnGeno==FALSE){

				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
				genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoValSimValues[selectedPhenoIndividualIndices]
				
				selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
				selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
			  }
			}
			if(selectionOnSimulated==FALSE){

				genoSelectionTable <- ApplySelection(genoValues,no_selected)
				phenoSelectionTable <- ApplySelection(phenoValues,no_selected)

				print(paste("genoSelectionTable-dim-",dim(genoSelectionTable)))
				print(paste("phenoSelectionTable-dim-",dim(phenoSelectionTable)))
				print(paste("selectedGenoIndividualIndices",is.na(phenoSelectionTable[,2])))

				if(selectionOnGeno==TRUE){

				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
				genoSelectedValues <- genoValues[selectedGenoIndividualIndices]
				phenoSelectedValues <- phenoValues[selectedGenoIndividualIndices]
				genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
				
				selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
				selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
			
				}else if(selectionOnGeno==FALSE){

				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
				genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
				phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
				genoSimSelectedValues <- genoValSimValues[selectedPhenoIndividualIndices]
				
				selectedNextGenGenoTable_GM <- nextGenGenoTable_GM[selectedPhenoIndividualIndices,]	
				selectedNextGenGenoTableMod_GM <- nextGenGenoTableMod_GM[selectedPhenoIndividualIndices,]
			  }
		
			}
	 
			selectedNextGenGenoTable_GM_List[[nParameter]] <- selectedNextGenGenoTable_GM
			selectedNextGenGenoTableMod_GM_List[[nParameter]] <- selectedNextGenGenoTableMod_GM
		}
	}


### Split Geno and Pheno values according to families  ########################
#### GM method for parent selection in cycle nCyc #####################################
       
	genoSelected_List <- foreach(k=1:nParameters) %dopar% (LTSGenoIMV2::getGM_Frontier_FS_selectedGenoList_Cyc_2_N(selectedNextGenGenoTable_GM_List[[k]],selectedNextGenGenoTableMod_GM_List[[k]],PredictionModel_Geno_List[[k]],PredictionModel_Pheno_List[[k]],modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List))

    
 
    save(genoSelected_List,simResults_List_Cycle,file=WorkSpaceName)
	rm(simResults_List_Cycle)
	gc()
## Reduce Parameter set to 30  
   
    for(nParameter in 1:nParameters){
    
      selectedParents_Geno_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[1]]
      CriterionValue_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nParameter]][[1]]
          
    }  
 
    count <- 1
	GainValues <- c()
	UsefulnessValues <- c()
	InbreedValues <- c()
	selectedParents_Geno_Par_List <- list()
	
    for (nPar1 in 1: nParameters){ 
	
	  nParameters2 <- length(selectedParents_Geno_List[[nCyc]][[nPar1]])
 
      for(nPar2 in 1:nParameters2){ 
	    GainValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][1]
	    UsefulnessValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][2]
	    InbreedValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][3]
	    selectedParents_Geno_Par_List[[count]]  <-  selectedParents_Geno_List[[nCyc]][[nPar1]][[nPar2]]
	 
	    count <- count+1
	  }
    }
 
	
	### 
	
    Quartile3 <- summary(GainValues)[5]
	Quartile1 <- summary(GainValues)[2]
	MedianG <- 	summary(GainValues)[3]	
		
	indices1 <- which(GainValues >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues <= MedianG) & (GainValues >= Quartile1))
	indices2B <- which((GainValues >= MedianG) & (GainValues <= Quartile3))
	indices2CG <- (which(GainValues <= MedianG+1 & GainValues >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
###########
				
	for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	   
			   selectedParents_Geno_List[[nCyc]][[nParameter]] <-  selectedParents_Geno_Par_List[[nParameterIndex]]
			   CriterionValue_List[[nCyc]][[nParameter]] <- c(GainValues[nParameterIndex],UsefulnessValues[nParameterIndex],InbreedValues[nParameterIndex])
	} 
   
  #  save.image(WorkspaceName)
    gc() 	
  
#####################################################################################################################
    for(nParameter in 1:nParameters){
	
	    genoValues <- genoValues_List[[nParameter]] 
        phenoValues <- phenoValues_List[[nParameter]] 
        genoValSimValues <- genoValSimValues_List[[nParameter]] 
	    phenoValSimValues <- phenoValSimValues_List[[nParameter]] 
	       
    	uniqueParentIDs <-  unique(as.vector(selectedParents_Geno_List[[nCyc]][[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs

### GM Output in selected individual ids and Num Progeny from selected parents
		
	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-   genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(genoValSimValues[GenoTableIndices_topN])
		
		topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
		topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- PhenoTableIndices_topN
		PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 

#########################
	
		
		no_selected <- numberSelected
		nCrosses <- noCrosses
		nProgeny <- no_Progeny
		cycle1_nProgeny <- nProgeny
		nFamilies <- 20
		nLines <- nFamilies*cycle1_nProgeny 
		nProgeny_per_Cross <- nLines/nCrosses			
	
				
		Cycle_GenoData  <- (array(0,c(1,nMarkers,2,nCrosses*nProgeny_per_Cross)))
		selectedCycle_GenoData  <- (array(0,c(1,nMarkers,2,no_selected)))	
		nCrsInit <- 1
		nCrsFinal <- nProgeny_per_Cross
		Cycle_Progeny_F5 <- Cycle_Progeny_F5_List[[nParameter]]
   
		for(nCrs in 1:nCrosses){

		  Cycle_GenoData[1,,,nCrsInit:nCrsFinal] <- Cycle_Progeny_F5[[nCrs]]
		 
		  nCrsInit <- nCrsFinal+1
		  nCrsFinal <-  nCrsFinal+ nProgeny_per_Cross 
		
		}
		
		if(selectionOnGeno == TRUE){
		   selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedGenoIndividualIndices]
		}else if(selectionOnGeno == FALSE){ 
		   selectedCycle_GenoData[1,,,] <- Cycle_GenoData[1,,,selectedPhenoIndividualIndices]
		}
			
			F5RILs_GM <- getF5RILs_FS_GM_Frontier(selectedCycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New,nProgeny_per_Cross)
		
			Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
			nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
		
	}
	
	#	save.image(WorkspaceName) 
	
	gc() 
		
    simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)
	
}

######
       Input_Data_List <- list(Cycle_Progeny_F5_List,nCrosses_List,nProgeny_per_Cross,PredictionModel_Geno_List,PredictionModel_Pheno_List,trainGeno_backingFileName_List,trainSimPhenoTablePreCycle_List,trainSimGenoValTablePreCycle_List,varE,CombinedTableCount)
	   
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,
	   attainedPhenoValues_List_Cycle,topNGenoValues_List_Cycle,GenoTableIndices_topN_List_Cycle,topNPhenoValues_List_Cycle,PhenoTableIndices_topN_List_Cycle,selectedParents_Geno_List,CriterionValue_List)	

	   return(list(Input_Data_List,simResults_List))

}

############################
######

 runSimulations20X_BD_Wght_GM_Frontier_V3A_Complete <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach){

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
	  
	   
	  #trainGenoNewTablePreCycle <- as.big.matrix(trainGenoNewTable)
	 
 	  trainGeno_backingFileName <- paste("trainTable_Cyc1",condition,"_",Rep,".bin",sep="") 
      trainGeno_descriptorFileName <- paste("trainTable_Cyc1",condition,"_",Rep,".desc",sep="")
		 			
	  trainGenoNewTablePreCycle <- deepcopy(as.big.matrix(trainGenoNewTable),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
	  
	  
	  
      trainSimPhenoTablePreCycle <- trainSimPhenoValTable
      trainSimGenoValTablePreCycle <- trainSimGenoValTable

	  trainTableWindowSize <- TrainTableWindowSize
	  
	  
## Weighting Parameters	  
	  Weighted <- weighted
	  WghtMethod <- wghtMethod

	  
## GM Parameters 
	  
	  GM_Param_List <- GM_Frontier_Param_List
	  no_cores <- GM_Param_List$no_cores
	  startCycle <- StartCycle
	  WorkspaceName <- WorkSpaceName
	  nCores_ForEach <- noCores_ForEach 
	  
	  save.image(WorkspaceName)
#################################################################################################
#################################################################################################
### Initialize variables	
	
   	  j<-1
 	  CombinedTableCount <- 1
################################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
################################################################################################	  	  
#### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	
## genotype table for cycle 1
	
	cycle1GenoTable_GM <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable_GM <- generateWeightedGenoTable(cycle1GenoTable_GM,no_QTL)
	cycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(cycle1GenoTable_GM)
	  
	genoValSimValues <- simulateGenotypicValues(cycle1GenoTable_GM,no_QTL) 
	phenoValSimValues <- simulatePhenotypicValues(cycle1GenoTable_GM,varE,no_QTL) 

	nIndividuals <- nrow(cycle1GenoTable_GM)
	
## Write Genotable DF 	
	cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 
	
	
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

 
	  selectedGenoIndividualIndices2 <- rep(0,no_selected) 
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)
	 

# #################################################################################
    nCrosses <- dim(cycle1_Geno_Data)[1]
    nProgeny <- dim(cycle1_Geno_Data)[4]
    no_selected <- 200

	
    Cycle_GenoData  <- (array(0,c(1,nMarkers,2,cycle1_nProgeny*nCrosses)))
   
    initIndex <- 1
    finalIndex <- nProgeny
   
    for(nCrs in 1:nCrosses){

	
	 Cycle_GenoData[1,,,initIndex:finalIndex] <- cycle1_Geno_Data[nCrs,,,] 
	 initIndex <- initIndex + nProgeny
	 finalIndex <- finalIndex + nProgeny
     }
	
	
	

#### GM method for parent selection in cycle 1 #####################################
 nParameters <- 30 
 selectedParents_Geno_List <- rep(list(rep(list(list()),nParameters)),nCycles)
 CriterionValue_List <- rep(list(rep(list(list()),nParameters)),nCycles) 

 
 genoSelected_List <- getGM_Frontier_FS_selectedGenoList_Cyc1(cycle1GenoTable_GM,cycle1GenoTableMod_GM,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List)

 selectedParents_Geno_List[[nCyc]] <- genoSelected_List[[1]]
 CriterionValue_List[[nCyc]] <- genoSelected_List[[2]]
 selectedParents_Geno <- genoSelected_List[[1]]
 nParameters <- genoSelected_List[[3]]
 


save(genoSelected_List,file=WorkspaceName)

 #save.image(WorkspaceName) 	
 gc()
 
################################################################################### 
## Initiate Output list 	  
	  
	  GenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  GenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	 
	  PhenoVal_NX_N_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  PhenoVal_Sim_NX_2k_3c_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  
      attainedPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      attainedGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	
	  topNGenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      GenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)

      topNPhenoValues_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
      PhenoTableIndices_topN_List_Cycle <- rep(list(rep(list(list()),nParameters)),nCycles)
	  
	  Cycle_Progeny_F5_List <- list()
      nCrosses_List <- list()
	  
	  
	 
##
 
   for(nParameter in 1:nParameters){
   
        uniqueParentIDs <-  unique(as.vector(selectedParents_Geno[[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
### GM Output in selected individual ids 

############################

	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]]<- genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <- max(genoValSimValues[GenoTableIndices_topN])
		
	    topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
	   
	    PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
	    topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <-  PhenoTableIndices_topN
	    PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		
	    
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
###################################################################################	
## Get F5 RILs for selected parents in cycle 1
		nProgeny_per_Cross <- 1
      
        F5RILs_GM <- getF5RILs_FS_GM_Frontier_Complete(Cycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New,nProgeny_per_Cross)

        Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
		
	}

	simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)
	
	save(simResults_List_Cycle,file=WorkspaceName)
	
    startCycle <- startCycle + 1 	 
	 
	gc() 
	
###########################################################################

    startCycle <- 2
    nCores <- nCores_ForEach 
    #registerDoParallel(nCores)

    nextGenGenoTable_GM_List <- list()
    nextGenGenoTableMod_GM_List <- list()
    newNextGenGenoTable_GM_List <- list()
	PredictionModel_Geno_List <- list()
	PredictionModel_Pheno_List <- list()
	genoValues_List <- list()
    phenoValues_List <- list()
	genoValSimValues_List <- list()
	phenoValSimValues_List <- list()
	
	PredictionModel_Geno_List <- list()
	PredictionModel_Pheno_List <- list()
		   
	trainGenoNewTablePreCycle_List <- list()
	trainGenoNewTablePreCycle_Mat_List <- list()
	
	trainSimPhenoTablePreCycle_List <- list()
	trainSimGenoValTablePreCycle_List <- list()
	
	trainGeno_backingFileName_List <- list()
    trainGeno_descriptorFileName_List <- list()	
     
for(nCyc in startCycle:nCycles){
      
    print(nCyc)

### Generate F5 RIL progeny from geno data for selected geno data
 
    for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <-  Cycle_Progeny_F5_List[[nParameter]]
		nCrosses <-  nCrosses_List[[nParameter]] 
		
		nextGenGenoTable_GM <- generateMinus101GenoFormat_GM_Frontier_V1_Complete(Cycle_Progeny_F5,nCrosses)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)
		
		nextGenGenoTable_GM_List[[nParameter]] <- nextGenGenoTable_GM
        nextGenGenoTableMod_GM_List[[nParameter]] <- nextGenGenoTableMod_GM
		newNextGenGenoTable_GM_List[[nParameter]] <- newNextGenGenoTable_GM 
		genoValSimValues_List[[nParameter]] <- genoValSimValues 
		phenoValSimValues_List[[nParameter]] <- phenoValSimValues 
		
## Write genotable df	
		
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
		
	}
	
############################################################################################### 
		if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
		
		  if(nCyc ==2){
		   
		    trainGenoNewTablePreCycle <- attach.big.matrix(trainGeno_descriptorFileName)
   		    trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)

			PredictionModels <- foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_List[[j]],genoValSimValues_List[[j]],phenoValSimValues_List[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,trainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount,nCyc))
		  }else if(nCyc >2){
		    # trainGenoNewTable_List[[nrep]] <- attach.big.matrix(trainGeno_descriptorFileName)
		    # trainGenoNewTable_Mat_List[[nrep]] <- bigmemory::as.matrix(trainGenoNewTable_List[[nrep]])
		    
			for(nParameter in 1:nParameters){
			  trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		      trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
			}
		    PredictionModels <- foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_List[[j]],genoValSimValues_List[[j]],phenoValSimValues_List[[j]],PredictionModel_Geno_List[[j]],PredictionModel_Pheno_List[[j]],modelType,updateType,trainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount,nCyc)) 
 
		  }
		}
		
		for(nParameter in 1:nParameters){
   
			PredictionModel_Geno_List[[nParameter]] <- PredictionModels[[nParameter]][[1]]
			PredictionModel_Pheno_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
			CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		   
		    PredictionModel_Geno <- PredictionModel_Geno_List[[nParameter]] 
		    PredictionModel_Pheno <- PredictionModel_Pheno_List[[nParameter]]  
			
			trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="") 
            trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
			
			trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
			trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		
			
			trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
			trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
			trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		
		    nextGenGenoTable_GM <- nextGenGenoTable_GM_List[[nParameter]] 
            nextGenGenoTableMod_GM <- nextGenGenoTableMod_GM_List[[nParameter]] 
		    newNextGenGenoTable_GM <- newNextGenGenoTable_GM_List[[nParameter]]
			
			Freq <- getFreq(nextGenGenoTable_GM,FavAllele)
		  
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

			genoValues_List[[nParameter]] <- genoValues
			phenoValues_List[[nParameter]] <- phenoValues
		
			genoValSimValues <- genoValSimValues_List[[nParameter]]  
			phenoValSimValues <- phenoValSimValues_List[[nParameter]]   


		}

#### GM method for parent selection in cycle nCyc #####################################
 
       
		genoSelected_List <- foreach(k=1:nParameters) %dopar% (LTSGenoIMV2::getGM_Frontier_FS_selectedGenoList_Cyc_2_N(nextGenGenoTable_GM_List[[k]],nextGenGenoTableMod_GM_List[[k]],PredictionModel_Geno_List[[k]],PredictionModel_Pheno_List[[k]],modelType,selectionOnGeno,no_cores,nFamilies,GM_Param_List))

		save(genoSelected_List,simResults_List_Cycle,file=WorkSpaceName)
		rm(simResults_List_Cycle)
## Reduce Parameter set to 30  
   
    for(nParameter in 1:nParameters){
    
      selectedParents_Geno_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[1]]
      CriterionValue_List[[nCyc]][[nParameter]] <- genoSelected_List[[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nParameter]][[1]]
          
    }  
 
    count <- 1
	GainValues <- c()
	UsefulnessValues <- c()
	InbreedValues <- c()
	selectedParents_Geno_Par_List <- list()
	
    for (nPar1 in 1: nParameters){ 
	
	  nParameters2 <- length(selectedParents_Geno_List[[nCyc]][[nPar1]])
 
      for(nPar2 in 1:nParameters2){ 
	    GainValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][1]
	    UsefulnessValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][2]
	    InbreedValues[count] <- CriterionValue_List[[nCyc]][[nPar1]][[nPar2]][3]
	    selectedParents_Geno_Par_List[[count]]  <-  selectedParents_Geno_List[[nCyc]][[nPar1]][[nPar2]]
	 
	    count <- count+1
	  }
    }
 
	
	### 
	
    Quartile3 <- summary(GainValues)[5]
	Quartile1 <- summary(GainValues)[2]
	MedianG <- 	summary(GainValues)[3]	
		
	indices1 <- which(GainValues >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues <= MedianG) & (GainValues >= Quartile1))
	indices2B <- which((GainValues >= MedianG) & (GainValues <= Quartile3))
	indices2CG <- (which(GainValues <= MedianG+1 & GainValues >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
###########
				
	for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	   
			   selectedParents_Geno_List[[nCyc]][[nParameter]] <-  selectedParents_Geno_Par_List[[nParameterIndex]]
			   CriterionValue_List[[nCyc]][[nParameter]] <- c(GainValues[nParameterIndex],UsefulnessValues[nParameterIndex],InbreedValues[nParameterIndex])
	} 
   
  #  save.image(WorkspaceName)
    gc() 	
  
#####################################################################################################################
    for(nParameter in 1:nParameters){
	
	    genoValues <- genoValues_List[[nParameter]] 
        phenoValues <- phenoValues_List[[nParameter]] 
        genoValSimValues <- genoValSimValues_List[[nParameter]] 
	    phenoValSimValues <- phenoValSimValues_List[[nParameter]] 
	       
    	uniqueParentIDs <-  unique(as.vector(selectedParents_Geno_List[[nCyc]][[nParameter]]))
   
 	    topNGenoSimValues <- genoValSimValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs 
			   
	    topNPhenoSimValues <- phenoValSimValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs
		
	    topNGenoValues <- genoValues[uniqueParentIDs]
	    GenoTableIndices_topN <- uniqueParentIDs
			   
	    topNPhenoValues <- phenoValues[uniqueParentIDs]
	    PhenoTableIndices_topN <- uniqueParentIDs

### GM Output in selected individual ids and Num Progeny from selected parents
		
	    GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-  genoValSimValues
	    GenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <-   genoValues
	    attainedGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(genoValSimValues[GenoTableIndices_topN])
		
		topNGenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNGenoValues
	    GenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- GenoTableIndices_topN
		GenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNGenoValues_List_Cycle[[nCyc]][[nParameter]]
		
		PhenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValSimValues
	    PhenoVal_NX_2k_3c_List_Cycle[[nCyc]][[nParameter]] <- phenoValues
	    attainedPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  max(phenoValSimValues[PhenoTableIndices_topN])
		
		topNPhenoValues_List_Cycle[[nCyc]][[nParameter]] <-  topNPhenoValues
	    PhenoTableIndices_topN_List_Cycle[[nCyc]][[nParameter]] <- PhenoTableIndices_topN
		PhenoVal_NX_N_3c_List_Cycle[[nCyc]][[nParameter]] <- topNPhenoValues_List_Cycle[[nCyc]][[nParameter]]
		
	    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 

#########################
	
		Cycle_Progeny_F5 <- Cycle_Progeny_F5_List[[nParameter]]
	
		nCrosses <- length(Cycle_Progeny_F5)
		nProgeny_per_Cross <- 1
		Cycle_GenoData  <- (array(0,c(1,nMarkers,2,nCrosses*nProgeny_per_Cross)))
		
		for(nCrs in 1:nCrosses){

			Cycle_GenoData[1,,,1] <- Cycle_Progeny_F5[[nCrs]]
	 		
		}
	
	     	
        F5RILs_GM <- getF5RILs_FS_GM_Frontier_Complete(Cycle_GenoData,selectedParents_Geno_List[[nCyc]][[nParameter]],selectionOnGeno,nMarkers,NAM_LinkMap_New,nProgeny_per_Cross)
	
	    Cycle_Progeny_F5_List[[nParameter]] <- F5RILs_GM[[1]]
        nCrosses_List[[nParameter]] <- F5RILs_GM[[2]]
	
	}
	

	gc() 
		
    simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

}

######
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	   return(simResults_List)

}

###### 

###### Errors : Fail to load packages on multiple nodes 
##### Don't: do not use '.packages' in foreach. Only use clusterEvalQ to load packages and donot reload... 
###### Errors : Cannot find function present in loaded packages
##### Dos : Refer package namespace for functions in foreach even after loading packages on clusterEvalQ


#Initial_Input_Variable_List <- list(noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainTableWindowSize,nFamilies,GM_Param_List,noCores_ForEach,WorkspaceName)

# Input_Data_List <- list(Cycle_Progeny_F5_List,nCrosses_List,nProgeny_per_Cross,PredictionModel_Geno_List,PredictionModel_Pheno_List,trainGeno_backingFileName_List,trainSimPhenoTablePreCycle_List,trainSimGenoValTablePreCycle_List,varE,CombinedTableCount)
   
# simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,attainedPhenoValues_List_Cycle,topNGenoValues_List_Cycle,GenoTableIndices_topN_List_Cycle,topNPhenoValues_List_Cycle,PhenoTableIndices_topN_List_Cycle,selectedParents_Geno_List,CriterionValue_List)	
