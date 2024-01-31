### Simulation Functions for Island selection with GM Frontier methods





runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3A_Old <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach){

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
	  
#################################################################################################	  	  	 
#### Get predicted geno and pheno values for cycle1 ###################################################
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

	  CombinedTableCount <- 1
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
################################################################################################
#### Get predicted geno and pheno values for cycle1 ###################################################
    
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

#### Split Geno and Pheno values according to families  ########################

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
    Cycle_Progeny_F5_List <- cycle1_GenoData_List

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
	
			
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
	
## Write Genotable DF 	
	cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 
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
  
  trainGenoNewTablePreCycle_List <- list()
  trainGenoNewTablePreCycle_Mat_List <- list()
  
  trainSimPhenoTablePreCycle_List <- list()
  trainSimGenoValTablePreCycle_List <- list()
	
  
  trainGeno_backingFileName_List <- list()
  trainGeno_descriptorFileName_List <- list()
 
    
 for(nCyc in startCycle:nCycles){
   
    #nCyc <- 2
    print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
    for(nFamily in 1:nFamilies){
    
	  for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]]
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
		
		
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
		
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
		
	}
   }

##	
	
	PredictionModel_Geno_Parameter_List <- list() 
	PredictionModel_Pheno_Parameter_List <- list()
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	nextGenGenoTable_GM_Fam_Par <- list()
	nextGenGenoTableMod_GM_Fam_Par <- list()
	newNextGenGenoTable_GM_Fam_Par <- list()
	genoValSimValues_Fam_Par <- list()
	phenoValSimValues_Fam_Par <- list()
	
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
    		
			nextGenGenoTable_GM_Fam_Par[[nParameter]] <- nextGenGenoTable_GM_Fam
			nextGenGenoTableMod_GM_Fam_Par[[nParameter]] <- nextGenGenoTableMod_GM_Fam
			newNextGenGenoTable_GM_Fam_Par[[nParameter]] <- newNextGenGenoTable_GM_Fam
			genoValSimValues_Fam_Par[[nParameter]] <- genoValSimValues_Fam
			phenoValSimValues_Fam_Par[[nParameter]] <- phenoValSimValues_Fam
			
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
		
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
	}	
############################################################################################### 
    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	    if(nCyc ==2){
		  
          trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		
		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount))
 		}else if(nCyc >2){ 
		   for(nParameter in 1:nParameters){
			  trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		          trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
			}
		 
 		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount))
		}		    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
		  trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="")
		  trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
			
		  trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
		  trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		
					
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
		 
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
                  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
	}
	 
     if(modelUpdate ==FALSE){ 
	    for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		  
		  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]
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
		 
			output_List_afterExchange <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::exchangeGenoData_GM(Cycle_Progeny_F5_Fam_CrossList[[k]],genoValSimValues_List_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
		  }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		  }
		}
		  
## Get F5 RILs for selected parents in cycle 1

		Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily
		familyInfo_List <- list()

		nProg <- nProgInit 
		nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	      for(nFamily in 1:length(families)){ 

          		Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
				  
			    familyInfo_List[[nParameter]] <- output_List_afterExchange[[nParameter]][[2]]
					
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
    
   ### 
	
    Quartile3 <- summary(GainValues[[nFamily]])[5]
	Quartile1 <- summary(GainValues[[nFamily]])[2]
	MedianG <- 	summary(GainValues[[nFamily]])[3]	
		
	indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues[[nFamily]] <= MedianG) & (GainValues[[nFamily]] >= Quartile1))
	indices2B <- which((GainValues[[nFamily]] >= MedianG) & (GainValues[[nFamily]] <= Quartile3))
	indices2CG <- (which(GainValues[[nFamily]] <= MedianG+1 & GainValues[[nFamily]] >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
######

    			
    for(nParameter in 1:nParameters){
               nParameterIndex <- indicesRange[nParameter]
	           selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <-  selectedParents_Geno_Par_List[[nFamily]][[nParameterIndex]]
			   CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- c(GainValues[[nFamily]][nParameterIndex],UsefulnessValues[[nFamily]][nParameterIndex],InbreedValues[[nFamily]][nParameterIndex])
		
	} 
   }

  save.image(WorkspaceName)
  gc() 	
 
######################################################################################################## 

    Cycle_GenoData_List <- Cycle_Progeny_F5_Fam_CrossList
		
       
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% (getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]][[j]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))

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

	print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]],mean))),sep=""))

		 
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
 
	simResults_List_Cycle <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)
	
	
	save(genoSelected_List,simResults_List_Cycle,file=WorkspaceName)

    rm(simResults_List_Cycle)
}

######
	simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	return(simResults_List)

}





runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3A <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach){

    
######## Get predicted geno and pheno values for cycle1 ###################################################
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

	  CombinedTableCount <- 1
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
# ###############################################################################################	  	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
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
    Cycle_Progeny_F5_List<- cycle1_GenoData_List

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
	
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
	
## Write Genotable DF 	
	cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 

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
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
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
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
 
  PredictionModel_Geno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Pheno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Geno_Parameter_List <- list() 
  PredictionModel_Pheno_Parameter_List <- list()
  
  genoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  genoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  
  trainGenoNewTablePreCycle_List <- list()
  trainGenoNewTablePreCycle_Mat_List <- list()
  
  trainSimPhenoTablePreCycle_List <- list()
  trainSimGenoValTablePreCycle_List <- list()
	
  trainGeno_backingFileName_List <- list()
  trainGeno_descriptorFileName_List <- list()
 
    
 for(nCyc in startCycle:nCycles){
   
    print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
    for(nFamily in 1:nFamilies){
    
	  for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]]
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
	
	
	
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	nextGenGenoTable_GM_Fam_Par <- list()
	nextGenGenoTableMod_GM_Fam_Par <- list()
	newNextGenGenoTable_GM_Fam_Par <- list()
	genoValSimValues_Fam_Par <- list()
	phenoValSimValues_Fam_Par <- list()
	
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
    		
			nextGenGenoTable_GM_Fam_Par[[nParameter]] <- nextGenGenoTable_GM_Fam
			nextGenGenoTableMod_GM_Fam_Par[[nParameter]] <- nextGenGenoTableMod_GM_Fam
			newNextGenGenoTable_GM_Fam_Par[[nParameter]] <- newNextGenGenoTable_GM_Fam
			genoValSimValues_Fam_Par[[nParameter]] <- genoValSimValues_Fam
			phenoValSimValues_Fam_Par[[nParameter]] <- phenoValSimValues_Fam

#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
		
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
	}
	
############################################################################################### 
    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	    if(nCyc ==2){
		  trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		
		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount))
 		}else if(nCyc >2){ 
		   for(nParameter in 1:nParameters){
			  trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		      trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
			}
		 
 		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount))
		}		    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
		  trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="")
		  trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
			
		  trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
		  trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		
					
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
		 
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
          newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
	}
	
	
	if(modelUpdate ==FALSE){ 
	    for(nParameter in 1:nParameters){
		
		  # PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModel_Geno
		  # PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModel_Pheno
		  
		  # PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  # PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		  
		  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]
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
		 
			output_List_afterExchange <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::exchangeGenoData_GM(Cycle_Progeny_F5_Fam_CrossList[[k]],genoValSimValues_List_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
		  }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		 }
		  
## Get F5 RILs for selected parents in cycle 1

		Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily

		nProg <- nProgInit 
		nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	      for(nFamily in 1:length(families)){ 

           nProg <- nProgInit 
	       nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					 
		        # Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nParameter]][[nFamily]][[1]][[1]][nCross,,nProgInit:nProgFinal] 
				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
				
		    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
			
			  familyInfo <- output_List_afterExchange[[nParameter]][[2]]
					
			
	      }


		  
	    }
	
		# Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
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
    
   ### 
	
    Quartile3 <- summary(GainValues[[nFamily]])[5]
	Quartile1 <- summary(GainValues[[nFamily]])[2]
	MedianG <- 	summary(GainValues[[nFamily]])[3]	
		
	indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues[[nFamily]] <= MedianG) & (GainValues[[nFamily]] >= Quartile1))
	indices2B <- which((GainValues[[nFamily]] >= MedianG) & (GainValues[[nFamily]] <= Quartile3))
	indices2CG <- (which(GainValues[[nFamily]] <= MedianG+1 & GainValues[[nFamily]] >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
######

    			
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
	nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
    nCrosses_inFamily <- nSel_inFamily

    nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
     for(nParameter in 1:nParameters){
	 		
	   for(nFamily in 1:length(families)){ 

           nProg <- nProgInit 
	       nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					 
		       				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
				
		    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
	    } 
	}

######################################################################################################## 

    Cycle_GenoData_List <- Cycle_Progeny_F5_Fam_CrossList
		
 
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% (getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]][[j]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))

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

	print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]],mean))),sep=""))

		 
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
	
    save.image(WorkspaceName)
	
 }
 
######
	simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	return(simResults_List)

}

####################



runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3A2 <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach){

    
######## Get predicted geno and pheno values for cycle1 ###################################################
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

	  CombinedTableCount <- 1
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
# ###############################################################################################	  	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
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
    Cycle_Progeny_F5_List<- cycle1_GenoData_List

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
	
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
	
## Write Genotable DF 	
	cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 

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
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
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
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
 
  PredictionModel_Geno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Pheno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Geno_Parameter_List <- list() 
  PredictionModel_Pheno_Parameter_List <- list()
  
  genoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  genoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  
  trainGenoNewTablePreCycle_List <- list()
  trainGenoNewTablePreCycle_Mat_List <- list()
  
  trainSimPhenoTablePreCycle_List <- list()
  trainSimGenoValTablePreCycle_List <- list()
	
  trainGeno_backingFileName_List <- list()
  trainGeno_descriptorFileName_List <- list()
 
    
 for(nCyc in startCycle:nCycles){
   
       print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
    
		for(nFamily in 1:nFamilies){
    
			for(nParameter in 1:nParameters){
	
				Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]]
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
	
	
	
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	nextGenGenoTable_GM_Fam_Par <- list()
	nextGenGenoTableMod_GM_Fam_Par <- list()
	newNextGenGenoTable_GM_Fam_Par <- list()
	genoValSimValues_Fam_Par <- list()
	phenoValSimValues_Fam_Par <- list()
	
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
    		
			nextGenGenoTable_GM_Fam_Par[[nParameter]] <- nextGenGenoTable_GM_Fam
			nextGenGenoTableMod_GM_Fam_Par[[nParameter]] <- nextGenGenoTableMod_GM_Fam
			newNextGenGenoTable_GM_Fam_Par[[nParameter]] <- newNextGenGenoTable_GM_Fam
			genoValSimValues_Fam_Par[[nParameter]] <- genoValSimValues_Fam
			phenoValSimValues_Fam_Par[[nParameter]] <- phenoValSimValues_Fam

#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
		
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
	}
	
############################################################################################### 
    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	    if(nCyc ==2){
		  trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		
		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount))
 		}else if(nCyc >2){ 
		   for(nParameter in 1:nParameters){
			  trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		      trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
			}
		 
 		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount))
		}		    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
		  trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="")
		  trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
			
		  trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
		  trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		
					
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
		 
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
          newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
	}
	
	
	if(modelUpdate ==FALSE){ 
	    for(nParameter in 1:nParameters){
		
		  # PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModel_Geno
		  # PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModel_Pheno
		  
		  # PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  # PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		  
		  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]
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
		 
			output_List_afterExchange <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::exchangeGenoData_GM(Cycle_Progeny_F5_Fam_CrossList[[k]],genoValSimValues_List_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
		  }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		 }
		  
## Get F5 RILs for selected parents in cycle 1

		Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily

		nProg <- nProgInit 
		nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	      for(nFamily in 1:length(families)){ 

           nProg <- nProgInit 
	       nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					 
		        # Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nParameter]][[nFamily]][[1]][[1]][nCross,,nProgInit:nProgFinal] 
				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
				
		    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
			
			  familyInfo <- output_List_afterExchange[[nParameter]][[2]]
					
			
	      }


		  
	    }
	
		# Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
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
    
   ### 
	
    Quartile3 <- summary(GainValues[[nFamily]])[5]
	Quartile1 <- summary(GainValues[[nFamily]])[2]
	MedianG <- 	summary(GainValues[[nFamily]])[3]	
		
	indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues[[nFamily]] <= MedianG) & (GainValues[[nFamily]] >= Quartile1))
	indices2B <- which((GainValues[[nFamily]] >= MedianG) & (GainValues[[nFamily]] <= Quartile3))
	indices2CG <- (which(GainValues[[nFamily]] <= MedianG+1 & GainValues[[nFamily]] >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
######

    			
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
	nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
    nCrosses_inFamily <- nSel_inFamily

    nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
     for(nParameter in 1:nParameters){
	 		
	   for(nFamily in 1:length(families)){ 

           nProg <- nProgInit 
	       nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					 
		       				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
				
		    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
	    } 
	}

######################################################################################################## 

    Cycle_GenoData_List <- Cycle_Progeny_F5_Fam_CrossList
		
 
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% (getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]][[j]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))

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

	print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]],mean))),sep=""))

		 
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
	
    save.image(WorkspaceName)
	
 }
 
######
	simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	return(simResults_List)

}

####################
############


runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3A_doMPI <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach){

    
######## Get predicted geno and pheno values for cycle1 ###################################################
     options(warn=-1)
	 library(Rmpi)
	 library(doMPI)
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

	  CombinedTableCount <- 1
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
# ###############################################################################################	  	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	#cl <- makeCluster(nCores_FE,outfile=paste("log_IM_GM",nCores_FE,sep=""))
	cl <- startMPIcluster(nCores_FE)
	registerDoMPI(cl)
	
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
    Cycle_Progeny_F5_List<- cycle1_GenoData_List

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
	
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
	
## Write Genotable DF 	
	cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 

######################################################################################   

#### GM method for parent selection in cycle 1 ################################################
# #### GM method for isolated families
    nParameters <- 30
 			 
 	selectedParents_Geno_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles) 
	CriterionValue_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

   exportDoMPI(cl,c('getIM_GM_Frontier_selectedGenoList_SingleFams_V2','cycle1GenoTable_GM_List','cycle1GenoTableMod_GM_List','PredictionModel_Geno','PredictionModel_Pheno','modelType','selectionOnGeno','no_cores','nFamilies','GM_Frontier_Param_List'))
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
   
    exportDoMPI(cl,c('getF5RILs_FS_IM_GM_Frontier','Cycle_GenoData_List','selectedParents_Geno_List','selectionOnGeno','nMarkers','NAM_LinkMap_New'))
    
	F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% (getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
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
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
 
  PredictionModel_Geno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Pheno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Geno_Parameter_List <- list() 
  PredictionModel_Pheno_Parameter_List <- list()
  
  genoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  genoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  
  trainGenoNewTablePreCycle_List <- list()
  trainGenoNewTablePreCycle_Mat_List <- list()
  
  trainSimPhenoTablePreCycle_List <- list()
  trainSimGenoValTablePreCycle_List <- list()
	
  trainGeno_backingFileName_List <- list()
  trainGeno_descriptorFileName_List <- list()
 
    
 for(nCyc in startCycle:nCycles){
   
    print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
 
    for(nFamily in 1:nFamilies){
    
	  for(nParameter in 1:nParameters){
	
	    Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]]
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
	
	
	
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	nextGenGenoTable_GM_Fam_Par <- list()
	nextGenGenoTableMod_GM_Fam_Par <- list()
	newNextGenGenoTable_GM_Fam_Par <- list()
	genoValSimValues_Fam_Par <- list()
	phenoValSimValues_Fam_Par <- list()
	
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
    		
			nextGenGenoTable_GM_Fam_Par[[nParameter]] <- nextGenGenoTable_GM_Fam
			nextGenGenoTableMod_GM_Fam_Par[[nParameter]] <- nextGenGenoTableMod_GM_Fam
			newNextGenGenoTable_GM_Fam_Par[[nParameter]] <- newNextGenGenoTable_GM_Fam
			genoValSimValues_Fam_Par[[nParameter]] <- genoValSimValues_Fam
			phenoValSimValues_Fam_Par[[nParameter]] <- phenoValSimValues_Fam

#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
		
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
	}
	
############################################################################################### 
    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	    if(nCyc ==2){
		  trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		  exportDoMPI(cl,c('getUpdatedGSModel','newNextGenGenoTable_GM_Fam_Par','genoValSimValues_Fam_Par','phenoValSimValues_Fam_Par','PredictionModel_Geno','PredictionModel_Pheno','modelType','updateType','TrainTableWindowSize','trainGenoNewTablePreCycle_Mat','trainSimPhenoTablePreCycle','trainSimGenoValTablePreCycle','CombinedTableCount'))
		 
		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount))
 		}else if(nCyc >2){ 
		   for(nParameter in 1:nParameters){
			  trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		      trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
			}
			exportDoMPI(cl,c('getUpdatedGSModel','newNextGenGenoTable_GM_Fam_Par','genoValSimValues_Fam_Par','phenoValSimValues_Fam_Par','PredictionModel_Geno_Parameter_List','PredictionModel_Pheno_Parameter_List','modelType','updateType','TrainTableWindowSize','trainGenoNewTablePreCycle_Mat_List','trainSimPhenoTablePreCycle_List','trainSimGenoValTablePreCycle_List','CombinedTableCount'))
		 
 		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount))
		}		    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
		  trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="")
		  trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
			
		  trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
		  trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		
					
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
		 
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
          newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
	}
	
	
	if(modelUpdate ==FALSE){ 
	    for(nParameter in 1:nParameters){
		
		  # PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModel_Geno
		  # PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModel_Pheno
		  
		  # PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  # PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		  
		  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]
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
		
			exportDoMPI(cl,c('exchangeGenoData_GM','Cycle_Progeny_F5_Fam_CrossList','genoValSimValues_List_Reverse','migrationSize','emigrantGroups','immigrantGroups','direction','Policy'))
		 
			output_List_afterExchange <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::exchangeGenoData_GM(Cycle_Progeny_F5_Fam_CrossList[[k]],genoValSimValues_List_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
		  }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		 }
		  
## Get F5 RILs for selected parents in cycle 1

		Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily

		nProg <- nProgInit 
		nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	      for(nFamily in 1:length(families)){ 

           nProg <- nProgInit 
	       nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					 
		        # Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nParameter]][[nFamily]][[1]][[1]][nCross,,nProgInit:nProgFinal] 
				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
				
		    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
			
			  familyInfo <- output_List_afterExchange[[nParameter]][[2]]
					
			
	      }


		  
	    }
	
		# Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
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
 
    exportDoMPI(cl,c('getIM_GM_Frontier_selectedGenoList_SingleFams_V2','nextGenGenoTable_GM_List','nextGenGenoTableMod_GM_List','PredictionModel_Geno','PredictionModel_Pheno','modelType','selectionOnGeno','no_cores','nFamilies','GM_Frontier_Param_List'))

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
    
   ### 
	
    Quartile3 <- summary(GainValues[[nFamily]])[5]
	Quartile1 <- summary(GainValues[[nFamily]])[2]
	MedianG <- 	summary(GainValues[[nFamily]])[3]	
		
	indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues[[nFamily]] <= MedianG) & (GainValues[[nFamily]] >= Quartile1))
	indices2B <- which((GainValues[[nFamily]] >= MedianG) & (GainValues[[nFamily]] <= Quartile3))
	indices2CG <- (which(GainValues[[nFamily]] <= MedianG+1 & GainValues[[nFamily]] >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
######

    			
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
	nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
    nCrosses_inFamily <- nSel_inFamily

    nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
     for(nParameter in 1:nParameters){
	 		
	   for(nFamily in 1:length(families)){ 

           nProg <- nProgInit 
	       nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					 
		       				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
				
		    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
	    } 
	}

######################################################################################################## 

    Cycle_GenoData_List <- Cycle_Progeny_F5_Fam_CrossList
		
    exportDoMPI(cl,c('getF5RILs_FS_IM_GM_Frontier','Cycle_GenoData_List','selectedParents_Geno_List','selectionOnGeno','nMarkers','NAM_LinkMap_New'))
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% (getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]][[j]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))

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

	print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]],mean))),sep=""))

		 
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
	
    save.image(WorkspaceName)
	
 }
 
######
	simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	return(simResults_List)

}

####################


runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3A_doSNOW <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach,GIndData){

    
######## Get predicted geno and pheno values for cycle1 ###################################################
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
 	
	  gIndData <- GIndData
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

	  CombinedTableCount <- 1
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
# ###############################################################################################	  	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)

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
## selected Line indices to select lines from cycle genotable

	
		genoSelectionTable <- ApplySelection_Discrete(SelCriterion_List,no_selected)
		genoSimSelectionTable <- genoSelectionTable
		selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
        selectedLineIndices <- selectedGenoIndividualIndices

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
    Cycle_Progeny_F5_List<- cycle1_GenoData_List

    if(migFreq !=0){
       
 	    if(nCyc==1 && nCyc%%migFreq==0){
	  
	    if(Policy != "GeneticClusters"){
		 
		 output_List_afterExchange <- exchangeGenoData_GM(Cycle_Progeny_F5_List,SelCriterion_List,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
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
	
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
	
## Write Genotable DF 	
	cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 

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
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
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
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
 
  PredictionModel_Geno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Pheno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Geno_Parameter_List <- list() 
  PredictionModel_Pheno_Parameter_List <- list()
  
  genoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  genoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  
  trainGenoNewTablePreCycle_List <- list()
  trainGenoNewTablePreCycle_Mat_List <- list()
  
  trainSimPhenoTablePreCycle_List <- list()
  trainSimGenoValTablePreCycle_List <- list()
	
  trainGeno_backingFileName_List <- list()
  trainGeno_descriptorFileName_List <- list()
 
    
 for(nCyc in startCycle:nCycles){
   
    print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
    
	
    nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
		genoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[4]]
		phenoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[5]]
	 }
	}
	
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	nextGenGenoTable_GM_Fam_Par <- list()
	nextGenGenoTableMod_GM_Fam_Par <- list()
	newNextGenGenoTable_GM_Fam_Par <- list()
	genoValSimValues_Fam_Par <- list()
	phenoValSimValues_Fam_Par <- list()
	
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
    		
			nextGenGenoTable_GM_Fam_Par[[nParameter]] <- nextGenGenoTable_GM_Fam
			nextGenGenoTableMod_GM_Fam_Par[[nParameter]] <- nextGenGenoTableMod_GM_Fam
			newNextGenGenoTable_GM_Fam_Par[[nParameter]] <- newNextGenGenoTable_GM_Fam
			genoValSimValues_Fam_Par[[nParameter]] <- genoValSimValues_Fam
			phenoValSimValues_Fam_Par[[nParameter]] <- phenoValSimValues_Fam

#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
	if(gIndData ==TRUE){
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
	}
  }
	
############################################################################################### 
    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	    if(nCyc ==2){
		  trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		
		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount))
 		}else if(nCyc >2){ 
		   for(nParameter in 1:nParameters){
			  trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		      trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
			}
		 
 		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount))
		}		    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
		  trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="")
		  trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
			
		  trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
		  trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		
					
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
		 
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
          newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
	}
	
	
	if(modelUpdate ==FALSE){ 
	    for(nParameter in 1:nParameters){
		
		  # PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModel_Geno
		  # PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModel_Pheno
		  
		  # PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  # PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		  
		  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]
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
	}

### Get Predicted Geno and Pheno Value List
 
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

#### Set Selection Criterion List

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


### Apply selection to cycle nCyc population 
### Genoselection table & Island migration for cycle nCyc

	
	genoSelectionTable <- ApplySelection_Discrete_Frontier(SelCriterion_List,no_selected,nFamilies,nParameters)
    selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) lapply(x,function(y) as.vector(y[,2]))) 
    selectedLineIndices <- selectedGenoIndividualIndices

######################################################################################   
## get correct table format for migration policy and exchange genodata functions

      Cycle_Progeny_F5_List_Fam_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
      nextGenGenoTable_GM_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters) 
      genoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      phenoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      SelCriterion_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)

      for(nFamily in 1:nFamilies){
    
	    for(nParameter in 1:nParameters){
	
	        Cycle_Progeny_F5 <-   Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]]
	 
	        Cycle_Progeny_F5_List_Fam_Reverse[[nParameter]][[nFamily]] <- Cycle_Progeny_F5
            SelCriterion_List_Reverse[[nParameter]][[nFamily]] <- SelCriterion_List[[nFamily]][[nParameter]]

	    }
      }
 
       
     for(nParameter in 1:nParameters){

	   for(nFamily in 1:nFamilies){

        	nextGenGenoTable_GM_List_Reverse[[nParameter]][[nFamily]] <- nextGenGenoTable_GM_List[[nFamily]][[nParameter]] 

            genoSelectionTable_Parameter[[nParameter]][[nFamily]] <- genoSelectionTable[[nFamily]][[nParameter]]

            }
     }

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

######################################################################################   
## get correct table format for exchange genodata functions


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
	

################

	
### Exchange selected geno data among pairs ##############
     
################################

    if(migFreq !=0){
   
 	  if(nCyc>=2 && nCyc%%migFreq==0){
	  
	         if(Policy != "GeneticClusters"){
		 
			output_List_afterExchange <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::exchangeGenoData_GM(Cycle_Progeny_F5_Fam_CrossList[[k]],SelCriterion_List_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
		 }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		 }
		  
###########################################################################
### Get CycleProgeny Genodata and genotype table after exchange

		Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily

		nProg <- nProgInit 
		nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	      for(nFamily in 1:length(families)){ 

           	nProg <- nProgInit 
	        nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					      				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
					    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
			
			  familyInfo <- output_List_afterExchange[[nParameter]][[2]]
	      }
	  
	    }
	
		# Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
	    Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 
	}   


### Get nextGenGenoTable_List after exchange
	
	nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
		genoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[4]]
		phenoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[5]]
	 }
	}
  
################################################################################################### GM method for isolated families

 system.time({

genoSelected_List <- foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% 
(getIM_GM_Frontier_selectedGenoList_SingleFams_V2(nextGenGenoTable_GM_List[[k]][[j]],nextGenGenoTableMod_GM_List[[k]][[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))

 })


## Reduce Parameter set to 30  

  #selectedParents_Geno_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
  #CriterionValue_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

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
    
   ### 
	
    Quartile3 <- summary(GainValues[[nFamily]])[5]
	Quartile1 <- summary(GainValues[[nFamily]])[2]
	MedianG <- 	summary(GainValues[[nFamily]])[3]	
		
	indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues[[nFamily]] <= MedianG) & (GainValues[[nFamily]] >= Quartile1))
	indices2B <- which((GainValues[[nFamily]] >= MedianG) & (GainValues[[nFamily]] <= Quartile3))
	indices2CG <- (which(GainValues[[nFamily]] <= MedianG+1 & GainValues[[nFamily]] >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
######

    			
    for(nParameter in 1:nParameters){
            nParameterIndex <- indicesRange[nParameter]
	        selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <-  selectedParents_Geno_Par_List[[nFamily]][[nParameterIndex]]
			CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- c(GainValues[[nFamily]][nParameterIndex],UsefulnessValues[[nFamily]][nParameterIndex],InbreedValues[[nFamily]][nParameterIndex])
		
	} 
   }

  save.image(WorkspaceName)
  gc() 	
   	
###########################################################################
## Get F5 RILs for selected parents in cycle nCyc

    Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
	nProgInit <- 1
	nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
    nCrosses_inFamily <- nSel_inFamily

    nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
     for(nParameter in 1:nParameters){
	 		
	   for(nFamily in 1:length(families)){ 

           nProg <- nProgInit 
	       nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					 
		       				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
				
		    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
	    } 
	}

######################################################################################################## 

    Cycle_GenoData_List <- Cycle_Progeny_F5_Fam_CrossList
		
    Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
    nCrosses_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
	
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% (getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]][[j]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))

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

print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]],mean))),sep=""))

		 
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
	
    save.image(WorkspaceName)
	
 }
 
######
	simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	return(simResults_List)

}

################ Selected pop 


runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3ASel_doSNOW <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach,GIndData){

    
######## Get predicted geno and pheno values for cycle1 ###################################################
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
 	
	  gIndData <- GIndData
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

  #if(startCycle ==1){
	
   	  Migration_List <- list()
	  PCA_List <- list()
	  j<-1
 	  selectedParents_Geno_List <- list() 
	  selectedParents_Pheno_List <- list()
	  selectedParents_NumProgeny_Geno_List <- list()
	  selectedParents_NumProgeny_Pheno_List <- list() 
	  
	  familyPairs_List <- list()
	  familyInd_List <- list()

	  CombinedTableCount <- 1
	  
	  no_selected_Fam <- 0.1 * cycle1_nProgeny
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
# ###############################################################################################	  	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	   
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
## selected Line indices to select lines from cycle genotable

	
		genoSelectionTable <- ApplySelection_Discrete(SelCriterion_List,no_selected)
		genoSimSelectionTable <- genoSelectionTable
		selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
        selectedLineIndices <- selectedGenoIndividualIndices

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
    Cycle_Progeny_F5_List<- cycle1_GenoData_List

    if(migFreq !=0){
       
 	   if(nCyc==1 && nCyc%%migFreq==0){
	  
	    if(Policy != "GeneticClusters"){
		 
		 output_List_afterExchange <- exchangeGenoData_GM_Sel(Cycle_Progeny_F5_List,SelCriterion_List,selectedLineIndices,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
		}else if(Policy == "GeneticClusters"){ 
		
		 
		 output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		}
		  
		Sel_Cycle_Progeny_F5_List <- (output_List_afterExchange)[[1]]
		
		familyInfo <- output_List_afterExchange[[2]]
		
		Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   }
	   
	    no_selected_Fam <- 10 
		
        Sel_Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,no_selected_Fam))
        nProgIndex <-  dim(Sel_Cycle_Progeny_F5_List[[nFamily]])[3]
	    nProg <-1
				
	    for(nFamily in 1:length(families)){ 
					 
		    Sel_Cycle_Progeny_F5[nFamily,,,nProg:nProgIndex] <- Sel_Cycle_Progeny_F5_List[[nFamily]][,,nProg:nProgIndex]
	    } 
		
	}
	
### selected genoTable data after exchange
	print(paste("Cyc-",nCyc,"before-genFn"))
    selCycle1GenoTable_GM <-  generateMinus101GenoFormat_V1(Sel_Cycle_Progeny_F5,nFamilies,no_selected_Fam)
	print(paste("Cyc-",nCyc,"after-genFn"))
    selNewCycle1GenoTable_GM <- generateWeightedGenoTable(selCycle1GenoTable_GM,no_QTL)
    selCycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(selCycle1GenoTable_GM)

## Get simulated geno and pheno values 	
 
    genoValSimValues <- simulateGenotypicValues(cycle1GenoTable_GM,no_QTL) 
    phenoValSimValues <- simulatePhenotypicValues(cycle1GenoTable_GM,varE,no_QTL) 

### i) Get cycle 1 genotype table and criterion value list for sets of singleislands and pairs of islands	
### Get list of cycle1 geno table for family ind and family pairs ######################################
	
	cycle1GenoTable_GM_List <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
	cycle1GenoTableMod_GM_List <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
    selCycle1GenoTable_GM_List <- rep(list(matrix(0,nrow=no_selected_Fam,ncol=nMarkers)),nFamilies)
	selCycle1GenoTableMod_GM_List <- rep(list(matrix(0,nrow=no_selected_Fam,ncol=nMarkers)),nFamilies)
	
	initIndex <- 1
	finalIndex <- cycle1_nProgeny_Fam
	
	initSelIndex <- 1
	finalSelIndex <- no_selected_Fam
	
	for(nFam in 1:nFamilies){ 
		
		cycle1GenoTable_GM_List[[nFam]] <- cycle1GenoTable_GM[initIndex:finalIndex,]
		cycle1GenoTableMod_GM_List[[nFam]] <- cycle1GenoTableMod_GM[initIndex:finalIndex,]
		
		selCycle1GenoTable_GM_List[[nFam]] <- selCycle1GenoTable_GM[initSelIndex:finalSelIndex,]
        selCycle1GenoTableMod_GM_List[[nFam]] <- selCycle1GenoTableMod_GM[initSelIndex:finalSelIndex,]
	    
		initIndex <- finalIndex+1
		finalIndex <- finalIndex + cycle1_nProgeny_Fam
		
	}
	
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
	
## Write Genotable DF 	
	#cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	#cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	#cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	#write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 

######################################################################################   

#### GM method for parent selection in cycle 1 ################################################
# #### GM method for isolated families
    nParameters <- 30
 			 
  selectedParents_Geno_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles) 
  CriterionValue_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  system.time({

    genoSelected_List <- foreach(k=1:nFamilies,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getIM_GM_Frontier_selectedGenoList_SingleFams_V2(selCycle1GenoTable_GM_List[[k]],selCycle1GenoTableMod_GM_List[[k]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))
  
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
    
    Sel_Cycle_GenoData_List <- Sel_Cycle_Progeny_F5_List

## Get F5 RILs for selected parents in cycle 1
   
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (getF5RILs_FS_IM_GM_Frontier_Sel(Sel_Cycle_GenoData_List[[k]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
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

#}

###########################################################################

  startCycle <- 2
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
 
  PredictionModel_Geno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Pheno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Geno_Parameter_List <- list() 
  PredictionModel_Pheno_Parameter_List <- list()
  
  genoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  genoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  
  trainGenoNewTablePreCycle_List <- list()
  trainGenoNewTablePreCycle_Mat_List <- list()
  
  trainSimPhenoTablePreCycle_List <- list()
  trainSimGenoValTablePreCycle_List <- list()
	
  trainGeno_backingFileName_List <- list()
  trainGeno_descriptorFileName_List <- list()
 
    
 for(nCyc in startCycle:nCycles){

    
    print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
    
	
    nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
		genoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[4]]
		phenoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[5]]
	 }
	}
	
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	nextGenGenoTable_GM_Fam_Par <- list()
	nextGenGenoTableMod_GM_Fam_Par <- list()
	newNextGenGenoTable_GM_Fam_Par <- list()
	genoValSimValues_Fam_Par <- list()
	phenoValSimValues_Fam_Par <- list()
	
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
    		
			nextGenGenoTable_GM_Fam_Par[[nParameter]] <- nextGenGenoTable_GM_Fam
			nextGenGenoTableMod_GM_Fam_Par[[nParameter]] <- nextGenGenoTableMod_GM_Fam
			newNextGenGenoTable_GM_Fam_Par[[nParameter]] <- newNextGenGenoTable_GM_Fam
			genoValSimValues_Fam_Par[[nParameter]] <- genoValSimValues_Fam
			phenoValSimValues_Fam_Par[[nParameter]] <- phenoValSimValues_Fam

#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
	if(gIndData ==TRUE){
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
	}
  }
	
############################################################################################### 
    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	    if(nCyc ==2){
		  trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		
		  PredictionModels <- foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount))
 		}else if(nCyc >2){ 
		   for(nParameter in 1:nParameters){
			  trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		      trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
			}
		 
 		  PredictionModels <- foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount))
		}		    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
		  trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="")
		  trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
			
		  trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
		  trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		
					
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
		 
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
          newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
	}
	
	
	if(modelUpdate ==FALSE){ 
	    for(nParameter in 1:nParameters){
		
		  # PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModel_Geno
		  # PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModel_Pheno
		  
		  # PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  # PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		  
		  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]
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
	}

### Get Predicted Geno and Pheno Value List
 
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

#### Set Selection Criterion List

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


### Apply selection to cycle nCyc population 
### Genoselection table & Island migration for cycle nCyc

	
	genoSelectionTable <- ApplySelection_Discrete_Frontier(SelCriterion_List,no_selected,nFamilies,nParameters)
    selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) lapply(x,function(y) as.vector(y[,2]))) 
    selectedLineIndices <- selectedGenoIndividualIndices

######################################################################################   
## get correct table format for migration policy and exchange genodata functions

    Cycle_Progeny_F5_List_Fam_V2  <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
    nCross_inFamily <- nCrosses/nFamilies

   for(nFamily in 1:nFamilies){ 
	for(nParameter in 1:nParameters){
            initIndex <- 1
    	    finalIndex <- 10 

            for(nCross in 1:nCross_inFamily){ 

		Cycle_Progeny_F5_List_Fam_V2[[nFamily]][[nParameter]][,,initIndex:finalIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,]
              
		initIndex <- finalIndex+1
		finalIndex <- finalIndex + no_selected_Fam
            }

     	}
     } 


      Cycle_Progeny_F5_List_Fam_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
      Sel_Cycle_Progeny_F5_List_Fam_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
      nextGenGenoTable_GM_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters) 
      genoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      phenoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      SelCriterion_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
	
	  selNextGenGenoTable_GM_List <- rep(list(rep(list(),nFamilies)),nParameters)
	  selNextGenGenoTableMod_GM_List <- rep(list(rep(list(),nFamilies)),nParameters)
          selNewNextGenGenoTable_GM_List <- rep(list(rep(list(),nFamilies)),nParameters)
	  
	  for(nFamily in 1:nFamilies){
    
	    for(nParameter in 1:nParameters){
	
	        Cycle_Progeny_F5 <-   Cycle_Progeny_F5_List_Fam_V2[[nFamily]][[nParameter]]
	        Cycle_Progeny_F5_List_Fam_Reverse[[nParameter]][[nFamily]] <- Cycle_Progeny_F5
			Sel_Cycle_Progeny_F5 <- Cycle_Progeny_F5[,,as.vector(selectedLineIndices[[nFamily]][[nParameter]])]
			Sel_Cycle_Progeny_F5_List_Fam_Reverse[[nParameter]][[nFamily]] <- Sel_Cycle_Progeny_F5
			SelCriterion_List_Reverse[[nParameter]][[nFamily]] <- SelCriterion_List[[nFamily]][[nParameter]]

	    }
     }
 
       
     for(nParameter in 1:nParameters){

	   for(nFamily in 1:nFamilies){

        	nextGenGenoTable_GM_List_Reverse[[nParameter]][[nFamily]] <- nextGenGenoTable_GM_List[[nFamily]][[nParameter]] 
			
            genoSelectionTable_Parameter[[nParameter]][[nFamily]] <- genoSelectionTable[[nFamily]][[nParameter]]
			
	        selNextGenGenoTable_GM_List[[nParameter]][[nFamily]] <- nextGenGenoTable_GM_List[[nFamily]][[nParameter]][(selectedLineIndices[[nFamily]][[nParameter]]),]
			
			selNextGenGenoTableMod_GM_List[[nParameter]][[nFamily]] <- nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]][(selectedLineIndices[[nFamily]][[nParameter]]),] 

        }
     }

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

######################################################################################   
## get correct table format for exchange genodata functions


        Cycle_Progeny_F5_Fam_ParamList <- rep(list(array(0,c(nFamilies,nMarkers,2,cycle1_nProgeny))),nParameters)

        Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)),nParameters)
		
	Sel_Cycle_Progeny_F5_Fam_ParamList <- rep(list(array(0,c(nFamilies,nMarkers,2,no_selected_Fam))),nParameters)

        Sel_Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,no_selected_Fam))),nFamilies)),nParameters)
        
	nProgInit <- 1
	nProgFinal <-  dim(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[1]])[4]
        nCrosses_inFamily <- nSel_inFamily

        nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	     for(nFamily in 1:length(families)){ 

              nProg <- nProgInit 
	      nProgIndex <- nProgFinal

               for(nCross in 1:nCross_inFamily){
					 
		Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,nProgInit:nProgFinal] 
				
                Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]][,,nProg:nProgIndex]<- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,nProgInit:nProgFinal] 
               
                			
		    
		nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal


	       }
	     } 
	}
	

        nProg <- 1 
	nProgIndex <- no_selected_Fam

        selectedLineIndices_Reverse <- rep(list(rep(list(c()),nFamilies)),nParameters)
	
        for(nParameter in 1:nParameters){
	 		
	     for(nFamily in 1:length(families)){ 

             
                Sel_Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,as.vector(selectedLineIndices[[nFamily]][[nParameter]])]
				
				Sel_Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]][,,nProg:nProgIndex] <-  Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]][,,as.vector(selectedLineIndices[[nFamily]][[nParameter]])]

               selectedLineIndices_Reverse[[nParameter]][[nFamily]] <-  selectedLineIndices[[nFamily]][[nParameter]]
          
	  }
	}
 

################
### Exchange selected geno data among pairs ##############
     
################################

    if(migFreq !=0){
   
 	  if(nCyc>=2 && nCyc%%migFreq==0){
	  
	     if(Policy != "GeneticClusters"){
		 
			output_List_afterExchange <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::exchangeGenoData_GM_Sel(Cycle_Progeny_F5_Fam_CrossList[[k]],SelCriterion_List_Reverse[[k]],selectedLineIndices_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
	      }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
	      }
		  
###########################################################################
### Get CycleProgeny Genodata and genotype table after exchange

		Sel_Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,no_selected_Fam))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily

		nProg <- nProgInit 
		nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	      for(nFamily in 1:length(families)){ 

           	nProg <- nProgInit 
	        nProgIndex <- nProgFinal

                for(nCross in 1:nCrosses_inFamily){
					      				
				Sel_Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
					    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
			
	        familyInfo <- output_List_afterExchange[[nParameter]][[2]]
	      }
	  
	    }
	
		# Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
	    Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 
	}  

    Sel_Cycle_Progeny_F5_Fam_CrossList_Inv <- rep(list(rep(list(array(0,c(nMarkers,2,no_selected_Fam))),nFamilies)),nParameters)
	for(nFamily in 1:nFamilies){
	 for(nParameter in 1:nParameters){
		Sel_Cycle_Progeny_F5_Fam_CrossList_Inv[[nParameter]][[nFamily]] <- Sel_Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]]   
	 }
	}


### Get selected nextGenGenoTable_List after exchange
	
	nextGenGenoTableDataList <-  foreach(nParameter=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getNextGenGenoTableData_Sel_V2(Sel_Cycle_Progeny_F5_Fam_CrossList_Inv[[nParameter]],nCrosses_List[[1]][[nParameter]],no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		selNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		selNextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		selNewNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
	 }
    }
  
#################### GM method for isolated families 

 system.time({

	genoSelected_List <- foreach(k=1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% 
	(getIM_GM_Frontier_selectedGenoList_SingleFams_V2(selNextGenGenoTable_GM_List[[k]][[j]],selNextGenGenoTableMod_GM_List[[k]][[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))

 })


## Reduce Parameter set to 30  

  #selectedParents_Geno_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
  #CriterionValue_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  for(nFamily in 1:nFamilies){
    for(nParameter in 1:nParameters){
   
      selectedParents_Geno_Fam_List[[nCyc]][[nFamily]][[nParameter]] <- genoSelected_List[[nFamily]][[nParameter]][[1]]
      CriterionValue_Fam_List[[nCyc]][[nFamily]][[nParameter]] <-genoSelected_List[[nFamily]][[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nFamily]][[nParameter]][[1]]
      
    }
  }  
 
 ###  
    
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
    
   ### 
	
    Quartile3 <- summary(GainValues[[nFamily]])[5]
	Quartile1 <- summary(GainValues[[nFamily]])[2]
	MedianG <- 	summary(GainValues[[nFamily]])[3]	
		
	indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues[[nFamily]] <= MedianG) & (GainValues[[nFamily]] >= Quartile1))
	indices2B <- which((GainValues[[nFamily]] >= MedianG) & (GainValues[[nFamily]] <= Quartile3))
	indices2CG <- (which(GainValues[[nFamily]] <= MedianG+1 & GainValues[[nFamily]] >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
######

    			
    for(nParameter in 1:nParameters){
            nParameterIndex <- indicesRange[nParameter]
	        selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <-  selectedParents_Geno_Par_List[[nFamily]][[nParameterIndex]]
			CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- c(GainValues[[nFamily]][nParameterIndex],UsefulnessValues[[nFamily]][nParameterIndex],InbreedValues[[nFamily]][nParameterIndex])
		
	} 
   
   }

  save.image(WorkspaceName)
  gc() 	
   	
####################################################################################################
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
	print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]],mean))),sep=""))
		 
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
	
    save.image(WorkspaceName)
	
 } 
 
######
	simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	return(simResults_List)

}

################ 




runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3ASel <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach,GIndData){

    
######## Get predicted geno and pheno values for cycle1 ###################################################
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
 	
	  gIndData <- GIndData
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

  #if(startCycle ==1){
	
   	  Migration_List <- list()
	  PCA_List <- list()
	  j<-1
 	  selectedParents_Geno_List <- list() 
	  selectedParents_Pheno_List <- list()
	  selectedParents_NumProgeny_Geno_List <- list()
	  selectedParents_NumProgeny_Pheno_List <- list() 
	  
	  familyPairs_List <- list()
	  familyInd_List <- list()

	  CombinedTableCount <- 1
	  
	  no_selected_Fam <- 0.1 * cycle1_nProgeny
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
# ###############################################################################################	  	  	  
# ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	   
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
## selected Line indices to select lines from cycle genotable

	
		genoSelectionTable <- ApplySelection_Discrete(SelCriterion_List,no_selected)
		genoSimSelectionTable <- genoSelectionTable
		selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
        selectedLineIndices <- selectedGenoIndividualIndices

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
    Cycle_Progeny_F5_List<- cycle1_GenoData_List

    if(migFreq !=0){
       
 	   if(nCyc==1 && nCyc%%migFreq==0){
	  
	    if(Policy != "GeneticClusters"){
		 
		 output_List_afterExchange <- exchangeGenoData_GM_Sel(Cycle_Progeny_F5_List,SelCriterion_List,selectedLineIndices,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
		}else if(Policy == "GeneticClusters"){ 
		
		 
		 output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		}
		  
		Sel_Cycle_Progeny_F5_List <- (output_List_afterExchange)[[1]]
		
		familyInfo <- output_List_afterExchange[[2]]
		
		Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   }
	   
	    no_selected_Fam <- 10 
		
        Sel_Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,no_selected_Fam))
        nProgIndex <-  dim(Sel_Cycle_Progeny_F5_List[[nFamily]])[3]
	    nProg <-1
				
	    for(nFamily in 1:length(families)){ 
					 
		    Sel_Cycle_Progeny_F5[nFamily,,,nProg:nProgIndex] <- Sel_Cycle_Progeny_F5_List[[nFamily]][,,nProg:nProgIndex]
	    } 
		
	}
	
### selected genoTable data after exchange
	print(paste("Cyc-",nCyc,"before-genFn"))
    selCycle1GenoTable_GM <-  generateMinus101GenoFormat_V1(Sel_Cycle_Progeny_F5,nFamilies,no_selected_Fam)
	print(paste("Cyc-",nCyc,"after-genFn"))
    selNewCycle1GenoTable_GM <- generateWeightedGenoTable(selCycle1GenoTable_GM,no_QTL)
    selCycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(selCycle1GenoTable_GM)

## Get simulated geno and pheno values 	
 
    genoValSimValues <- simulateGenotypicValues(cycle1GenoTable_GM,no_QTL) 
    phenoValSimValues <- simulatePhenotypicValues(cycle1GenoTable_GM,varE,no_QTL) 

### i) Get cycle 1 genotype table and criterion value list for sets of singleislands and pairs of islands	
### Get list of cycle1 geno table for family ind and family pairs ######################################
	
	cycle1GenoTable_GM_List <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
	cycle1GenoTableMod_GM_List <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
    selCycle1GenoTable_GM_List <- rep(list(matrix(0,nrow=no_selected_Fam,ncol=nMarkers)),nFamilies)
	selCycle1GenoTableMod_GM_List <- rep(list(matrix(0,nrow=no_selected_Fam,ncol=nMarkers)),nFamilies)
	
	initIndex <- 1
	finalIndex <- cycle1_nProgeny_Fam
	
	initSelIndex <- 1
	finalSelIndex <- no_selected_Fam
	
	for(nFam in 1:nFamilies){ 
		
		cycle1GenoTable_GM_List[[nFam]] <- cycle1GenoTable_GM[initIndex:finalIndex,]
		cycle1GenoTableMod_GM_List[[nFam]] <- cycle1GenoTableMod_GM[initIndex:finalIndex,]
		
		selCycle1GenoTable_GM_List[[nFam]] <- selCycle1GenoTable_GM[initSelIndex:finalSelIndex,]
        selCycle1GenoTableMod_GM_List[[nFam]] <- selCycle1GenoTableMod_GM[initSelIndex:finalSelIndex,]
	    
		initIndex <- finalIndex+1
		finalIndex <- finalIndex + cycle1_nProgeny_Fam
		
	}
	
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
	
## Write Genotable DF 	
	#cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	#cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	#cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	#write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 

######################################################################################   

#### GM method for parent selection in cycle 1 ################################################
# #### GM method for isolated families
    nParameters <- 30
 			 
  selectedParents_Geno_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles) 
  CriterionValue_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  system.time({

    genoSelected_List <- foreach(k=1:nFamilies,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getIM_GM_Frontier_selectedGenoList_SingleFams_V2(selCycle1GenoTable_GM_List[[k]],selCycle1GenoTableMod_GM_List[[k]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))
  
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
    
    Sel_Cycle_GenoData_List <- Sel_Cycle_Progeny_F5_List

## Get F5 RILs for selected parents in cycle 1
   
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getF5RILs_FS_IM_GM_Frontier_Sel(Sel_Cycle_GenoData_List[[k]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
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

#}

###########################################################################

  startCycle <- 2
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
 
  PredictionModel_Geno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Pheno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Geno_Parameter_List <- list() 
  PredictionModel_Pheno_Parameter_List <- list()
  
  genoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  genoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  
  trainGenoNewTablePreCycle_List <- list()
  trainGenoNewTablePreCycle_Mat_List <- list()
  
  trainSimPhenoTablePreCycle_List <- list()
  trainSimGenoValTablePreCycle_List <- list()
	
  trainGeno_backingFileName_List <- list()
  trainGeno_descriptorFileName_List <- list()
 
    
 for(nCyc in startCycle:nCycles){

    
    print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
    
	
    nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
		genoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[4]]
		phenoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[5]]
	 }
	}
	
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	nextGenGenoTable_GM_Fam_Par <- list()
	nextGenGenoTableMod_GM_Fam_Par <- list()
	newNextGenGenoTable_GM_Fam_Par <- list()
	genoValSimValues_Fam_Par <- list()
	phenoValSimValues_Fam_Par <- list()
	
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
    		
			nextGenGenoTable_GM_Fam_Par[[nParameter]] <- nextGenGenoTable_GM_Fam
			nextGenGenoTableMod_GM_Fam_Par[[nParameter]] <- nextGenGenoTableMod_GM_Fam
			newNextGenGenoTable_GM_Fam_Par[[nParameter]] <- newNextGenGenoTable_GM_Fam
			genoValSimValues_Fam_Par[[nParameter]] <- genoValSimValues_Fam
			phenoValSimValues_Fam_Par[[nParameter]] <- phenoValSimValues_Fam

#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
	if(gIndData ==TRUE){
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
	}
  }
	
############################################################################################### 
    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	    if(nCyc ==2){
		  trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		
		  PredictionModels <- foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount))
 		}else if(nCyc >2){ 
		   for(nParameter in 1:nParameters){
			  trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		      trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
			}
		 
 		  PredictionModels <- foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount))
		}		    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
		  trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="")
		  trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
			
		  trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
		  trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		
					
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
		 
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
          newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
	}
	
	
	if(modelUpdate ==FALSE){ 
	    
		for(nParameter in 1:nParameters){
		
		  # PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModel_Geno
		  # PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModel_Pheno
		  
		  # PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  # PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		  
		  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]
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
   }

### Get Predicted Geno and Pheno Value List
 
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

#### Set Selection Criterion List

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


### Apply selection to cycle nCyc population 
### Genoselection table & Island migration for cycle nCyc

	
	genoSelectionTable <- ApplySelection_Discrete_Frontier(SelCriterion_List,no_selected,nFamilies,nParameters)
    selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) lapply(x,function(y) as.vector(y[,2]))) 
    selectedLineIndices <- selectedGenoIndividualIndices

######################################################################################   
## get correct table format for migration policy and exchange genodata functions

    Cycle_Progeny_F5_List_Fam_V2  <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
    nCross_inFamily <- nCrosses/nFamilies

   for(nFamily in 1:nFamilies){ 
	for(nParameter in 1:nParameters){
            initIndex <- 1
    	    finalIndex <- 10 

            for(nCross in 1:nCross_inFamily){ 

		Cycle_Progeny_F5_List_Fam_V2[[nFamily]][[nParameter]][,,initIndex:finalIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,]
              
		initIndex <- finalIndex+1
		finalIndex <- finalIndex + no_selected_Fam
            }

     	}
     } 


      Cycle_Progeny_F5_List_Fam_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
      Sel_Cycle_Progeny_F5_List_Fam_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
      nextGenGenoTable_GM_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters) 
      genoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      phenoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      SelCriterion_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
	
	  selNextGenGenoTable_GM_List <- rep(list(rep(list(),nFamilies)),nParameters)
	  selNextGenGenoTableMod_GM_List <- rep(list(rep(list(),nFamilies)),nParameters)
          selNewNextGenGenoTable_GM_List <- rep(list(rep(list(),nFamilies)),nParameters)
	  
	  for(nFamily in 1:nFamilies){
    
	    for(nParameter in 1:nParameters){
	
	        Cycle_Progeny_F5 <-   Cycle_Progeny_F5_List_Fam_V2[[nFamily]][[nParameter]]
	        Cycle_Progeny_F5_List_Fam_Reverse[[nParameter]][[nFamily]] <- Cycle_Progeny_F5
			Sel_Cycle_Progeny_F5 <- Cycle_Progeny_F5[,,as.vector(selectedLineIndices[[nFamily]][[nParameter]])]
			Sel_Cycle_Progeny_F5_List_Fam_Reverse[[nParameter]][[nFamily]] <- Sel_Cycle_Progeny_F5
			SelCriterion_List_Reverse[[nParameter]][[nFamily]] <- SelCriterion_List[[nFamily]][[nParameter]]

	    }
     }
 
       
     for(nParameter in 1:nParameters){

	   for(nFamily in 1:nFamilies){

        	nextGenGenoTable_GM_List_Reverse[[nParameter]][[nFamily]] <- nextGenGenoTable_GM_List[[nFamily]][[nParameter]] 
			
            genoSelectionTable_Parameter[[nParameter]][[nFamily]] <- genoSelectionTable[[nFamily]][[nParameter]]
			
	        selNextGenGenoTable_GM_List[[nParameter]][[nFamily]] <- nextGenGenoTable_GM_List[[nFamily]][[nParameter]][(selectedLineIndices[[nFamily]][[nParameter]]),]
			
			selNextGenGenoTableMod_GM_List[[nParameter]][[nFamily]] <- nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]][(selectedLineIndices[[nFamily]][[nParameter]]),] 

        }
     }

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

######################################################################################   
## get correct table format for exchange genodata functions


        Cycle_Progeny_F5_Fam_ParamList <- rep(list(array(0,c(nFamilies,nMarkers,2,cycle1_nProgeny))),nParameters)

        Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)),nParameters)
		
	Sel_Cycle_Progeny_F5_Fam_ParamList <- rep(list(array(0,c(nFamilies,nMarkers,2,no_selected_Fam))),nParameters)

        Sel_Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,no_selected_Fam))),nFamilies)),nParameters)
        
	nProgInit <- 1
	nProgFinal <-  dim(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[1]])[4]
        nCrosses_inFamily <- nSel_inFamily

        nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	     for(nFamily in 1:length(families)){ 

              nProg <- nProgInit 
	      nProgIndex <- nProgFinal

               for(nCross in 1:nCross_inFamily){
					 
		Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,nProgInit:nProgFinal] 
				
                Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]][,,nProg:nProgIndex]<- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,nProgInit:nProgFinal] 
               
                			
		    
		nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal


	       }
	     } 
	}
	

        nProg <- 1 
	nProgIndex <- no_selected_Fam

        selectedLineIndices_Reverse <- rep(list(rep(list(c()),nFamilies)),nParameters)
	
        for(nParameter in 1:nParameters){
	 		
	     for(nFamily in 1:length(families)){ 

             
                Sel_Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,as.vector(selectedLineIndices[[nFamily]][[nParameter]])]
				
				Sel_Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]][,,nProg:nProgIndex] <-  Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]][,,as.vector(selectedLineIndices[[nFamily]][[nParameter]])]

               selectedLineIndices_Reverse[[nParameter]][[nFamily]] <-  selectedLineIndices[[nFamily]][[nParameter]]
          
	  }
	}
 

################
### Exchange selected geno data among pairs ##############
     
################################

    if(migFreq !=0){
   
 	  if(nCyc>=2 && nCyc%%migFreq==0){
	  
	     if(Policy != "GeneticClusters"){
		 
			output_List_afterExchange <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::exchangeGenoData_GM_Sel(Cycle_Progeny_F5_Fam_CrossList[[k]],SelCriterion_List_Reverse[[k]],selectedLineIndices_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
	      }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
	      }
		  
###########################################################################
### Get CycleProgeny Genodata and genotype table after exchange

		Sel_Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,no_selected_Fam))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily

		nProg <- nProgInit 
		nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	      for(nFamily in 1:length(families)){ 

           	nProg <- nProgInit 
	        nProgIndex <- nProgFinal

                for(nCross in 1:nCrosses_inFamily){
					      				
				Sel_Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
					    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
			
	        familyInfo <- output_List_afterExchange[[nParameter]][[2]]
	      }
	  
	    }
	
		# Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
	    Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 
	}   


### Get selected nextGenGenoTable_List after exchange
	
	nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getNextGenGenoTableData_Sel_V2(Sel_Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		selNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		selNextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		selNewNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
	 }
    }
  
#################### GM method for isolated families 

 system.time({

	genoSelected_List <- foreach(k=1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% 
	(getIM_GM_Frontier_selectedGenoList_SingleFams_V2(selNextGenGenoTable_GM_List[[k]][[j]],selNextGenGenoTableMod_GM_List[[k]][[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))

 })


## Reduce Parameter set to 30  

  #selectedParents_Geno_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
  #CriterionValue_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  for(nFamily in 1:nFamilies){
    for(nParameter in 1:nParameters){
   
      selectedParents_Geno_Fam_List[[nCyc]][[nFamily]][[nParameter]] <- genoSelected_List[[nFamily]][[nParameter]][[1]]
      CriterionValue_Fam_List[[nCyc]][[nFamily]][[nParameter]] <-genoSelected_List[[nFamily]][[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nFamily]][[nParameter]][[1]]
      
    }
  }  
 
 ###  
    
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
    
   ### 
	
    Quartile3 <- summary(GainValues[[nFamily]])[5]
	Quartile1 <- summary(GainValues[[nFamily]])[2]
	MedianG <- 	summary(GainValues[[nFamily]])[3]	
		
	indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues[[nFamily]] <= MedianG) & (GainValues[[nFamily]] >= Quartile1))
	indices2B <- which((GainValues[[nFamily]] >= MedianG) & (GainValues[[nFamily]] <= Quartile3))
	indices2CG <- (which(GainValues[[nFamily]] <= MedianG+1 & GainValues[[nFamily]] >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	
######

    			
    for(nParameter in 1:nParameters){
            nParameterIndex <- indicesRange[nParameter]
	        selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <-  selectedParents_Geno_Par_List[[nFamily]][[nParameterIndex]]
			CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- c(GainValues[[nFamily]][nParameterIndex],UsefulnessValues[[nFamily]][nParameterIndex],InbreedValues[[nFamily]][nParameterIndex])
		
	} 
   
   }

  save.image(WorkspaceName)
  gc() 	
   	
####################################################################################################
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
	print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]],mean))),sep=""))
		 
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
	
    save.image(WorkspaceName)
	
 } 
 
######
	simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	return(simResults_List)

}

################ 




runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3ASel_Reduced <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach,GIndData){

    
######## Get predicted geno and pheno values for cycle1 ###################################################
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
 	
	  gIndData <- GIndData
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

  #if(startCycle ==1){
	
   	  Migration_List <- list()
	  PCA_List <- list()
	  j<-1
 	  selectedParents_Geno_List <- list() 
	  selectedParents_Pheno_List <- list()
	  selectedParents_NumProgeny_Geno_List <- list()
	  selectedParents_NumProgeny_Pheno_List <- list() 
	  
	  familyPairs_List <- list()
	  familyInd_List <- list()

	  CombinedTableCount <- 1
	  
	  no_selected_Fam <- 0.1 * cycle1_nProgeny
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
# ###############################################################################################	  	  	  
# ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)
	   
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
## selected Line indices to select lines from cycle genotable

	
		genoSelectionTable <- ApplySelection_Discrete(SelCriterion_List,no_selected)
		genoSimSelectionTable <- genoSelectionTable
		selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
        selectedLineIndices <- selectedGenoIndividualIndices

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
    Cycle_Progeny_F5_List<- cycle1_GenoData_List

    if(migFreq !=0){
       
 	   if(nCyc==1 && nCyc%%migFreq==0){
	  
	    if(Policy != "GeneticClusters"){
		 
		 output_List_afterExchange <- exchangeGenoData_GM_Sel(Cycle_Progeny_F5_List,SelCriterion_List,selectedLineIndices,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
		}else if(Policy == "GeneticClusters"){ 
		
		 
		 output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		}
		  
		Sel_Cycle_Progeny_F5_List <- (output_List_afterExchange)[[1]]
		
		familyInfo <- output_List_afterExchange[[2]]
		
		Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   }
	   
	    no_selected_Fam <- 10 
		
        Sel_Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,no_selected_Fam))
        nProgIndex <-  dim(Sel_Cycle_Progeny_F5_List[[nFamily]])[3]
	    nProg <-1
				
	    for(nFamily in 1:length(families)){ 
					 
		    Sel_Cycle_Progeny_F5[nFamily,,,nProg:nProgIndex] <- Sel_Cycle_Progeny_F5_List[[nFamily]][,,nProg:nProgIndex]
	    } 
		
	}
	
### selected genoTable data after exchange
	print(paste("Cyc-",nCyc,"before-genFn"))
    selCycle1GenoTable_GM <-  generateMinus101GenoFormat_V1(Sel_Cycle_Progeny_F5,nFamilies,no_selected_Fam)
	print(paste("Cyc-",nCyc,"after-genFn"))
    selNewCycle1GenoTable_GM <- generateWeightedGenoTable(selCycle1GenoTable_GM,no_QTL)
    selCycle1GenoTableMod_GM <- generate012GenoFormat_Inverse(selCycle1GenoTable_GM)

## Get simulated geno and pheno values 	
 
    genoValSimValues <- simulateGenotypicValues(cycle1GenoTable_GM,no_QTL) 
    phenoValSimValues <- simulatePhenotypicValues(cycle1GenoTable_GM,varE,no_QTL) 

### i) Get cycle 1 genotype table and criterion value list for sets of singleislands and pairs of islands	
### Get list of cycle1 geno table for family ind and family pairs ######################################
	
	cycle1GenoTable_GM_List <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
	cycle1GenoTableMod_GM_List <- rep(list(matrix(0,nrow=cycle1_nProgeny_Fam,ncol=nMarkers)),nFamilies)
    selCycle1GenoTable_GM_List <- rep(list(matrix(0,nrow=no_selected_Fam,ncol=nMarkers)),nFamilies)
	selCycle1GenoTableMod_GM_List <- rep(list(matrix(0,nrow=no_selected_Fam,ncol=nMarkers)),nFamilies)
	
	initIndex <- 1
	finalIndex <- cycle1_nProgeny_Fam
	
	initSelIndex <- 1
	finalSelIndex <- no_selected_Fam
	
	for(nFam in 1:nFamilies){ 
		
		cycle1GenoTable_GM_List[[nFam]] <- cycle1GenoTable_GM[initIndex:finalIndex,]
		cycle1GenoTableMod_GM_List[[nFam]] <- cycle1GenoTableMod_GM[initIndex:finalIndex,]
		
		selCycle1GenoTable_GM_List[[nFam]] <- selCycle1GenoTable_GM[initSelIndex:finalSelIndex,]
        selCycle1GenoTableMod_GM_List[[nFam]] <- selCycle1GenoTableMod_GM[initSelIndex:finalSelIndex,]
	    
		initIndex <- finalIndex+1
		finalIndex <- finalIndex + cycle1_nProgeny_Fam
		
	}
	
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
	
## Write Genotable DF 	
	#cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	#cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	#cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	#write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 

######################################################################################   

#### GM method for parent selection in cycle 1 ################################################
# #### GM method for isolated families
    nParameters <- 10
 			 
  selectedParents_Geno_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles) 
  CriterionValue_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  system.time({

    genoSelected_List <- foreach(k=1:nFamilies,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getIM_GM_Frontier_selectedGenoList_SingleFams_V2_Reduced(selCycle1GenoTable_GM_List[[k]],selCycle1GenoTableMod_GM_List[[k]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))
  
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
    
    Sel_Cycle_GenoData_List <- Sel_Cycle_Progeny_F5_List

## Get F5 RILs for selected parents in cycle 1
   
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getF5RILs_FS_IM_GM_Frontier_Sel(Sel_Cycle_GenoData_List[[k]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
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

#}

###########################################################################

  startCycle <- 2
  nParameters <- 10
  nFamilies <- 20
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
 
  PredictionModel_Geno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Pheno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Geno_Parameter_List <- list() 
  PredictionModel_Pheno_Parameter_List <- list()
  
  genoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  genoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  
  trainGenoNewTablePreCycle_List <- list()
  trainGenoNewTablePreCycle_Mat_List <- list()
  
  trainSimPhenoTablePreCycle_List <- list()
  trainSimGenoValTablePreCycle_List <- list()
	
  trainGeno_backingFileName_List <- list()
  trainGeno_descriptorFileName_List <- list()
 
    
 for(nCyc in startCycle:nCycles){

    
    print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
    
	
    nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
		genoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[4]]
		phenoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[5]]
	 }
	}
	
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	nextGenGenoTable_GM_Fam_Par <- list()
	nextGenGenoTableMod_GM_Fam_Par <- list()
	newNextGenGenoTable_GM_Fam_Par <- list()
	genoValSimValues_Fam_Par <- list()
	phenoValSimValues_Fam_Par <- list()
	
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
    		
			nextGenGenoTable_GM_Fam_Par[[nParameter]] <- nextGenGenoTable_GM_Fam
			nextGenGenoTableMod_GM_Fam_Par[[nParameter]] <- nextGenGenoTableMod_GM_Fam
			newNextGenGenoTable_GM_Fam_Par[[nParameter]] <- newNextGenGenoTable_GM_Fam
			genoValSimValues_Fam_Par[[nParameter]] <- genoValSimValues_Fam
			phenoValSimValues_Fam_Par[[nParameter]] <- phenoValSimValues_Fam

#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
	if(gIndData ==TRUE){
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
	}
  }
	
############################################################################################### 
    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	    if(nCyc ==2){
		  trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		
		  PredictionModels <- foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount))
 		}else if(nCyc >2){ 
		   for(nParameter in 1:nParameters){
			  trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		      trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
			}
		 
 		  PredictionModels <- foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount))
		}		    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
		  trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="")
		  trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
			
		  trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
		  trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		
					
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
		 
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
          newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
	}
	
	
	if(modelUpdate ==FALSE){ 
	    
		for(nParameter in 1:nParameters){
		
		  # PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModel_Geno
		  # PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModel_Pheno
		  
		  # PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  # PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		  
		  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]
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
   }

### Get Predicted Geno and Pheno Value List
 
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

#### Set Selection Criterion List

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


### Apply selection to cycle nCyc population 
### Genoselection table & Island migration for cycle nCyc

	
	genoSelectionTable <- ApplySelection_Discrete_Frontier(SelCriterion_List,no_selected,nFamilies,nParameters)
    selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) lapply(x,function(y) as.vector(y[,2]))) 
    selectedLineIndices <- selectedGenoIndividualIndices

######################################################################################   
## get correct table format for migration policy and exchange genodata functions

    Cycle_Progeny_F5_List_Fam_V2  <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
    nCross_inFamily <- nCrosses/nFamilies

   for(nFamily in 1:nFamilies){ 
	for(nParameter in 1:nParameters){
            initIndex <- 1
    	    finalIndex <- 10 

            for(nCross in 1:nCross_inFamily){ 

		Cycle_Progeny_F5_List_Fam_V2[[nFamily]][[nParameter]][,,initIndex:finalIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,]
              
		initIndex <- finalIndex+1
		finalIndex <- finalIndex + no_selected_Fam
            }

     	}
     } 


      Cycle_Progeny_F5_List_Fam_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
      Sel_Cycle_Progeny_F5_List_Fam_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
      nextGenGenoTable_GM_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters) 
      genoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      phenoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      SelCriterion_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
	
	  selNextGenGenoTable_GM_List <- rep(list(rep(list(),nFamilies)),nParameters)
	  selNextGenGenoTableMod_GM_List <- rep(list(rep(list(),nFamilies)),nParameters)
          selNewNextGenGenoTable_GM_List <- rep(list(rep(list(),nFamilies)),nParameters)
	  
	  for(nFamily in 1:nFamilies){
    
	    for(nParameter in 1:nParameters){
	
	        Cycle_Progeny_F5 <-   Cycle_Progeny_F5_List_Fam_V2[[nFamily]][[nParameter]]
	        Cycle_Progeny_F5_List_Fam_Reverse[[nParameter]][[nFamily]] <- Cycle_Progeny_F5
			Sel_Cycle_Progeny_F5 <- Cycle_Progeny_F5[,,as.vector(selectedLineIndices[[nFamily]][[nParameter]])]
			Sel_Cycle_Progeny_F5_List_Fam_Reverse[[nParameter]][[nFamily]] <- Sel_Cycle_Progeny_F5
			SelCriterion_List_Reverse[[nParameter]][[nFamily]] <- SelCriterion_List[[nFamily]][[nParameter]]

	    }
     }
 
       
     for(nParameter in 1:nParameters){

	   for(nFamily in 1:nFamilies){

        	nextGenGenoTable_GM_List_Reverse[[nParameter]][[nFamily]] <- nextGenGenoTable_GM_List[[nFamily]][[nParameter]] 
			
            genoSelectionTable_Parameter[[nParameter]][[nFamily]] <- genoSelectionTable[[nFamily]][[nParameter]]
			
	        selNextGenGenoTable_GM_List[[nParameter]][[nFamily]] <- nextGenGenoTable_GM_List[[nFamily]][[nParameter]][(selectedLineIndices[[nFamily]][[nParameter]]),]
			
			selNextGenGenoTableMod_GM_List[[nParameter]][[nFamily]] <- nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]][(selectedLineIndices[[nFamily]][[nParameter]]),] 

        }
     }

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

######################################################################################   
## get correct table format for exchange genodata functions


        Cycle_Progeny_F5_Fam_ParamList <- rep(list(array(0,c(nFamilies,nMarkers,2,cycle1_nProgeny))),nParameters)

        Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)),nParameters)
		
		Sel_Cycle_Progeny_F5_Fam_ParamList <- rep(list(array(0,c(nFamilies,nMarkers,2,no_selected_Fam))),nParameters)

        Sel_Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,no_selected_Fam))),nFamilies)),nParameters)
        
	nProgInit <- 1
	nProgFinal <-  dim(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[1]])[4]
        nCrosses_inFamily <- nSel_inFamily

        nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	     for(nFamily in 1:length(families)){ 

              nProg <- nProgInit 
	      nProgIndex <- nProgFinal

               for(nCross in 1:nCross_inFamily){
					 
				Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,nProgInit:nProgFinal] 
				
                Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]][,,nProg:nProgIndex]<- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nCross]][1,,,nProgInit:nProgFinal] 
           
				nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	       }
	     } 
	}
	

        nProg <- 1 
		nProgIndex <- no_selected_Fam

        selectedLineIndices_Reverse <- rep(list(rep(list(c()),nFamilies)),nParameters)
	
        for(nParameter in 1:nParameters){
	 		
	     for(nFamily in 1:length(families)){ 

             
                Sel_Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,as.vector(selectedLineIndices[[nFamily]][[nParameter]])]
				
				Sel_Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]][,,nProg:nProgIndex] <-  Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]][,,as.vector(selectedLineIndices[[nFamily]][[nParameter]])]

				selectedLineIndices_Reverse[[nParameter]][[nFamily]] <-  selectedLineIndices[[nFamily]][[nParameter]]
          
		  }
		}
 

################
### Exchange selected geno data among pairs ##############
     
################################

    if(migFreq !=0){
   
 	  if(nCyc>=2 && nCyc%%migFreq==0){
	  
	     if(Policy != "GeneticClusters"){
		 
			output_List_afterExchange <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::exchangeGenoData_GM_Sel(Cycle_Progeny_F5_Fam_CrossList[[k]],SelCriterion_List_Reverse[[k]],selectedLineIndices_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
	      }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
	      }
		  
###########################################################################
### Get CycleProgeny Genodata and genotype table after exchange

		Sel_Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,no_selected_Fam))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily

		nProg <- nProgInit 
		nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	      for(nFamily in 1:length(families)){ 

           	nProg <- nProgInit 
	        nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					      				
				Sel_Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
					    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
			
	        familyInfo <- output_List_afterExchange[[nParameter]][[2]]
	      }
	  
	    }
	
		# Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
	    Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 
	}   


### Get selected nextGenGenoTable_List after exchange
	
	nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getNextGenGenoTableData_Sel_V2(Sel_Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		selNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		selNextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		selNewNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
	 }
    }
  
#################### GM method for isolated families 

 system.time({

	genoSelected_List <- foreach(k=1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (getIM_GM_Frontier_selectedGenoList_SingleFams_V2_Reduced(selNextGenGenoTable_GM_List[[k]][[j]],selNextGenGenoTableMod_GM_List[[k]][[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))

 })


## Reduce Parameter set to 10  
  nParameters <- 10
 
  #selectedParents_Geno_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
  #CriterionValue_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  for(nFamily in 1:nFamilies){
    for(nParameter in 1:nParameters){
   
      selectedParents_Geno_Fam_List[[nCyc]][[nFamily]][[nParameter]] <- genoSelected_List[[nFamily]][[nParameter]][[1]]
      CriterionValue_Fam_List[[nCyc]][[nFamily]][[nParameter]] <-genoSelected_List[[nFamily]][[nParameter]][[2]]
      selectedParents_Geno <- genoSelected_List[[nFamily]][[nParameter]][[1]]
      
    }
  }  
 
 ###  
    
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
    
   ### 
	
    Quartile3 <- summary(GainValues[[nFamily]])[5]
	Quartile1 <- summary(GainValues[[nFamily]])[2]
	MedianG <- 	summary(GainValues[[nFamily]])[3]	
		
	indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues[[nFamily]] <= MedianG) & (GainValues[[nFamily]] >= Quartile1))
	indices2B <- which((GainValues[[nFamily]] >= MedianG) & (GainValues[[nFamily]] <= Quartile3))
	indices2CG <- (which(GainValues[[nFamily]] <= MedianG+1 & GainValues[[nFamily]] >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	#indicesRange <- c(indicesRange1G,indicesRange2G,indicesRange3G)
	indicesRange <- indicesRange3G
######

    			
    for(nParameter in 1:nParameters){
            nParameterIndex <- indicesRange[nParameter]
	        selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <-  selectedParents_Geno_Par_List[[nFamily]][[nParameterIndex]]
			CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- c(GainValues[[nFamily]][nParameterIndex],UsefulnessValues[[nFamily]][nParameterIndex],InbreedValues[[nFamily]][nParameterIndex])
	} 
  }

  save.image(WorkspaceName)
  gc() 	
   	
  
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(j=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::getF5RILs_FS_IM_GM_Frontier_Sel(Sel_Cycle_GenoData_List[[k]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
    nCrosses_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
	
  	for(nFamily in 1:nFamilies){
          for(nParameter in 1:nParameters){
          	Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]] <- F5RILs_GM_Fam[[nFamily]][[nParameter]][[1]]
          	nCrosses_List[[nFamily]][[nParameter]] <- F5RILs_GM_Fam[[nFamily]][[nParameter]][[2]]
 	
         }
    }

    gc()
 
	
	
####################################################################################################
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
	print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]],mean))),sep=""))
		 
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
	
    save.image(WorkspaceName)
	
 } 
 
######
	simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	return(simResults_List)

}

################ 


runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3ASel_Reduced_doSNOW <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach,GIndData,SelCriteria){

    
######## Get predicted geno and pheno values for cycle1 ###################################################
     options(warn=-1)
### Assign Variables #########################################################################
## Selection Parameters 

  	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues
	  selCriteria <- SelCriteria
	
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
 	
	  gIndData <- GIndData
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

	  CombinedTableCount <- 1
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
# ###############################################################################################	  	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)

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
## selected Line indices to select lines from cycle genotable

	
		genoSelectionTable <- ApplySelection_Discrete(SelCriterion_List,no_selected)
		genoSimSelectionTable <- genoSelectionTable
		selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
        selectedLineIndices <- selectedGenoIndividualIndices

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
    Cycle_Progeny_F5_List<- cycle1_GenoData_List

    if(migFreq !=0){
       
 	    if(nCyc==1 && nCyc%%migFreq==0){
	  
	    if(Policy != "GeneticClusters"){
		 
		 output_List_afterExchange <- exchangeGenoData_GM(Cycle_Progeny_F5_List,SelCriterion_List,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
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
	
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
	
## Write Genotable DF 	
	cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 

######################################################################################   

#### GM method for parent selection in cycle 1 ################################################
# #### GM method for isolated families
    nParameters <- 10
 			 
 	selectedParents_Geno_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles) 
	CriterionValue_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  system.time({

    genoSelected_List <- foreach(k=1:nFamilies) %dopar% (getIM_GM_Frontier_selectedGenoList_SingleFams_V2_Reduced(cycle1GenoTable_GM_List[[k]],cycle1GenoTableMod_GM_List[[k]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))
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
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
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
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
 
  PredictionModel_Geno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Pheno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Geno_Parameter_List <- list() 
  PredictionModel_Pheno_Parameter_List <- list()
  
  genoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  genoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  
  trainGenoNewTablePreCycle_List <- list()
  trainGenoNewTablePreCycle_Mat_List <- list()
  
  trainSimPhenoTablePreCycle_List <- list()
  trainSimGenoValTablePreCycle_List <- list()
	
  trainGeno_backingFileName_List <- list()
  trainGeno_descriptorFileName_List <- list()
 
    
 for(nCyc in startCycle:nCycles){
   
    print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
   
    nProgeny_in_Cross <- 1 
	
    nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],nProgeny_in_Cross,no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
		genoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[4]]
		phenoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[5]]
	 }
	}
	
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	nextGenGenoTable_GM_Fam_Par <- list()
	nextGenGenoTableMod_GM_Fam_Par <- list()
	newNextGenGenoTable_GM_Fam_Par <- list()
	genoValSimValues_Fam_Par <- list()
	phenoValSimValues_Fam_Par <- list()
	
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
    		
			nextGenGenoTable_GM_Fam_Par[[nParameter]] <- nextGenGenoTable_GM_Fam
			nextGenGenoTableMod_GM_Fam_Par[[nParameter]] <- nextGenGenoTableMod_GM_Fam
			newNextGenGenoTable_GM_Fam_Par[[nParameter]] <- newNextGenGenoTable_GM_Fam
			genoValSimValues_Fam_Par[[nParameter]] <- genoValSimValues_Fam
			phenoValSimValues_Fam_Par[[nParameter]] <- phenoValSimValues_Fam

#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
	if(gIndData ==TRUE){
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
	}
  }
	
############################################################################################### 
    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	    if(nCyc ==2){
		  trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		
		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount))
 		}else if(nCyc >2){ 
		   for(nParameter in 1:nParameters){
			  #trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGeno_FileName<- trainGeno_FileName_List[[nParameter]]
			  trainGenoNewTablePreCycle <- read.table(trainGeno_FileName,sep="\t")
			}
		
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		      trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
		}
		 
 		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount))
				    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
          trainGeno_FileName<- trainGeno_FileName_List[[nParameter]]
		  trainGenoNewTablePreCycle <- read.table(trainGeno_FileName,sep="\t")		  
		  trainGeno_FileName <- paste("trainTable_",selCriteria,"_",BD,"_",Policy,"_",condition,"_",Rep,"_",nParameter,".txt",sep="") 
		  write.table(trainGenoNewTablePreCycle,trainGeno_FileName,sep="\t") 
		  rm(trainGenoNewTablePreCycle)
		   
		  # trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="")
		  # trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
		  # trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
		  # trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		  # trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
					
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		 		 
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
          newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
	
 }
	
	if(modelUpdate ==FALSE){ 
	    for(nParameter in 1:nParameters){
		
		  # PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModel_Geno
		  # PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModel_Pheno
		  
		  # PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  # PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		  
		  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]
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
	}

### Get Predicted Geno and Pheno Value List
 
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

#### Set Selection Criterion List

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


### Apply selection to cycle nCyc population 
### Genoselection table & Island migration for cycle nCyc

	
	genoSelectionTable <- ApplySelection_Discrete_Frontier(SelCriterion_List,no_selected,nFamilies,nParameters)
    selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) lapply(x,function(y) as.vector(y[,2]))) 
    selectedLineIndices <- selectedGenoIndividualIndices

######################################################################################   
## get correct table format for migration policy and exchange genodata functions

      Cycle_Progeny_F5_List_Fam_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
      nextGenGenoTable_GM_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters) 
      genoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      phenoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      SelCriterion_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)

      for(nFamily in 1:nFamilies){
    
	    for(nParameter in 1:nParameters){
	
	        Cycle_Progeny_F5 <-   Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]]
	 
	        Cycle_Progeny_F5_List_Fam_Reverse[[nParameter]][[nFamily]] <- Cycle_Progeny_F5
            SelCriterion_List_Reverse[[nParameter]][[nFamily]] <- SelCriterion_List[[nFamily]][[nParameter]]

	    }
      }
 
       
     for(nParameter in 1:nParameters){

	   for(nFamily in 1:nFamilies){

        	nextGenGenoTable_GM_List_Reverse[[nParameter]][[nFamily]] <- nextGenGenoTable_GM_List[[nFamily]][[nParameter]] 

            genoSelectionTable_Parameter[[nParameter]][[nFamily]] <- genoSelectionTable[[nFamily]][[nParameter]]

            }
     }

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

######################################################################################   
## get correct table format for exchange genodata functions


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
	

################

	
### Exchange selected geno data among pairs ##############
     
################################

    if(migFreq !=0){
   
 	  if(nCyc>=2 && nCyc%%migFreq==0){
	  
	         if(Policy != "GeneticClusters"){
		 
			output_List_afterExchange <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::exchangeGenoData_GM(Cycle_Progeny_F5_Fam_CrossList[[k]],SelCriterion_List_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
		 }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		 }
		  
###########################################################################
### Get CycleProgeny Genodata and genotype table after exchange

		Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily

		nProg <- nProgInit 
		nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	      for(nFamily in 1:length(families)){ 

           	nProg <- nProgInit 
	        nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					      				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
					    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
			
			  familyInfo <- output_List_afterExchange[[nParameter]][[2]]
	      }
	  
	    }
	
		# Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
	    Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 
	}   


### Get nextGenGenoTable_List after exchange
	
	nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
		genoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[4]]
		phenoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[5]]
	 }
	}
  
################################################################################################### GM method for isolated families

 system.time({

genoSelected_List <- foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% 
(getIM_GM_Frontier_selectedGenoList_SingleFams_V2_Reduced(nextGenGenoTable_GM_List[[k]][[j]],nextGenGenoTableMod_GM_List[[k]][[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))

 })


## Reduce Parameter set to 30  

  #selectedParents_Geno_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
  #CriterionValue_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

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
    
   ### 
	
    Quartile3 <- summary(GainValues[[nFamily]])[5]
	Quartile1 <- summary(GainValues[[nFamily]])[2]
	MedianG <- 	summary(GainValues[[nFamily]])[3]	
		
	indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues[[nFamily]] <= MedianG) & (GainValues[[nFamily]] >= Quartile1))
	indices2B <- which((GainValues[[nFamily]] >= MedianG) & (GainValues[[nFamily]] <= Quartile3))
	indices2CG <- (which(GainValues[[nFamily]] <= MedianG+1 & GainValues[[nFamily]] >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange3G)
	
######
    			
    for(nParameter in 1:nParameters){
            nParameterIndex <- indicesRange[nParameter]
	        selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <-  selectedParents_Geno_Par_List[[nFamily]][[nParameterIndex]]
			CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- c(GainValues[[nFamily]][nParameterIndex],UsefulnessValues[[nFamily]][nParameterIndex],InbreedValues[[nFamily]][nParameterIndex])
		
	} 
  }

  save.image(WorkspaceName)
  gc() 	
   	
###########################################################################
## Get F5 RILs for selected parents in cycle nCyc

    Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
	nProgInit <- 1
	nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
    nCrosses_inFamily <- nSel_inFamily

    nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
     for(nParameter in 1:nParameters){
	 		
	   for(nFamily in 1:length(families)){ 

           nProg <- nProgInit 
	       nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					 
		       				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
				
		    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
	    } 
	}

######################################################################################################## 

    Cycle_GenoData_List <- Cycle_Progeny_F5_Fam_CrossList
		
    Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
    nCrosses_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
	
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% (getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]][[j]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))

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

print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]],mean))),sep=""))

		 
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
	
    save.image(WorkspaceName)
	
 
 
 }
 
######
	simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	return(simResults_List)

}


########

runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3A_Reduced <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach,GIndData){

    
######## Get predicted geno and pheno values for cycle1 ###################################################
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
 	
	  gIndData <- GIndData
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

	  CombinedTableCount <- 1
###########################################################################################
	
	  trainGenoNewTable_List <- list() 
	  trainSimGenoValTable_List <- list() 
	  trainSimPhenoTable_List <- list()
	  
# ###############################################################################################	  	  	  # ### Get predicted geno and pheno values for cycle1 ###################################################
    
	nCyc <-1
	print(nCyc)

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
## selected Line indices to select lines from cycle genotable

	
		genoSelectionTable <- ApplySelection_Discrete(SelCriterion_List,no_selected)
		genoSimSelectionTable <- genoSelectionTable
		selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
        selectedLineIndices <- selectedGenoIndividualIndices

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
    Cycle_Progeny_F5_List<- cycle1_GenoData_List

    if(migFreq !=0){
       
 	    if(nCyc==1 && nCyc%%migFreq==0){
	  
	    if(Policy != "GeneticClusters"){
		 
		 output_List_afterExchange <- exchangeGenoData_GM(Cycle_Progeny_F5_List,SelCriterion_List,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
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
	
#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
	
## Write Genotable DF 	
	cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable_GM,AlleleConversionTable_Combined)
 	cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,".txt",sep=""),sep="\t") 

######################################################################################   

#### GM method for parent selection in cycle 1 ################################################
# #### GM method for isolated families
    nParameters <- 10
 			 
 	selectedParents_Geno_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles) 
	CriterionValue_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  system.time({

    genoSelected_List <- foreach(k=1:nFamilies) %dopar% (getIM_GM_Frontier_selectedGenoList_SingleFams_V2_Reduced(cycle1GenoTable_GM_List[[k]],cycle1GenoTableMod_GM_List[[k]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))
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
        
	Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
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
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
 
  PredictionModel_Geno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Pheno_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  PredictionModel_Geno_Parameter_List <- list() 
  PredictionModel_Pheno_Parameter_List <- list()
  
  genoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  genoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  phenoValSimValues_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  
  trainGenoNewTablePreCycle_List <- list()
  trainGenoNewTablePreCycle_Mat_List <- list()
  
  trainSimPhenoTablePreCycle_List <- list()
  trainSimGenoValTablePreCycle_List <- list()
	
  trainGeno_backingFileName_List <- list()
  trainGeno_descriptorFileName_List <- list()

  Cycle_Progeny_F5_Fam_Cyc1 <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies) 

  for(nFamily in 1:nFamilies){ 
    for(nParameter in 1:nParameters){ 
         for(nProgeny in 1:cycle1_nProgeny){
  		Cycle_Progeny_F5_Fam_Cyc1[[nFamily]][[nParameter]][,,nProgeny] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][[nProgeny]][1,,,1] 

            } 
    }
 } 
  
  Cycle_Progeny_F5_List_Fam <- Cycle_Progeny_F5_Fam_Cyc1

    
 for(nCyc in startCycle:nCycles){
   
    print(nCyc)
     
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families,
### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
    nProgeny_in_Cross <- 1  
	
    nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],nProgeny_in_Cross,no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
		genoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[4]]
		phenoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[5]]
	 }
	}
	
	genoValues_Parameter_List<- list()
	phenoValues_Parameter_List <-list()
	
	nextGenGenoTable_GM_Fam_Par <- list()
	nextGenGenoTableMod_GM_Fam_Par <- list()
	newNextGenGenoTable_GM_Fam_Par <- list()
	genoValSimValues_Fam_Par <- list()
	phenoValSimValues_Fam_Par <- list()
	
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
    		
			nextGenGenoTable_GM_Fam_Par[[nParameter]] <- nextGenGenoTable_GM_Fam
			nextGenGenoTableMod_GM_Fam_Par[[nParameter]] <- nextGenGenoTableMod_GM_Fam
			newNextGenGenoTable_GM_Fam_Par[[nParameter]] <- newNextGenGenoTable_GM_Fam
			genoValSimValues_Fam_Par[[nParameter]] <- genoValSimValues_Fam
			phenoValSimValues_Fam_Par[[nParameter]] <- phenoValSimValues_Fam

#####  Get population vector

		populations <- rep(0,nFamilies*100)
		init <- 1
		for(nPop in 1:nFamilies){

			final <- nPop*100
			populations[init:final]<- rep(nPop,100)

			init <- final+1

		}
        		
## Write genotable df	
	if(gIndData ==TRUE){
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("PM_GM_GenoDataInd_Con-",condition,"_Rep-",Rep,"_Cyc_",nCyc,"_Par-",nParameter,".txt",sep=""),sep="\t")
	}
  }
	
############################################################################################### 
    if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	    if(nCyc ==2){
		  trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		
		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount))
 		}else if(nCyc >2){ 
		   for(nParameter in 1:nParameters){
			  trainGeno_descriptorFileName <- trainGeno_descriptorFileName_List[[nParameter]]
			  trainGenoNewTablePreCycle_List[[nParameter]] <- attach.big.matrix(trainGeno_descriptorFileName)
		      trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- bigmemory::as.matrix(trainGenoNewTablePreCycle_List[[nParameter]])
			}
		 
 		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount))
		}		    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
		  trainGeno_backingFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".bin",sep="")
		  trainGeno_descriptorFileName <- paste("trainTable_",condition,"_",Rep,"_",nParameter,".desc",sep="")
			
		  trainGeno_backingFileName_List[[nParameter]] <- trainGeno_backingFileName 
		  trainGeno_descriptorFileName_List[[nParameter]] <- trainGeno_descriptorFileName
		
					
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  trainGenoNewTablePreCycle_List[[nParameter]] <- deepcopy(as.big.matrix(PredictionModels[[nParameter]][[4]]),backingfile=trainGeno_backingFileName,descriptorfile=trainGeno_descriptorFileName)
		 
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
		  
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
          newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
	}
	
	
	if(modelUpdate ==FALSE){ 
	    for(nParameter in 1:nParameters){
		
		  # PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModel_Geno
		  # PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModel_Pheno
		  
		  # PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  # PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		  
		  newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]
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
	}

### Get Predicted Geno and Pheno Value List
 
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

#### Set Selection Criterion List

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


### Apply selection to cycle nCyc population 
### Genoselection table & Island migration for cycle nCyc

	
	genoSelectionTable <- ApplySelection_Discrete_Frontier(SelCriterion_List,no_selected,nFamilies,nParameters)
    selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) lapply(x,function(y) as.vector(y[,2]))) 
    selectedLineIndices <- selectedGenoIndividualIndices

######################################################################################   
## get correct table format for migration policy and exchange genodata functions

      Cycle_Progeny_F5_List_Fam_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)
      nextGenGenoTable_GM_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters) 
      genoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      phenoSelectionTable_Parameter <- rep(list(rep(list(),nFamilies)),nParameters)
      SelCriterion_List_Reverse <- rep(list(rep(list(),nFamilies)),nParameters)

      for(nFamily in 1:nFamilies){
    
	    for(nParameter in 1:nParameters){
	
	        Cycle_Progeny_F5 <-   Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]]
	 
	        Cycle_Progeny_F5_List_Fam_Reverse[[nParameter]][[nFamily]] <- Cycle_Progeny_F5
            SelCriterion_List_Reverse[[nParameter]][[nFamily]] <- SelCriterion_List[[nFamily]][[nParameter]]

	    }
      }
 
       
     for(nParameter in 1:nParameters){

	   for(nFamily in 1:nFamilies){

        	nextGenGenoTable_GM_List_Reverse[[nParameter]][[nFamily]] <- nextGenGenoTable_GM_List[[nFamily]][[nParameter]] 

            genoSelectionTable_Parameter[[nParameter]][[nFamily]] <- genoSelectionTable[[nFamily]][[nParameter]]

            }
     }

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

######################################################################################   
## get correct table format for exchange genodata functions


        Cycle_Progeny_F5_Fam_ParamList <- rep(list(array(0,c(nFamilies,nMarkers,2,cycle1_nProgeny))),nParameters)
        Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)),nParameters)
        
	nProgInit <- 1
	nProgFinal <-  dim(Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]])[3]
    nCrosses_inFamily <- nSel_inFamily

        nProg <- nProgInit 
	nProgIndex <- nProgFinal
      
    
        for(nParameter in 1:nParameters){
	 		
	    for(nFamily in 1:length(families)){ 

              nProg <- nProgInit 
	      nProgIndex <- nProgFinal

            			 
		Cycle_Progeny_F5_Fam_ParamList[[nParameter]][nFamily,,,nProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] 
                Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]] <- Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]]
		
		#nProg <- nProgIndex+1 
                #nProgIndex <- nProgIndex + nProgFinal
            }
	} 
		

################
	
### Exchange selected geno data among pairs ##############
     
################################

    if(migFreq !=0){
   
 	  if(nCyc>=2 && nCyc%%migFreq==0){
	  
	         if(Policy != "GeneticClusters"){
		 
			output_List_afterExchange <- foreach(k=1:nParameters,.packages="LTSGenoIMV2") %dopar% (LTSGenoIMV2::exchangeGenoData_GM(Cycle_Progeny_F5_Fam_CrossList[[k]],SelCriterion_List_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    
		 }else if(Policy == "GeneticClusters"){ 

			output_List_afterExchange <- exchangeGenoData_GeneticClusters(Cycle_Progeny_F5_List_Fam_Reverse,migrationSize,emigrantGroups,immigrantGroups,emigrantGroup_Prob,direction,Policy)
		 
		 }
		  
###########################################################################
### Get CycleProgeny Genodata and genotype table after exchange

		Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily

		nProg <- nProgInit 
		nProgIndex <- nProgFinal
    
        for(nParameter in 1:nParameters){
	 		
	      for(nFamily in 1:length(families)){ 

           	nProg <- nProgInit 
	        nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					      				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
					    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
			
			  familyInfo <- output_List_afterExchange[[nParameter]][[2]]
	      }
	  
	    }
	
		# Cycle_Progeny_F5_List_Fam <- (output_List_afterExchange)[[1]]
	    Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 
	}   


### Get nextGenGenoTable_List after exchange
	
	nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies,.packages="LTSGenoIMV2") %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]],nCrosses_List[[nFamily]][[nParameter]],nProgeny_in_Cross,no_QTL,varE))
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
		genoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[4]]
		phenoValSimValues_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[5]]
	 }
	}
  
################################################################################################### GM method for isolated families

 system.time({

genoSelected_List <- foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% 
(getIM_GM_Frontier_selectedGenoList_SingleFams_V2_Reduced(nextGenGenoTable_GM_List[[k]][[j]],nextGenGenoTableMod_GM_List[[k]][[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))

 })


## Reduce Parameter set to 30  

  #selectedParents_Geno_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)
  #CriterionValue_Fam_List <-  rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

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
    
   ### 
	
    Quartile3 <- summary(GainValues[[nFamily]])[5]
	Quartile1 <- summary(GainValues[[nFamily]])[2]
	MedianG <- 	summary(GainValues[[nFamily]])[3]	
		
	indices1 <- which(GainValues[[nFamily]] >= Quartile3)
    indicesSortRange1 <- sort.int(GainValues[[nFamily]][indices1],decreasing=TRUE,index.return=TRUE)
	indicesRange1G <- indicesSortRange1[[2]][1:10]
			
	indices2A <- which((GainValues[[nFamily]] <= MedianG) & (GainValues[[nFamily]] >= Quartile1))
	indices2B <- which((GainValues[[nFamily]] >= MedianG) & (GainValues[[nFamily]] <= Quartile3))
	indices2CG <- (which(GainValues[[nFamily]] <= MedianG+1 & GainValues[[nFamily]] >= MedianG-1))[1]
	indicesRange2G <- c(sample(indices2A,4),indices2CG,sample(indices2B,5))

	indices3 <- which(GainValues[[nFamily]] <= Quartile1)
    indicesSortRange3 <- sort.int(GainValues[[nFamily]][indices3],decreasing=FALSE,index.return=TRUE) 
	indicesRange3G <- indicesSortRange3[[2]][1:10]
			
	indicesRange <- c(indicesRange3G)
	
######
    			
    for(nParameter in 1:nParameters){
            nParameterIndex <- indicesRange[nParameter]
	        selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <-  selectedParents_Geno_Par_List[[nFamily]][[nParameterIndex]]
			CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- c(GainValues[[nFamily]][nParameterIndex],UsefulnessValues[[nFamily]][nParameterIndex],InbreedValues[[nFamily]][nParameterIndex])
		
	} 
  }

  save.image(WorkspaceName)
  gc() 	
   	
###########################################################################
## Get F5 RILs for selected parents in cycle nCyc

    Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
	nProgInit <- 1
	nProgFinal <-  dim(output_List_afterExchange[[nParameter]][[1]][[nFamily]])[3]
    nCrosses_inFamily <- nSel_inFamily

    nProg <- nProgInit 
	nProgIndex <- nProgFinal
    
     for(nParameter in 1:nParameters){
	 		
	   for(nFamily in 1:length(families)){ 

           nProg <- nProgInit 
	       nProgIndex <- nProgFinal

            for(nCross in 1:nCrosses_inFamily){
					 
		       				
				Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
				
		    
		        nProg <- nProgIndex+1 
                nProgIndex <- nProgIndex + nProgFinal
	        }
	    } 
	}

######################################################################################################## 

    Cycle_GenoData_List <- Cycle_Progeny_F5_Fam_CrossList
		
    Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
    nCrosses_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
	
      
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParameters) %dopar% (getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]][[j]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))

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

print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]],mean))),sep=""))

		 
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
	
    save.image(WorkspaceName)
	
 }
 
######
	simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	return(simResults_List)

}


######## 
