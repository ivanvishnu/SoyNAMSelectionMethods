
######## GM_Frontier Method on Islands 





#####

runSimulations20X_IslandSelection_BD_Wght_GM_Frontier_V3A_Complete_nParam1 <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,TrainTableWindowSize,SClassPlot,WorkSpaceName,StartCycle,GM_Frontier_Param_List,noCores_ForEach,GIndData,SelCriteria){ 

    
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

	  nProgeny_in_Family <- 100
	  nCrosses_in_Family <- nSel_inFamily
      nProgeny_per_Cross <- nProgeny_in_Family/nCrosses_in_Family
      
	  nProgeny_per_Family <- nCrosses_in_Family*nProgeny_per_Cross

    
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
		selectedGenoIndividualIndices <- lapply(genoSelectionTable,function(x) as.vector(x[,2])) 
        selectedLineIndices <- selectedGenoIndividualIndices
		
		
	  genoSelectionTable_List <- rep(list(),nFamilies)
      phenoSelectionTable_List <- rep(list(),nFamilies)
    
      for(nFamily in 1:nFamilies){

        				
            genoSelectionTable_List[[nFamily]] <- genoSelectionTable[[nFamily]]
				        
       }
    

### Migration Policy to identify immigration and emigration groups ##############################################


	    migrationGroups <- migrationPolicy_GM_V2(Cycle1GenoTable_GM_List,genoSelectionTable_List,Policy,nFamilies,no_QTL)
	
	    emigrantGroups <- migrationGroups[[1]]
	    immigrantGroups <- migrationGroups[[2]]
	
	

######################################################################################   
### Exchange selected geno data among pairs ##############
		Cycle_Progeny_F5_List_Fam <- cycle1_GenoData_List
    	 
		output_List_afterExchange <- LTSGenoIMV2::exchangeGenoData_GM(Cycle_Progeny_F5_List_Fam,SelCriterion_List,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    
		  
###########################################################################
### Get CycleProgeny Genodata and genotype table after exchange

		Cycle_Progeny_F5_Fam_CrossList <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[1]][[nFamily]])[3]
		nCrosses_inFamily <- nSel_inFamily

		nProg <- 1 
		nProgIndex <-cycle1_nProgeny 
	 		
	    for(nFamily in 1:length(families)){ 
                   		      				
			Cycle_Progeny_F5_Fam_CrossList[[nFamily]][,,nProgInit:nProgIndex] <- output_List_afterExchange[[1]][[nFamily]][,,nProgInit:nProgFinal]
					    
	    }
			
	    familyInfo <- output_List_afterExchange[[2]]
	    Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
### Get selected cycle1GenoTable_List after exchange

        nProgeny_per_Cross <- 1
	
		cycle1GenoTableDataList <- foreach(nFamily =1:nFamilies) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_Fam_CrossList[[nFamily]],nProgeny_in_Family,nProgeny_per_Cross,no_QTL,varE))

        Cycle1GenoTable_GM_List <- rep(list(list()),nFamilies)
        Cycle1GenoTableMod_GM_List <- rep(list(list()),nFamilies)
        Newcycle1GenGenoTable_GM_List <- rep(list(list()),nFamilies)

        for(nFamily in 1:nFamilies){
	

			Cycle1GenoTable_GM_List[[nFamily]]<-  cycle1GenoTableDataList[[nFamily]][[1]]
			Cycle1GenoTableMod_GM_List[[nFamily]] <- cycle1GenoTableDataList[[nFamily]][[2]]
			Newcycle1GenGenoTable_GM_List[[nFamily]] <- cycle1GenoTableDataList[[nFamily]][[3]]
		
		}
      
### i) Get cycle 1 genotype table and criterion value list for sets of singleislands and pairs of islands	
####  Get population vector
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
	nParameters <- 1
	cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	write.table(cycle1GenoData_DF,paste("IM_GM_GenoDataInd_Con-",condition,"_",selCriteria,"_",Policy,"_Rep-",Rep,"_Par-",nParameters,"_Cyc_",nCyc,".txt",sep=""),sep="\t")

######################################################################################   
#### GM method for parent selection in cycle 1 ################################################

    nParameters <- 1
			 
    selectedParents_Geno_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles) 
    CriterionValue_List <- rep(list(rep(list(rep(list(list()),nParameters)),nFamilies)),nCycles)

  system.time({

    genoSelected_List <- foreach(k=1:nFamilies) %dopar% (LTSGenoIMV2::getIM_GM_Frontier_selectedGenoList_SingleFams_V2(Cycle1GenoTable_GM_List[[k]],Cycle1GenoTableMod_GM_List[[k]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k))
  })
 
   for(nFamily in 1:nFamilies){
   
      nParameters <- length(genoSelected_List[[nFamily]][[1]])
	  for(nParameter in 1:nParameters){
		selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <- genoSelected_List[[nFamily]][[1]][[nParameter]]
		CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- genoSelected_List[[nFamily]][[2]][[nParameter]]
		
      }
   }


#######################################################################################

    print(paste(nCyc,"-",nParameter,"-",mean(genoValSimValues),sep="")) 
	
		
       nParams_List <- list()
       nParams_List_Cyc <- rep(list(list()),nCycles)
	for(k in 1:nFamilies){
	   nParams <- length(genoSelected_List[[k]][[1]])
	   nParams_List[[k]] <- nParams
	}
   
       nParams_List_Cyc[[nCyc]] <- nParams_List
	
## Get F5 RILs for selected parents in cycle 1
        
   F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:(nParams_List[[k]])) %dopar% (LTSGenoIMV2::getF5RILs_FS_IM_GM_Frontier(Cycle_Progeny_F5_Fam_CrossList[[k]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))

  
        
   Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
    nCrosses_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
	
  	for(nFamily in 1:nFamilies){
	      nParameters <- length(F5RILs_GM_Fam[[nFamily]])
          for(nParameter in 1:nParameters){
          	Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]] <- F5RILs_GM_Fam[[nFamily]][[nParameter]][[1]]
          	nCrosses_List[[nFamily]][[nParameter]] <- F5RILs_GM_Fam[[nFamily]][[nParameter]][[2]]
 	
         }
    }

    gc()
 
      nParameters <- 1
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
 
   nParameters <- length(F5RILs_GM_Fam[[nFamily]])
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
  nParameters <- 1
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
	
  trainGeno_FileName_List <- list()
  trainGeno_descriptorFileName_List <- list()

 Cycle_Progeny_F5_List <- Cycle_Progeny_F5_List_Fam

Cycle_Progeny_F5_List_Fam_V2 <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)

    nProgeny_in_Family <- 100
    nCrosses_in_Family <- nSel_inFamily
    nProgeny_per_Cross <- nProgeny_in_Family/nCrosses_in_Family
    nProgeny_per_Family <- nCrosses_in_Family*nProgeny_per_Cross


  for(nFamily in 1:nFamilies) {
        nParameters <- length(F5RILs_GM_Fam[[nFamily]])

 	for(nParameter in 1:nParameters){ 
             init <- 1
             final <- nProgeny_per_Cross
            
            for(nCrs in 1:nCrosses_in_Family ){
                    
                         

                    Cycle_Progeny_F5_List_Fam_V2[[nFamily]][[nParameter]][,,init:final] <- Cycle_Progeny_F5_List[[nFamily]][[nParameter]][[nCrs]][1,,,1:nProgeny_per_Cross]

                    init <- final+1
                    final <- final+nProgeny_per_Cross
            } 
      } 

 }

 for(nCyc in startCycle:nCycles){
 
     
    print(nCyc)
     
	 
	
#####  Get population vector

	populations <- rep(0,nFamilies*100)
	init <- 1
	for(nPop in 1:nFamilies){

		final <- nPop*100
		populations[init:final]<- rep(nPop,100)

		init <- final+1
	}
        		 
	 
### Given cycle geno data,selected parents and number of progeny list for pairs of families and isolated families, ### returns Cycle_Progeny_F5_FamPair,Cycle_Progeny_F5_List_FamPair,Cycle_Progeny_F5_Fam,Cycle_Progeny_F5_List_Fam 
### Generate F5 RIL progeny from geno data for selected geno data
    nextGenGenoTableDataList <- list()
	for(nFamily in 1:nFamilies){  
		nextGenGenoTableDataList[[nFamily]] <-  foreach(nParameter=1:(nParams_List[[nFamily]])) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_List_Fam_V2[[nFamily]][[nParameter]],nProgeny_in_Family,nProgeny_per_Cross,no_QTL,varE))
	}
	
    for(nFamily in 1:nFamilies){

      nParameters <- length(nextGenGenoTableDataList[[nFamily]])
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

          
  

        nParameters_Vec <- unlist(lapply(nextGenGenoTable_GM_List,length))

        nParameters <- length(nextGenGenoTableDataList[[nFamily]])

      	for(nParameter in 1:(max(nParameters_Vec))){
		
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


## Write genotable df	
	if(gIndData ==TRUE){
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable_GM_Fam,AlleleConversionTable_Combined)
        nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	    write.table(nextGenGenoData_DF,paste("IM_GM_GenoDataInd_Con-",condition,"_",selCriteria,"_",Policy,"_Rep-",Rep,"_Par-",nParameter,"_Cyc_",nCyc,".txt",sep=""),sep="\t")
	}
   
  }
  
  
 
	
############################################################################################### 
  
   if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	
	    if(nCyc ==2){
		  
		  trainGenoNewTablePreCycle_Mat <- bigmemory::as.matrix(trainGenoNewTablePreCycle)
		  nParameters <- length(newNextGenGenoTable_GM_Fam_Par)
		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat,trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount,nCyc))
		  
 		}else if(nCyc >2){
		
		    nParameters <- length(trainGeno_FileName_List)
            for(nParameter in 1:nParameters){
			  
			  trainGeno_FileName <- trainGeno_FileName_List[[nParameter]] 
			  trainGenoNewTablePreCycle_Mat_List[[nParameter]] <- read.table(trainGeno_FileName,sep="\t")
			}
		 
 		  PredictionModels <- foreach(j=1:nParameters) %dopar% (LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable_GM_Fam_Par[[j]],genoValSimValues_Fam_Par[[j]],phenoValSimValues_Fam_Par[[j]],PredictionModel_Geno_Parameter_List[[j]],PredictionModel_Pheno_Parameter_List[[j]],modelType,updateType,TrainTableWindowSize,trainGenoNewTablePreCycle_Mat_List[[j]],trainSimPhenoTablePreCycle_List[[j]],trainSimGenoValTablePreCycle_List[[j]],CombinedTableCount,nCyc))
		}	    	
				
		for(nParameter in 1:nParameters){
		
		  PredictionModel_Geno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[1]] 
		  PredictionModel_Pheno_Parameter_List[[nParameter]] <- PredictionModels[[nParameter]][[2]]
		  CombinedTableCount <- PredictionModels[[nParameter]][[3]]
		  
		  trainGeno_FileName <- paste("trainTable_",BD,"_",selCriteria,"_",Policy,"_",condition,"_",Rep,"_",nParameter,".txt",sep="") 
		 		  
		  trainGeno_FileName_List[[nParameter]] <- trainGeno_FileName
			
		  trainGenoNewTablePreCycle <- PredictionModels[[nParameter]][[4]]
		  trainSimPhenoTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[5]]
		  trainSimGenoValTablePreCycle_List[[nParameter]] <- PredictionModels[[nParameter]][[6]]
			
		  write.table(trainGenoNewTablePreCycle,trainGeno_FileName,sep="\t")
		 
		  PredictionModel_Geno <- PredictionModel_Geno_Parameter_List[[nParameter]] 
		  PredictionModel_Pheno <- PredictionModel_Pheno_Parameter_List[[nParameter]]
		   
          newNextGenGenoTable_GM_Fam <- newNextGenGenoTable_GM_Fam_Par[[nParameter]]

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
		
		    genoValues_FamPair <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_FamPair <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
			genoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues_Fam <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)
			
			genoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Geno)
			phenoValues <- PredictBayesPhenoValues(newNextGenGenoTable_GM_Fam,PredictionModel_Pheno)

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

               
		nParameter <- 1
	     newNextGenGenoTable_GM_Fam <- apply(newNextGenGenoTable_GM_Fam_Par[[nParameter]],2,function(x) as.numeric(as.character(x)))
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
 
    #nParameters <- 1
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

	nParameters <- 1
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
	
	        Cycle_Progeny_F5 <-   Cycle_Progeny_F5_List_Fam_V2[[nFamily]][[nParameter]]
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


	   migrationGroups <- migrationPolicy_GM_V2(nextGenGenoTable_GM_List_Reverse[[1]],genoSelectionTable_Parameter[[1]],Policy,nFamilies,no_QTL)
	
	    emigrantGroups <- migrationGroups[[1]]
	    immigrantGroups <- migrationGroups[[2]]
	
	
######################################################################################   
## get correct table format for exchange genodata functions
       

        Cycle_Progeny_F5_Fam_ParamList <- rep(list(array(0,c(nFamilies,nMarkers,2,cycle1_nProgeny))),nParameters)

        Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)),nParameters)
		
		selectedLineIndices_Reverse <- rep(list(rep(list(c()),nFamilies)),nParameters)
             

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
				selectedLineIndices_Reverse[[nParameter]][[nFamily]] <-  selectedLineIndices[[nFamily]][[nParameter]]
			} 
		}
	
################
### Exchange selected geno data among pairs ##############
###############################

    if(migFreq !=0){
   
 	  if(nCyc>=2 && nCyc%%migFreq==0){
	  
			output_List_afterExchange <- foreach(k=1:nParameters) %dopar% (LTSGenoIMV2::exchangeGenoData_GM(Cycle_Progeny_F5_Fam_CrossList[[k]],SelCriterion_List_Reverse[[k]],migrationSize,emigrantGroups,immigrantGroups,direction,Policy))
	    		migEvent <- TRUE
	      
###########################################################################
### Get CycleProgeny Genodata and genotype table after exchange


		Cycle_Progeny_F5_Fam_CrossList <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)
        
		nProgInit <- 1
		nProgFinal <-  dim(output_List_afterExchange[[1]][[1]][[1]])[3]
				
		for(nFamily in 1:length(families)){ 
		
			Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]][,,nProgInit:nProgFinal] <- output_List_afterExchange[[nParameter]][[1]][[nFamily]][,,nProgInit:nProgFinal]
		}
	  
	        familyInfo <- output_List_afterExchange[[1]][[2]]
			
	    	Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	   } 

          else if(nCyc>=2 && nCyc%%migFreq!=0){ 
              migEvent <- FALSE

          }
   
    }


 if(migEvent==TRUE){
### Get selected nextGenGenoTable_List after exchange
   
    nProgeny_per_Cross <- 1
	
    nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies) %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_Fam_CrossList[[nFamily]][[nParameter]],nProgeny_in_Family,nProgeny_per_Cross,no_QTL,varE))

   
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
	  
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
	 }
    }
  
 } else if(migEvent==FALSE){ 

    nProgeny_per_Cross <- 1
	
    nextGenGenoTableDataList <- foreach(nFamily =1:nFamilies) %:% foreach(nParameter=1:nParameters) %dopar% (LTSGenoIMV2::getNextGenGenoTableData(Cycle_Progeny_F5_Fam_CrossList[[nParameter]][[nFamily]],nProgeny_in_Family,nProgeny_per_Cross,no_QTL,varE))

   
  nextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  nextGenGenoTableMod_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
  newNextGenGenoTable_GM_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
	  
	
    for(nFamily in 1:nFamilies){
	
	 for(nParameter in 1:nParameters){
    
		nextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-  nextGenGenoTableDataList[[nFamily]][[nParameter]][[1]]
		nextGenGenoTableMod_GM_List[[nFamily]][[nParameter]] <- nextGenGenoTableDataList[[nFamily]][[nParameter]][[2]]
		newNextGenGenoTable_GM_List[[nFamily]][[nParameter]] <-nextGenGenoTableDataList[[nFamily]][[nParameter]][[3]]
		
	 }
    }
  
}
  
#################### GM method for isolated families  %:% foreach(j=1:nParameters) 

 j <- 1

 system.time({

	genoSelected_List <- foreach(k=1:nFamilies,.export=c("selectionOnGeno","nFamilies","modelType","GM_Frontier_Param_List")) %dopar% tryCatch(LTSGenoIMV2::getIM_GM_Frontier_selectedGenoList_SingleFams_V2(nextGenGenoTable_GM_List[[k]][[j]],nextGenGenoTableMod_GM_List[[k]][[j]],PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,k),error= function (e) print(paste("Error-Family-",k,e)))

        
 })
#########################################3
     

   save.image(WorkspaceName)
 
### 
## Reduce Parameter set to 30   

   for(nFamily in 1:nFamilies){
   
     nParameters <- length(genoSelected_List[[nFamily]][[1]])
	 
     for(nParameter in 1:nParameters){
		selectedParents_Geno_List[[nCyc]][[nFamily]][[nParameter]] <- genoSelected_List[[nFamily]][[1]][[nParameter]]
		CriterionValue_List[[nCyc]][[nFamily]][[nParameter]] <- genoSelected_List[[nFamily]][[2]][[nParameter]]
		
      }
   }


####


    Cycle_GenoData_List <- Cycle_Progeny_F5_Fam_CrossList

    nParams_List <- list()
	for(k in 1:nFamilies){
	   nParams <- length(genoSelected_List[[k]][[1]])
	   nParams_List[[k]] <- nParams
	}
	
	nParams_List_Cyc[[nCyc]] <- nParams_List

## Get F5 RILs for selected parents in cycle 1
        
    F5RILs_GM_Fam <-  foreach(k=1:nFamilies) %:% foreach(j=1:nParams_List[[k]]) %dopar% (LTSGenoIMV2::getF5RILs_FS_IM_GM_Frontier(Cycle_GenoData_List[[k]][[j]],selectedParents_Geno_List[[nCyc]][[k]][[j]],selectionOnGeno,nMarkers,NAM_LinkMap_New))

       
    Cycle_Progeny_F5_List_Fam <- rep(list(rep(list(list()),nParameters)),nFamilies)	
    nCrosses_List <- rep(list(rep(list(list()),nParameters)),nFamilies)
	
    for(nFamily in 1:nFamilies){
	      nParameters <- length(F5RILs_GM_Fam[[nFamily]])
          for(nParameter in 1:nParameters){
          	Cycle_Progeny_F5_List_Fam[[nFamily]][[nParameter]] <- F5RILs_GM_Fam[[nFamily]][[nParameter]][[1]]
          	nCrosses_List[[nFamily]][[nParameter]] <- F5RILs_GM_Fam[[nFamily]][[nParameter]][[2]]
 	
         }
    }

    gc()
   nParameters <- 1
   Cycle_Progeny_F5_List <- Cycle_Progeny_F5_List_Fam
   nProgeny_in_Family <- 100
   nCrosses_in_Family <- nSel_inFamily
   nProgeny_per_Cross <- nProgeny_in_Family/nCrosses_in_Family
   nProgeny_per_Family <- nCrosses_in_Family*nProgeny_per_Cross
   
    
   Cycle_Progeny_F5_List_Fam_V2 <- rep(list(rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nParameters)),nFamilies)

   nProgeny_per_Family <- nCrosses_in_Family*nProgeny_per_Cross

   for(nFamily in 1:nFamilies) {

 	for(nParameter in 1:nParameters){ 
             init <- 1
             final <- nProgeny_per_Cross
            
            for(nCrs in 1:nCrosses_in_Family ){
                    
                         

                    Cycle_Progeny_F5_List_Fam_V2[[nFamily]][[nParameter]][,,init:final] <- Cycle_Progeny_F5_List[[nFamily]][[nParameter]][[nCrs]][1,,,1:nProgeny_per_Cross]

                    init <- final+1
                    final <- final+nProgeny_per_Cross
            } 
      } 
   }
   	
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
	print(paste(nCyc,"-Fam20-",mean(unlist(lapply(GenoVal_Sim_NX_2k_3c_List_Cycle[[nCyc]][[20]][[1]],mean))),sep=""))
		 
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
# ######
	 simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	 return(simResults_List)

 }


######## GM_Frontier Method on Admixed Populations


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

####### GM_Frontier Method on Admixed Populations - Two part code to get output in early cycles

### Part 1

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
	nParameters <- 30
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
### Part 2

runSimulations20X_BD_Wght_GM_Frontier_V3A_Resume_Split <- function(Initial_Input_Variable_List,Input_Data_List,simResults_List,StartCycle,FinishCycle){

	  options(warn=-1) 	
### Assign Variables #########################################################################
### Selection Parameters 

      no_QTL <- Initial_Input_Variable_List[[1]]
	  h2 <- Initial_Input_Variable_List[[2]]
  	  no_selected <- Initial_Input_Variable_List[[3]]
	  nCrosses <- Initial_Input_Variable_List[[4]]
	  nProgeny <- Initial_Input_Variable_List[[5]]
	  nCycles <- Initial_Input_Variable_List[[6]]
	  
	  cycle1_nProgeny <- nProgeny
  
  	  selectionOnGeno <- Initial_Input_Variable_List[[7]]
	  selectionOnSimulated <- Initial_Input_Variable_List[[8]]
	   
	  NAM_LinkMap_New <- Initial_Input_Variable_List[[9]]
  
##  Model Update Parameters
	  
      modelRetrain <- Initial_Input_Variable_List[[10]]
	  retrainFrequency <-  Initial_Input_Variable_List[[11]]
	  modelType <- Initial_Input_Variable_List[[12]]
	  modelUpdate <- Initial_Input_Variable_List[[13]]
      updateFrequency <- Initial_Input_Variable_List[[14]]
	  updateType <- Initial_Input_Variable_List[[15]]

      nMarkers <- 4289
	  condition <- Initial_Input_Variable_List[[16]]
	  Rep <- Initial_Input_Variable_List[[17]] 
 
## Weighting Parameters	  
	  Weighted <- Initial_Input_Variable_List[[18]] 
	  WghtMethod <- Initial_Input_Variable_List[[19]] 
   
## BD Parameter	  
	  BD <- Initial_Input_Variable_List[[20]] 

## FormatConversionTable and trainTableWindowSize
	  AlleleConversionTable_Combined <- Initial_Input_Variable_List[[21]]
	  trainTableWindowSize <- Initial_Input_Variable_List[[22]]
	  
	  nFamilies <- Initial_Input_Variable_List[[23]]
  
## GM Parameters 
	  
	  GM_Param_List <- Initial_Input_Variable_List[[24]]
	  nCores_ForEach <- Initial_Input_Variable_List[[25]]
	  
	  no_cores <- GM_Param_List$no_cores
	  noCores_ForEach
	  WorkspaceName <- Initial_Input_Variable_List[[26]]
	  nCores_ForEach <- noCores_ForEach 
	  nCores <- nCores_ForEach 
	  
### Input Data 
  
	  
	  Cycle_Progeny_F5_List <- Input_Data_List[[1]]
	  nCrosses_List <- Input_Data_List[[2]]
	  nProgeny_per_Cross <- Input_Data_List[[3]]
	  PredictionModel_Geno_List <- Input_Data_List[[4]]
	  PredictionModel_Pheno_List <- Input_Data_List[[5]]
	  trainGeno_backingFileName_List <- Input_Data_List[[6]]
	  trainSimPhenoTablePreCycle_List <- Input_Data_List[[7]]
	  trainSimGenoValTablePreCycle_List <- Input_Data_List[[8]]
	  varE <- Input_Data_List[[9]]
	  CombinedTableCount <- Input_Data_List[[10]]
  

	  nLines <- nFamilies * cycle1_nProgeny 
      nProgeny_per_Cross <- nLines/nCrosses		
  
	 
## Initiate Output list 

		
      GenoVal_Sim_NX_2k_3c_List_Cycle <-  simResults_List[[1]]	
	  GenoVal_NX_2k_3c_List_Cycle <- simResults_List[[2]]
	  GenoVal_NX_N_3c_List_Cycle <- simResults_List[[3]]
	  
	  PhenoVal_Sim_NX_2k_3c_List_Cycle <- simResults_List[[4]]
	  PhenoVal_NX_2k_3c_List_Cycle <- simResults_List[[5]] 
	  PhenoVal_NX_N_3c_List_Cycle <- simResults_List[[6]]
	       
      attainedGenoValues_List_Cycle <- simResults_List[[7]]
	  attainedPhenoValues_List_Cycle <- simResults_List[[8]]
	  
	  topNGenoValues_List_Cycle <- simResults_List[[9]]
      GenoTableIndices_topN_List_Cycle <- simResults_List[[10]]

      topNPhenoValues_List_Cycle <- simResults_List[[11]]
      PhenoTableIndices_topN_List_Cycle <- simResults_List[[12]]
	  
	  selectedParents_Geno_List <- simResults_List[[13]]
	  CriterionValue_List <- simResults_List[[14]]
	
	 
#################################################################################################
#################################################################################################
	  

    startCycle <- StartCycle
	finishCycle <- FinishCycle
    
    
    nextGenGenoTable_GM_List <- list()
    nextGenGenoTableMod_GM_List <- list()
    newNextGenGenoTable_GM_List <- list()
	
	genoValues_List <- list()
    phenoValues_List <- list()
	genoValSimValues_List <- list()
	phenoValSimValues_List <- list()
				   
	trainGenoNewTablePreCycle_List <- list()
	trainGenoNewTablePreCycle_Mat_List <- list()
	
	trainGeno_descriptorFileName_List <- list()	
 
    selectedNextGenGenoTable_GM_List <- list() 
    selectedNextGenGenoTableMod_GM_List <- list() 

for(nCyc in startCycle:finishCycle){
      
	print(nCyc)
	
	nParameters <- length(Cycle_Progeny_F5_List)

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

    
  if(nCyc > startCycle){
     save(genoSelected_List,simResults_List_Cycle,file=WorkSpaceName)
	 rm(simResults_List_Cycle)
	 gc()
  }
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
       Input_Data_List <- list(Cycle_Progeny_F5_List,nCrosses_List,PredictionModel_Geno_List,PredictionModel_Pheno_List,trainGeno_backingFileName_List,trainSimPhenoTablePreCycle_List,trainSimGenoValTablePreCycle_List)
	   
	  
	   simResults_List <- list(GenoVal_Sim_NX_2k_3c_List_Cycle,GenoVal_NX_2k_3c_List_Cycle,GenoVal_NX_N_3c_List_Cycle,PhenoVal_Sim_NX_2k_3c_List_Cycle,PhenoVal_NX_2k_3c_List_Cycle,PhenoVal_NX_N_3c_List_Cycle,attainedGenoValues_List_Cycle,
	   attainedPhenoValues_List_Cycle,topNGenoValues_List_Cycle,GenoTableIndices_topN_List_Cycle,topNPhenoValues_List_Cycle,PhenoTableIndices_topN_List_Cycle,selectedParents_Geno_List,CriterionValue_List)

	   return(list(Input_Data_List,simResults_List))

}


######






