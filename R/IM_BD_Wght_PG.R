
runSimulations20X_IslandSelection_BD_Wght_PG <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,MigrationSize,policy,migrationFrequency,policyChangeFrequency,Direction,NAM_LinkMap,ModelRetrain,RetrainFrequency,ModelType,ModelUpdate,UpdateFrequency,UpdateType,i,k,weighted,wghtMethod,BreedDesign,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,SclassPlot,SelCriteria,gIndStats){


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

      selCriteria <- SelCriteria

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
	  
	  j<-1
 	  CombinedTableCount <- 1
	  trainTableWindowSize <- 10
	  
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
	  GIndStats <- gIndStats

      AlleleConversionTable_Combined <-  alleleConversionTable_Combined
	 
### Initialize variables
	 
	  Migration_List <- list()
	  PCA_List <- list()
	  j<-1


	  
	  nCyc <-1
	 	  
	  trainGenoNewTablePreCycle <- (trainGenoNewTable)
      trainSimPhenoTablePreCycle <- trainSimPhenoValTable
      trainSimGenoValTablePreCycle <- trainSimGenoValTable
	  
### PG_ Variables list 

      percentHeterozygous_List<- rep(list(rep(0,nCycles)),nFamilies)
	  percentHeterozygous_QTL_List<- rep(list(rep(0,nCycles)),nFamilies)
	  FavorableAlleleFreq <- rep(list(rep(list(list()),nCycles)),nFamilies)
	  FavorableAllele <- rep(list(rep(list(list()),nCycles)),nFamilies)
	  FavorableQTLAlleleFreq <- rep(list(rep(list(list()),nCycles)),nFamilies)
	  FavorableQTLAllele <- rep(list(rep(list(list()),nCycles)),nFamilies)
	  percentHet <- rep(0,nCycles)

	  QTL_GenoTable_List <- list()


	  Diff_Stats_List_Global <- list()
	  Diff_Stats_List_Local <- rep(list(rep(list(list()),nFamilies)),nCycles)

 
	  
		
#### /nSel_inFamily
###############################################################################################
## Get predicted geno and pheno values for cycle1 for combined population set

	cycle1GenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
	newCycle1GenoTable <- generateWeightedGenoTable(cycle1GenoTable,no_QTL)
	nIndividuals <- nrow(cycle1GenoTable)
	
	percentHet[nCyc] <- getPercentHeterozygous(cycle1GenoTable) 
	
	
	
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

		F <- list() 
		ExpHet <- list()
		Freq_GInd <- list() 
		MAF <- list() 
### get island PG parameters

        cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable,AlleleConversionTable_Combined)
 	    cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")


		# nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		# nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
     if(GIndStats == TRUE){
		Diff_Stats <- diff_stats(cycle1GenoTable_GInd)
	    Diff_Stats_List_Local[[nCyc]] <- as.big.matrix(Diff_Stats[[1]])
	    Diff_Stats_List_Global[[nCyc]] <- Diff_Stats[[2]]
		
		F[[nCyc]] <- inbreeding(cycle1GenoTable_GInd,res.type="estimate")
        ExpHet[[nCyc]] <- Hs(cycle1GenoTable_GInd)
		Freq_GInd[[nCyc]] <- makefreq(cycle1GenoTable_GInd)
		MAF[[nCyc]] <- minorAllele(cycle1GenoTable_GInd)
	  }
		
		cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
		write.table(cycle1GenoData_DF,paste("IM_",BD,"_",selCriteria,"_",migFreq,"_",Policy,"_GenoDataInd_Rep",Rep,"_Cyc_",nCyc,sep=""),sep="\t")
		
		rm(cycle1GenoTable_GInd)
		rm(Diff_Stats)
		gc()

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
    
	  selectedGenoIndividualIndices_List_Cyc <- rep(list(rep(list(rep(0,no_selected)),nFamilies)),nCycles)


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
	
	
#,no_QTL
	if(Policy != "GeneticClusters"){
	migrationGroups <- migrationPolicy(cycle1GenoTable,genoSelectionTable,Policy,nFamilies,no_QTL)
	
	emigrantGroups <- migrationGroups[[1]]
	immigrantGroups <- migrationGroups[[2]]
	
	}else if (Policy == "GeneticClusters"){ 
	
	migrationGroups <- migrationPolicy(cycle1GenoTable,genoSelectionTable,Policy,nFamilies,no_QTL)
	
	emigrantGroups <- migrationGroups[[1]]
	immigrantGroups <- migrationGroups[[2]]
	emigrantGroup_Prob <- migrationGroups[[3]]
	
	}
	

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

				
#		gind2genalex(cycle1GenoTable_GInd,"cycle1GenoTable_GInd.csv",geo=FALSE)	
###########################################################################
# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################
   
    for(nCyc in 2:nCycles){
	   
	   #nCyc <- 2
		print(nCyc)

        selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		for(nFamily in 1:nFamilies){

			if(nCyc==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
					
					selectedGenoIndividualIndices_List_Cyc[[(nCyc-1)]] <- selectedGenoIndividualIndices_List
					
					
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
					
					selectedGenoIndividualIndices_List_Cyc[[(nCyc-1)]] <- selectedPhenoIndividualIndices_List
					
				}
	
			}

			if(nCyc>2){

		
				selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				
				selectedGenoIndividualIndices_List_Cyc[[(nCyc-1)]] <- selectedGenoIndividualIndices_List
				
						 
		    }
	    }


	
### Exchange selected geno data among pairs
    if(migFreq !=0){
      if(nCyc>=2 && nCyc%%migFreq==0){
	   
	   	output_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,migrationSize,emigrantGroups,immigrantGroups,direction,Policy)
	    		
		selectedGenoData_List_afterExchange <- (output_List_afterExchange)[[1]]
		familyInfo <- output_List_afterExchange[[2]]
		Migration_List[[nCyc]] <- list(emigrantGroups,immigrantGroups,familyInfo)
		
	  }else if(nCyc>=2 && nCyc%% migFreq !=0){
	  
		selectedGenoData_List_afterExchange <- selectedGenoData_List 
		
	  }
    }
	if(migFreq ==0){
		selectedGenoData_List_afterExchange <- selectedGenoData_List 
	}
### Generate F5 RIL progeny for selected geno data
 
	
	F5RILs <- getF5RILs_BD(BD,selectedGenoData_List_afterExchange,nFamilies,nSel_inFamily,nProgeny,nMarkers,NAM_LinkMap_New,Rep)
		
	Cycle_Progeny_F5 <- F5RILs[[1]]
	Cycle_Progeny_F5_List <- F5RILs[[2]]
	
	

 # dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

####################################################################################################################################################################################################
### This section works with full data table

 
		nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nFamilies,nCrosses_inFamily*nProgeny)
		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
		
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
		
			
		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

		initIndex <- 1
		finalIndex <- cycle1_nProgeny
		for(nFamily in 1:nFamilies){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)
			
			QTL_GenoTable_List[[nFamily]] <- getQTLTable(nextGenGenoTable_List[[nFamily]],no_QTL)
			nextGenGenoTable_AlleleFormat_Local <- getAlleleFormat(nextGenGenoTable_List[[nFamily]],AlleleConversionTable_Combined)
 		    		
			
		    Freq_Allele_List <-	getFreq_BasePoln(nextGenGenoTable_List[[nFamily]],MarkerEffects)
		    FavorableAlleleFreq[[nFamily]][[nCyc]] <- Freq_Allele_List[[1]]
		    FavorableAllele[[nFamily]][[nCyc]]  <- Freq_Allele_List[[2]]
		    percentHeterozygous_List[[nFamily]][nCyc] <- getPercentHeterozygous(nextGenGenoTable_List[[nFamily]])


		    Freq_QTLAllele_List <-	getFreq_BasePoln(QTL_GenoTable_List[[nFamily]],MarkerEffects)
		    FavorableQTLAlleleFreq[[nFamily]][[nCyc]] <- Freq_QTLAllele_List[[1]]
		    FavorableQTLAllele[[nFamily]][[nCyc]]  <- Freq_QTLAllele_List[[2]]
		    percentHeterozygous_QTL_List[[nFamily]][nCyc] <- getPercentHeterozygous(QTL_GenoTable_List[[nFamily]])
			
			initIndex <- finalIndex + 1
		    finalIndex <- finalIndex + cycle1_nProgeny
		
		}

		genoSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)
	
		
### get island PG parameters

      	nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")
      
	   if(GIndStats == TRUE){

		Diff_Stats <- diff_stats(nextGenGenoTable_GInd)
	    Diff_Stats_List_Local[[nCyc]] <- as.big.matrix(Diff_Stats[[1]])
	    Diff_Stats_List_Global[[nCyc]] <- Diff_Stats[[2]]
		F[[nCyc]] <- inbreeding(nextGenGenoTable_GInd,res.type="estimate")
        ExpHet[[nCyc]] <- Hs(nextGenGenoTable_GInd)
		Freq_GInd[[nCyc]] <- makefreq(nextGenGenoTable_GInd)
		MAF[[nCyc]] <- minorAllele(nextGenGenoTable_GInd)
	   }
		
		nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
		write.table(nextGenGenoData_DF,paste("IM_",BD,"_",selCriteria,"_",migFreq,"_",Policy,"_GenoDataInd_Rep",Rep,"_Cyc_",nCyc,sep=""),sep="\t")
				
		rm(nextGenGenoTable_GInd)
		rm(Diff_Stats)
		gc()
		
if(sclassPlot==TRUE){		
		
	if(nCyc%%5 == 0 && sd(genoSimValues !=0)){	

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

if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	        if(nCyc >2){
			   trainGenoNewTablePreCycle <- read.table(trainGeno_FileName,sep="\t")
			}

			PredictionModels <- LTSGenoIMV2::getUpdatedGSModel(newNextGenGenoTable,genoSimValues,phenoSimValues,PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,trainTableWindowSize,as.matrix(trainGenoNewTablePreCycle),trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount,nCyc)
 
			PredictionModel_Geno <- PredictionModels[[1]]
			PredictionModel_Pheno <- PredictionModels[[2]]
			CombinedTableCount <- PredictionModels[[3]]
			trainGenoNewTablePreCycle <- PredictionModels[[4]]
			trainSimPhenoTablePreCycle <- PredictionModels[[5]]
			trainSimGenoValTablePreCycle <- PredictionModels[[6]]
			
			
			
			trainGeno_FileName <- paste("trainTable_",selCriteria,"_",BD,"_",Policy,"_",condition,"_",Rep,"_", ".txt",sep="") 
			write.table(trainGenoNewTablePreCycle,trainGeno_FileName,sep="\t") 
			rm(trainGenoNewTablePreCycle)
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

	   
	Criterion_List <- getCriterionValueList(BD,genoValues,phenoValues,genoSimValues,phenoSimValues,nCrosses_inFamily,nProgeny,selectionOnGeno,selectionOnSimulated,nFamilies)

	genoValues_List <- Criterion_List[[1]]
	phenoValues_List <- Criterion_List[[2]]
	genoSimValues_List <- Criterion_List[[3]]
	phenoSimValues_List <- Criterion_List[[4]] 
	 
### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
		    if(selectionOnGeno ==TRUE){
		
				genoSelectionTable <- ApplySelection_Discrete(genoSimValues_List,no_selected)
		    
			}else if(selectionOnGeno ==FALSE){ 

				genoSelectionTable <- ApplySelection_Discrete(phenoSimValues_List,no_selected)
            }			
			
		}
		if (selectionOnSimulated==FALSE){
			
			if(selectionOnGeno ==TRUE){
			
				genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			}else if(selectionOnGeno ==FALSE){
			
				 genoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
			}
		}

		if(nCyc>=2 && nCyc%%migPolicyFreq==0 && sd(genoSimValues !=0)){
			
	
			if(Policy != "GeneticClusters"){
	           migrationGroups <- migrationPolicy(nextGenGenoTable,genoSelectionTable,Policy,nFamilies,no_QTL)
	           emigrantGroups <- migrationGroups[[1]]
	           immigrantGroups <- migrationGroups[[2]]
	
	        }else if(Policy == "GeneticClusters"){ 
	         
			   migrationGroups_Prev <- migrationGroups
	           migrationGroups <- migrationPolicy(nextGenGenoTable,genoSelectionTable,Policy,nFamilies,no_QTL)
	 	       emigrantGroups <- migrationGroups[[1]]
	           immigrantGroups <- migrationGroups[[2]]
	           emigrantGroup_Prob <- migrationGroups[[3]]
	           RetainMigPolicy <- migrationGroups[[4]] 
			 
			   if(RetainMigPolicy == TRUE){ 
			     migrationGroups <- migrationGroups_Prev
			   }
	 		}
				
		}
		
#####################################################################################################################
## Assign output variables for cycle n

	for(nFamily in 1:nFamilies){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,nCyc] <- genoSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,nCyc] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,nCyc]<- phenoSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,nCyc] <- phenoValues_List[[nFamily]]

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
			genoSimSelectedValues <- genoSimValues_List[[nFamily]][selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoSimValues_List[[nFamily]][selectedGenoIndividualIndices]
				
			selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
			attainedGenoValues_List[[nFamily]][nCyc] <- max(genoSimSelectedValues)
			
			
			if(selectionOnSimulated ==TRUE){
		        if(selectionOnGeno ==TRUE){
		
				GenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- as.vector(genoSelectionTable[[nFamily]][,1])
			    PhenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- phenoSimValues_List[[nFamily]][selectedGenoIndividualIndices]
				
			   } else if(selectionOnGeno ==FALSE){ 

				GenoVal_NX_N_3c_List[[nFamily]][,nCyc] <- genoSimValues_List[[nFamily]][selectedGenoIndividualIndices]
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
    
	if(GIndStats ==TRUE){
	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List,PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List,selectedGenoIndividualIndices_List_Cyc,QTL_GenoTable_List,FavorableAllele,FavorableAlleleFreq,percentHeterozygous_List,Freq_QTLAllele_List,FavorableQTLAlleleFreq,percentHeterozygous_QTL_List,Diff_Stats_List_Local,Diff_Stats_List_Global,F,ExpHet,Freq_GInd,MAF)
    }else if(GIndStats ==FALSE){
	
	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List,PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List,selectedGenoIndividualIndices_List_Cyc,QTL_GenoTable_List,FavorableAllele,FavorableAlleleFreq,percentHeterozygous_List,FavorableQTLAllele,FavorableQTLAlleleFreq,percentHeterozygous_QTL_List)
	
	}
	 

	  return(simResults_List)
  
 }
 
			
####


runSimulations20X_RRBLUP_WGS_V3_PG <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,ModelRetrain,RetrainFrequency,ModelUpdate,UpdateType,UpdateFrequency,i,k,ModelType,weighted,weightingMethod,AlphaShape,BetaShape,TimeHorizon,alleleConversionTable_Combined,trainGenoNewTable,trainSimGenoValTable,trainSimPhenoValTable,BreedingDesign,NAM_LinkMap,SelCriteria,gIndStats){


### Assign Variables

	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues

      selCriteria <- SelCriteria 
	  GIndStats <- gIndStats
####################################################

	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList

	  varE <- varEList

	  condition <- i
	  Rep <- k
	  modelType <- ModelType
	  modelRetrain <- ModelRetrain
	  updateFrequency <- UpdateFrequency
	  retrainFrequency <- RetrainFrequency
	  updateType <- UpdateType
      modelUpdate <- ModelUpdate

	  nFamilies <- 20
      cycle1_nProgeny <- 100
	  nMarkers <- 4289
      nIndividuals <- 2000

      nProgeny <- nIndividuals/no_selected 

      Weighted <- weighted
	  WghtMethod <- weightingMethod  
	  alphaShape <- AlphaShape
	  betaShape <- BetaShape
	  timeHorizon <- TimeHorizon
	  alphaPar_Vector <- rep(0,nCycles)
      BD <- BreedingDesign
         
	  NAM_LinkMap_New <- NAM_LinkMap 

       trainGenoNewTablePreCycle <- (trainGenoNewTable)
       trainSimPhenoTablePreCycle <- trainSimPhenoValTable
       trainSimGenoValTablePreCycle <- trainSimGenoValTable
	   trainTableWindowSize <- 10
	   AlleleConversionTable_Combined <-  alleleConversionTable_Combined 
	   
	   j<-1
 	   CombinedTableCount <- 1
			  
### PG_ Variables list 

      percentHeterozygous_List<- rep(0,nCycles)
	  percentHeterozygous_QTL_List<- rep(0,nCycles)
	  
	  FavorableAlleleFreq <- rep(list(list()),nCycles)
	  FavorableAllele <- rep(list(list()),nCycles)
	  FavorableQTLAlleleFreq <-rep(list(list()),nCycles)
	  FavorableQTLAllele <- rep(list(list()),nCycles)
	  percentHet <- rep(0,nCycles)

	  QTL_GenoTable_List <- list()


	  Diff_Stats_List_Global <- list()
	  Diff_Stats_List_Local <- rep(list(list()),nCycles)

#########################################################################################################
	  nCyc <- 1
      print(paste("cycle_No",nCyc))
	
####### Predict geno and pheno values ###################################################################

	  cycle1GenoTable <- generateMinus101GenoFormat(cycle1_Geno_Data,nFamilies,cycle1_nProgeny,nFamilies,no_selected)
	  newCycle1GenoTable <- generateWeightedGenoTable(cycle1GenoTable,no_QTL)

	  if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[1])
	  }else if(selectionOnGeno==FALSE){
	    	MarkerEffects <- unlist(PredictionModel_Pheno[1])
	  }

      Freq_List <- getFreq_BasePoln(cycle1GenoTable,MarkerEffects)
	  Freq <- Freq_List[[1]]
	  FavAllele <- Freq_List[[2]]
	  
	  nextGenGenoTable <- cycle1GenoTable
	  QTL_GenoTable <- getQTLTable(nextGenGenoTable,no_QTL)

	  Freq_Allele <-	getFreq_BasePoln(nextGenGenoTable,MarkerEffects)
	  FavorableAlleleFreq[[nCyc]] <- Freq_Allele[[1]]
	  FavorableAllele[[nCyc]]  <- Freq_Allele[[2]]
	  percentHeterozygous_List[nCyc] <- getPercentHeterozygous(nextGenGenoTable)

	  Freq_QTLAllele <-	getFreq_BasePoln(QTL_GenoTable,MarkerEffects)
	  FavorableQTLAlleleFreq[[nCyc]] <- Freq_QTLAllele[[1]]
	  FavorableQTLAllele[[nCyc]]  <- Freq_QTLAllele[[2]]
	  percentHeterozygous_QTL_List[nCyc] <- getPercentHeterozygous(QTL_GenoTable)
	 
###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(pop in 1:20){

			final <- pop*100
			populations[init:final]<- rep(pop,100)
			init <- final+1
		}
		
      cycle1GenoTable_AlleleFormat <- getAlleleFormat(cycle1GenoTable,AlleleConversionTable_Combined)
 	  cycle1GenoTable_GInd <- df2genind(cycle1GenoTable_AlleleFormat,pop=populations,sep="")
		
	  cycle1GenoData_DF <- genind2df(cycle1GenoTable_GInd)
	  write.table(cycle1GenoData_DF,paste("PM_",selCriteria,"_",BD,"_GenoDataInd_Con_",condition,"_Rep_",Rep,"_Cyc_",nCyc,sep=""),sep="\t")
	  	
##### gIndStats
###### Predict geno and pheno values ###################################################################

        if(Weighted==TRUE && WghtMethod=="JM"){
			genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newCycle1GenoTable,PredictionModel_Geno,Freq)
			phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newCycle1GenoTable,PredictionModel_Pheno,Freq)
        
		}else if (Weighted==TRUE && WghtMethod =="DW"){
		
		    cycleNumber <- nCyc
		    genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newCycle1GenoTable,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)
			phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newCycle1GenoTable,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

			genoValues <- genoValues_PredList[[1]]
			phenoValues <- phenoValues_PredList[[1]]
			alphaPar_Vector[i] <- genoValues_PredList[[2]]
		
		}else if(Weighted==FALSE){
            genoValues  <- PredictRRBLUPPhenoValues(newCycle1GenoTable,PredictionModel_Geno)
            phenoValues <- PredictRRBLUPPhenoValues(newCycle1GenoTable,PredictionModel_Pheno)
        }

	
### Initiate Arrays and Variables  #######################################################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nFamilies*cycle1_nProgeny)*nCycles),nrow=(nFamilies*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_NX_N_3c <- matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nFamilies*cycle1_nProgeny)*nCycles),nrow=(nFamilies*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nFamilies*cycle1_nProgeny)*nCycles),nrow=(nFamilies*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nFamilies*cycle1_nProgeny)*nCycles),nrow=(nFamilies*cycle1_nProgeny),ncol=nCycles)

      attainedGenoValues <- rep(0,nCycles)
	 
	  selectedGenoIndividualIndices_List <- list()
      selectedPhenoIndividualIndices_List <- list()

######## Functions on Cycle 1 Geno Data ################################################################################

	if(selectionOnSimulated ==FALSE){
	
		if(selectionOnGeno==TRUE){
	
		  sortedGenoValues <- sort(genoValues,decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:no_selected]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:no_selected]
          topNPhenoValues <- phenoValues[GenoTableIndices_topN] 
		  topNGenoSimValues <- genoSimValues[GenoTableIndices_topN]
		  PhenoTableIndices_topN <- GenoTableIndices_topN
		}else if (selectionOnGeno==FALSE) {
		  	  
		  sortedPhenoValues <- sort(phenoValues,decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:no_selected]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:no_selected]
		  topNGenoValues <- genoValues[PhenoTableIndices_topN]
		  topNGenoSimValues <- genoSimValues[PhenoTableIndices_topN]
		  GenoTableIndices_topN <- PhenoTableIndices_topN
		} 
  
    }
	if(selectionOnSimulated ==TRUE){
	
		if(selectionOnGeno==TRUE){
		 
   		  sortedGenoValues <- sort(genoSimValues,decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:no_selected]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:no_selected]
		  topNPhenoValues <- phenoValues[GenoTableIndices_topN] 
		  topNGenoSimValues <- genoSimValues[GenoTableIndices_topN]
		  PhenoTableIndices_topN <- GenoTableIndices_topN
		
		}else if(selectionOnGeno==FALSE){

		  sortedPhenoValues <- sort(phenoSimValues,decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:no_selected]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:no_selected]
		  topNGenoValues <- genoValues[PhenoTableIndices_topN]
		  topNGenoSimValues <- genoSimValues[PhenoTableIndices_topN]
		  GenoTableIndices_topN <- PhenoTableIndices_topN
	    }
	}
		        
##################################################################################################################
## Assign output variables for cycle 1

	  GenoVal_Sim_NX_2k_3c[,1] <- genoSimValues
	  GenoVal_NX_2k_3c[,1] <- genoValues
	  GenoVal_NX_N_3c[,1] <- topNGenoValues
    
	  attainedGenoValues[1] <- max(topNGenoSimValues)
	  
	  selectedGenoIndividualIndices <- GenoTableIndices_topN
      selectedPhenoIndividualIndices <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c[,1] <- phenoSimValues
	  PhenoVal_NX_2k_3c[,1] <- phenoValues
	  PhenoVal_NX_N_3c[,1] <- topNPhenoValues

      selectedGenoIndividualIndices_List[[nCyc]] <- selectedGenoIndividualIndices
      selectedPhenoIndividualIndices_List[[nCyc]] <- selectedPhenoIndividualIndices
	  

      rm(cycle1GenoTable)
	  rm(newCycle1GenoTable)

############### Cycles 2*29 ###################################################################################
### extract genotype data of selected individuals
  
  for(nCyc in 2:nCycles){

    print(paste("cycleNo-",nCyc))
    
    if(nCyc==2){

		if(selectionOnGeno == TRUE){
			selectedGenoData <- extractGenoData(selectedGenoIndividualIndices,cycle1_Geno_Data,no_selected)
		}else if(selectionOnGeno == FALSE){
			selectedGenoData <- extractGenoData(selectedPhenoIndividualIndices,cycle1_Geno_Data,no_selected)
		}

    }else if(nCyc>2){

		if(selectionOnGeno == TRUE){

			selectedGenoData <- extractGenoData(selectedGenoIndividualIndices,Cycle_Progeny_F5,no_selected)
		}else if(selectionOnGeno == FALSE){
			selectedGenoData <- extractGenoData(selectedPhenoIndividualIndices,Cycle_Progeny_F5,no_selected)
		}

    }

	
	Cycle_Progeny_F5 <- getF5RILs_BD_WGS(BD,selectedGenoData,no_selected,nProgeny,nMarkers,NAM_LinkMap_New)
		
########################
## nextGenGenoTable <- generate012GenoFormat(Cycle_Progeny_F5,no_selected)
## generateMinus101GenoFormat <- function(Cycle_Progeny_Data,number_of_Crosses,n_Progeny,NFamilies,no_Selected)
    
    nextGenGenoTable <- generateMinus101GenoFormat_V1(Cycle_Progeny_F5,nCrosses,nProgeny) 
   
    Freq <- getFreq(nextGenGenoTable,FavAllele)
    genoSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
    phenoSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

    nextGenGenoData <- Cycle_Progeny_F5
    QTL_GenoTable <- getQTLTable(nextGenGenoTable,no_QTL)

    newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
	Freq_Allele <-	getFreq_BasePoln(nextGenGenoTable,MarkerEffects)
	FavorableAlleleFreq[[nCyc]] <- Freq_Allele[[1]]
	FavorableAllele[[nCyc]]  <- Freq_Allele[[2]]
	percentHeterozygous_List[nCyc] <- getPercentHeterozygous(nextGenGenoTable)

	Freq_QTLAllele <-	getFreq_BasePoln(QTL_GenoTable,MarkerEffects)
    FavorableQTLAlleleFreq[[nCyc]] <- Freq_QTLAllele[[1]]
	FavorableQTLAllele[[nCyc]]  <- Freq_QTLAllele[[2]]
	percentHeterozygous_QTL_List[nCyc] <- getPercentHeterozygous(QTL_GenoTable)
			
    nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
	print(dim(nextGenGenoTable)) 
	print(dim(nextGenGenoTable_AlleleFormat))
	
	 
###  Get population vector

		populations <- rep(0,20*100)
		init <- 1
		for(pop in 1:20){

			final <- pop*100
			populations[init:final]<- rep(pop,100)
			init <- final+1
		}
		
### Get GInd objects 	
	
	print(length(populations))
    nextGenGenoTable_GInd <- df2genind(as.data.frame(nextGenGenoTable_AlleleFormat),pop=populations,sep="")
    nextGenGenoData_DF <- genind2df(nextGenGenoTable_GInd)
	write.table(nextGenGenoData_DF,paste("PM_",selCriteria,"_",BD,"_GenoDataInd_Con_",condition,"_Rep_",Rep,"_Cyc_",nCyc,sep=""),sep="\t")

   
################################################
  gc()

	if(modelUpdate ==TRUE && nCyc%%updateFrequency==0){
	
	        if(nCyc >2){
			   trainGenoNewTablePreCycle <- read.table(trainGeno_FileName,sep="\t")
			}

			PredictionModels <- getUpdatedGSModel(newNextGenGenoTable,genoSimValues,phenoSimValues,PredictionModel_Geno,PredictionModel_Pheno,modelType,updateType,trainTableWindowSize,as.matrix(trainGenoNewTablePreCycle),trainSimPhenoTablePreCycle,trainSimGenoValTablePreCycle,CombinedTableCount,nCyc)
 
			PredictionModel_Geno <- PredictionModels[[1]]
			PredictionModel_Pheno <- PredictionModels[[2]]
			CombinedTableCount <- PredictionModels[[3]]
			trainGenoNewTablePreCycle <- PredictionModels[[4]]
			trainSimPhenoTablePreCycle <- PredictionModels[[5]]
			trainSimGenoValTablePreCycle <- PredictionModels[[6]]
			
			
			
			trainGeno_FileName <- paste("trainTable_",selCriteria,"_",BD,"_",condition,"_",Rep,"_", ".txt",sep="") 
			write.table(trainGenoNewTablePreCycle,trainGeno_FileName,sep="\t") 
			rm(trainGenoNewTablePreCycle)
	}


#################get predicted geno and pheno values with JM and unweighted GP models
  	
	
	if(Weighted==TRUE && WghtMethod =="JM"){
		genoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable,PredictionModel_Geno,Freq)
		phenoValues <- PredictRRBLUPPhenoValues_WGS_Jannink1(newNextGenGenoTable,PredictionModel_Pheno,Freq)
    }else if(Weighted ==TRUE && WghtMethod =="DW"){
	
		cycleNumber <- nCyc
		genoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable,PredictionModel_Geno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)
		phenoValues_PredList <- PredictRRBLUPPhenoValues_WGS_DWModel(newNextGenGenoTable,PredictionModel_Pheno,Freq,alphaShape,betaShape,cycleNumber,timeHorizon)

		genoValues <- genoValues_PredList[[1]]
		phenoValues <- phenoValues_PredList[[1]]

		alphaPar_Vector[i] <- genoValues_PredList[[2]]
		
	}else if(Weighted==FALSE){
		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)
    }
	
#############
    print(paste("GenoVal_Length",length(genoValues)))
    print(paste("PhenoVal_Length",length(phenoValues)))

    GenoVal_Sim_NX_2k_3c[,nCyc] <- genoSimValues
    GenoVal_NX_2k_3c[,nCyc] <- genoValues
   
    PhenoVal_Sim_NX_2k_3c[,nCyc]<- phenoSimValues
    PhenoVal_NX_2k_3c[,nCyc] <- phenoValues
	
###### 
   
### Selection based on simulated or predicted values

    if(selectionOnSimulated ==TRUE){

		genoSelectionTable <- ApplySelection(genoSimValues,no_selected)
		phenoSelectionTable <- ApplySelection(phenoSimValues,no_selected)


		if(selectionOnGeno==TRUE){

			selectedGenoIndividualIndices <- as.vector(genoSelectionTable[,2])
			genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
			phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]

		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
			phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
		}

		GenoVal_NX_N_3c[,nCyc] <- genoSimSelectedValues
	    PhenoVal_NX_N_3c[,nCyc] <- phenoSimSelectedValues
		attainedGenoValues[nCyc] <- max(genoSimSelectedValues)

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
		        genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
		
		}else if(selectionOnGeno==FALSE){

			selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[,2])
			genoSelectedValues <- genoValues[selectedPhenoIndividualIndices]
			phenoSelectedValues <- phenoValues[selectedPhenoIndividualIndices]
		    genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
		}
		    GenoVal_NX_N_3c[,nCyc]<- genoSelectedValues
	        PhenoVal_NX_N_3c[,nCyc]<- phenoSelectedValues
        	attainedGenoValues[nCyc] <- max(genoSimSelectedValues)
	}
	
	selectedGenoIndividualIndices_List[[nCyc]] <- selectedGenoIndividualIndices
	selectedPhenoIndividualIndices_List[[nCyc]] <- selectedPhenoIndividualIndices

   	rm(nextGenGenoTable)
	rm(newNextGenGenoTable)
    gc()
  }
  
  
   
  simResults_List <- list(GenoVal_Sim_NX_2k_3c,GenoVal_NX_2k_3c,GenoVal_NX_N_3c,PhenoVal_Sim_NX_2k_3c,PhenoVal_NX_2k_3c,PhenoVal_NX_N_3c,attainedGenoValues,selectedGenoIndividualIndices_List,selectedPhenoIndividualIndices_List,FavorableAlleleFreq,FavorableAllele,FavorableQTLAlleleFreq,FavorableQTLAllele,percentHeterozygous_List,percentHeterozygous_QTL_List)

  return(simResults_List)

}

#############################################################


generateProgeny <- function(Parent1,Parent2,Number_Progeny,NAM_LinkMap){

  NAM_LinkMap_New <- NAM_LinkMap
  RF_Vector <- as.vector(NAM_LinkMap_New[,3])
  rr <- as.vector(NAM_LinkMap_New[,3])
  number_progeny <- Number_Progeny

  x1<-Parent1
  x2<- Parent2


  a <- Number_Progeny
  n <- dim(x1)[1]

  y<-array(0,dim=c(n,2,a))


  for (i in 1:a) {

    p1<-cumsum(runif(n)<=rr)%%2
    p2<-cumsum(runif(n)<=rr)%%2 #creating gamete from p2


    y[p1==0,1,i]=x1[p1==0,1]

    y[p1==1,1,i] = x1[p1==1,2]

    y[p2==0,2,i] = x2[p2==0,1]

    y[p2==1,2,i] = x2[p2==1,2]
  }

  progeny <- y

  return(progeny)

}


   
getPGStats <- function(CycleGenoTable_AF,Populations){ 

        GenTable_AF <- CycleGenoTable_AF
		populations <- Populations
		
### get island PG parameters
        CycleGenoTable_GInd <- df2genind(GenTable_AF,pop=populations,sep="")
 
		Diff_Stats <- diff_stats(CycleGenoTable_GInd)
	    # Diff_Stats_List_Local <- as.big.matrix(Diff_Stats[[1]])
	    # Diff_Stats_List_Global <- Diff_Stats[[2]]
		
		F <- inbreeding(CycleGenoTable_GInd,res.type="estimate")
        ExpHet <- Hs(CycleGenoTable_GInd)
		Freq_GInd <- makefreq(CycleGenoTable_GInd)
		MAF <- minorAllele(CycleGenoTable_GInd)
		
		return(list(Diff_Stats,F,ExpHet,Freq_GInd,MAF))

}


 
getPGStatsGM_Cyc1 <- function(CycleGenoTable_AF,nextGenGenoTable_NumF,Populations,MarkerEffects,noQTL){ 

        GenTable_AF <- CycleGenoTable_AF
		nextGenGenoTable <- nextGenGenoTable_NumF
		populations <- Populations
		no_QTL <- noQTL
					
### get PG stats from genotable

        # nextGenGenoTable_Numeric <- getNumericFormat(GenTable_AF,AlleleConversionTable_Combined)
		# nextGenGenoTable <- apply(nextGenGenoTable_Numeric,2,as.numeric)
		
		QTL_GenoTable <- getQTLTable(nextGenGenoTable,no_QTL)
		Freq_Allele <-	getFreq_BasePoln(nextGenGenoTable,MarkerEffects)
		FavorableAlleleFreq <- Freq_Allele[[1]]
		FavorableAllele <- Freq_Allele[[2]]
		percentHeterozygous_List <- getPercentHeterozygous(nextGenGenoTable)

		Freq_QTLAllele <- getFreq_BasePoln(QTL_GenoTable,MarkerEffects)
		FavorableQTLAlleleFreq <- Freq_QTLAllele[[1]]
		FavorableQTLAllele <- Freq_QTLAllele[[2]]
		percentHeterozygous_QTL <- getPercentHeterozygous(QTL_GenoTable)
		
### get PG parameters from GenInd

        CycleGenoTable_GInd <- df2genind(GenTable_AF,pop=populations,sep="")
 
		Diff_Stats <- diff_stats(CycleGenoTable_GInd)
	    
		# Diff_Stats_List_Local <- as.big.matrix(Diff_Stats[[1]])
	    # Diff_Stats_List_Global <- Diff_Stats[[2]]
		
		F <- inbreeding(CycleGenoTable_GInd,res.type="estimate")
        ExpHet <- Hs(CycleGenoTable_GInd)
		Freq_GInd <- makefreq(CycleGenoTable_GInd)
		MAF <- minorAllele(CycleGenoTable_GInd)
		
		return(list(FavorableAlleleFreq,FavorableAllele,percentHeterozygous_List,FavorableQTLAlleleFreq,FavorableQTLAllele,percentHeterozygous_QTL,Diff_Stats,F,ExpHet,Freq_GInd,MAF))

}

###

 
getPGStatsGM <- function(CycleGenoTable_AF,nextGenGenoTable_NumF,Populations,FavAllele,MarkerEffects,noQTL){ 

        GenTable_AF <- CycleGenoTable_AF
		populations <- Populations
		no_QTL <- noQTL
		nextGenGenoTable <- nextGenGenoTable_NumF
### get PG stats from genotable

        #nextGenGenoTable_Numeric <- getNumericFormat(GenTable_AF,AlleleConversionTable_Combined)
		#nextGenGenoTable <- apply(nextGenGenoTable_Numeric,2,as.numeric)
		
		QTL_GenoTable <- getQTLTable(nextGenGenoTable,no_QTL)
		
		Freq_Allele <-	getFreq(nextGenGenoTable,FavAllele)
		FavorableAlleleFreq <- Freq_Allele[[1]]
		FavorableAllele  <- Freq_Allele[[2]]
		percentHeterozygous_List <- getPercentHeterozygous(nextGenGenoTable)

		Freq_QTLAllele <- getFreq(QTL_GenoTable,FavAllele)
		FavorableQTLAlleleFreq <- Freq_QTLAllele[[1]]
		FavorableQTLAllele  <- Freq_QTLAllele[[2]]
		percentHeterozygous_QTL <- getPercentHeterozygous(QTL_GenoTable)
		
### get PG parameters from GenInd

        CycleGenoTable_GInd <- df2genind(GenTable_AF,pop=populations,sep="")
 		Diff_Stats <- diff_stats(CycleGenoTable_GInd)
	    
				
		F <- inbreeding(CycleGenoTable_GInd,res.type="estimate")
        ExpHet <- Hs(CycleGenoTable_GInd)
		Freq_GInd <- makefreq(CycleGenoTable_GInd)
		MAF <- minorAllele(CycleGenoTable_GInd)
		
		return(list(FavorableAlleleFreq,FavorableAllele,FavorableQTLAlleleFreq,FavorableQTLAllele,percentHeterozygous_QTL,Diff_Stats,F,ExpHet,Freq_GInd,MAF))

}


###########################################


 # getNumericFormat <- function(GenoTable_AlleleFormat,AlleleConversionTable_Combined){

    # options(warn=-1)
	# genoTable <- GenoTable_AlleleFormat
	# nIndividuals <- nrow(genoTable)
	
	# IA3023_alleles <- as.character(AlleleConversionTable_Combined[,8])
	# Alt_Homozygous_alleles<- as.character(AlleleConversionTable_Combined[,10])
	# Het_alleles <- as.character(AlleleConversionTable_Combined[,22]) 
	
	# IA3023_alleles_Indices <- apply(genoTable,2,function(x) which(as.character(x) == (IA3023_alleles))) 
	# Het_alleles_Indices <- apply(genoTable,2,function(x) which(as.character(x) == (Het_alleles)))
	# Alt_Homozygous_alleles_Indices <- apply(genoTable,2,function(x) which(as.character(x) == (Alt_Homozygous_alleles)))
	
    
	# nMarkers <- length(IA3023_alleles) 
	
	
	# for(i in 1:nMarkers){
		
		# if(length(IA3023_alleles_Indices) !=0){genoTable[,(IA3023_alleles_Indices[[i]])] <- "1"}
		# if(length(Het_alleles_Indices) !=0){ genoTable[,(Het_alleles_Indices[[i]])] <- "0"}
		# if(length(Alt_Homozygous_alleles_Indices)!=0){genoTable[,(Alt_Homozygous_alleles_Indices[[i]])] <- "-1"}
	# }

	# return(genoTable)
 # } 
 

