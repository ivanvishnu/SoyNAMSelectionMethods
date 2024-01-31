

#### PM Model Update Methods
  # if(nCyc%%updateFrequency ==0){
			 
		# if(nCyc >2){
			   # trainGenoNewTablePreCycle <- read.table(trainGeno_FileName,sep="\t")
		# }
# ## Extract training and test data for model build
		# IndividualNames <-rep(0,(nCrosses*nProgeny))

		# IndividualNames <- paste("Ind",c(1:(nCrosses*nProgeny)),sep="")
		
# ##################################################################################################

		# names(phenoSimValues)<-IndividualNames

		# Mean_Fixed<- rep(1,nIndividuals)

		# phenotypicValuesSimTable<- cbind(phenoSimValues,Mean_Fixed)


# #### Table for Simulated Genotypic Values ####################################

		# names(genoSimValues)<-IndividualNames

		# Mean_Fixed<- rep(1,nIndividuals)

		# genotypicValuesSimTable <- cbind(genoSimValues,Mean_Fixed)


# ### Model Retrain 

	# if(modelRetrain==TRUE){ 
			# PredDefined <- FALSE

			# while(PredDefined ==FALSE){

			# indices<-c(1:nIndividuals)

			# trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))
			# testIndices<- indices[which(!indices %in% trainIndices)]

	# #### Training Set table for Genotype and Simulated Phenotype ###################

			# # trainGenoTable <- genoTable_Mod[trainIndices,]
			# trainGenoNewTable <- newNextGenGenoTable[trainIndices,]
			# trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]
			# trainSimGenoTable <-  genotypicValuesSimTable[trainIndices,]

	# #### Validation Set table for Genotype and Simulated Phenotype ###################

			# # testGenoTable<- genoTable_Mod[testIndices,]
			# testGenoNewTable <- newNextGenGenoTable[testIndices,]
			# testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]
			# testSimGenoTable <-  genotypicValuesSimTable[testIndices,]

			# print(length(trainSimGenoTable[,1]))
			# print(sd(trainSimGenoTable[,1]))	
			# print(sd(trainSimPhenoTable[,1]))
			# print(i%%updateFrequency) 
			# print(nCyc%%updateFrequency)
			# print(summary(trainSimPhenoTable[,1]))		
	# ############### Change RRBLUP to Bayes
		  # if(!is.null(sd(trainSimGenoTable[,1])) && !is.null(sd(trainSimPhenoTable[,1])) && !is.na(sd(trainSimGenoTable[,1])) && !is.na(sd(trainSimPhenoTable[,1]))){
			# if(sd(trainSimGenoTable[,1])!=0 && sd(trainSimPhenoTable[,1])!=0 && nCyc%%retrainFrequency==0){

				# if(modelType=="RRBLUP"){
					# PredictionModel_Geno <- (buildRRBLUPModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable))
					# PredictionModel_Pheno <-(buildRRBLUPModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))
				# }else if(modelType=="RRBLUP_REML"){
					# PredictionModel_Geno <- (buildRRBLUPModel_REML(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable))
					# PredictionModel_Pheno <-(buildRRBLUPModel_REML(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))
				# }
			# }
		  # }
		
			# if(length(levels(factor(PredictionModel_Pheno[[1]])))>1 && mean(PredictionModel_Pheno[[1]])>0) {
							 # PredDefined <- TRUE }      
			# }
	   # }

	# gc()	
# ### 
		
    # if(modelUpdate ==TRUE && updateType == "FullSet"){
	
	   # PredictionModel_Geno_PreCycle <- PredictionModel_Geno
       # PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno

        # set.seed(25+Rep)
	    # indices<-c(1:nIndividuals)

        # trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        # testIndices<- indices[which(!indices %in% trainIndices)]

# #### Training Set table for Genotype and Simulated Phenotype ###################

        # trainGenoNewTable <- newNextGenGenoTable[trainIndices,]

        # trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]

        # trainSimGenoValTable <-  genotypicValuesSimTable[trainIndices,]

# #### Validation Set table for Genotype and Simulated Phenotype ###################

        # testGenoNewTable <- newNextGenGenoTable[testIndices,]

        # testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]

        # testSimGenoValTable <-  genotypicValuesSimTable[testIndices,]


# ##### Build RRBLUP prediction model every cycle

       # # trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))
	    # trainGenoNewTableComb <- (rbind(as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))

        # trainSimGenoValTableComb <- rbind((trainSimGenoValTablePreCycle),trainSimGenoValTable)

        # trainSimPhenoValTableComb <- rbind((trainSimPhenoTablePreCycle),trainSimPhenoTable)

        # print(dim(trainGenoNewTableComb))

        # if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){
            # if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0){ 

				# if(modelType=="BayesB"){
					# PredictionModel_Geno <- (buildBayesBModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					# PredictionModel_Pheno <-(buildBayesBModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
				# }else if(modelType=="BL"){
			
					# PredictionModel_Geno <- (buildBLModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					# PredictionModel_Pheno <-(buildBLModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
			    # }else if(modelType=="RRBLUP"){

                        # PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable))

                        # PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable))

                # }else if(modelType=="RRBLUP_REML"){

                        # PredictionModel_Geno <- buildRRBLUPModel_REML(as.matrix(trainGenoNewTableComb), trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        # PredictionModel_Pheno <- buildRRBLUPModel_REML(as.matrix(trainGenoNewTableComb), trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

                # }
				
		    # }
		 
		    # trainGenoNewTablePreCycle <-  trainGenoNewTableComb
            # trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)
            # trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
			
			# if((length(levels(factor(PredictionModel_Pheno[[1]]))) == 1)||(length(levels(factor(PredictionModel_Geno[[1]]))) ==1)){

                         # PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         # PredictionModel_Geno <- PredictionModel_Geno_PreCycle

			# }
		    # rm(trainGenoNewTableComb)
            # rm(trainSimGenoValTableComb)
            # rm(trainSimPhenoValTableComb)
            # gc()

       # }
	# }
	
	# if(modelUpdate ==TRUE && updateType=="TrainingCycleWindow"){
			
		# PredictionModel_Geno_PreCycle <- PredictionModel_Geno
		# PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno
        
	    # indices<-c(1:nIndividuals)

        # trainIndices <-sample(c(1:nIndividuals),(0.8*nIndividuals))

        # testIndices <- indices[which(!indices %in% trainIndices)]
		
		# nTrainIndices <- length(trainIndices)

# #### Training Set table for Genotype and Simulated Phenotype ###################

        # trainGenoNewTable <- newNextGenGenoTable[trainIndices,]
        # trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]
        # trainSimGenoValTable <- genotypicValuesSimTable[trainIndices,]

# #### Validation Set table for Genotype and Simulated Phenotype ###################
        # testGenoNewTable <- newNextGenGenoTable[testIndices,]
        # testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]
        # testSimGenoValTable <- genotypicValuesSimTable[testIndices,]

# ##### Build RRBLUP prediction model every cycle


        # #trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))
		
		# trainGenoNewTableComb <- (rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))
		# trainSimGenoValTableComb <- rbind((trainSimGenoValTablePreCycle),trainSimGenoValTable)
        # trainSimPhenoValTableComb <- rbind((trainSimPhenoTablePreCycle),trainSimPhenoTable)
        # print(dim(trainGenoNewTableComb))
		# Avg_GenoSim_Prev <- mean(trainSimGenoValTableComb[,1])
		# Avg_GenoSim_Current <- mean(trainSimGenoValTable[,1])

		# change_GenoSim <- Avg_GenoSim_Prev - Avg_GenoSim_Current
		# combinedTableCount <- combinedTableCount+1
         # print(change_GenoSim)

        # if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){

           # if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0 && change_GenoSim !=0 ){ 

           # if(modelType=="RRBLUP"){

                        # PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable))

                        # PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable))

           # }else if(modelType=="RRBLUP_REML"){


                        # PredictionModel_Geno <- buildRRBLUPModel_REML((trainGenoNewTableComb),trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        # PredictionModel_Pheno <- buildRRBLUPModel_REML((trainGenoNewTableComb),trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

           # }

	     # } 

           # print(sd(PredictionModel_Pheno[[1]])) 
           # print(sd(PredictionModel_Geno[[1]])) 
           # print(length(levels(factor(PredictionModel_Pheno[[1]]))))
           # print(length(levels(factor(PredictionModel_Geno[[1]]))))	 

           # genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
           # phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)
           # print(sd(phenoValues))
           # print(sd(genoValues))

 
           # if(length(levels(factor(PredictionModel_Pheno[[1]]))) ==1 || length(levels(factor(PredictionModel_Geno[[1]]))) ==1 || (is.na(sd(phenoValues)) || sd(phenoValues)==0) || (is.na(sd(genoValues)) || sd(genoValues)==0) || change_GenoSim==0){

                         # PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         # PredictionModel_Geno <- PredictionModel_Geno_PreCycle

           # }
	 # }
		
		# if(combinedTableCount < trainTableWindowSize){
			
			    # trainGenoNewTablePreCycle <-  trainGenoNewTableComb

                # trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)

                # trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
					
		# }else if(combinedTableCount >= trainTableWindowSize){
			
			    # trainGenoNewTablePreCycle <-  trainGenoNewTableComb[-c(1:nTrainIndices),]

                # trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb[-c(1:nTrainIndices),])

                # trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb[-c(1:nTrainIndices),])
				
				# combinedTableCount <- (trainTableWindowSize-1)
					
		# }
		
	
                # rm(trainGenoNewTableComb)
                # rm(trainSimGenoValTableComb)
                # rm(trainSimPhenoValTableComb)
                # gc()
       # }

    # nLines_trainGeno <- nrow(trainGenoNewTablePreCycle)
	# rownames(trainGenoNewTablePreCycle) <- paste("Ind",c(1:nLines_trainGeno),sep="")
	# rownames(trainSimPhenoTablePreCycle) <- paste("Ind",c(1:nLines_trainGeno),sep="")
	# rownames(trainSimGenoValTablePreCycle) <- paste("Ind",c(1:nLines_trainGeno),sep="")


    # trainGeno_FileName <- paste("trainTable_",selCriteria,"_",BD,"_",condition,"_",Rep,"_", ".txt",sep="") 
    # write.table(trainGenoNewTablePreCycle,trainGeno_FileName,sep="\t") 
	# rm(trainGenoNewTablePreCycle)
  # }
   

#### Island Model GP model update method



# if(modelUpdate == TRUE && nCyc%%updateFrequency==0){

	# IndividualNames <-rep(0,(nCrosses*nProgeny))

    # for(nInd in 1:(nCrosses*nProgeny)){

        # IndividualNames[nInd]<- paste("Ind",nInd,sep="")

    # }

# ################################################################################################
	# names(phenoValSimValues)<- IndividualNames
    # Mean_Fixed<- rep(1,nIndividuals)

    # phenotypicValuesSimTable<- cbind(phenoValSimValues,Mean_Fixed)

# ## Table for Simulated Genotypic Values ####################################


    # names(genoValSimValues)<-IndividualNames
    # Mean_Fixed<- rep(1,nIndividuals)

    # genotypicValuesSimTable <- cbind(genoValSimValues,Mean_Fixed)

# ###################

	# if(modelRetrain==TRUE){

		# indices<-c(1:nIndividuals)

		# trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))
		# testIndices<- indices[which(!indices %in% trainIndices)]

	# ## Training Set table for Genotype and Simulated Phenotype ###################

		# trainGenoTable <- genoTable_Mod[trainIndices,]
		# trainGenoNewTable <- newNextGenGenoTable[trainIndices,]
		# trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]
		# trainSimGenoTable <-  genotypicValuesSimTable[trainIndices,]

	# ## Validation Set table for Genotype and Simulated Phenotype ###################

		# testGenoTable<- genoTable_Mod[testIndices,]
		# testGenoNewTable <- newNextGenGenoTable[testIndices,]
		# testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]
		# testSimGenoTable <-  genotypicValuesSimTable[testIndices,]

# ############# Changed RRBLUP to Bayes

	   # if(sd(trainSimGenoTable[,1])!=0 && sd(trainSimPhenoTable[,1])!=0 && nCyc%%updateFrequency==0){
	  
	    # if(modelType=="BayesB"){ 
	   		# PredictionModel_Geno <- (buildBayesBModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable,h2))
			# PredictionModel_Pheno <-(buildBayesBModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
		# }else if(modelType=="BL"){ 
		
			# PredictionModel_Geno <- (buildBLModel(trainGenoNewTable,trainSimGenoTable,testGenoNewTable,testSimGenoTable,h2))
			# PredictionModel_Pheno <-(buildBLModel(trainGenoNewTable,trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
		# }else if(modelType=="RRBLUP"){

                        # PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable))

                        # PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableiComb),trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))

        # }else if(modelType=="RRBLUP_REML"){


                        # PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        # PredictionModel_Pheno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

        # }
		
	  # }
	  # gc()
	  
	# }
	
	# if(modelUpdate ==TRUE && updateType == "EqPartition"){	
	
	
	    # PredictionModel_Geno_PreCycle <- PredictionModel_Geno
        # PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno
        
	    # indices<-c(1:nIndividuals)

        # trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        # testIndices<- indices[which(!indices %in% trainIndices)]

# ## Training Set table for Genotype and Simulated Phenotype ###################

        # trainGenoNewTable <- newNextGenGenoTable[trainIndices,]

        # trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]

        # trainSimGenoValTable <-  genotypicValuesSimTable[trainIndices,]

# ## Validation Set table for Genotype and Simulated Phenotype ###################

        # testGenoNewTable <- newNextGenGenoTable[testIndices,]

        # testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]

        # testSimGenoValTable <-  genotypicValuesSimTable[testIndices,]

		
# Training Set List with equal partitioning among cycles 		
		
	# trainGenoNewTable_List[[nCyc]] <- as.big.matrix(trainGenoNewTable)
	# trainSimGenoValTable_List[[nCyc]] <- as.big.matrix(trainSimGenoValTable)
	# trainSimPhenoTable_List[[nCyc]] <- as.big.matrix(trainSimPhenoTable)
      
	# nTrainingSets <- length(trainGenoNewTable_List)
	# nTrainIndividuals <- (0.8*nIndividuals) %/% nCyc
        # nTrainIndividualsPlusRemainder <- nTrainIndividuals + ((0.8*nIndividuals) %% nCyc)
		
# ### Build prediction model 

        # trainGenoNewTableComb <- c() 
		# trainSimGenoValTableComb <- c() 
		# trainSimPhenoValTableComb <- c()
		 
		
        # trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))
 
        # for(nTrainSets in 1:nTrainingSets){
		   
		    # if(nTrainSets < nTrainingSets){
			   # trainingIndices <- sample(c(1:(0.8*nIndividuals)),nTrainIndividuals)
		    # }else if(nTrainSets == nTrainingSets){
			   # trainingIndices <- sample(c(1:(0.8*nIndividuals)),nTrainIndividualsPlusRemainder)
		    # }   
			   
			   # trainGenoTableInCycle <- trainGenoNewTable_List[[nTrainSets]][trainingIndices,]
			   # trainSimGenoValTableInCycle <- trainSimGenoValTable_List[[nTrainSets]][trainingIndices,]
			   # trainSimPhenoTableInCycle <- trainSimPhenoTable_List[[nTrainSets]][trainingIndices,]
			   
			   # trainGenoNewTableComb <- rbind((trainGenoNewTableComb),bigmemory::as.matrix(trainGenoTableInCycle))
			   # trainSimGenoValTableComb <- rbind((trainSimGenoValTableComb),bigmemory::as.matrix(trainSimGenoValTableInCycle))
                           # trainSimPhenoValTableComb <- rbind((trainSimPhenoValTableComb),bigmemory::as.matrix(trainSimPhenoTableInCycle))
		# }
		
        # print(dim(trainGenoNewTableComb))

        # if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){
            # if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0 && nCyc%%updateFrequency==0){ 

				# if(modelType=="BayesB"){
					# PredictionModel_Geno <- (buildBayesBModel((trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					# PredictionModel_Pheno <-(buildBayesBModel((trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
				# }else if(modelType=="BL"){
			
					# PredictionModel_Geno <- (buildBLModel((trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					# PredictionModel_Pheno <-(buildBLModel((trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
			    # }else if(modelType=="RRBLUP"){

                        # PredictionModel_Geno <- (buildRRBLUPModel((trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable))

                        # PredictionModel_Pheno <-(buildRRBLUPModel((trainGenoNewTableiComb),trainSimPhenoTable,testGenoNewTable,testSimPhenoTable))

                # }else if(modelType=="RRBLUP_REML"){



                        # PredictionModel_Geno <- buildRRBLUPModel_REML((trainGenoNewTableComb), trainSimGenoValTableComb,(testGenoNewTable),testSimGenoValTable)

                        # PredictionModel_Pheno <- buildRRBLUPModel_REML((trainGenoNewTableComb), trainSimPhenoValTableComb,(testGenoNewTable),testSimPhenoTable)

                # }
		    # }
	# }
		 
		    # if((length(levels(factor(PredictionModel_Pheno[[1]]))) == 1 && sd(PredictionModel_Pheno[[1]])==0) ||(length(levels(factor(PredictionModel_Geno[[1]]))) ==1 && sd(PredictionModel_Geno[[1]])==0)){

                         # PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         # PredictionModel_Geno <- PredictionModel_Geno_PreCycle

                     # }
		   

            # rm(trainGenoNewTableComb)
            # rm(trainSimGenoValTableComb)
            # rm(trainSimPhenoValTableComb)
            # gc()

        # }
		 
	# if(modelUpdate ==TRUE && updateType == "FullSet"){
	
	   # PredictionModel_Geno_PreCycle <- PredictionModel_Geno
       # PredictionModel_Pheno_PreCycle <- PredictionModel_Pheno

	    # indices<-c(1:nIndividuals)

        # trainIndices<-sample(c(1:nIndividuals),(0.8*nIndividuals))

        # testIndices<- indices[which(!indices %in% trainIndices)]

# ## Training Set table for Genotype and Simulated Phenotype ###################

        # trainGenoNewTable <- newNextGenGenoTable[trainIndices,]

        # trainSimPhenoTable <- phenotypicValuesSimTable[trainIndices,]

        # trainSimGenoValTable <-  genotypicValuesSimTable[trainIndices,]

# ## Validation Set table for Genotype and Simulated Phenotype ###################

        # testGenoNewTable <- newNextGenGenoTable[testIndices,]

        # testSimPhenoTable <- phenotypicValuesSimTable[testIndices,]

        # testSimGenoValTable <-  genotypicValuesSimTable[testIndices,]


# ### Build RRBLUP prediction model every cycle


        # trainGenoNewTableComb <- as.big.matrix(rbind(bigmemory::as.matrix(trainGenoNewTablePreCycle),as.matrix(trainGenoNewTable)))

        # trainSimGenoValTableComb <- rbind((trainSimGenoValTablePreCycle),trainSimGenoValTable)

        # trainSimPhenoValTableComb <- rbind((trainSimPhenoTablePreCycle),trainSimPhenoTable)

        # print(dim(trainGenoNewTableComb))

        # if(!is.null(sd(trainSimGenoValTableComb[,1])) && !is.null(sd(trainSimPhenoValTableComb[,1])) && !is.na(sd(trainSimGenoValTableComb[,1])) && !is.na(sd(trainSimPhenoValTableComb[,1]))){
            # if(sd(trainSimGenoValTableComb[,1])!=0 && sd(trainSimPhenoValTableComb[,1])!=0 && nCyc%%updateFrequency==0){ 

				# if(modelType=="BayesB"){
					# PredictionModel_Geno <- (buildBayesBModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					# PredictionModel_Pheno <-(buildBayesBModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoTable,testGenoNewTable,testSimPhenoTable,h2))
				# }else if(modelType=="BL"){
			
					# PredictionModel_Geno <- (buildBLModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable,h2))
					# PredictionModel_Pheno <-(buildBLModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable,h2))
			    # }else if(modelType=="RRBLUP"){

                        # PredictionModel_Geno <- (buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableComb),trainSimGenoValTableComb,testGenoNewTable,testSimGenoValTable))

                        # PredictionModel_Pheno <-(buildRRBLUPModel(bigmemory::as.matrix(trainGenoNewTableiComb),trainSimPhenoValTableComb,testGenoNewTable,testSimPhenoTable))

                # }else if(modelType=="RRBLUP_REML"){



                        # PredictionModel_Geno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimGenoValTableComb,as.matrix(testGenoNewTable),testSimGenoValTable)

                        # PredictionModel_Pheno <- buildRRBLUPModel_REML(bigmemory::as.matrix(trainGenoNewTableComb), trainSimPhenoValTableComb,as.matrix(testGenoNewTable),testSimPhenoTable)

                # }
				
		    # }
		 
		    # trainGenoNewTablePreCycle <-  trainGenoNewTableComb
            # trainSimPhenoTablePreCycle <-   (trainSimPhenoValTableComb)
            # trainSimGenoValTablePreCycle <- (trainSimGenoValTableComb)
			
	# if((length(levels(factor(PredictionModel_Pheno[[1]]))) == 1 && sd(PredictionModel_Pheno[[1]])==0) ||(length(levels(factor(PredictionModel_Geno[[1]]))) ==1 && sd(PredictionModel_Geno[[1]])==0)){

                         # PredictionModel_Pheno <- PredictionModel_Pheno_PreCycle
                         # PredictionModel_Geno <- PredictionModel_Geno_PreCycle

       # }
			
			

            # rm(trainGenoNewTableComb)
            # rm(trainSimGenoValTableComb)
            # rm(trainSimPhenoValTableComb)
            # gc()

        # }
# }
# }	
