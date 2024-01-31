### Selection Methods


###### IM GM Frontier method for optimal cross selection

getIM_GM_Frontier_selectedGenoList_SingleFams_V2 <- function(NextGenGenoTable_GM_List,NextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,nFamily){


	nextGenGenoTable_GM_List <- NextGenGenoTable_GM_List 
	nextGenGenoTableMod_GM_List <- NextGenGenoTableMod_GM_List 
	nSnglFams <- nFamilies

##  GM Parameters       
        GM_Param_List <- GM_Frontier_Param_List
        selectedParents_Geno_List<- list()
        selectedCriterion_List <- list()

### GM_Parameters

        Mc.cores= GM_Param_List$no_cores
        No_cores = GM_Param_List$no_cores
        Mutprob= GM_Param_List$mutprob
        NpopGA= GM_Param_List$npopGA
        NitGA= GM_Param_List$nitGA
        Method= GM_Param_List$method
        no_cores <- No_cores
        #Plotiters= GM_Param_List$plotiters
        #Nelite= GM_Param_List$nelite


         rNames_Family1 <- c()

                  for(indPairs in 1:nrow(nextGenGenoTable_GM_List)){

                    rNames_Family1 <-c(rNames_Family1,paste(nFamily," Ind",indPairs))
                  }

                  rownames(nextGenGenoTable_GM_List) <- rNames_Family1

				  M <- (nextGenGenoTable_GM_List)
                  K <- GenomicMatingV2::Amat.pieces(M,pieces=15,mc.cores=no_cores)

                  Markers <- nextGenGenoTableMod_GM_List

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

                        baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutionsFrontier(Markers,Markers2=NULL,K,markerEffects_Geno,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,method=Method)
                        # ,nmates=NMates ; ,plotiters=Plotiters
                        N <- nrow(M)
                        ebvs = Markers %*% markerEffects_Geno
                        names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")

                        if(nrow(baseGASolns_GenoModel[[1]]) > 1){

                        Quartile3 <- summary(baseGASolns_GenoModel[[1]][,1])[5]
                        Quartile1 <- summary(baseGASolns_GenoModel[[1]][,1])[2]


                        indices1 <- which(baseGASolns_GenoModel[[1]][,1]>= Quartile3)
			indicesSortRange1 <- sort.int(baseGASolns_GenoModel[[1]][indices1,1],decreasing=TRUE,index.return=TRUE)
                        indicesRange1 <- indicesSortRange1[[2]][1:10]

                        indices2 <- which((baseGASolns_GenoModel[[1]][,1]<= Quartile3) & (baseGASolns_GenoModel[[1]][,1]>= Quartile1))
                        indicesRange2 <- sample(indices2,10,replace=TRUE)

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
 }else if(selectionOnGeno==FALSE){
## selection on phenotypic values 
                       
                    baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutionsFrontier(Markers,Markers2=NULL,K,markerEffects_Pheno,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,method=Method)

#nmates=NMates plotiters=Plotiters

       ###   length(baseGASolns_GenoModel) - 2 ; 
       # case1 : returns only one set of criterion values and baseGASolns is not null 
       # case2: returns > 1 set of criterion values 

       # case 2a:  more than 30 sets and there is a range of genovalues. So requires reduction of parameter set to 30 sets
       # case 2b: less than 30 sets and there is a range of genovalues. Requires no reduction of parameter set
       # case 2c: number of parameter sets could be less than or greater than 30, but there is no range of genovalues (with defined genovalues) 
       # case 2d: number of parameter sets could be less than or greater than 30, with undefined (NAN genovalues) 

           
                        N <- nrow(M)
                        ebvs = Markers %*% markerEffects_Pheno
                        names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")

               		                  
				if(length(baseGASolns_GenoModel[[1]])%%3 ==0){
			 
                        # case 1	

                                if(length(baseGASolns_GenoModel[[1]])/3 == 1){
					                nParameters <- 1
                			        selectedParents_Geno_List[[nParameters]] <-  baseGASolns_GenoModel[[2]][[1]]
                                    selectedCriterion_List[[nParameters]] <- baseGASolns_GenoModel[[1]]
                                          
                                } 
 # case 2	    
                                if(length(baseGASolns_GenoModel[[1]])/3 > 1){  
                                     
			# case 2a:                              
                            		if(length(baseGASolns_GenoModel[[1]]) >=30 && !anyNA(baseGASolns_GenoModel[[1]][,1])){
             
                                            if((summary(baseGASolns_GenoModel[[1]][,1])[5]-summary(baseGASolns_GenoModel[[1]][,1])[2]) !=0){

											Quartile3 <- summary(baseGASolns_GenoModel[[1]][,1])[5]
											Quartile1 <- summary(baseGASolns_GenoModel[[1]][,1])[2]


											indices1 <- which(baseGASolns_GenoModel[[1]][,1]>= Quartile3)
											indicesSortRange1 <- sort.int(baseGASolns_GenoModel[[1]][indices1,1],decreasing=TRUE,index.return=TRUE)
											indicesRange1 <- indicesSortRange1[[2]][1:10]

											indices2 <- which((baseGASolns_GenoModel[[1]][,1]<= Quartile3) & (baseGASolns_GenoModel[[1]][,1]>= Quartile1))

											indicesRange2 <- sample(indices2,10,replace=TRUE)

											indices3 <- which(baseGASolns_GenoModel[[1]][,1]<= Quartile1)
											indicesSortRange3 <- sort.int(baseGASolns_GenoModel[[1]][indices3,1],decreasing=FALSE,index.return=TRUE)
											indicesRange3 <- indicesSortRange3[[2]][1:10]

											indicesRange <- c(indicesRange1,indicesRange2,indicesRange3)
											nParameters <- length(indicesRange)
											#for(nParameter in 1:nParameters){

													nParameterIndex <- indicesRange[nParameters]

													selectedParents_Geno_List[[1]] <-  baseGASolns_GenoModel[[2]][[nParameterIndex]]
													selectedCriterion_List[[1]] <- baseGASolns_GenoModel[[1]][nParameterIndex,]
											#}
										}
										
									}
                                                                 
 # case 2b: 

                        		if(length(baseGASolns_GenoModel[[1]][,1])<=30 && length(baseGASolns_GenoModel[[1]][,1])> 1 && !anyNA(baseGASolns_GenoModel[[1]][,1])){

									if(summary(baseGASolns_GenoModel[[1]][,1])[5] -summary(baseGASolns_GenoModel[[1]][,1])[2] !=0){

										nParameters <- length(baseGASolns_GenoModel[[1]][,1])
										#for(nParameter in 1:nParameters){

											selectedParents_Geno_List[[1]] <-  baseGASolns_GenoModel[[2]][[nParameters]]
											selectedCriterion_List[[1]] <- baseGASolns_GenoModel[[1]][nParameters,]

										#}
									}
								    # selectedParents_Geno_List[[1]] <-  baseGASolns_GenoModel[[2]][[nParameters]]
								    # selectedCriterion_List[[1]] <- baseGASolns_GenoModel[[1]][nParameters,]
								
								} 

			#case 2c:   
								if(length(baseGASolns_GenoModel[[1]][,1]) >=1 && !anyNA(baseGASolns_GenoModel[[1]][,1])){

                                        if((summary(baseGASolns_GenoModel[[1]][,1])[5] -summary(baseGASolns_GenoModel[[1]][,1])[2]) ==0){
										   nParameters <- 1
										   selectedParents_Geno_List[[nParameters]] <-  baseGASolns_GenoModel[[2]][[1]]
										   selectedCriterion_List[[nParameters]] <- baseGASolns_GenoModel[[1]][1,]
				        
                                        }
                                }

			#case 2d: 

								if(length(baseGASolns_GenoModel[[1]][,1])>=1 && anyNA(baseGASolns_GenoModel[[1]][,1])){ 
				           
									nParameters <- 1
									selectedParents_Geno_List[[nParameters]] <-  baseGASolns_GenoModel[[2]][[1]]
									selectedCriterion_List[[nParameters]] <- baseGASolns_GenoModel[[1]][1,]
				          
								}
                    
							}
                } 
                              

            }
    
            return(list(selectedParents_Geno_List,selectedCriterion_List))
    }

########

getIM_GM_Frontier_selectedGenoList_SingleFams_V2A <- function(NextGenGenoTable_GM_List,NextGenGenoTableMod_GM_List,
    PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,nFamily)
{ 
    nextGenGenoTable_GM_List <- NextGenGenoTable_GM_List
    nextGenGenoTableMod_GM_List <- NextGenGenoTableMod_GM_List
    nSnglFams <- nFamilies
    GM_Param_List <- GM_Frontier_Param_List
    selectedParents_Geno_List <- list()
    selectedCriterion_List <- list()
    Mc.cores = GM_Param_List$no_cores
    No_cores = GM_Param_List$no_cores
    Mutprob = GM_Param_List$mutprob
    NpopGA = GM_Param_List$npopGA
    NitGA = GM_Param_List$nitGA
    Method = GM_Param_List$method
    no_cores <- No_cores
    rNames_Family1 <- c()
    for (indPairs in 1:nrow(nextGenGenoTable_GM_List)) {
        rNames_Family1 <- c(rNames_Family1, paste(nFamily, " Ind",
            indPairs))
    }
    rownames(nextGenGenoTable_GM_List) <- rNames_Family1
    M <- (nextGenGenoTable_GM_List)
    K <- GenomicMatingV2::Amat.pieces(M, pieces = 15, mc.cores = no_cores)
    Markers <- nextGenGenoTableMod_GM_List
    if (modelType == "BayesB" || modelType == "BL") {
        markerEffects_Geno <- as.vector(unlist(PredictionModel_Geno[[3]]))
        markerEffects_Pheno <- as.vector(unlist(PredictionModel_Pheno[[3]]))
    }
    if (modelType == "RRBLUP" || modelType == "RRBLUP_REML") {
        markerEffects_Geno <- as.vector(unlist(PredictionModel_Geno[[1]]))
        markerEffects_Pheno <- as.vector(unlist(PredictionModel_Pheno[[1]]))
    }
    
     if (selectionOnGeno == TRUE) {
        baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutionsFrontier(Markers,
            Markers2 = NULL, K, markerEffects_Geno, mc.cores = Mc.cores,
            mutprob = Mutprob, npopGA = NpopGA, nitGA = NitGA,
            method = Method)
        N <- nrow(M)
        ebvs = Markers %*% markerEffects_Geno
        names(ebvs) <- colnames(K) <- rownames(K) <- paste("l",
            1:N, sep = "")
        if (nrow(baseGASolns_GenoModel[[1]]) > 1) {
            Quartile3 <- summary(baseGASolns_GenoModel[[1]][,
                1])[5]
            Quartile1 <- summary(baseGASolns_GenoModel[[1]][,
                1])[2]
            indices1 <- which(baseGASolns_GenoModel[[1]][, 1] >=
                Quartile3)
            indicesSortRange1 <- sort.int(baseGASolns_GenoModel[[1]][indices1,
                1], decreasing = TRUE, index.return = TRUE)
            indicesRange1 <- indicesSortRange1[[2]][1:10]
            indices2 <- which((baseGASolns_GenoModel[[1]][, 1] <=
                Quartile3) & (baseGASolns_GenoModel[[1]][, 1] >=
                Quartile1))
            indicesRange2 <- sample(indices2,10,replace=TRUE)
            indices3 <- which(baseGASolns_GenoModel[[1]][, 1] <=
                Quartile1)
            indicesSortRange3 <- sort.int(baseGASolns_GenoModel[[1]][indices3,
                1], decreasing = FALSE, index.return = TRUE)
            indicesRange3 <- indicesSortRange3[[2]][1:10]
            indicesRange <- c(indicesRange1, indicesRange2, indicesRange3)
            nParameters <- length(indicesRange)
            for (nParameter in 1:nParameters) {
                nParameterIndex <- indicesRange[nParameter]
                selectedParents_Geno_List[[nParameter]] <- baseGASolns_GenoModel[[2]][[nParameterIndex]]
                selectedCriterion_List[[nParameter]] <- baseGASolns_GenoModel[[1]][nParameterIndex,]
            }
        }
    }
    else if (selectionOnGeno == FALSE) {
        baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutionsFrontier(Markers,
            Markers2 = NULL, K, markerEffects_Pheno, mc.cores = Mc.cores,
            mutprob = Mutprob, npopGA = NpopGA, nitGA = NitGA,
            method = Method)
        N <- nrow(M)
        ebvs = Markers %*% markerEffects_Pheno
        names(ebvs) <- colnames(K) <- rownames(K) <- paste("l",
            1:N, sep = "")
        if (length(baseGASolns_GenoModel[[1]])%%3 == 0) {
            if (length(baseGASolns_GenoModel[[1]])/3 == 1) {
                nParameters <- 1
                selectedParents_Geno_List[[nParameters]] <- baseGASolns_GenoModel[[2]][[1]]
                selectedCriterion_List[[nParameters]] <- baseGASolns_GenoModel[[1]]
            }
            if (length(baseGASolns_GenoModel[[1]])/3 > 1){
                if (length(baseGASolns_GenoModel[[1]]) > 30 && !anyNA(baseGASolns_GenoModel[[1]][, 1])){
                  if ((summary(baseGASolns_GenoModel[[1]][, 1])[5] - summary(baseGASolns_GenoModel[[1]][,1])[2]) !=0){
                    Quartile3 <- summary(baseGASolns_GenoModel[[1]][,1])[5]
                    Quartile1 <- summary(baseGASolns_GenoModel[[1]][,1])[2]
                        indices1 <- which(baseGASolns_GenoModel[[1]][,1] >= Quartile3)
                    indicesSortRange1 <- sort.int(baseGASolns_GenoModel[[1]][indices1,
                      1], decreasing = TRUE, index.return = TRUE)
                    indicesRange1 <- indicesSortRange1[[2]][1:10]
                    indices2 <- which((baseGASolns_GenoModel[[1]][,
                      1] <= Quartile3) & (baseGASolns_GenoModel[[1]][,
                      1] >= Quartile1))
                    indicesRange2 <- sample(indices2,10,replace=TRUE)
                    indices3 <- which(baseGASolns_GenoModel[[1]][,
                      1] <= Quartile1)
                    indicesSortRange3 <- sort.int(baseGASolns_GenoModel[[1]][indices3,1], decreasing = FALSE, index.return = TRUE)
                    indicesRange3 <- indicesSortRange3[[2]][1:10]
                    indicesRange <- c(indicesRange1, indicesRange2,indicesRange3)
                    nParameters <- length(indicesRange)
                    for(nParameter in 1:nParameters){
                    	nParameterIndex <- indicesRange[nParameter]
                    	selectedParents_Geno_List[[nParameter]] <- baseGASolns_GenoModel[[2]][[nParameterIndex]]
                    	selectedCriterion_List[[nParameter]] <- baseGASolns_GenoModel[[1]][nParameterIndex,]
                    }
                  }
                }
                if(length(baseGASolns_GenoModel[[1]][, 1]) <= 30 && length(baseGASolns_GenoModel[[1]][, 1]) > 1 && !anyNA(baseGASolns_GenoModel[[1]][, 1])){
                  if (summary(baseGASolns_GenoModel[[1]][, 1])[5] - summary(baseGASolns_GenoModel[[1]][, 1])[2] != 0){
                    nParameters <- length(baseGASolns_GenoModel[[1]][,1])
                    for(nParameter in 1:nParameters){
                       nParameterIndex <- nParameter
                       selectedParents_Geno_List[[nParameter]] <- baseGASolns_GenoModel[[2]][[nParameterIndex]]
                       selectedCriterion_List[[nParameter]] <- baseGASolns_GenoModel[[1]][nParameterIndex,]
                    }
                  }
                }
                if (length(baseGASolns_GenoModel[[1]][, 1]) >= 1 && !anyNA(baseGASolns_GenoModel[[1]][, 1])) {
                  if ((summary(baseGASolns_GenoModel[[1]][, 1])[5] - summary(baseGASolns_GenoModel[[1]][, 1])[2]) ==   0) {
                    nParameters <- 1
                    selectedParents_Geno_List[[nParameters]] <- baseGASolns_GenoModel[[2]][[1]]
                    selectedCriterion_List[[nParameters]] <- baseGASolns_GenoModel[[1]][1,]
                  }
                }
                if (length(baseGASolns_GenoModel[[1]][, 1]) >= 1 && anyNA(baseGASolns_GenoModel[[1]][, 1])){
                  nParameters <- 1
                  selectedParents_Geno_List[[nParameters]] <- baseGASolns_GenoModel[[2]][[1]]
                  selectedCriterion_List[[nParameters]] <- baseGASolns_GenoModel[[1]][1,]
                }
            }
        }
    }
    return(list(selectedParents_Geno_List, selectedCriterion_List))
} 


####### 


getIM_GM_Frontier_selectedGenoList_SingleFams_V2_Reduced <- function(nextGenGenoTable_GM_List,nextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,GM_Frontier_Param_List,nFamily){
  	
### GM method for nextGen Geno data set
    
	nSnglFams <- nFamilies

##  GM Parameters	
	GM_Param_List <- GM_Frontier_Param_List
	selectedParents_Geno_List<- list()
	selectedCriterion_List <- list()
	
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
	no_cores <- No_cores
#######
	
	# for(nFamily in 1:nSnglFams){
                 		  
		  rNames_Family1 <- c() 
		 		 
		  for(indPairs in 1:nrow(nextGenGenoTable_GM_List)){ 
		
		    rNames_Family1 <-c(rNames_Family1,paste(nFamily," Ind",indPairs))
		  }
		
		  rownames(nextGenGenoTable_GM_List) <- rNames_Family1
		 		
         
		  M <- (nextGenGenoTable_GM_List)
		  K <- GenomicMatingV2::Amat.pieces(M,pieces=15,mc.cores=no_cores) 
	    
		  Markers <- nextGenGenoTableMod_GM_List
		
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
			
			if(nrow(baseGASolns_GenoModel[[1]]) > 1)

			Quartile3 <- summary(baseGASolns_GenoModel[[1]][,1])[5]
			Quartile1 <- summary(baseGASolns_GenoModel[[1]][,1])[2]
			
			
			indices1 <- which(baseGASolns_GenoModel[[1]][,1]>= Quartile3)
            indicesSortRange1 <- sort.int(baseGASolns_GenoModel[[1]][indices1,1],decreasing=TRUE,index.return=TRUE)
			indicesRange1 <- indicesSortRange1[[2]][1:10]
			
			indices2 <- which((baseGASolns_GenoModel[[1]][,1]<= Quartile3) & (baseGASolns_GenoModel[[1]][,1]>= Quartile1))
			indicesRange2 <- sample(indices2,10,replace=TRUE) 
			
            indices3 <- which(baseGASolns_GenoModel[[1]][,1]<= Quartile1)
            indicesSortRange3 <- sort.int(baseGASolns_GenoModel[[1]][indices3,1],decreasing=FALSE,index.return=TRUE) 
			indicesRange3 <- indicesSortRange3[[2]][1:10]
			
			indicesRange <- indicesRange3
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
			
			if(!is.null(dim(baseGASolns_GenoModel[[1]]))){
			Quartile3 <- summary(baseGASolns_GenoModel[[1]][,1])[5]
			Quartile1 <- summary(baseGASolns_GenoModel[[1]][,1])[2]
			
			
			indices1 <- which(baseGASolns_GenoModel[[1]][,1]>= Quartile3)
            indicesSortRange1 <- sort.int(baseGASolns_GenoModel[[1]][indices1,1],decreasing=TRUE,index.return=TRUE)
			indicesRange1 <- indicesSortRange1[[2]][1:10]
			
			indices2 <- which((baseGASolns_GenoModel[[1]][,1]<= Quartile3) & (baseGASolns_GenoModel[[1]][,1]>= Quartile1))
			indicesRange2 <- sample(indices2,10,replace=TRUE) 
			
            indices3 <- which(baseGASolns_GenoModel[[1]][,1]<= Quartile1)
            indicesSortRange3 <- sort.int(baseGASolns_GenoModel[[1]][indices3,1],decreasing=FALSE,index.return=TRUE) 
			indicesRange3 <- indicesSortRange3[[2]][1:10]
			
			indicesRange <- indicesRange3
			nParameters <- length(indicesRange)
			for(nParameter in 1:nParameters){
	        
			   nParameterIndex <- indicesRange[nParameter]
	        	       
			   selectedParents_Geno_List[[nParameter]] <-  baseGASolns_GenoModel[[2]][[nParameterIndex]]
			   selectedCriterion_List[[nParameter]] <- baseGASolns_GenoModel[[1]][nParameterIndex,]
		
			  } 
			}else if(is.null(dim(baseGASolns_GenoModel[[1]]))){ 
                nParameters <- 1 

              selectedParents_Geno_List[[nParameters]] <-  baseGASolns_GenoModel[[2]]
			  selectedCriterion_List[[nParameters]] <- baseGASolns_GenoModel[[1]]
            }
        }
	  
	   	
       		return(list(selectedParents_Geno_List,selectedCriterion_List))
    }



###########################

####### GM method for optimal cross selection
	

 getGM_selectedGenoList_FamilyPairs<- function(nextGenGenoTable_GM_List,nextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,familyPairs,GM_Param_List){

	nextGenGenoTable_GM_List <- nextGenGenoTable_GM_List 
	nextGenGenoTableMod_GM_List <- nextGenGenoTableMod_GM_List 
	
	selectedParents_Geno_List <- list()
	selectedParents_NumProgeny_Geno_List <- list()
	
	nFamilyPairs <- dim(familyPairs)[1]
	
	selectedParents_Geno_FamPair_List <- list()
	
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
		
	
	for(nFamPairs in 1:nFamilyPairs){
         
          family1<- familyPairs[nFamPairs,1]
	      family2<- familyPairs[nFamPairs,2]
		  rNames_Family <- c() 
		  # rNames_Family2 <- c()
		
		  for(famPairs in 1:nrow(nextGenGenoTable_GM_List[[nFamPairs]])){ 
		
		     rNames_Family <-c(rNames_Family,paste(nFamPairs," Ind",famPairs))
		     # rNames_Family2 <-c(rNames_Family2,paste(family2," Ind",famPairs))
		   }
		
		  rownames(nextGenGenoTable_GM_List[[nFamPairs]]) <- rNames_Family
		 
 		  # rownames(nextGenGenoTable_GM_List[[family2]]) <- rNames_Family2
			
		  # M <- rbind(nextGenGenoTable_GM_List[[family1]],nextGenGenoTable_GM_List[[family2]])
		 
 		  M <- nextGenGenoTable_GM_List[[nFamPairs]]
          K <- GenomicMatingV2::Amat.pieces(M,pieces=15,mc.cores=No_cores) 
			
		  # Markers <- rbind(nextGenGenoTableMod_GM_List[[family1]],nextGenGenoTableMod_GM_List[[family2]])
		  
		  Markers <- nextGenGenoTableMod_GM_List[[nFamPairs]]
		  
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

			selectedParents_Geno_List[[nFamPairs]] <- selectedParents_Geno
		
			selectedParents_NumProgeny_Geno_List[[nFamPairs]] <- as.vector(selectedParents_NumProgeny_Geno)
			
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
		    selectedParents_Geno_List[[nFamPairs]] <- selectedParents_Geno
			selectedParents_NumProgeny_Geno_List[[nFamPairs]] <- as.vector(selectedParents_NumProgeny_Geno)
			
			
		}
	
	    
####  Exchange size list among family pairs	
        
			parent1 <- rep(0,length(selectedParents_Geno_List[[nFamPairs]])/2)
			parent2 <- rep(0,length(selectedParents_Geno_List[[nFamPairs]])/2)
			
			# selectedParents_Geno_fam1 <- rep(0,length(selectedParents_Geno[[nFamPairs]])/2)
			# selectedParents_Geno_fam2 <- rep(0,length(selectedParents_Geno[[nFamPairs]])/2)
			
			selectedParents_Geno_fam1 <- c() 
			selectedParents_Geno_fam2  <- c()
						
			nParPairs <-1 
	   
			while(nParPairs <length(selectedParents_Geno_List[[nFamPairs]])){ 
				parent1[nParPairs] <-selectedParents_Geno_List[[nFamPairs]][nParPairs]
				parent2[nParPairs] <- selectedParents_Geno_List[[nFamPairs]][nParPairs+1]
				nParPairs <- nParPairs+2
			}
			parent1 <- parent1[seq(1,length(selectedParents_Geno_List[[nFamPairs]]),by=2)]
			parent2 <- parent2[seq(1,length(selectedParents_Geno_List[[nFamPairs]]),by=2)] 
				
			## sort families into family pairs
			for(nParents in 1:length(parent1)){

			   if(!is.na(parent1[nParents])){
			
				if(parent1[nParents]<=100){
						selectedParents_Geno_fam1 <- c(selectedParents_Geno_fam1,parent1[nParents])
				} else if (parent1[nParents]>=100){ 
						selectedParents_Geno_fam2 <- c(selectedParents_Geno_fam2,parent1[nParents])
				}
			}	
			for(nParents in 1:length(parent2)){
			
			  if(!is.na(parent2[nParents])){
			
				if(parent2[nParents]<=100){
						selectedParents_Geno_fam1 <- c(selectedParents_Geno_fam1,parent2[nParents])
				} else if (parent2[nParents]>=100){ 
						selectedParents_Geno_fam2 <- c(selectedParents_Geno_fam2,parent2[nParents])
				}
			  }

			}
		 }
			
			selectedParents_Geno_Fam1 <- unique(selectedParents_Geno_fam1)
			selectedParents_Geno_Fam2 <- (unique(selectedParents_Geno_fam2)-100)
		   
						
			selectedParents_Geno_FamPair_List[[nFamPairs]] <- list(selectedParents_Geno_Fam1,selectedParents_Geno_Fam2)
			
	    }
	
	
	
		return(list(selectedParents_Geno_List,selectedParents_NumProgeny_Geno_List,selectedParents_Geno_FamPair_List))
	
	
}
  

 
 getGM_selectedGenoList_SingleFams <- function(nextGenGenoTable_GM_List,nextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,families,GM_Param_List){
  	
### GM method for nextGen Geno data set
    
	nSnglFams <- nFamilyInd

##  GM Parameters	
	
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
	
	
	for(nFamily in 1:nSnglFams){
                 		  
		  rNames_Family1 <- c() 
		  fam1 <- families[nFamily]
		 
		  for(indPairs in 1:nrow(nextGenGenoTable_GM_List[[nFamily]])){ 
		
		    rNames_Family1 <-c(rNames_Family1,paste(nFamily," Ind",indPairs))
		  }
		
		  rownames(nextGenGenoTable_GM_List[[nFamily]]) <- rNames_Family1
		 		
         
		  M <- (nextGenGenoTable_GM_List[[nFamily]])
		  K <- GenomicMatingV2::Amat.pieces(M,pieces=15,mc.cores=no_cores) 
	    
		  Markers <- nextGenGenoTableMod_GM_List[[nFamily]]
		
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

			selectedParents_Geno_List[[nFamily]] <- selectedParents_Geno
		
			selectedParents_NumProgeny_Geno_List[[nFamily]] <- as.vector(selectedParents_NumProgeny_Geno)
			
		}else if(selectionOnGeno==FALSE){
		
## selection on phenotypic values 
			
			baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutions(Markers,Markers2=NULL,K,markerEffects_Pheno,minparents=Minparents,impinbreedstepsize=Impinbreedstepsize,impvar=Impvar,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,impforinbreed=Impforinbreed,plotiters=Plotiters,nelite=Nelite,method=Method,plotMates=PlotMates)
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Pheno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
			
			tableselected_Geno <- table(factor(unlist(c(baseGASolns_GenoModel[[1]])), levels=1:N))
		    tableselected2_Geno<- table(factor(paste(baseGASolns_GenoModel[[1]][,1],baseGASolns_GenoModel[[1]][,2], sep="_")))
		    tableselected2_Geno<- sort(tableselected2_Geno[tableselected2_Geno>0], decreasing=T)
	 
			selectedParents_Geno <- as.numeric(unlist(strsplit(names(tableselected2_Geno),"_")))
			selectedParents_NumProgeny_Geno <- c() 
		
			for(nGeno in 1:length(tableselected2_Geno)){
				selectedParents_NumProgeny_Geno <- c(selectedParents_NumProgeny_Geno,rep(tableselected2_Geno[nGeno],2))
			} 
		    selectedParents_Geno_List[[nFamily]] <- selectedParents_Geno
			selectedParents_NumProgeny_Geno_List[[nFamily]] <- as.vector(selectedParents_NumProgeny_Geno)
		}
		
      }
		
        return(list(selectedParents_Geno_List,selectedParents_NumProgeny_Geno_List))
  
   }
	
#####################################################################


 getGM_selectedGenoList_SingleFams_V2 <- function(nextGenGenoTable_GM_List,nextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,modelType,selectionOnGeno,no_cores,nFamilies,families,GM_Param_List){
  	
### GM method for nextGen Geno data set
    
	nSnglFams <- nFamilies

##  GM Parameters	
	
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
	
	
	for(nFamily in 1:nSnglFams){
                 		  
		  rNames_Family1 <- c() 
		  fam1 <- families[nFamily]
		 
		  for(indPairs in 1:nrow(nextGenGenoTable_GM_List[[nFamily]])){ 
		
		    rNames_Family1 <-c(rNames_Family1,paste(nFamily," Ind",indPairs))
		  }
		
		  rownames(nextGenGenoTable_GM_List[[nFamily]]) <- rNames_Family1
		 		
         
		  M <- (nextGenGenoTable_GM_List[[nFamily]])
		  K <- Amat.pieces(M,pieces=15,mc.cores=no_cores) 
	    
		  Markers <- nextGenGenoTableMod_GM_List[[nFamily]]
		
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
			baseGASolns_GenoModel <- getGaSolutions(Markers,Markers2=NULL,K,markerEffects_Geno,minparents=Minparents,impinbreedstepsize=Impinbreedstepsize,impvar=Impvar,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,impforinbreed=Impforinbreed,plotiters=Plotiters,nelite=Nelite,method=Method,plotMates=PlotMates)
			
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

			selectedParents_Geno_List[[nFamily]] <- selectedParents_Geno
		
			selectedParents_NumProgeny_Geno_List[[nFamily]] <- as.vector(selectedParents_NumProgeny_Geno)
			
		}else if(selectionOnGeno==FALSE){
		
## selection on phenotypic values 
			
			baseGASolns_GenoModel <- getGaSolutions(Markers,Markers2=NULL,K,markerEffects_Pheno,minparents=Minparents,impinbreedstepsize=Impinbreedstepsize,impvar=Impvar,mc.cores=Mc.cores,mutprob=Mutprob,npopGA=NpopGA,nitGA=NitGA,impforinbreed=Impforinbreed,plotiters=Plotiters,nelite=Nelite,method=Method,plotMates=PlotMates)
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Pheno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
			
			tableselected_Geno <- table(factor(unlist(c(baseGASolns_GenoModel[[1]])), levels=1:N))
		    tableselected2_Geno<- table(factor(paste(baseGASolns_GenoModel[[1]][,1],baseGASolns_GenoModel[[1]][,2], sep="_")))
		    tableselected2_Geno<- sort(tableselected2_Geno[tableselected2_Geno>0], decreasing=T)
	 
			selectedParents_Geno <- as.numeric(unlist(strsplit(names(tableselected2_Geno),"_")))
			selectedParents_NumProgeny_Geno <- c() 
		
			for(nGeno in 1:length(tableselected2_Geno)){
				selectedParents_NumProgeny_Geno <- c(selectedParents_NumProgeny_Geno,rep(tableselected2_Geno[nGeno],2))
			} 
		    selectedParents_Geno_List[[nFamily]] <- selectedParents_Geno
			selectedParents_NumProgeny_Geno_List[[nFamily]] <- as.vector(selectedParents_NumProgeny_Geno)
		}
		
      }
		
        return(list(selectedParents_Geno_List,selectedParents_NumProgeny_Geno_List))
  
   }
	
#####################################################################
	 
 getGM_selectedGenoList <- function(nextGenGenoTable_GM_List,nextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,selectionOnGeno,no_cores,nFamilies){

	nextGenGenoTable_GM_List <- nextGenGenoTable_GM_List 
	nextGenGenoTableMod_GM_List <- nextGenGenoTableMod_GM_List 
	
	selectedParents_Geno_List <- list()
	selectedParents_Pheno_List <- list()
	  
	selectedParents_NumProgeny_Geno_List <- list()
	selectedParents_NumProgeny_Pheno_List <- list() 
		
	
	for(nFamily in 1:nFamilies){

        M <- (nextGenGenoTable_GM_List[[nFamily]])
        K <- GenomicMatingV2::Amat.pieces(M,pieces=15,mc.cores=no_cores) 
	    
		
		
		Markers <- (nextGenGenoTableMod_GM_List[[nFamily]])
		
		markerEffects_Geno <- as.vector(unlist(PredictionModel_Geno[[3]]))
		markerEffects_Pheno <- as.vector(unlist(PredictionModel_Pheno[[3]]))
			 
				
## selection on genotypic values 
		
		if(selectionOnGeno ==TRUE){
			baseGASolns_GenoModel <- GenomicMatingV2::getGaSolutions(Markers,K,markerEffects_Geno,minparents=1,impinbreedstepsize=.06,impvar=.001,mc.cores=no_cores,mutprob=0.01,npopGA=500,nitGA=40,impforinbreed=10,plotiters=FALSE,nelite=5) 
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Geno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
		
			tableselected_Geno <-table(factor(unlist(c(baseGASolns_GenoModel[[1]])), levels=1:N))
			tableselected2_Geno<-table(factor(paste(baseGASolns_GenoModel[[1]][,1],baseGASolns_GenoModel[[1]][,2], sep="_")))
			tableselected2_Geno <- sort(tableselected2_Geno[tableselected2_Geno>0], decreasing=T)
	 	 	selectedParents_Geno <- as.numeric(unlist(strsplit(names(tableselected2_Geno),"_")))
			selectedParents_NumProgeny_Geno <- c() 
		
			for(nGeno in 1:length(tableselected2_Geno)){
				selectedParents_NumProgeny_Geno <- c(selectedParents_NumProgeny_Geno,rep(tableselected2_Geno[nGeno],2))
			} 

			selectedParents_Geno_List[[nFamily]] <- selectedParents_Geno
		
			selectedParents_NumProgeny_Geno_List[[nFamily]] <- as.vector(selectedParents_NumProgeny_Geno)
			
		}else if(selectionOnGeno==FALSE){
		
		## selection on phenotypic values 
			baseGASolns_PhenoModel <- GenomicMatingV2::getGaSolutions(Markers,K,markerEffects_Pheno,minparents=1,impinbreedstepsize=.06,impvar=.001,mc.cores=no_cores,mutprob=0.01,npopGA=500,nitGA=40,impforinbreed=10,plotiters=FALSE,nelite=5) 
			
			
			N <- nrow(M)
			ebvs = Markers %*% markerEffects_Geno
			names(ebvs)<-colnames(K)<-rownames(K)<-paste("l",1:N,sep="")
			
			## gets the list of pairs 
			tableselected_Pheno <- table(factor(unlist(c(baseGASolns_PhenoModel[[1]])), levels=1:N))
		    tableselected2_Pheno<- table(factor(paste(baseGASolns_PhenoModel[[1]][,1],baseGASolns_PhenoModel[[1]][,2], sep="_")))
		    tableselected2_Pheno<- sort(tableselected2_Pheno[tableselected2_Pheno>0], decreasing=T)
	 
			selectedParents_Pheno <- as.numeric(unlist(strsplit(names(tableselected2_Pheno),"_")))
			selectedParents_NumProgeny_Pheno <- c() 
		
		    ## gets the list of num_progeny per pair
			for(nPheno in 1:length(tableselected2_Pheno)){
				selectedParents_NumProgeny_Pheno <- c(selectedParents_NumProgeny_Pheno,rep(tableselected2_Pheno[nPheno],2))
			} 
		    selectedParents_Pheno_List[[nFamily]] <- selectedParents_Pheno
			selectedParents_NumProgeny_Pheno_List[[nFamily]] <- as.vector(selectedParents_NumProgeny_Pheno)
			
			
		}
		
	}
	
	
		return(list(selectedParents_Pheno_List,selectedParents_NumProgeny_Pheno_List))
	
}
	
################################################################################################################### 


####################

ApplySelection_Discrete_Frontier <- function(PhenoValues_List,no_selected,n_Families,n_Parameters){
    nParameters <- n_Parameters
    nSel <- no_selected
    phenoValues_List <- PhenoValues_List
	nFamilies <- n_Families 
	selected_PG_table_list <- rep(list(rep(list(list()),nParameters)),nFamilies)
 
    if(nSel ==20){
		
		nSel_inGroup <- 2

	    for(nFamily in 1:nFamilies){
		
		  for(nParameter in 1:nParameters){

			sortedPhenoValues <- sort(phenoValues_List[[nFamily]][[nParameter]],decreasing=TRUE,index.return=TRUE)
			topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inGroup]
			GenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inGroup]

			selected_PG_table<- cbind(topNPhenoValues,GenoTableIndices_topN)
			selected_PG_table_list[[nFamily]][[nParameter]] <- selected_PG_table
		  }
		}
		
	}else if(nSel!=20){
		nSel_inGroup <- nSel/20
		for(nFamily in 1:nFamilies){
	  	  for(nParameter in 1:nParameters){

			sortedPhenoValues <- sort(phenoValues_List[[nFamily]][[nParameter]],decreasing=TRUE,index.return=TRUE)
	   	    topNPhenoValues<- sortedPhenoValues[[1]][1:nSel_inGroup]
		    GenoTableIndices_topN <-  sortedPhenoValues[[2]][1:nSel_inGroup]

		    selected_PG_table<- cbind(topNPhenoValues,GenoTableIndices_topN)
		    selected_PG_table_list[[nFamily]][[nParameter]] <- selected_PG_table
		   } 
		}
	}
 

	return(selected_PG_table_list)

 }
 



######################################################################

####### GM Frontier method for optimal cross selection

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
			indicesRange2 <- sample(indices2,10,replace=TRUE) 
			
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
			indicesRange2 <- sample(indices2,10,replace=TRUE) 
			
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
			indicesRange2 <- sample(indices2,10,replace=TRUE) 
			
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
			indicesRange2 <- sample(indices2,10,replace=TRUE) 
			
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
		
	
	return(list(selectedParents_Geno_List,selectedCriterion_List,nParameters))
	
	
} 



############################

getCriterionValueList_GM_V2 <- function(genoValues_Fam,phenoValues_Fam,genoValSimValues_Fam,phenoValSimValues_Fam,cycle_nProgeny_Fam,nFamilies)
 {

####       
	    initIndex <- 1
	    finalIndex <- cycle_nProgeny_Fam
	    genoValues_Fam_List <- list()
	    phenoValues_Fam_List <- list()
	    genoValSimValues_Fam_List <- list()
	    phenoValSimValues_Fam_List <- list()

	    for(nFamily in 1:nFamilies){

			genoValues_Fam_List[[nFamily]] <- (genoValues_Fam)[initIndex:finalIndex]
			phenoValues_Fam_List[[nFamily]] <- (phenoValues_Fam)[initIndex:finalIndex]
			genoValSimValues_Fam_List[[nFamily]] <- genoValSimValues_Fam[initIndex:finalIndex]
			phenoValSimValues_Fam_List[[nFamily]] <- phenoValSimValues_Fam[initIndex:finalIndex]
			initIndex <- finalIndex+1
			finalIndex <- finalIndex+cycle_nProgeny_Fam

	    }
	
	     

		return(list(genoValues_Fam_List,phenoValues_Fam_List,genoValSimValues_Fam_List,phenoValSimValues_Fam_List))

}

############


##########################################################################
#### Function to select parents for next gen


ApplySelection <- function(PhenoValues,no_selected){


   nSel<-no_selected

   phenoValues<- PhenoValues

   sortedPhenoValues <- sort(phenoValues,decreasing=TRUE,index.return=TRUE)

    topNPhenoValues<- sortedPhenoValues[[1]][1:nSel]
    GenoTableIndices_topN <-  sortedPhenoValues[[2]][1:nSel]

    selected_PG_table<- cbind(topNPhenoValues,GenoTableIndices_topN)

  
    return(selected_PG_table)

 }

# testSelection <- ApplySelection(genoValues,20)

######################################################################################################

#########################################################################
#### Function to select parents for next gen for discrete selection

######## Fn 11

ApplySelection_Discrete <- function(PhenoValues_List,no_selected){

    nSel <- no_selected
    selected_PG_table_list <- list()
    phenoValues_List<- PhenoValues_List


    if(nSel ==20) {

		# nSel_inGroup <- nSel/20
		nSel_inGroup <- 2

	    for( i in 1:20){

			sortedPhenoValues <- sort(phenoValues_List[[i]],decreasing=TRUE,index.return=TRUE)
			topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inGroup]
			GenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inGroup]

			selected_PG_table<- cbind(topNPhenoValues,GenoTableIndices_topN)
			selected_PG_table_list[[i]] <- selected_PG_table
		}
	}

    else if(nSel!=20){
    nSel_inGroup <- nSel/20
    for(i in 1:20){

	    sortedPhenoValues <- sort(phenoValues_List[[i]],decreasing=TRUE,index.return=TRUE)
	   	topNPhenoValues<- sortedPhenoValues[[1]][1:nSel_inGroup]
		GenoTableIndices_topN <-  sortedPhenoValues[[2]][1:nSel_inGroup]

		selected_PG_table<- cbind(topNPhenoValues,GenoTableIndices_topN)
		selected_PG_table_list[[i]] <- selected_PG_table

	}
 }

	return(selected_PG_table_list)

 }

 
#########################################################################
## Fn 12:
#### Function to select parents for next gen for discrete selection

ApplySelection_IslandModel <- function(PhenoValues_List,no_selected,migrate=FALSE){


   nSel<-no_selected

   selected_PG_table_list <- list()

   phenoValues<- PhenoValues_List

   if(nSel ==20) {

		# nSel_inGroup <- nSel/20
		nSel_inGroup <- 2

	    for( i in 1:20){

			sortedPhenoValues <- sort(phenoValues_List[[i]],decreasing=TRUE,index.return=TRUE)
			topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inGroup]
			GenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inGroup]

			selected_PG_table<- cbind(topNPhenoValues,GenoTableIndices_topN)
			selected_PG_table_list[[i]] <- selected_PG_table
		}

	}

  if(nSel ==60){

	nSel_inGroup <- nSel/20

    for(i in 1:20){

		sortedPhenoValues <- sort(phenoValues_List[[i]],decreasing=TRUE,index.return=TRUE)
	   	topNPhenoValues<- sortedPhenoValues[[1]][1:nSel_inGroup]
		GenoTableIndices_topN <-  sortedPhenoValues[[2]][1:nSel_inGroup]

		selected_PG_table<- cbind(topNPhenoValues,GenoTableIndices_topN)
		selected_PG_table_list[[i]] <- selected_PG_table
    }
  }

  if(nSel == 200){

	nSel_inGroup <- nSel/20

	for(i in 1:20){

	sortedPhenoValues <- sort(phenoValues,decreasing=TRUE,index.return=TRUE)

		topNPhenoValues<- sortedPhenoValues[[1]][1:nSel_inGroup]
		GenoTableIndices_topN <-  sortedPhenoValues[[2]][1:nSel_inGroup]

		selected_PG_table<- cbind(topNPhenoValues,GenoTableIndices_topN)
		selected_PG_table_list[[i]] <- selected_PG_table

	}

  }


  return(selected_PG_table_List)

 }

# testSelection <- ApplySelection(genoValues,20)
