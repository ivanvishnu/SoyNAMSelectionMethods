### Exchange and Migration Methods

#########################
 
extractGenoData_FamilyWise_GM <- function(GenoTableIndices1,ProgenyGenoData){
  
    GenoTableIndices <- GenoTableIndices1
	progenyGenoData <- ProgenyGenoData
	
	
	
	nSelectedIndividuals <- length(GenoTableIndices) 
  
	genoData <- array(0,c(nSelectedIndividuals,nMarkers,2))
	
	
				
	for(j in 1:nSelectedIndividuals){
		
			individualIndex <- GenoTableIndices[j]
			
			genoData[j,,] <- progenyGenoData[,,individualIndex] 
		
	}
		
	 dimnames(genoData)[[1]] <- as.vector(GenoTableIndices)
		
	  
	return(genoData) 
	
}
        
# ###############################################################################################3


## Fn 

  # exchangeGenoData_GM <- function(Cycle_Progeny_F5_List_Fam,genoValSimValues_Fam_List,MigrationSize,EmigrantGroups,ImmigrantGroups,Direction,Policy){

		# nMarkers <-4289
		# selected_Geno_table_List <- Cycle_Progeny_F5_List_Fam
				
		# # selected_PG_table_List <- Selected_PG_table_List
		# migrationSize <- MigrationSize

	  	# emigrantGroups<- EmigrantGroups
		# immigrantGroups <- ImmigrantGroups
 
		# direction <- Direction
		# nInFamilies <- length(immigrantGroups)
		
		# FamilyInfo <- c()

		# if(Policy=="BestIsland" || Policy=="RandomBest"){ 
				# nOutFamilies <- 1
						
		# }else if(Policy=="GeneticClusters"){
				# nOutFamilies <- length(emigrantGroups)
				# if(nOutFamilies >= nSel_inFamily){ 
					# nOutFamilies <- nSel_inFamily -1
				# }
		# }else if(Policy=="FullyConnected"){ 
				# nOutFamilies <- (nSel_inFamily/2)
		# }
			
			
		# immigrationSize_per_Family <- migrationSize/nOutFamilies
		# emigrationSize_per_Family <- migrationSize/nOutFamilies
				
		# for(i in 1:nInFamilies){
			
			   					
				# inFamily <- as.numeric(immigrantGroups[i]) 
				
				# inFamily_replaceIndices <- sort(genoValSimValues_Fam_List[[inFamily]],decreasing=TRUE,index.return=TRUE)[[2]][11:(emigrationSize_per_Family+immigrationSize_per_Family)]
				# inFamily_emigrateIndices <- sort(genoValSimValues_Fam_List[[inFamily]],decreasing=TRUE,index.return=TRUE)[[2]][1:emigrationSize_per_Family]
				
								
				# for(j in 1:nOutFamilies){ 
						
										
					# outFamily <- as.numeric(sample(emigrantGroups,1)) 
					# outFamily_emigrateIndices <- sort(genoValSimValues_Fam_List[[outFamily]],decreasing=TRUE,index.return=TRUE)[[2]][1:emigrationSize_per_Family]
					# outFamily_replaceIndices <- sort(genoValSimValues_Fam_List[[outFamily]],decreasing=TRUE,index.return=TRUE)[[2]][11:(emigrationSize_per_Family+immigrationSize_per_Family)]
						
					# FamilyInfo <- rbind(FamilyInfo,c(inFamily,outFamily))
			
				    # #selected_Geno_table_List
				
				    # #Cycle_Progeny_F5_List_Fam[[inFamily]][,,nProgIndex]
			
# ## RM includes only the first and second individual in most cases
				

					# if(direction ==2){
					
						# selected_Geno_table_List[[outFamily]][,,outFamily_replaceIndices] <- selected_Geno_table_List[[inFamily]][,,inFamily_emigrateIndices]
						# selected_Geno_table_List[[inFamily]][,,inFamily_replaceIndices] <- selected_Geno_table_List[[outFamily]][,,outFamily_emigrateIndices]
					# }
					# if(direction==1){
					    # selected_Geno_table_List[[inFamily]][,,inFamily_replaceIndices] <- selected_Geno_table_List[[outFamily]][,,outFamily_emigrateIndices]
					
					# }
				  			
			    # }
	 	# }
		# colnames(FamilyInfo) <- c("In-Family","Out-Family")

		# return(list((selected_Geno_table_List),as.matrix(FamilyInfo)))

# }


###### 06/03/2020

exchangeGenoData_GM <- function(Cycle_Progeny_F5_List_Fam, genoValSimValues_Fam_List,
    MigrationSize, EmigrantGroups, ImmigrantGroups, Direction,Policy)
{   
    nMarkers <- 4289
    selected_Geno_Table_List <- Cycle_Progeny_F5_List_Fam
    sorted_Selected_Geno_Table_List <-  selected_Geno_Table_List
    nFamilies <- length(Cycle_Progeny_F5_List_Fam)
	 
    migrationSize <- MigrationSize
	nSel_inFamily <- nFamilies/(migrationSize)
	
    emigrantGroups <- EmigrantGroups
    immigrantGroups <- ImmigrantGroups
    direction <- Direction
    nInFamilies <- length(immigrantGroups)
    selCriterion_List <- list()

    FamilyInfo <- c()
    for(nFamily in 1:nFamilies){
    	sortedIndices <- sort(genoValSimValues_Fam_List[[nFamily]],decreasing = TRUE, index.return = TRUE)[[2]]
   	    sorted_Selected_Geno_Table_List[[nFamily]] <- selected_Geno_Table_List[[nFamily]][,,sortedIndices]
        selCriterion_List[[nFamily]] <- genoValSimValues_Fam_List[[nFamily]][sortedIndices]
    }


    if(migrationSize == 1){
        if (Policy == "BestIsland" || Policy == "RandomBest") {
            nOutFamilies <- 1
        }
        else if (Policy == "FullyConnected") {
            nOutFamilies <- (nSel_inFamily/2)
        }

        for (i in 1:nInFamilies) {
            inFamily <- as.numeric(immigrantGroups[i])
            for (j in 1:nOutFamilies) {
                replaceSlots <- j + 1
                finalIndex <- replaceSlots
                outFamily <- as.numeric(sample(emigrantGroups,
                  1))
                inIndex <- 1
                outIndex <- 1
                FamilyInfo <- rbind(FamilyInfo, c(inFamily, outFamily))
                if ((selCriterion_List[[inFamily]][finalIndex] <= selCriterion_List[[outFamily]][1])) {
                  if (direction == 2) {
                    replaceSlot1_OutFamily_GenoData <- sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex]
                    sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex] <- sorted_Selected_Geno_Table_List[[inFamily]][,,inIndex ]
                    sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex] <- sorted_Selected_Geno_Table_List[[outFamily]][,,outIndex]
                  }
                  if (direction == 1) {
                    sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex] <- sorted_Selected_Geno_Table_List[[outFamily]][,,outIndex]
                  }
                }
                else if (selCriterion_List[[inFamily]][finalIndex] > selCriterion_List[[outFamily]][1]) {
                  sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex] <- sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex]
                }
            }
        }
    }
    if (migrationSize == 2) {
        if (Policy == "BestIsland" || Policy == "RandomBest") {
            nOutFamilies <- 1
        }
        else if (Policy == "FullyConnected") {
            nOutFamilies <- (nSel_inFamily/2)
        }
        for (i in 1:nInFamilies) {
            inFamily <- as.numeric(immigrantGroups[i])
            for (j in 1:nOutFamilies) {
                replaceSlots <- j + 1
                finalIndex1 <- replaceSlots
                finalIndex2 <- finalIndex1 + 1
                outFamily <- as.numeric(sample(emigrantGroups,1))
                 inIndex1 <- 1
                inIndex2 <- 2
                outIndex1 <- 1
                outIndex2 <- 2
                FamilyInfo <- rbind(FamilyInfo, c(inFamily, outFamily))
                if ((selCriterion_List[[inFamily]][finalIndex1] <= selCriterion_List[[outFamily]][1])) {
                  if (direction == 2){
                    replaceSlot1_OutFamily_GenoData <- sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex1]
                    sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex1] <- sorted_Selected_Geno_Table_List[[inFamily]][,,inIndex1]
                    replaceSlot2_OutFamily_GenoData <- sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex2]
                    sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex2] <- sorted_Selected_Geno_Table_List[[inFamily]][,,inIndex2]
                    sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex1] <- replaceSlot1_OutFamily_GenoData   
          	        sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex2] <- replaceSlot2_OutFamily_GenoData
                   }
                  if (direction == 1) {
                    sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex1] <- sorted_Selected_Geno_Table_List[[outFamily]][,,outIndex1]
                    sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex2] <- sorted_Selected_Geno_Table_List[[outFamily]][,,outIndex2]
                  }
                }
				else if (selCriterion_List[[inFamily]][finalIndex1] > selCriterion_List[[outFamily]][1]) {
                  sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex1] <- sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex1]
                  sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex2] <- sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex2]
                }
            }
        }
    }
    colnames(FamilyInfo) <- c("In-Family", "Out-Family")
    return(list((sorted_Selected_Geno_Table_List), as.matrix(FamilyInfo)))
} 


# #################################

exchangeGenoData_GM_Sel <-function (Selected_Cycle_Progeny_F5_List_Fam,Selected_SelCriterion_List,
    SelectedGenoIndividualIndices,MigrationSize,EmigrantGroups,
    ImmigrantGroups,Direction,Policy)
{ 
    nMarkers <- 4289
    selected_Geno_Table_List <- Cycle_Progeny_F5_List_Fam
    sorted_Selected_Geno_Table_List <-  selected_Geno_Table_List
    nFamilies <- length(Cycle_Progeny_F5_List_Fam)
	 
    migrationSize <- MigrationSize
	nSel_inFamily <- nFamilies/(migrationSize)
	
    emigrantGroups <- EmigrantGroups
    immigrantGroups <- ImmigrantGroups
    direction <- Direction
    nInFamilies <- length(immigrantGroups)
    selCriterion_List <- list()

    FamilyInfo <- c()
    for(nFamily in 1:nFamilies){
    	sortedIndices <- sort(genoValSimValues_Fam_List[[nFamily]],decreasing = TRUE, index.return = TRUE)[[2]]
   	    sorted_Selected_Geno_Table_List[[nFamily]] <- selected_Geno_Table_List[[nFamily]][,,sortedIndices]
        selCriterion_List[[nFamily]] <- genoValSimValues_Fam_List[[nFamily]][sortedIndices]
    }


    if(migrationSize == 1){
        if (Policy == "BestIsland" || Policy == "RandomBest") {
            nOutFamilies <- 1
        }
        else if (Policy == "FullyConnected") {
            nOutFamilies <- (nSel_inFamily/2)
        }

        for (i in 1:nInFamilies) {
            inFamily <- as.numeric(immigrantGroups[i])
            for (j in 1:nOutFamilies) {
                replaceSlots <- j + 1
                finalIndex <- replaceSlots
                outFamily <- as.numeric(sample(emigrantGroups,
                  1))
                inIndex <- 1
                outIndex <- 1
                FamilyInfo <- rbind(FamilyInfo, c(inFamily, outFamily))
                if ((selCriterion_List[[inFamily]][finalIndex] <= selCriterion_List[[outFamily]][1])) {
                  if (direction == 2) {
                    replaceSlot1_OutFamily_GenoData <- sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex]
                    sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex] <- sorted_Selected_Geno_Table_List[[inFamily]][,,inIndex ]
                    sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex] <- sorted_Selected_Geno_Table_List[[outFamily]][,,outIndex]
                  }
                  if (direction == 1) {
                    sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex] <- sorted_Selected_Geno_Table_List[[outFamily]][,,outIndex]
                  }
                }
                else if (selCriterion_List[[inFamily]][finalIndex] > selCriterion_List[[outFamily]][1]) {
                  sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex] <- sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex]
                }
            }
        }
    }
    if (migrationSize == 2) {
        if (Policy == "BestIsland" || Policy == "RandomBest") {
            nOutFamilies <- 1
        }
        else if (Policy == "FullyConnected") {
            nOutFamilies <- (nSel_inFamily/2)
        }
        for (i in 1:nInFamilies) {
            inFamily <- as.numeric(immigrantGroups[i])
            for (j in 1:nOutFamilies) {
                replaceSlots <- j + 1
                finalIndex1 <- replaceSlots
                finalIndex2 <- finalIndex1 + 1
                outFamily <- as.numeric(sample(emigrantGroups,
                  1))
                 inIndex1 <- 1
                inIndex2 <- 2
                outIndex1 <- 1
                outIndex2 <- 2
                FamilyInfo <- rbind(FamilyInfo, c(inFamily, outFamily))
                if ((selCriterion_List[[inFamily]][finalIndex1] <= selCriterion_List[[outFamily]][1])) {
                  if (direction == 2){
                    replaceSlot1_OutFamily_GenoData <- sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex1]
                    sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex1] <- sorted_Selected_Geno_Table_List[[inFamily]][,,inIndex1]
                    replaceSlot2_OutFamily_GenoData <- sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex2]
                    sorted_Selected_Geno_Table_List[[outFamily]][,,finalIndex2] <- sorted_Selected_Geno_Table_List[[inFamily]][,,inIndex2]
                    sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex1] <- replaceSlot1_OutFamily_GenoData   
          	        sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex2] <- replaceSlot2_OutFamily_GenoData
                   }
                  if (direction == 1) {
                    sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex1] <- sorted_Selected_Geno_Table_List[[outFamily]][,,outIndex1]
                    sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex2] <- sorted_Selected_Geno_Table_List[[outFamily]][,,outIndex2]
                  }
                }
				else if (selCriterion_List[[inFamily]][finalIndex1] > selCriterion_List[[outFamily]][1]) {
                  sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex1] <- sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex1]
                  sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex2] <- sorted_Selected_Geno_Table_List[[inFamily]][,,finalIndex2]
                }
            }
        }
    }
    colnames(FamilyInfo) <- c("In-Family", "Out-Family")
    return(list((sorted_Selected_Geno_Table_List), as.matrix(FamilyInfo)))
}  





 migrationPolicy_GM <- function(GenoTable,SelCriterion_List,Policy,nFamilies){ 
 
	selected_PG_List <-SelCriterion_List
	nSel_inFamily <- length(selected_PG_List[[1]])
	
	if(Policy=="BestIsland"){ 
	 	 
     	maxPhenoValues <- rep(0,nFamilies)

		for(nFamily in 1:nFamilies){

			maxPhenoValues[nFamily] <- max(selected_PG_List[[nFamily]])

		}

		sortedGroups<- (sort.int(maxPhenoValues,decreasing=TRUE,index.return=TRUE))[[2]]

		emigrantGroups <- sortedGroups[1:10]
		immigrantGroups <- sortedGroups[11:20]
		
		OutGroup <- emigrantGroups[1] 
		InGroup <- immigrantGroups

    
	}else if(Policy=="RandomBest"){  
	
			
		maxPhenoValues <- rep(0,nFamilies)

		for(nFamily in 1:nFamilies){

			maxPhenoValues[nFamily] <- max(selected_PG_List[[nFamily]])

		}

		sortedGroups<- (sort.int(maxPhenoValues,decreasing=TRUE,index.return=TRUE))[[2]]

		emigrantGroups <- sortedGroups[1:10]
		immigrantGroups <- sortedGroups[11:20]
		
		OutGroup <- emigrantGroups
		InGroup <- immigrantGroups
	
	} else if(Policy=="GeneticClusters"){
  
	  populations <- rep(0,20*100)
	  init <- 1
		for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1

		}

       
	  d_EucDist<- as.matrix(dist(GenoTable))
	  d_EucMat <- as.matrix(d_EucDist)
	  clusterGen_Euc <- hclust(as.dist(d_EucMat)) 
	  c.Tree <- c.Tree <- cutree(clusterGen_Euc,k=5) 
	  
	 	  
	  Group <- c() 
	  for(i in 1:nFamilies){ 
 
		Group <- c(Group,rep(i,100))
	  }
	   
	  l <- length(levels(factor(c.Tree))) 
	  groups <- list()
	  for(i in 1:l){ 
	 
		indices <- which(c.Tree %in% i) 
		
		groups[[i]] <- levels(factor(Group[indices]))
	 
	  }
	 
	  nGroups<- length(groups)
	  groupIndices <- c(1:nGroups)
	  
	  
	  outGroupIndex <- sample(groupIndices,1)
	  inGroupIndex <- sample(setdiff(groupIndices,outGroupIndex),1)
	  
	  OutGroup<- groups[[outGroupIndex]]
	  InGroup <- groups[[inGroupIndex]]
	 
	  	  
	 }else if(Policy=="FullyConnected"){ 
	 
		FCGroups<- combn(c(1:20),2)
		indices<- sample(c(1:dim(FCGroups)[2]),10)
		OutGroup <- FCGroups[1,indices]
		InGroup <-  FCGroups[2,indices]
		
	}
	  
	 
	  emigrantGroups<- OutGroup
	  immigrantGroups <- InGroup
 
	  direction <- Direction
	  nInFamilies <- length(immigrantGroups)
		
	  FamilyInfo <- c()


	    if(Policy=="BestIsland"){
			nOutFamilies <- 1
					
			for(i in 1:nInFamilies){

				inFamily <- as.numeric(immigrantGroups[i])
				outFamily <- as.numeric(emigrantGroups[1])
				FamilyInfo <- rbind(FamilyInfo,c(inFamily,outFamily))
			}
		}else if(Policy=="RandomBest"){ 
		
			nOutFamilies <- 1
					
			for(i in 1:nInFamilies){

			    for(j in 1:nOutFamilies){
				inFamily <- as.numeric(immigrantGroups[i])
				outFamily <- as.numeric(sample(emigrantGroups,1))
				FamilyInfo <- rbind(FamilyInfo,c(inFamily,outFamily))
				}
			}
				
		}else if(Policy=="GeneticClusters"){
			nOutFamilies <- length(emigrantGroups)
		
		
		    for(i in 1:nInFamilies){

			inFamily <- as.numeric(immigrantGroups[i])
			
				for(j in 1:nOutFamilies){ 
					outFamily <- as.numeric(sample(emigrantGroups,1))
					FamilyInfo <- rbind(FamilyInfo,c(inFamily,outFamily))
				}
		    }
		}else if(Policy=="FullyConnected"){ 
			nOutFamilies <- length(emigrantGroups)
			
			FamilyInfo <- cbind(immigrantGroups,emigrantGroups)
				
		}
		
				
		colnames(FamilyInfo) <- c("In-Family","Out-Family")

		return(as.matrix(FamilyInfo))

}
 
 
######

 migrationPolicy_GM_V2 <- function(GenoTable,Selected_PG_table_List,Policy,nFamilies,noQTL){ 
 
        no_QTL <- noQTL
	if(Policy=="BestIsland"){ 
	 	 
        
		selected_PG_table_List <- Selected_PG_table_List
	
		maxPhenoValues <- rep(0,nFamilies)

		for(nFamily in 1:nFamilies){

			maxPhenoValues[nFamily] <- max(selected_PG_table_List[[nFamily]][,1])

		}

		sortedGroups<- (sort.int(maxPhenoValues,decreasing=TRUE,index.return=TRUE))[[2]]

		emigrantGroups <- sortedGroups[1:10]
		immigrantGroups <- sortedGroups[11:20]
		
		OutGroup <- emigrantGroups[1] 
		InGroup <- immigrantGroups

      
	
	}else if(Policy=="RandomBest"){  
	
		selected_PG_table_List <- Selected_PG_table_List
	
		maxPhenoValues <- rep(0,nFamilies)

		for(nFamily in 1:nFamilies){

			maxPhenoValues[nFamily] <- max(selected_PG_table_List[[nFamily]][,1])

		}

		sortedGroups<- (sort.int(maxPhenoValues,decreasing=TRUE,index.return=TRUE))[[2]]

		emigrantGroups <- sortedGroups[1:10]
		immigrantGroups <- sortedGroups[11:20]
		
		OutGroup <- emigrantGroups
		InGroup <- immigrantGroups

 
	} else if(Policy=="GeneticClusters"){
	
###  Get population vector

	  populations <- rep(0,20*100)
	  init <- 1
	  for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1
      }

	  genoValSimValues <- simulateGenotypicValues(GenoTable,no_QTL)   
	    
		
	  d_EucDist<- as.matrix(dist(GenoTable))
	  d_EucMat <- as.matrix(d_EucDist)
	  clusterGen_Euc <- hclust(as.dist(d_EucMat)) 
	  c.Tree <- cutree(clusterGen_Euc,k=3) 
	  
	  # d <- pairwise_Gst_Nei(GenoTable_GInd)
      # dMat <- as.matrix(d)
	  # clusterGen <- hclust(as.dist(dMat))
	  # c.Tree <- cutree(clusterGen,k=5) 
	  
	  Group <- c() 
	  for(i in 1:nFamilies){ 
 
		Group <- c(Group,rep(i,100))
	  }
	   
	  l <- length(levels(factor(c.Tree))) 
	  groups <- list()
	  group_Avg_GenoSim <- c()
	  group_SD_GenoSim <- c()
	  Std_group_Avg_GenoSim <- c()
	  for(i in 1:l){ 
	 
		indices <- which(c.Tree %in% i) 
		groups[[i]] <- levels(factor(Group[indices]))
	    group_Avg_GenoSim[i]   <- mean(genoValSimValues[indices])
        group_SD_GenoSim[i] <- sd(genoValSimValues[indices])
		
		Std_group_Avg_GenoSim[i] <- group_Avg_GenoSim[i]/length(groups)
		
	  }
	  
	  groups_sorted <- sort.int(Std_group_Avg_GenoSim,decreasing=TRUE,index.return=TRUE)
	 
	  nGroups<- length(groups)
	  groupIndices <- c(1:nGroups)
	  
	  if(length(which(group_SD_GenoSim %in% 0)) >= 1){ 
		RetainMigPolicy <- TRUE
	  }else if (length(which(group_SD_GenoSim %in% 0)) == 0){ 
	    RetainMigPolicy <- FALSE
	  }
	  
## Trial1 Outgroup and In group are randomly sampled from distance based genetic clusters
	  # outGroupIndex <- sample(groupIndices,1)
	  # inGroupIndex <- sample(setdiff(groupIndices,outGroupIndex),1)
	  
## Trial 2 Outgroup and Ingroup are groups with maximum geentic distance 
	  # Avg cluster genetic distance
      
      outGroupIndex <- groups_sorted[[2]][1]
	  inGroupIndex1 <-  groups_sorted[[2]][2]
	  inGroupIndex2 <-  groups_sorted[[2]][3]
	  OutGroup <- groups[[outGroupIndex]] 
	  
	 	
	  families_in_group <- levels(factor(OutGroup))
	  family_in_group_Avg_GenoSim <- c()
		   
	  for (k in families_in_group){ 
		   
		    familyIndices <- which(Group %in% as.numeric(k)) 
			
			family_in_group_Avg_GenoSim[k] <- mean(genoValSimValues[familyIndices])
				    
	  } 
	  
	  fraction_Avg_GenoSim <-  family_in_group_Avg_GenoSim/group_Avg_GenoSim[outGroupIndex]
	  OutGroupSampleProb <- abs(fraction_Avg_GenoSim/sum(fraction_Avg_GenoSim))
		
	  InGroup1 <- groups[[inGroupIndex1]]
	  InGroup2 <- groups[[inGroupIndex2]]
	  InGroup <- c(InGroup1,InGroup2)
		
	}else if(Policy=="FullyConnected"){ 
	 
		OutGroup <- c(1:nFamilies)
		InGroup <- c(1:nFamilies)
	  
	}
	if(Policy != "GeneticClusters"){
	  output <- list(OutGroup,InGroup) 
	} else if (Policy == "GeneticClusters"){ 
	   output <- list(OutGroup,InGroup,abs(OutGroupSampleProb),RetainMigPolicy)
	}
	  
	return(output)
 
 }


 
 

## Migration Policy Function 

## Given policy topology and number of families, returns immigrant and emigrant groups

## User defined migration size and frequency compared to a dynamic policy, where migration is triggered when the variance of an island drops below a threshold
 
 # Policy <- "BestIsland"	
 
 # Policy <-"GeneticClusters"

 # Policy <- "RandomBest"
 
 # Policy <- "FullyConnected" 
 
 # (All islands are connected to each other (island is randomly selected, but only half of the selected in any island is replaced)

    # Policy <-"GeneticClusters"
    # GenoTable <- nextGenGenoTable
    # migPol <- migrationPolicy(GenoTable,genoSelectionTable,Policy,nFamilies)

	

## Fn 45: 

 migrationPolicy <- function(GenoTable,Selected_PG_table_List,Policy,nFamilies,noQTL){ 
 
    no_QTL <- noQTL
	if(Policy=="BestIsland"){ 
	 	 
        
		selected_PG_table_List <- Selected_PG_table_List
	
		maxPhenoValues <- rep(0,nFamilies)

		for(nFamily in 1:nFamilies){

			maxPhenoValues[nFamily] <- max(selected_PG_table_List[[nFamily]][,1])

		}

		sortedGroups<- (sort.int(maxPhenoValues,decreasing=TRUE,index.return=TRUE))[[2]]

		emigrantGroups <- sortedGroups[1:10]
		immigrantGroups <- sortedGroups[11:20]
		
		OutGroup <- emigrantGroups[1] 
		InGroup <- immigrantGroups

      
	
	}else if(Policy=="RandomBest"){  
	
		selected_PG_table_List <- Selected_PG_table_List
	
		maxPhenoValues <- rep(0,nFamilies)

		for(nFamily in 1:nFamilies){

			maxPhenoValues[nFamily] <- max(selected_PG_table_List[[nFamily]][,1])

		}

		sortedGroups<- (sort.int(maxPhenoValues,decreasing=TRUE,index.return=TRUE))[[2]]

		emigrantGroups <- sortedGroups[1:10]
		immigrantGroups <- sortedGroups[11:20]
		
		OutGroup <- emigrantGroups
		InGroup <- immigrantGroups

       
	
	} else if(Policy=="GeneticClusters"){
	
	  
###  Get population vector

	  populations <- rep(0,20*100)
	  init <- 1
	  for(i in 1:20){

			final <- i*100
			populations[init:final]<- rep(i,100)

			init <- final+1
      }

	  genoValSimValues <- simulateGenotypicValues(GenoTable,no_QTL)   
	    
		
	  d_EucDist<- as.matrix(dist(GenoTable))
	  d_EucMat <- as.matrix(d_EucDist)
	  clusterGen_Euc <- hclust(as.dist(d_EucMat)) 
	  c.Tree <- cutree(clusterGen_Euc,k=3) 
	  
	  # d <- pairwise_Gst_Nei(GenoTable_GInd)
      # dMat <- as.matrix(d)
	  # clusterGen <- hclust(as.dist(dMat))
	  # c.Tree <- cutree(clusterGen,k=5) 
	  
	  Group <- c() 
	  for(i in 1:nFamilies){ 
 
		Group <- c(Group,rep(i,100))
	  }
	   
	  l <- length(levels(factor(c.Tree))) 
	  groups <- list()
	  group_Avg_GenoSim <- c()
	  Std_group_Avg_GenoSim <- c()
	  for(i in 1:l){ 
	 
		indices <- which(c.Tree %in% i) 
		groups[[i]] <- levels(factor(Group[indices]))
	    group_Avg_GenoSim[i]   <- mean(genoValSimValues[indices]) 
		Std_group_Avg_GenoSim[i] <- group_Avg_GenoSim/length(groups)
	  }
	  
	  groups_sorted <- sort.int(Std_group_Avg_GenoSim,decreasing=TRUE,index.return=TRUE)
	 
	  nGroups<- length(groups)
	  groupIndices <- c(1:nGroups)
	  
	## Trial1 Outgroup and In group are randomly sampled from distance based genetic clusters
	  # outGroupIndex <- sample(groupIndices,1)
	  # inGroupIndex <- sample(setdiff(groupIndices,outGroupIndex),1)
	  
	## Trial 2 Outgroup and Ingroup are groups with maximum geentic distance 
	  # Avg cluster genetic distance
      
      outGroupIndex <- groups_sorted[[2]][1]
	  inGroupIndex1 <-  groups_sorted[[2]][2]
	  inGroupIndex2 <-  groups_sorted[[2]][3]
	  
	  OutGroup<- groups[[outGroupIndex]]
	  
	  InGroup1 <- groups[[inGroupIndex1]]
	  InGroup2 <- groups[[inGroupIndex2]]
		  	
	  InGroup <- c(InGroup1,InGroup2)
		
	 }else if(Policy=="FullyConnected"){ 
	 
		OutGroup <- c(1:nFamilies)
		InGroup <- c(1:nFamilies)
	  
	 }
	  
	  return(list(OutGroup,InGroup))
 
 }
 
 
# Selected_GenoData_List
### Exchange selected geno data among pairs
# selectedGenoData_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction,Policy)

## Best Island - Emigrants only out of the best island in to immigrant group, which is randomly sampled for nInFamilies (10)
## Random Best- emigrants out of best island randomly selected from the top 10, with immigrant islands being selected at random 
## Fully Connected - All islands are potentially connected with all other islands (both are sampled randomly for half of the immigrant to be replaced)
## Genetic clusters 

# selectedGenoData_List_afterExchange <- exchangeGenoData(selectedGenoData_List,genoSelectionTable,1,emigrantGroups,immigrantGroups,direction,Policy)

# selected_PG_table_List <- genoSelectionTable
# migrationSize <- 1

## Fn 46: 

    exchangeGenoData <- function(Selected_GenoData_List,Selected_PG_table_List,MigrationSize,EmigrantGroups,ImmigrantGroups,Direction,Policy){

		nMarkers <-4289
		selected_Geno_table_List <- Selected_GenoData_List
		selected_PG_table_List <- Selected_PG_table_List
		migrationSize <- MigrationSize

	    nSel_inFamily <- dim(selected_PG_table_List[[1]])[1]
		emigrantGroups<- EmigrantGroups
		immigrantGroups <- ImmigrantGroups
 
		direction <- Direction
		nInFamilies <- length(immigrantGroups)
		
		FamilyInfo <- c()

		if(migrationSize==1){

			if(Policy=="BestIsland" || Policy=="RandomBest"){ 
				nOutFamilies <- 1
						
			}else if(Policy=="GeneticClusters"){
				nOutFamilies <- length(emigrantGroups)
				if(nOutFamilies >= nSel_inFamily){ 
					nOutFamilies <- nSel_inFamily -1
				}
			}else if(Policy=="FullyConnected"){ 
				nOutFamilies <- (nSel_inFamily/2)
			}
			
			
			for(i in 1:nInFamilies){
			
			   	# outFamily <- sample(emigrantGroups,1)
				# inFamily <- sample(immigrantGroups,1)
				
				inFamily <- as.numeric(immigrantGroups[i])
				
				for(j in 1:nOutFamilies){ 
					
					
					replaceSlots <- j+1
					# finalIndex <- sample(replaceSlots,1)
					finalIndex <- replaceSlots
					
					outFamily <- as.numeric(sample(emigrantGroups,1))
					
					inIndex <- 1
					outIndex <- 1
				   
					FamilyInfo <- rbind(FamilyInfo,c(inFamily,outFamily))
				
				
				## GM includes only the first and second individual in most cases
				

				if((selected_PG_table_List[[inFamily]][finalIndex,1] <= selected_PG_table_List[[outFamily]][1,1])){

					if(direction ==2){
					
						replaceSlot1_OutFamily_GenoData<-  selected_Geno_table_List[[outFamily]][finalIndex,,]
						selected_Geno_table_List[[outFamily]][finalIndex,,] <- selected_Geno_table_List[[inFamily]][inIndex,,]
						selected_Geno_table_List[[inFamily]][finalIndex,,] <- selected_Geno_table_List[[outFamily]][outIndex,,]
					}
					if(direction==1){
					   selected_Geno_table_List[[inFamily]][finalIndex,,] <- selected_Geno_table_List[[outFamily]][outIndex,,]
					}
				}else if(selected_PG_table_List[[inFamily]][finalIndex,1] > selected_PG_table_List[[outFamily]][1,1]) {

					 selected_Geno_table_List[[inFamily]][finalIndex,,] <- selected_Geno_table_List[[inFamily]][finalIndex,,]

				}
			   }
				
			}
	    }
		
		if(migrationSize==2){
		
           if(Policy=="BestIsland" || Policy=="RandomBest"){ 
				nOutFamilies <- 1
						
			}else if(Policy=="GeneticClusters"){
				nOutFamilies <- length(emigrantGroups)/2
				if(nOutFamilies >= nSel_inFamily){ 
					nOutFamilies <- (nSel_inFamily -1)/2
				}
			}else if(Policy=="FullyConnected"){ 
				nOutFamilies <- (nSel_inFamily/2)
			}
	
		for(i in 1:nInFamilies){
		
		    

			# outFamily <- sample(emigrantGroups,1)
			# inFamily <- sample(immigrantGroups,1)
			
			inFamily <- as.numeric(immigrantGroups[i])
			
			for(j in 1:nOutFamilies){ 
                
				
				replaceSlots <- j+1
				finalIndex1 <- replaceSlots
				finalIndex2 <- finalIndex1+1 
				
				outFamily <- as.numeric(sample(emigrantGroups,1))
				
				inIndex1 <- 1
				inIndex2 <- 2
				outIndex1 <- 1
			    outIndex2 <- 2
			    FamilyInfo <- rbind(FamilyInfo,c(inFamily,outFamily))
				
				
			
## GM includes only the first and second individual in most cases
			

		    if((selected_PG_table_List[[inFamily]][finalIndex1,1] <= selected_PG_table_List[[outFamily]][1,1])){

				if(direction ==2){
				
				    replaceSlot1_OutFamily_GenoData<- selected_Geno_table_List[[outFamily]][finalIndex1,,]
					selected_Geno_table_List[[outFamily]][finalIndex1,,] <- selected_Geno_table_List[[inFamily]][inIndex1,,]
					
					replaceSlot2_OutFamily_GenoData <- selected_Geno_table_List[[outFamily]][finalIndex2,,] 
					selected_Geno_table_List[[outFamily]][finalIndex2,,] <- selected_Geno_table_List[[inFamily]][inIndex2,,]
					
					selected_Geno_table_List[[inFamily]][finalIndex1,,] <- replaceSlot1_OutFamily_GenoData
				    selected_Geno_table_List[[inFamily]][finalIndex2,,] <- replaceSlot2_OutFamily_GenoData
							
				}
				if(direction==1){

				  selected_Geno_table_List[[inFamily]][finalIndex1,,] <- selected_Geno_table_List[[outFamily]][outIndex1,,]
				  selected_Geno_table_List[[inFamily]][finalIndex2,,] <- selected_Geno_table_List[[outFamily]][outIndex2,,]
				}
				
				# selected_Geno_table_List[[inFamily]][finalIndex1,,] <- replaceSlot1_OutFamily_GenoData
				# selected_Geno_table_List[[inFamily]][finalIndex2,,] <- replaceSlot2_OutFamily_GenoData
				

			}else if(selected_PG_table_List[[inFamily]][finalIndex1,1] > selected_PG_table_List[[outFamily]][1,1]) {

				 selected_Geno_table_List[[inFamily]][finalIndex1,,] <- selected_Geno_table_List[[inFamily]][finalIndex1,,]
				 selected_Geno_table_List[[inFamily]][finalIndex2,,] <- selected_Geno_table_List[[inFamily]][finalIndex2,,]

			 }
            }
			
		 }
	    }
		
			colnames(FamilyInfo) <- c("In-Family","Out-Family")

		return(list((selected_Geno_table_List),as.matrix(FamilyInfo)))

}
 
 
####

################################################################################################################## Fn 25:
## Migration of individuals
##  Migration size

    # selected_PG_table_List <- genoSelectionTable
	# MigrationSize <- 1
	# Direction <- 2

## Fn 33: 
migrate2X <- function(Selected_PG_table_List,MigrationSize,cycleNo,Direction){

		nMarkers <- 4289

		selected_PG_table_List <- Selected_PG_table_List

		migrationSize <- MigrationSize

		direction <- Direction



		maxPhenoValues <- rep(0,20)

		for(nFamily in 1:20){

			maxPhenoValues[nFamily] <- max(selected_PG_table_List[[nFamily]][,1])

		}

		sortedGroups<- (sort.int(maxPhenoValues,decreasing=TRUE,index.return=TRUE))[[2]]

		# a <- rep(0,20)

		# for(nFamily in 1:20){

			# a[nFamily] <- which(sortedGroups %in% nFamily)

		# }

		emigrantGroups <- sortedGroups[1:10]
		immigrantGroups <- sortedGroups[11:20]


		for(i in 1:10){


			outFamily <- emigrantGroups[1]
			inFamily <- immigrantGroups[i]

			finalIndex <- 2


			if(selected_PG_table_List[[inFamily]][finalIndex,1] <= selected_PG_table_List[[outFamily]][1,1]){

				if(direction==2){

				   selected_PG_table_List[[outFamily]][finalIndex,1]<- selected_PG_table_List[[inFamily]][finalIndex,1]
				   selected_PG_table_List[[outFamily]][finalIndex,2]<- selected_PG_table_List[[inFamily]][finalIndex,2]
				}


				selected_PG_table_List[[inFamily]][finalIndex,1]<- selected_PG_table_List[[outFamily]][1,1]
				selected_PG_table_List[[inFamily]][finalIndex,2]<- selected_PG_table_List[[outFamily]][1,2]

			} else {

				selected_PG_table_List[[inFamily]][finalIndex,1]<- selected_PG_table_List[[inFamily]][finalIndex,1]
				selected_PG_table_List[[inFamily]][finalIndex,2]<- selected_PG_table_List[[inFamily]][finalIndex,2]

			}

		}

		  return(list(selected_PG_table_List,emigrantGroups,immigrantGroups))

  }

 # testMigrate <- migrate2X(genoSelectionTable,1)

  # selected_PG_table_List<- genoSelectionTable
  # migrationSize <-1
  # EmigrantGroups <- emigrantGroups
  # ImmigrantGroups <- immigrantGroups
  # direction <-2


## Fn 34:
migrate <- function(Selected_PG_table_List,MigrationSize,EmigrantGroups,ImmigrantGroups,Direction){

		selected_PG_table_List <- Selected_PG_table_List
		migrationSize <- MigrationSize

		emigrantGroups<- EmigrantGroups
		immigrantGroups <- ImmigrantGroups
		direction <- Direction
		nMarkers <- 4289

		for(i in 1:10){

			# outFamily <- sample(emigrantGroups,1)
			# inFamily <- sample(immigrantGroups,1)

			outFamily <- emigrantGroups[1]
			inFamily <- immigrantGroups[i]


			finalIndex <- 2

			# if(selected_PG_table_List[[inFamily]][finalIndex,1] >= selected_PG_table_List[[outFamily]][1,1]){

				# inFamilyTemp <-inFamily
				# inFamily <- outFamily
				# outFamily <- inFamilyTemp
			# }

		    if(selected_PG_table_List[[inFamily]][finalIndex,1] <= selected_PG_table_List[[outFamily]][1,1]){

				if(direction==2){

				   selected_PG_table_List[[outFamily]][finalIndex,1]<- selected_PG_table_List[[inFamily]][finalIndex,1]
				   selected_PG_table_List[[outFamily]][finalIndex,2]<- selected_PG_table_List[[inFamily]][finalIndex,2]
				}

				selected_PG_table_List[[inFamily]][finalIndex,1]<- selected_PG_table_List[[outFamily]][1,1]
				selected_PG_table_List[[inFamily]][finalIndex,2]<- selected_PG_table_List[[outFamily]][1,2]

			}else {

				selected_PG_table_List[[inFamily]][finalIndex,1]<- selected_PG_table_List[[inFamily]][finalIndex,1]
				selected_PG_table_List[[inFamily]][finalIndex,2]<- selected_PG_table_List[[inFamily]][finalIndex,2]
			}

		}

		  return(selected_PG_table_List)

}

#####

####
# Selected_PG_table_List <- genoSelectionTable
# MigrationSize <- 1
# EmigrantGroups <- emigrantGroups
# ImmigrantGroups <- immigrantGroups
# Direction <-2

## Fn 35:
migrateN_v1 <- function(Selected_PG_table_List,MigrationSize,EmigrantGroups,ImmigrantGroups,Direction){

		nMarkers <-4289
		selected_PG_table_List <- Selected_PG_table_List
		migrationSize <- MigrationSize

		emigrantGroups<- EmigrantGroups
		immigrantGroups <- ImmigrantGroups
		direction <- Direction
		restImmigrantGroups <- c(emigrantGroups[-1],immigrantGroups)
		nIG <- length(restImmigrantGroups)

		for(i in 1:nIG){

			# outFamily <- sample(emigrantGroups,1)
			# inFamily <- sample(immigrantGroups,1)

			outFamily <- emigrantGroups[1]
			inFamily <- restImmigrantGroups[i]


			finalIndex <- 2

			# if(selected_PG_table_List[[inFamily]][finalIndex,1] >= selected_PG_table_List[[outFamily]][1,1]){

				# inFamilyTemp <-inFamily
				# inFamily <- outFamily
				# outFamily <- inFamilyTemp
			# }

		    if(selected_PG_table_List[[inFamily]][finalIndex,1] <= selected_PG_table_List[[outFamily]][1,1]){

				if(direction==2){

				   selected_PG_table_List[[outFamily]][finalIndex,1]<- selected_PG_table_List[[inFamily]][finalIndex,1]
				   selected_PG_table_List[[outFamily]][finalIndex,2]<- selected_PG_table_List[[inFamily]][finalIndex,2]
				}

				selected_PG_table_List[[inFamily]][finalIndex,1]<- selected_PG_table_List[[outFamily]][1,1]
				selected_PG_table_List[[inFamily]][finalIndex,2]<- selected_PG_table_List[[outFamily]][1,2]

			}else {

				selected_PG_table_List[[inFamily]][finalIndex,1]<- selected_PG_table_List[[inFamily]][finalIndex,1]
				selected_PG_table_List[[inFamily]][finalIndex,2]<- selected_PG_table_List[[inFamily]][finalIndex,2]
			}

		}

		  return(selected_PG_table_List)

}

####################



#####################################################################################################
## Four issues with the current settings
## No_of Progeny (bulk- 100/cross, discrete 10/cross- 10/family)
##  Single-round robin design
## (first crossed to the rest)
## (migrate samples ingroup and outgroup every generation)
