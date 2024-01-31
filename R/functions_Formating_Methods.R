### Functions to generate table and format
###
   	
  getNextGenGenoTableData <- function(Cycle_Progeny_F5_List_Fam,nCrosses_List,NProgeny_in_Cross,no_QTL,varE){

	    Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Fam
		nCrosses <-  nCrosses_List
		nProgeny_in_Cross <- NProgeny_in_Cross 
		
		
		nextGenGenoTable_GM <- generateMinus101GenoFormat_IM_GM_Frontier_V2(Cycle_Progeny_F5,nCrosses,nProgeny_in_Cross)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)

        return(list(nextGenGenoTable_GM,nextGenGenoTableMod_GM,newNextGenGenoTable_GM,genoValSimValues,phenoValSimValues))
  }
	
	
	#nextGenGenoTable_GM <- generateMinus101GenoFormat_IM_GM_Frontier_V2(Cycle_Progeny_F5,nCrosses,nProgeny_in_Cross)

################ 


  getNextGenGenoTableData_Sel_V2 <- function(Cycle_Progeny_F5_List_Fam,nCrosses_List,nProgeny_per_Cross,noQTL,VarE){

	    Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Fam
		nProgeny_in_Fam <- dim(Cycle_Progeny_F5)[3]
		nCrosses <-  nCrosses_List
		nProgeny <- nProgeny_in_Fam
		nProgeny_in_Cross <- nProgeny_per_Cross
		no_QTL <- noQTL 
		varE <- VarE
					
		nextGenGenoTable_GM <- generateMinus101GenoFormat_IM_GM_Frontier_V2(Cycle_Progeny_F5,nCrosses,nProgeny_in_Cross)
		nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
		newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
	
		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE,no_QTL)

        return(list(nextGenGenoTable_GM,nextGenGenoTableMod_GM,newNextGenGenoTable_GM,genoValSimValues,phenoValSimValues))
	}	
	
################ 


getNextGenGenoTableData_Complete_V2 <- 
function (Cycle_Progeny_F5_List_Fam,NProgeny_in_Family,no_QTL, varE) 
{
    Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Fam
    nProgeny_in_Family <- NProgeny_in_Family
    nextGenGenoTable_GM <- generateMinus101GenoFormat_IM_GM_Frontier_V3_Complete(Cycle_Progeny_F5,nProgeny_in_Family)
    nextGenGenoTableMod_GM <- generate012GenoFormat_Inverse(nextGenGenoTable_GM)
    newNextGenGenoTable_GM <- generateWeightedGenoTable(nextGenGenoTable_GM,no_QTL)
    genoValSimValues <- simulateGenotypicValues(nextGenGenoTable_GM,no_QTL)
    phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable_GM,varE, no_QTL)
    return(list(nextGenGenoTable_GM, nextGenGenoTableMod_GM,newNextGenGenoTable_GM, genoValSimValues, phenoValSimValues))
}



####### 

generateMinus101GenoFormat_IM_GM_Frontier_V3_Complete <- 
function (Cycle_Progeny_Data,noProgeny_in_Family) 
{
    Cycle_Progeny <- Cycle_Progeny_Data
    nMarkers <- 4289
    nProgeny <- noProgeny_in_Family
    Cycle_Progeny_table <- c()
    
    cycle_family_table <- array(0, c(1,nMarkers,nProgeny))
    for (k in 1:nProgeny) {
        cycle_family_table[1,, k] <- (Cycle_Progeny[, 1, k]+Cycle_Progeny[, 2, k]) - 1
    
     }
     Cycle_Progeny_table <- cbind(Cycle_Progeny_table, cycle_family_table[1,, ])
   
    genoTable <- t(Cycle_Progeny_table)
    genoTable_Mod <- apply(genoTable, 2, as.numeric)
    return(genoTable_Mod)
}


	
### Generate Minus101 format genotype table for frontier progeny set

# generateMinus101GenoFormat_IM_GM_Frontier_V2 <- function(Cycle_Progeny_Data,nFamilies,nSel_in_Family){

  # Cycle_Progeny <- Cycle_Progeny_Data
  # #nCrosses <- number_of_Crosses
  
  # nMarkers <- 4289
  # nProgeny <- nSel_in_Family
 
  # Cycle_Progeny_table <-c()
  
  # for(i in 1:nFamilies){
   # cycle_family_table<- array(0,c(1,nMarkers,nProgeny))
  
    # for(k in 1:nProgeny){
   
	 # cycle_family_table[1,,k] <- (Cycle_Progeny[,1,k] + Cycle_Progeny[,2,k]) -1

	# }
   
    # Cycle_Progeny_table <-cbind(Cycle_Progeny_table,cycle_family_table[1,,])

  # }
  
 # IndividualNames <-rep(0,(nFamilies*nProgeny))

  # for(i in 1:(nFamilies*nProgeny)){
    # IndividualNames[i]<- paste("Ind",i,sep="")
  # }

  # # IndividualNames <- c("Animal_ID",IndividualNames)
  # #### RowNames

  # markerNames <-rep(0,nMarkers)


  # for(i in 1:(nMarkers)){
    # markerNames[i]<- paste("m",i,sep="")
  # }
  
  # colnames(Cycle_Progeny_table)<- IndividualNames
  # row.names(Cycle_Progeny_table)<- markerNames

  # genoTable <- t(Cycle_Progeny_table)
  
  # genoTable_Mod <- apply(genoTable,2,as.numeric)

  # ########## Translate 0-1-2 code to -1-0-1 code in genotable

  # #genoTable_Mod <- apply(genoTable,2,function(x)(gsub("\\b0\\b","minus1",as.matrix(x))))
  
  # # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b1\\b","0",as.matrix(x))))
  # # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b2\\b","1",as.matrix(x))))
  # # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\bminus1\\b","-1",as.matrix(x))))

  # return(genoTable_Mod)
# }



generateMinus101GenoFormat_IM_GM_Frontier_V2 <- function(Cycle_Progeny_Data,nFamilies,nSel_in_Family){   
   
    Cycle_Progeny <- Cycle_Progeny_Data
      
    nMarkers <- 4289
    nProgeny <- nSel_in_Family
	
              
    Cycle_Progeny_table <- (Cycle_Progeny[,1,] + Cycle_Progeny[,2,])- 1
       
    nLines <- dim(Cycle_Progeny_table)[2]
    IndividualNames <- paste("Ind", c(1:nLines), sep = "")
    markerNames <- paste("m",c(1:nMarkers), sep = "")
    
    colnames(Cycle_Progeny_table) <- IndividualNames
    row.names(Cycle_Progeny_table) <- markerNames
    genoTable <- t(Cycle_Progeny_table)
    genoTable_Mod <- apply(genoTable, 2, as.numeric)
    return(genoTable_Mod)
} 





### Generate Minus101 format genotype table for frontier progeny set

generateMinus101GenoFormat_IM_GM_Frontier_V1<- function(Cycle_Progeny_Data,number_of_Crosses,noProgeny_in_Cross){

  Cycle_Progeny <- Cycle_Progeny_Data
  nCrosses <- number_of_Crosses
  
  nMarkers <- 4289
  nProgeny <- noProgeny_in_Cross
 
  Cycle_Progeny_table <-c()
  
  for(i in 1:nCrosses){
   cycle_family_table<- array(0,c(1,nMarkers,nProgeny))
  
    for(k in 1:nProgeny){
   
	 cycle_family_table[1,,k] <- (Cycle_Progeny[[i]][1,,1,k] + Cycle_Progeny[[i]][1,,2,k]) -1

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

  genoTable <- t(Cycle_Progeny_table)
  
  genoTable_Mod <- apply(genoTable,2,as.numeric)

  ########## Translate 0-1-2 code to -1-0-1 code in genotable

  #genoTable_Mod <- apply(genoTable,2,function(x)(gsub("\\b0\\b","minus1",as.matrix(x))))
  
  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b1\\b","0",as.matrix(x))))
  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b2\\b","1",as.matrix(x))))
  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\bminus1\\b","-1",as.matrix(x))))

  return(genoTable_Mod)
}


 