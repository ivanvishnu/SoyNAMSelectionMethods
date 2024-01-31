###############


getF5RILs_FS_IM_GM_Frontier <- 
function(GenoData_List,SelectedParents_Geno_List,SelectionOnGeno,nMarkers,NAM_LinkMap){
  
  BD <- "GM"
  NAM_LinkMap_New <- NAM_LinkMap
  
  GenoData_List <- GenoData_List
  
  selectedParents_Geno_List <- SelectedParents_Geno_List
  selectionOnGeno <- SelectionOnGeno
##### Generate F5 RIL Progeny from selected parent set for all families
  
#####Initialize lists for progeny data for all families every cylce

nParents <- dim(selectedParents_Geno_List)[1]
   
nUniqueParents <-  length(unique(as.vector(selectedParents_Geno_List)))
   
nProgeny <- 1
   
nPairs <- 1
  
  
### Loop througn selected parent pairs to generate progeny for each family
   
 if(nParents > 1){
  
     Cycle_Progeny_F5_temp1_List <- rep(list(array(0,c(1,nMarkers,2,nProgeny))),nParents)
 } 

 

while(nPairs <= (nParents)) {
 
Cycle_Progeny_F1 <- array(0,c(1,nMarkers,2,1))
   
Cycle_Progeny_F2 <- array(0,c(1,nMarkers,2,nProgeny))
   
Cycle_Progeny_F3 <- array(0,c(1,nMarkers,2,nProgeny))
   
Cycle_Progeny_F4 <- array(0,c(1,nMarkers,2,nProgeny))
   
Cycle_Progeny_F5_temp1 <- array(0,c(1,nMarkers,2,nProgeny))
   
  
   
  
###################################################################################

  
### Get Parent pair geno data
   
parent1Index <- selectedParents_Geno_List[nPairs,1]
   
parent2Index <- selectedParents_Geno_List[nPairs,2]
   
Parent1<-  GenoData_List[,,parent1Index]
   
Parent2<-  GenoData_List[,,parent2Index]
   
  
  
j<-1
   
###  F1progeny
   
  
  
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

### Get Cycle_Progeny F5 data in correct format
 
  Cycle_Progeny_F5_List_Temp <- Cycle_Progeny_F5_temp1_List

  Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Temp
   
  return(list(Cycle_Progeny_F5,nParents))
   
  
  
} 


###############

getF5RILs_FS_IM_GM_Frontier_Sel <- 
function(GenoData_List,SelectedParents_Geno_List,SelectionOnGeno,nMarkers,NAM_LinkMap){
  
  BD <- "GM"
  NAM_LinkMap_New <- NAM_LinkMap
  
  GenoData_List <- GenoData_List
  
  selectedParents_Geno_List <- SelectedParents_Geno_List
  selectionOnGeno <- SelectionOnGeno
  
##### Generate F5 RIL Progeny from selected parent set for all families
#####Initialize lists for progeny data for all families every cylce

nParents <- dim(selectedParents_Geno_List)[1]
nProgeny_Fam <- 100
   
nUniqueParents <-  length(unique(as.vector(selectedParents_Geno_List)))
  
nProgeny <- nProgeny_Fam/nParents
  
nPairs <- 1
  
  
### Loop througn selected parent pairs to generate progeny for each family
   
 if(nParents > 1){
  
     Cycle_Progeny_F5_temp1_List <- rep(list(array(0,c(1,nMarkers,2,nProgeny))),nParents)
 } 

 

while(nPairs <= (nParents)) {
 
Cycle_Progeny_F1 <- array(0,c(1,nMarkers,2,1))
   
Cycle_Progeny_F2 <- array(0,c(1,nMarkers,2,nProgeny))
   
Cycle_Progeny_F3 <- array(0,c(1,nMarkers,2,nProgeny))
   
Cycle_Progeny_F4 <- array(0,c(1,nMarkers,2,nProgeny))
   
Cycle_Progeny_F5_temp1 <- array(0,c(1,nMarkers,2,nProgeny))
   
    
###################################################################################
### Get Parent pair geno data
   
parent1Index <- selectedParents_Geno_List[nPairs,1]
   
parent2Index <- selectedParents_Geno_List[nPairs,2]
   
Parent1<-  GenoData_List[,,parent1Index]
   
Parent2<-  GenoData_List[,,parent2Index]
  
j<-1
   
###  F1progeny
    
  
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

### Get Cycle_Progeny F5 data in correct format
 
  Cycle_Progeny_F5_List_Temp <- Cycle_Progeny_F5_temp1_List

  Cycle_Progeny_F5 <- Cycle_Progeny_F5_List_Temp
   
  return(list(Cycle_Progeny_F5,nParents))
   
  
  
} 


###############


	
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

############
	
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


############


getF5RILs_GM_Fam_V2 <- function(Cycle_GenoData_List,selectedParents_Geno_Fam_List,selectedParents_NumProgeny_Geno_Fam_List,nFamilies,nMarkers,NAM_LinkMap,families,cycle1_nProgeny){
		   
		NAM_LinkMap_New <- NAM_LinkMap
		
##### Generate F5 RIL Progeny from selected parent set for all families 
##### Initialize lists for progeny data for all families every cylce 
			
 		GenoData_List <- Cycle_GenoData_List	
				
		selectedParents_Geno_Fam <- selectedParents_Geno_Fam_List
		selectedParents_NumProgeny_Geno_Fam_List <- selectedParents_NumProgeny_Geno_Fam_List
			
		cycle1_nProgeny_Fam <- cycle1_nProgeny
								
		nProgeny_list <- list()
		nSel_inFamily_list <- list()
		
		Cycle_Progeny_F5_List_Fam_Temp<- rep(list(list()),nFamilies)
		
		cycle1_nProgeny <- 100 
		cycle_nProgeny <- cycle1_nProgeny 
		selectedGenoData_List <- list()
	
                Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,cycle_nProgeny))),nFamilies)
	        Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,cycle_nProgeny))	
#######################################################################################################################

## Generate F5 progeny for individual families
		

	    for(nfamily in 1:nFamilies){
	
	      			
			familyGenoData <- GenoData_List[[nfamily]] 

					
			nProgenyVec <- (selectedParents_NumProgeny_Geno_Fam_List[[nfamily]])
			
		        nParents <- length(nProgenyVec)
			nSelIndividuals_inFamily <- nrow(familyGenoData)

			
			nProgeny_Factor <- cycle1_nProgeny_Fam/nSelIndividuals_inFamily
						
			nProgeny_list[[nfamily]] <- (selectedParents_NumProgeny_Geno_Fam_List[[nfamily]])

			nSel_inFamily_list[[nfamily]] <-  (length(selectedParents_NumProgeny_Geno_Fam_List[[nfamily]]))/2
				
			nSel_inFamily <- nSel_inFamily_list[[nfamily]]		
		
                        max_nProgenyVec <- max(as.vector(nProgenyVec))				
			Cycle_Progeny_F5_temp2_List <- rep(list(array(0,c(1,nMarkers,2,max_nProgenyVec))),nSel_inFamily) 
			
			nPairs <- 1
			
			genoIndividualIndices <- selectedParents_Geno_Fam[[nfamily]]
									
			while(nPairs < (nParents)) {
			
				nProgeny <- nProgeny_Factor * nProgenyVec[nPairs]
								 
				Cycle_Progeny_F1 <- array(0,c(1,nMarkers,2,1)) 
				
				Cycle_Progeny_F2 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F3 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F4 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F5_temp2 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				parent1Ind <- selectedParents_Geno_Fam[[nfamily]][nPairs]
				parent2Ind <- selectedParents_Geno_Fam[[nfamily]][nPairs+1]
			
				parent1Index <- which(genoIndividualIndices %in% parent1Ind)[1]
				parent2Index <- which(genoIndividualIndices %in% parent2Ind)[1]
				
				#Parent1 <-  familyGenoData[,,parent1Index]
			        #Parent2 <-  familyGenoData[,,parent2Index]

				Parent1 <-  familyGenoData[parent1Index,,]
			       Parent2 <-  familyGenoData[parent2Index,,]
###  F1	progeny

				j <- 1
		
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
				Cycle_Progeny_F5_temp2[j,,,m] <- progeny1 
			  }
			  
			    Cycle_Progeny_F5_temp2_List[[(nPairs+1)/2]] <- Cycle_Progeny_F5_temp2
				
				nPairs <- nPairs + 2
			}
			
			   Cycle_Progeny_F5_List_Fam_Temp[[nfamily]] <- Cycle_Progeny_F5_temp2_List
	    
	    }

	
##################################################################		
		
      	 
	    for(nFamily in 1:length(families)){ 
		
		    		
		    nSel_inFamily <- length(Cycle_Progeny_F5_List_Fam_Temp[[nFamily]])
                    nProgenyVec_List <- nProgeny_list[[nFamily]]
		 
		    nProgIndex <- 0
                    nCrs <- 1
                    nProg <-1
			
		    for(nSel in 1:nSel_inFamily){

						 
			 	nProgeny <-  nProgenyVec_List[nCrs]
				initProg <-  nProgIndex+1
			        nProgIndex <- nProgIndex+nProgeny
									 		   
				Cycle_Progeny_F5[nFamily,,,initProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam_Temp[[nFamily]][[nSel]][1,,,nProg:nProgeny] 
                                Cycle_Progeny_F5_List[[nFamily]][,,initProg:nProgIndex] <- Cycle_Progeny_F5_List_Fam_Temp[[nFamily]][[nSel]][1,,,nProg:nProgeny]

			        nCrs <- nCrs +2
		    } 
	       }

            			
		return(list(Cycle_Progeny_F5,Cycle_Progeny_F5_List)) 
  
    }




getF5RILs_BD_WGS <- function(BD,SelectedGenoData,nSelected,nProgeny,nMarkers,NAM_LinkMap){

	NAM_LinkMap_New <- NAM_LinkMap
    selectedGenoData <- SelectedGenoData

    if(BD== "BD1"){ 
        
	    if(nSelected==2){
		    nCrosses <- 1
		}else{
		    nCrosses <- nSelected
		}	   
		
		Cycle_Progeny_F1 <- array(0,c(nSelected,nMarkers,2,nProgeny))

		Cycle_Progeny_F2 <- array(0,c(nSelected,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nSelected,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nSelected,nMarkers,2,nProgeny))

		Cycle_Progeny_F5 <- array(0,c(nSelected,nMarkers,2,nProgeny))

###################################################################################
	   
		for( j in 1:(nCrosses)){

## Get parents from selected GenoData List

		  if(j ==1 && nCrosses ==1){
			  Parent1<-  selectedGenoData[j,,]
			  Parent2<-  selectedGenoData[j,,]
		  }else if ((nCrosses >1)&&(j<=(nCrosses -1))) {
		       Parent1<-  selectedGenoData[j,,]
			   Parent2<-  selectedGenoData[j+1,,]
		  }else if(j==nCrosses){
		       Parent1<-  selectedGenoData[j,,]
			   Parent2<-  selectedGenoData[1,,]
		  }
		  
## Get F1- F5 RILs	  
		  
		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		 
		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5[j,,,m] <- progeny1
		  }

		}
    
} else if (BD =="BD2"){ 

    if(nSelected ==2){
			nCrosses <- 1
	}else{nCrosses <- nSelected}


		Cycle_Progeny_F1 <- array(0,c(nCrosses,nMarkers,2,1))

		Cycle_Progeny_F2 <- array(0,c(nCrosses,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nCrosses,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nCrosses,nMarkers,2,nProgeny))

		Cycle_Progeny_F5 <- array(0,c(nCrosses,nMarkers,2,nProgeny))

###################################################################################

		for( j in 1:(nCrosses)){

## Get parents from selectedGenoData_List 
		
		  if(j ==1 && nCrosses ==1){
			  Parent1<-  selectedGenoData[j,,]
			  Parent2<-  selectedGenoData[j+1,,]
		  }else if ((nCrosses>1)&&(j<=(nCrosses%/%2))) {
		       Parent1<-  selectedGenoData[1,,]
			   Parent2<-  selectedGenoData[j+1,,]
		  }else if ((nCrosses >1)&&(j>(nCrosses %/%2))&&(j<=(nCrosses-1))) {
		       Parent1<-  selectedGenoData[2,,]
			   Parent2<-  selectedGenoData[(j-(nCrosses %/%2)+1),,]
		  }else if(j==nCrosses){
		       Parent1<-  selectedGenoData[j,,]
			   Parent2<-  selectedGenoData[1,,]
		  }
	
## Get F1- F5 RILs
	
		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1

		 
		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5[j,,,m] <- progeny1
		  }

		}

		
} else if(BD == "RM"){ 

## RM 
        if(nSelected ==2){
			nCrosses <- 1
	    }else{nCrosses <- nSelected}

       
	    selectedGenoData <- selectedGenoData

		Cycle_Progeny_F1 <- array(0,c(nCrosses,nMarkers,2,1))

		Cycle_Progeny_F2 <- array(0,c(nCrosses,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nCrosses,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nCrosses,nMarkers,2,nProgeny))

		Cycle_Progeny_F5  <- array(0,c(nCrosses,nMarkers,2,nProgeny))

###################################################################################

			    
		for( j in 1:(nCrosses)){
		
		  parentIndices <- c(1:nCrosses)
	      parent1Index <-  sample(parentIndices,1)
		  parent2Index <- sample(parentIndices[-(parent1Index)],1)

## Get parents from selected GenoData List

		  if(j ==1 && nCrosses ==1){
			  Parent1 <-  selectedGenoData[parent1Index,,]
			  Parent2 <-  selectedGenoData[parent1Index,,]
		  }else if ((nCrosses>1)&&(j<=(nCrosses))){
		       Parent1 <-  selectedGenoData[parent1Index,,]
			   Parent2 <-  selectedGenoData[parent2Index,,]
		  }


## Get F1- F5 RILs	  
		  
		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5[j,,,m] <- progeny1
		  }

		}

     }	     


    return(Cycle_Progeny_F5)


}

#############################################################




### Fn 39: 
##################################################################################################################################

getF5RILs_BD <- function(BD,selectedGenoData_List,nFamilies,nSel_inFamily,nProgeny,nMarkers,NAM_LinkMap,Rep){

   NAM_LinkMap_New <- NAM_LinkMap

   selectedGenoData_List_afterExchange <- selectedGenoData_List
   
   initSeed <- 25+Rep
   set.seed(initSeed)

  if(BD == "BD1"){ 
        
	  selectedGenoData_List_afterExchange <- selectedGenoData_List

	  Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nSel_inFamily,nMarkers,2,nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nSel_inFamily))

	 for(nFamily in 1:nFamilies) {

		Cycle_Progeny_F1 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F2 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

###################################################################################

		if(nSel_inFamily==2){
				nCrosses_inFamily <- 1
		}else{
				nCrosses_inFamily <- nSel_inFamily
		}
			
	    # nCrosses <- nSel_inFamily

		for( j in 1:(nCrosses_inFamily)){

## Get parents from selected GenoData List

		  if(j ==1 && nCrosses_inFamily ==1){
			  Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=(nCrosses_inFamily-1))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if(j==nCrosses_inFamily){
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
		  }
		  
		  
	
		  

## Get F1- F5 RILs	  
		  
		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  # for(m in 1:nProgeny){

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){
		
			indCount <- 1

			for(nSel in 1:nSel_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,indCount] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
					
					indCount <- indCount+1
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){
		
			indCount <-1	

		    for(nSelFamily in 1:nSel_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,indCount] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]

					indCount <- indCount+1

				}
			}
		}

}
  if(BD == "BD2"){ 

    if(nSel_inFamily==2){
			nCrosses_inFamily <- 1
	}else{nCrosses_inFamily <- nSel_inFamily}


	Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))),nFamilies)

	Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nCrosses_inFamily*nProgeny))),nFamilies)

	Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nCrosses_inFamily))

	for(nFamily in 1:nFamilies) {

		Cycle_Progeny_F1 <- array(0,c(nCrosses_inFamily,nMarkers,2,1))

		Cycle_Progeny_F2 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

###################################################################################

		for( j in 1:(nCrosses_inFamily)){

## Get parents from selectedGenoData_List 
		
		if(j ==1 && nCrosses_inFamily ==1){
			  Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=(nCrosses_inFamily%/%2))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j>(nCrosses_inFamily%/%2))&&(j<=(nCrosses_inFamily-1))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][2,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][(j-(nCrosses_inFamily%/%2)+1),,]
		  }else if(j==nCrosses_inFamily){
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
		  }

		  
	
## Get F1- F5 RILs
	
		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  # for(m in 1:nProgeny){

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1

		  # }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){
		
		    indCount <-1 

			for(nSel in 1:nCrosses_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,indCount] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
					
					indCount <- indCount +1
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){
		
			indCount <-1 

		    for(nSelFamily in 1:nCrosses_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,indCount] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]

					indCount <- indCount +1 
				}
			}
		}
    
		
		
} 
  if(BD == "RM_SP"){

### Generate F5 RIL progeny for selected geno data
    if(nSel_inFamily==2){
			nCrosses_inFamily <- 1
	}else{nCrosses_inFamily <- nSel_inFamily}


	Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))),nFamilies)

	Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nCrosses_inFamily*nProgeny))),nFamilies)

	Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nCrosses_inFamily))

	for(nFamily in 1:nFamilies) {



		Cycle_Progeny_F1 <- array(0,c(nCrosses_inFamily,nMarkers,2,1))

		Cycle_Progeny_F2 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nCrosses_inFamily,nMarkers,2,nProgeny))


###################################################################################


		for( j in 1:(nCrosses_inFamily)){

		  if(j ==1 && nCrosses_inFamily ==1){
			  Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=(nCrosses_inFamily%/%2))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j>(nCrosses_inFamily%/%2))&&(j<=(nCrosses_inFamily-1))) {
		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][2,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][(j-(nCrosses_inFamily%/%2)+1),,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List_afterExchange[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List_afterExchange[[nFamily]][1,,]
		  }

		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  # for(m in 1:nProgeny){

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1

		  # }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){
		
		    indCount <-1

			for(nSel in 1:nCrosses_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,indCount] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
					
					indCount <- indCount+1
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nCrosses_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}

# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

} 
  if(BD == "RM"){
## RM
    
	  selectedGenoData_List_afterExchange <- selectedGenoData_List

	  Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nSel_inFamily,nMarkers,2,nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nSel_inFamily))

	 for(nFamily in 1:nFamilies){

		Cycle_Progeny_F1 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F2 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

###################################################################################

		if(nSel_inFamily==2){
				nCrosses_inFamily <- 1
		}else{
				nCrosses_inFamily <- nSel_inFamily
		}
			
	    # nCrosses <- nSel_inFamily

		for( j in 1:(nCrosses_inFamily)){
		
		 parentIndices <- c(1:nCrosses_inFamily)
	     parent1Index <-  sample(parentIndices,1)
		 parent2Index <- sample(parentIndices[-(parent1Index)],1)

## Get parents from selected GenoData List

		  if(j ==1 && nCrosses_inFamily ==1){
			  Parent1 <-  selectedGenoData_List_afterExchange[[nFamily]][parent1Index,,]
			  Parent2 <-  selectedGenoData_List_afterExchange[[nFamily]][parent1Index,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=(nCrosses_inFamily-1))){
		       Parent1 <-  selectedGenoData_List_afterExchange[[nFamily]][parent1Index,,]
			   Parent2 <-  selectedGenoData_List_afterExchange[[nFamily]][parent2Index,,]
		  }


## Get F1- F5 RILs	  
		  
		  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

		  # for(m in 1:nProgeny){

		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,]  <- progeny1


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F2[j,,,m]
			Parent2<- Cycle_Progeny_F2[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F3[j,,,m] <- progeny1
		  }


		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F3[j,,,m]
			Parent2<- Cycle_Progeny_F3[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F4[j,,,m] <- progeny1
		  }

		  for(m in 1:nProgeny){

			Parent1<- Cycle_Progeny_F4[j,,,m]
			Parent2<- Cycle_Progeny_F4[j,,,m]
			progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
			Cycle_Progeny_F5_temp[j,,,m] <- progeny1
		  }

		}

	     Cycle_Progeny_F5_List_temp[[nFamily]]  <- Cycle_Progeny_F5_temp

    }

### Get Cycle_Progeny F5 data in correct format

	    for(nFamily in 1:nFamilies){
		
			indCount <- 1

			for(nSel in 1:nSel_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,indCount] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
					
					indCount <- indCount+1
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format

		for(nFamily in 1:nFamilies){
		
			indCount <-1	

		    for(nSelFamily in 1:nSel_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,indCount] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]

					indCount <- indCount+1

				}
			}
		}

}
	
	return(list(Cycle_Progeny_F5,Cycle_Progeny_F5_List))


}
############################################################# 
########################################################################################################################################################



getF5RILs_GM <- function(BD,selectedGenoData_List,selectedParents_NumProgeny_Geno_List,nFamilies,nMarkers,NAM_LinkMap,familyPairs,families){

	    BD <- "GM"
		NAM_LinkMap_New <- NAM_LinkMap
		
##### Generate F5 RIL Progeny from selected parent set for all families 
#####	Initialize lists for progeny data for all families every cylce 
	
		Cycle_Progeny_F5_List_temp<- rep(list(list()),nFamilies)
		
		Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,100))),nFamilies) 
		
		Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,100))
		
		nProgeny_list <- list()
		nSel_inFamily_list <- list()
		
				
		for (nFamily in 1:nFamilies){ 
		
				
            if(selectionOnGeno==TRUE){

				nProgeny_list[[nFamily]] <- (selectedParents_NumProgeny_Geno_List[[nFamily]])

				nSel_inFamily_list[[nFamily]] <-  (length(selectedParents_NumProgeny_Geno_List[[nFamily]]))/2
				
			}else if(selectionOnGeno==FALSE){

				nProgeny_list[[nFamily]] <- (selectedParents_NumProgeny_Pheno_List[[nFamily]])

				nSel_inFamily_list[[nFamily]] <-  (length(selectedParents_NumProgeny_Pheno_List[[nFamily]]))/2
			}	

		}
			
	
		for(nFamily in 1:nFamilies){
		
					
			nSel_inFamily <- nSel_inFamily_list[[nFamily]]
		 
		    nParents <- nSel_inFamily*2
		 
			nProgenyVec <- nProgeny_list[[nFamily]]
			nProgenyVecPerPair <- nProgenyVec[seq(1,length(nProgenyVec),by=2)]
		
            Cycle_Progeny_F5_temp1_List <- rep(list(array(0,c(1,nMarkers,2,nProgeny))),nSel_inFamily) 
				  
			
			nPairs <- 1
			
	### Loop througn selected parent pairs to generate progeny for each family
			
			while(nPairs < (nParents)) {
			
				nProgeny <- nProgenyVec[nPairs]

									 
				Cycle_Progeny_F1 <- array(0,c(1,nMarkers,2,1)) 
				
				Cycle_Progeny_F2 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F3 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F4 <- array(0,c(1,nMarkers,2,nProgeny)) 
				
				Cycle_Progeny_F5_temp1 <- array(0,c(1,nMarkers,2,nProgeny)) 

			
	###################################################################################
		 
		### Get Parent pair geno data
				Parent1<-  selectedGenoData_List[[nFamily]][nPairs,,]
				Parent2<-  selectedGenoData_List[[nFamily]][nPairs+1,,]
			  
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
			  
			   Cycle_Progeny_F5_temp1_List[[(nPairs+1)/2]]<- Cycle_Progeny_F5_temp1
			
				nPairs <- nPairs+2
			 
			}
			
			   Cycle_Progeny_F5_List_temp[[nFamily]] <- Cycle_Progeny_F5_temp1_List

		}	
 
	   
### Get Cycle_Progeny F5 data in correct format 
   
	   
       	 
	    for(nFamily in 1:nFamilies){ 
		
		    nSel_inFamily <- length(Cycle_Progeny_F5_List_temp[[nFamily]])
		 
		    nProgIndex <- 1
			
			for(nSel in 1:nSel_inFamily){
			 
			 	nProgeny <-  dim(Cycle_Progeny_F5_List_temp[[nFamily]][[nSel]])[4]
			    
				nProg <-1
				#nProgIndex <- nProgIndex
				while(nProg <= nProgeny){
				 		   
				    Cycle_Progeny_F5[nFamily,,,nProgIndex] <- Cycle_Progeny_F5_List_temp[[nFamily]][[nSel]][1,,,nProg] 

					nProg <- nProg+1
					nProgIndex <- nProgIndex+1
				}
			} 
		}


			

        for(nFamily in 1:nFamilies){
								  

                nSel_inFamily <- length(Cycle_Progeny_F5_List_temp[[nFamily]])

                nProgIndex <- 1

                for(nSelFamily in 1:nSel_inFamily){

                        nProgeny <-  dim(Cycle_Progeny_F5_List_temp[[nFamily]][[nSelFamily]])[4]

                        nProg <-1

                            while(nProg <= nProgeny){

                                        Cycle_Progeny_F5_List[[nFamily]][,,nProgIndex] <- Cycle_Progeny_F5_List_temp[[nFamily]][[nSelFamily]][1,,,nProg]

                                        nProg <- nProg+1
                                        nProgIndex <-  nProgIndex+1

                            }
                    }
        }




	return(list(Cycle_Progeny_F5,Cycle_Progeny_F5_List)) 
  
  }


###########################################################################################################


