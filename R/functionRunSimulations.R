## Module 4: Functions to run the simulation
## Apply selection to choose parents for the next generation ##

## Fn 1:
### Extract Genotype data in long array format ########################################################
### Meant for genotype data with 20 crosses and 100 progeny with 2 columns of markers and 4000 loci

extractGenoData <- function(GenoTableIndices1,ProgenyGenoData,no_selected){

  GenoTableIndices <- GenoTableIndices1
  nFamilies <- dim(ProgenyGenoData)[1]
  nMarkers <-4289
  if(nFamilies == 20){

	n<- length(GenoTableIndices)

	genoData <- array(0,c(n,nMarkers,2))

	for(j in 1:n){

		if(GenoTableIndices[j]<100){
			i<-1
		}

		if(GenoTableIndices[j]>=100 && GenoTableIndices[j]!=5000 ){
			i <- (GenoTableIndices[j] %/% 100)+1
		}

		if(GenoTableIndices[j]>=100 && GenoTableIndices[j]!=2000 ){
			i <- (GenoTableIndices[j] %/% 100)+1
		}


		else if(GenoTableIndices[j]==5000 ){
			i<-50
			l<-100
		}

		else if(GenoTableIndices[j]==2000 ){
			i<-20
			l<-100
		}

		l <- GenoTableIndices[j] %% 100

		if(l==0){

			l=100
		}

		genoData[j,,] <- ProgenyGenoData[i,,,l]
    }
  }

  if(nFamilies ==50){

	n<- length(GenoTableIndices)

	genoData <- array(0,c(n,nMarkers,2))

    for(j in 1:n){

		if(GenoTableIndices[j]<40){
			i<-1
		}

		if(GenoTableIndices[j]>=40 && GenoTableIndices[j]!= 5000 ){
			i <- (GenoTableIndices[j] %/% 40)+1
		}

		if(GenoTableIndices[j]>=40 && GenoTableIndices[j]!= 2000 ){
			i <- (GenoTableIndices[j] %/% 40)+1
		}


		# else if(GenoTableIndices[j]==5000 ){
			# i<-50
			# l<-40
		# }

		else if(GenoTableIndices[j]==2000 ){
			i<-50
			l<-40
		}

		l <- GenoTableIndices[j] %% 40

		if(l==0){

			l=40
		}

		genoData[j,,] <- ProgenyGenoData[i,,,l]

	}
  }

  if(nFamilies == 200){

	n<- length(GenoTableIndices)

	genoData <- array(0,c(n,nMarkers,2))

    for(j in 1:n){

		if(GenoTableIndices[j]< 10){
			i<-1
		}

		if(GenoTableIndices[j]>=10 && GenoTableIndices[j]!= 5000 ){
			i <- (GenoTableIndices[j] %/% 10)+1
		}

		if(GenoTableIndices[j]>=10 && GenoTableIndices[j]!= 2000 ){
			i <- (GenoTableIndices[j] %/% 10)+1
		}


		# else if(GenoTableIndices[j]==5000 ){
			# i<-50
			# l<-40
		# }

		else if(GenoTableIndices[j]==2000 ){
			i<-200
			l<-10
		}

		l <- GenoTableIndices[j] %% 10

		if(l==0){

			l=10
		}

		genoData[j,,] <- ProgenyGenoData[i,,,l]

	}
  }

 	return(genoData)
}

## Fn 2:
## Function Test #######################################################################
# gen1_top20_genoData <- extractGenoData(GenoTableIndices_top20,F5_Progeny)
# gen1_top50_genoData <- extractGenoData(selectedGenoIndividualIndices,Cycle_Progeny)
########### Generate Doubled Haploid ####################################################

generateDoubledHaploid <- function(selectedGenoData){

  no_of_lines <- dim(selectedGenoData)[1]
  nMarkers <- dim(selectedGenoData)[2]
  nChr <- dim(selectedGenoData)[3]

  DH_Lines<- array(0,c(2*no_of_lines,nMarkers,nChr))

  j <-1

  for( i in 1: no_of_lines){

    DH_Lines[j,,] <- cbind(selectedGenoData[i,,1],selectedGenoData[i,,1])
    DH_Lines[j+1,,]<-  cbind(selectedGenoData[i,,2],selectedGenoData[i,,2])
    j<- j+2

  }

  return(DH_Lines)
}


### test generateDoubledHaploid
# DH_test <- generateDoubledHaploid(gen1_top20_genoData)

## Fn 2a:
########### Generate Doubled Haploid for (2a) ####################################################################

generateDoubledHaploid_2a <- function(selectedGenoData){

  no_of_lines <- dim(selectedGenoData)[1]
  nMarkers <- dim(selectedGenoData)[2]
  nChr <- dim(selectedGenoData)[3]
  nProgeny <- dim(selectedGenoData)[4]

  DH_Lines<- array(0,c(no_of_lines,nMarkers,nChr,nProgeny,2))
  DH_Lines_New <- array(0,c(no_of_lines,nMarkers,nChr,nProgeny))


  for( i in 1: no_of_lines){
    for( k in 1: nProgeny){

      index <-sample(c(1,2),1)

      DH_Lines_New[i,,,k] <- cbind(selectedGenoData[i,,index,k],selectedGenoData[i,,index,k])
      #DH_Lines[j,,,k,2] <-  cbind(selectedGenoData[i,,2,k],selectedGenoData[i,,2,k])


      }
  }

     return(DH_Lines_New)
  }


##  test generateDoubledHaploid_2a

#DH_test <- generateDoubledHaploid_2a(Cycle_Progeny)


#################################################################################################################
## Fn 3:
### Function to generate progeny given two parents and the number of Progeny
### generateProgeny(Parent1,Parent2,Number_of_progeny)

# generateProgeny <- function(Parent1,Parent2,Number_Progeny,NAM_LinkMap){

  # NAM_LinkMap_New <- NAM_LinkMap
  # RF_Vector <- as.vector(NAM_LinkMap_New[,3])
  # rr <- as.vector(NAM_LinkMap_New[,3])
  # number_progeny <- Number_Progeny

  # x1<-Parent1
  # x2<- Parent2


  # a <- Number_Progeny
  # n <- dim(x1)[1]

  # y<-array(0,dim=c(n,2,a))


  # for (i in 1:a) {

    # p1<-cumsum(runif(n)<=rr)%%2
    # p2<-cumsum(runif(n)<=rr)%%2 #creating gamete from p2


    # y[p1==0,1,i]=x1[p1==0,1]

    # y[p1==1,1,i] = x1[p1==1,2]

    # y[p2==0,2,i] = x2[p2==0,1]

    # y[p2==1,2,i] = x2[p2==1,2]
  # }

  # progeny <- y

  # return(progeny)


# }

## Fn 4: Function to generate genotype table in 012 format
######### Function to generate GenoTable in 012 format ########################

generate012GenoFormat <- function(Cycle_Progeny_Data,number_of_Crosses,n_Progeny){

  Cycle_Progeny <- Cycle_Progeny_Data


  nCrosses <- number_of_Crosses
  nProgeny<- n_Progeny
  nMarkers <- 4289

  Cycle_Progeny_table <-c()
  cycle_family_table<- array(0,c(1,nMarkers,nProgeny))

  for(i in 1:nCrosses){

    for(k in 1:nProgeny){


      cycle_family_table[1,,k] <- (Cycle_Progeny[i,,1,k] + Cycle_Progeny[i,,2,k])

    }

    Cycle_Progeny_table <-cbind(Cycle_Progeny_table,cycle_family_table[1,,1:nProgeny])

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


  return(genoTable)
}

## Test Function generate012GenoFormat
#Cycle_test <- generate012GenoFormat(Cycle_Progeny,50)

#############################################################################################################
## Fn 5:

generateMinus101GenoFormat <- function(Cycle_Progeny_Data,number_of_Crosses,n_Progeny,NFamilies,no_Selected){


  Cycle_Progeny <- Cycle_Progeny_Data

  nFamilies <- NFamilies
  nProgeny <- n_Progeny
  nMarkers <- 4289
  nSelected <- no_Selected

   # if(nSelected ==20){
	# cycle_family_table<- array(0,c(1,nMarkers,nProgeny))
   # }
  # cycle_family_table<- array(0,c(nCrosses,nMarkers,nProgeny))

   Cycle_Progeny_table <- c()
   cycle_family_table <- array(0,c(nFamilies,nMarkers,nProgeny))



  for(m in 1:nFamilies){

		for(k in 1:nProgeny){


			cycle_family_table[m,,k] <- (Cycle_Progeny[m,,1,k] + Cycle_Progeny[m,,2,k])

		}

        Cycle_Progeny_table <- cbind(Cycle_Progeny_table,cycle_family_table[m,,])
  }


  IndividualNames <-rep(0,(nFamilies*nProgeny))

  for(i in 1:(nFamilies*nProgeny)){
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

## Fn 6:
## Test Function generate012GenoFormat
## Cycle_test <- generateMinus101GenoFormat(F5_Progeny,20)

  generate012GenoFormat_Inverse<- function(genoTable_minus101){

      genoTable_Mod <- genoTable_minus101
	  genoTable_Mod <- apply(genoTable_Mod,2,function(x)(x+1))
	  genoTable_Mod <- apply(genoTable_Mod,2,as.numeric)
	

	  return(genoTable_Mod)

   }

	 

	  # genoTable_Mod <- Cycle_test
	  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b-1\\b","minus1",as.matrix(x))))
	  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b1\\b","2",as.matrix(x))))
	  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b0\\b","1",as.matrix(x))))
	  # # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b-1\\b","0",as.matrix(x))))
	  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\bminus1\\b","0",as.matrix(x))))
	  # genoTable_Mod <- apply(genoTable_Mod,2,as.numeric)

	  
	  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b1\\b","2",as.matrix(x))))
	  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b0\\b","1",as.matrix(x))))
	  # # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b-1\\b","0",as.matrix(x))))
	  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\bminus1\\b","0",as.matrix(x))))
	 



## Test Function generate012GenoFormat
# Cycle_test <- generateMinus101GenoFormat_V1(F5_Progeny,20,100)
# Cycle_012Test <- generate012GenoFormat_Inverse(Cycle_test)

## Fn 7
##################################################################################################

generateMinus101GenoFormat_V1 <- function(Cycle_Progeny_Data,number_of_Crosses,n_Progeny){

  Cycle_Progeny <- Cycle_Progeny_Data

  nCrosses <- number_of_Crosses
  nProgeny <- n_Progeny
  nMarkers <- 4289
  # nSelected <- noSelected
  Cycle_Progeny_table <-c()
  cycle_family_table<- array(0,c(nCrosses,nMarkers,nProgeny))


  for(i in 1:nCrosses){

    for(k in 1:nProgeny){


		cycle_family_table[i,,k] <- (Cycle_Progeny[i,,1,k] + Cycle_Progeny[i,,2,k]) -1 

	 }

    Cycle_Progeny_table <-cbind(Cycle_Progeny_table,cycle_family_table[i,,])

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

  ########## Translate 0-1-2 code to -1-0-1 code in genotable

  # genoTable_Mod <- apply(genoTable,2,function(x)(gsub("\\b0\\b","minus1",as.matrix(x))))
  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b1\\b","0",as.matrix(x))))
  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b2\\b","1",as.matrix(x))))
  # genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\bminus1\\b","-1",as.matrix(x))))
 
  genoTable_Mod <- apply(genoTable,2,as.numeric)
  return(genoTable_Mod)
}


## Test Function generate012GenoFormat
# Cycle_test <- generateMinus101GenoFormat(F5_Progeny,20)

## Fn 8
##################################################################################################

generateMinus101GenoFormat_V2 <- function(Cycle_Progeny_Data,number_of_Crosses,n_Progeny){

  Cycle_Progeny <- Cycle_Progeny_Data

  nCrosses <- 1
  nProgeny <- n_Progeny
  nMarkers <- 4289

  cycle_family_table<- array(0,c(nMarkers,nProgeny))


  # for(i in 1:nCrosses){

    for(k in 1:nProgeny){


		cycle_family_table[,k] <- (Cycle_Progeny[,1,k] + Cycle_Progeny[,2,k])

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


  colnames(cycle_family_table)<- IndividualNames
  row.names(cycle_family_table)<- markerNames

  genoTable<- t(cycle_family_table)

  ########## Translate 0-1-2 code to -1-0-1 code in genotable

  genoTable_Mod <- apply(genoTable,2,function(x)(gsub("\\b0\\b","minus1",as.matrix(x))))
  genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b1\\b","0",as.matrix(x))))
  genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\b2\\b","1",as.matrix(x))))
  genoTable_Mod <- apply(genoTable_Mod,2,function(x)(gsub("\\bminus1\\b","-1",as.matrix(x))))
  genoTable_Mod <- apply(genoTable_Mod,2,as.numeric)


  return(genoTable_Mod)
}

## Test Function generate012GenoFormat

######################################################################################################
## Fn 13
#### Function to simulate genotypic values give Genotype table from each cycle of crossing
######################################################################################################

simulateGenotypicValues <- function(nextGenGenoTable,noLoci) {

    genoTable <- nextGenGenoTable
    nIndividuals <-nrow(nextGenGenoTable)
    genotypicValuesSim <-rep(0,nIndividuals)
    nLoci <- noLoci
	nMarkers <- ncol(genoTable)
    ##### Extract QTL table ######

    if(nLoci== 40){
      QTL_indices<- seq(50,(nMarkers/2),by=53)
    }


    if(nLoci == 400){
      QTL_indices_1<- seq(8,nMarkers,by=12)
      QTL_indices_2 <- seq(50,nMarkers,by=107)
      n_3 <- 400-(length(QTL_indices_1))-(length(QTL_indices_2))
      QTL_indices_3<- sample(c(1:nMarkers),n_3)

      QTL_indices <- c(QTL_indices_1,QTL_indices_2,QTL_indices_3)
    }


    if(nLoci == nMarkers){
      QTL_indices<- c(1:nMarkers)
    }


    QTL_table <- genoTable[,c(QTL_indices)]

############### Simulate Genotypic Values ######################

    if(nLoci == 40){

		evenIndices<- seq(2,40,by=2)
		oddIndices <- seq(1,39,by=2)

		for(i in 1:nIndividuals){

			genotypicValuesSim[i] <- sum(QTL_table[i,evenIndices]*5)+sum(QTL_table[i,oddIndices]*-5)

		}
}

	if(nLoci == 400){
		evenIndices<- seq(2,400,by=2)
		oddIndices <- seq(1,399,by=2)

		for(i in 1:nIndividuals){

			genotypicValuesSim[i] <- sum(QTL_table[i,evenIndices]*0.5)+sum(QTL_table[i,oddIndices]*-0.5)

		}
}

	if(nLoci == nMarkers){

		if(nMarkers %% 2 ==1){
			oddIndices <- seq(1,nMarkers,by=2)
			evenIndices<- seq(2,(nMarkers-1),by=2)
		}

		else if(nMarkers %% 2 ==0){
			oddIndices <- seq(1,(nMarkers-1),by=2)
			evenIndices<- seq(2,nMarkers,by=2)
		}

		for(i in 1:nIndividuals){

			genotypicValuesSim[i] <- sum(QTL_table[i,evenIndices]*0.05)+sum(QTL_table[i,oddIndices]*-0.05)

		}
    }

    return(genotypicValuesSim)
}

###########################################################################################
## Fn 14:
######### Function to simulate Phenotypic Values ##########################################

simulatePhenotypicValues <- function(nextGenGenoTable,varE,noLoci){


    genoTable <- nextGenGenoTable
    nIndividuals <-nrow(nextGenGenoTable)
    genotypicValuesSim <-rep(0,nIndividuals)
    phenotypicValuesSim <-rep(0,nIndividuals)
    nMarkers <- ncol(genoTable)
    errorVar <- varE
    errorSD <- sqrt(errorVar)
    nLoci <- noLoci
#######################################################################
##### Extract QTL table ###### ###### ###### ###### ###### ###### ######

    if(nLoci== 40){
      QTL_indices<- seq(50,(nMarkers/2),by=53)
    }


    if(nLoci == 400){
      QTL_indices_1<- seq(8,nMarkers,by=12)
      QTL_indices_2 <- seq(50,nMarkers,by=107)
      n_3 <- 400-(length(QTL_indices_1))-(length(QTL_indices_2))
      QTL_indices_3<- sample(c(1:nMarkers),n_3)

      QTL_indices <- c(QTL_indices_1,QTL_indices_2,QTL_indices_3)
    }


    if(nLoci == nMarkers){
      QTL_indices<- c(1:nMarkers)
    }


    QTL_table <- genoTable[,c(QTL_indices)]

############### Simulate Genotypic Values #############################

    if(nLoci == 40){

      evenIndices<- seq(2,40,by=2)
      oddIndices <- seq(1,39,by=2)


      for(i in 1:nIndividuals){

        genotypicValuesSim[i] <- sum(QTL_table[i,evenIndices]*5)+sum(QTL_table[i,oddIndices]*-5)
        phenotypicValuesSim[i] <- genotypicValuesSim[i] + rnorm(1,mean=0,sd=errorSD)
      }
    }

    if(nLoci == 400){

      evenIndices<- seq(2,400,by=2)
      oddIndices <- seq(1,399,by=2)

      for(i in 1:nIndividuals){

        genotypicValuesSim[i] <- sum(QTL_table[i,evenIndices]*0.5)+sum(QTL_table[i,oddIndices]*-0.5)
        phenotypicValuesSim[i] <- genotypicValuesSim[i] + rnorm(1,mean=0,sd=errorSD)

      }
    }

    if(nLoci == nMarkers){

      if(nMarkers %% 2 ==1){
        oddIndices <- seq(1,nMarkers,by=2)
        evenIndices<- seq(2,(nMarkers-1),by=2)
      }

      else if(nMarkers %% 2 ==0){
        oddIndices <- seq(1,(nMarkers-1),by=2)
        evenIndices<- seq(2,nMarkers,by=2)
      }

      for(i in 1:nIndividuals){

        genotypicValuesSim[i] <- sum(QTL_table[i,evenIndices]*0.05)+sum(QTL_table[i,oddIndices]*-0.05)
        phenotypicValuesSim[i] <- genotypicValuesSim[i] + rnorm(1,mean=0,sd=errorSD)
      }
    }

    return(phenotypicValuesSim)
}

## Fn 15:
### Test Simulated Phenotypic Values
#########################################################################
getQTLTable <- function(GenoTable,noLoci){


	genoTable <- GenoTable
    nIndividuals <-nrow(genoTable)
    genotypicValuesSim <-rep(0,nIndividuals)
    nLoci <- noLoci
	nMarkers <- ncol(genoTable)
    ##### Extract QTL table ######

    if(nLoci== 40){
      QTL_indices<- seq(50,(nMarkers/2),by=53)
    }


    if(nLoci == 400){
      QTL_indices_1<- seq(8,nMarkers,by=12)
      QTL_indices_2 <- seq(50,nMarkers,by=107)
      n_3 <- 400-(length(QTL_indices_1))-(length(QTL_indices_2))
      QTL_indices_3<- sample(c(1:nMarkers),n_3)

      QTL_indices <- c(QTL_indices_1,QTL_indices_2,QTL_indices_3)
    }


    if(nLoci == nMarkers){
      QTL_indices<- c(1:nMarkers)
    }


    QTL_table <- genoTable[,c(QTL_indices)]


	return(QTL_table)

}

###

## Fn 16:
#### Generate Weighted Genotypic Table ##################################

generateWeightedGenoTable <- function(nextGenGenoTable,noLoci) {

	  genoTable <- nextGenGenoTable
	  nLoci <- noLoci
	  nMarkers <- ncol(genoTable)

	  if(nLoci== 40){
		  QTL_indices<- seq(50,(nMarkers/2),by=53)
		}


	  if(nLoci == 400){
		QTL_indices_1<- seq(8,nMarkers,by=12)
		QTL_indices_2 <- seq(50,nMarkers,by=107)
		n_3 <- 400-(length(QTL_indices_1))-(length(QTL_indices_2))
		QTL_indices_3<- sample(c(1:nMarkers),n_3)

		QTL_indices <- c(QTL_indices_1,QTL_indices_2,QTL_indices_3)
	  }


	  if(nLoci == nMarkers){
		QTL_indices<- c(1:nMarkers)
	  }



	  QTL_table <- genoTable[,c(QTL_indices)]
	  QTL_table_New <- QTL_table

	############ Weighting QTL based on QTL effect  #########################################


	  for(i in 1:nrow(genoTable)){



		if(nLoci==40) {
		  evenIndices<- seq(2,40,by=2)
		  oddIndices <- seq(1,39,by=2)

		  QTL_table_New[i,evenIndices] <- QTL_table[i,evenIndices]*2
		  QTL_table_New[i,oddIndices] <- QTL_table[i,oddIndices]*-2

		}

		if(nLoci==400) {
		  evenIndices<- seq(2,400,by=2)
		  oddIndices <- seq(1,399,by=2)

		  QTL_table_New[i,evenIndices] <- QTL_table[i,evenIndices]*2
		  QTL_table_New[i,oddIndices] <- QTL_table[i,oddIndices]*-2

		}


		if(nLoci==nMarkers){

		  if(nMarkers %% 2 ==1){
			oddIndices <- seq(1,nMarkers,by=2)
			evenIndices<- seq(2,(nMarkers-1),by=2)
		  }

		  else if(nMarkers %% 2 ==0){
			oddIndices <- seq(1,(nMarkers-1),by=2)
			evenIndices<- seq(2,nMarkers,by=2)
		  }

			  QTL_table_New[i,evenIndices] <- QTL_table[i,evenIndices]*2
			  QTL_table_New[i,oddIndices] <- QTL_table[i,oddIndices]*-2
		}
	  }

	  genoTable_New <- genoTable
	  genoTable_New[,c(QTL_indices)]<- QTL_table_New


	  return(genoTable_New)
}

# testweightGeno <- generateWeightedGenoTable(nextGenGenoTable,400)





### Function to estimate average percent heterozygous across all sampled loci...
### GenoTable in minus101 format ################################################################################
## Fn 18

  getPercentHeterozygous<- function(nextgenoTable){

    genoTable_Mod <- nextgenoTable
    nMarkers <- dim(genoTable_Mod)[2]
	nIndividuals <- dim(genoTable_Mod)[1]
    markerTest <- rep(0,nMarkers)
    no_genotype_1<- rep(0,nMarkers)
    no_genotype_0 <- rep(0,nMarkers)
    no_genotype_minus1 <- rep(0,nMarkers)
    percent_heterozygous_Mod <-rep(0,nMarkers)


    for(i in 1:nMarkers){

      markerTest[i]<- length(levels(factor(genoTable_Mod[,i])))
      no_genotype_1[i] <-length(genoTable_Mod[genoTable_Mod[,i]==1,i])
      no_genotype_0[i] <- length(genoTable_Mod[genoTable_Mod[,i]==0,i])
      no_genotype_minus1[i] <- length(genoTable_Mod[genoTable_Mod[,i]==-1,i])


      percent_heterozygous_Mod[i] <- no_genotype_0[i]/(no_genotype_1[i] + no_genotype_0[i] + no_genotype_minus1[i])

    }

    Avg_percent_heterozygous_Mod <- mean(percent_heterozygous_Mod)

    return(Avg_percent_heterozygous_Mod)
 }

# testpercentHet <- getPercentHeterozygous(nextGenGenoTable)

## Fn 19:
###############################################################################
## Function to get parent combination indices

getParentCombnIndices <- function(no_selected) {

	Parent_Combn_indices <- c()

	if(no_selected ==20) {

		combn_ParentVectorIndices <- c(2,3)
		for(l in 2:20){
			combn_ParentVectorIndices <- cbind(combn_ParentVectorIndices,c(1,l))
		}

		Parent_Combn_indices <- as.matrix(combn_ParentVectorIndices)


	}

	if(no_selected ==50) {

		combn_ParentVectorIndices <- c(2,3)

		for(l in 2:25){
			combn_ParentVectorIndices <- cbind(combn_ParentVectorIndices,c(1,l))
		}

		for(m in 25:50){
			combn_ParentVectorIndices <- cbind(combn_ParentVectorIndices,c(2,m))
		}

		Parent_Combn_indices <- as.matrix(combn_ParentVectorIndices)
	}


	if(no_selected ==100) {

		combn_ParentVectorIndices <- c(2,3)
		for(l in 2:20){
			combn_ParentVectorIndices <- cbind(combn_ParentVectorIndices,c(1,l))
		}

		Parent_Combn_indices <- as.matrix(combn_ParentVectorIndices)

	}

	if(no_selected ==200) {

		combn_ParentVectorIndices <- c(9,10)
		init  <- 2
		final  <- 25

		for(k in 1:8){

			for(m in init:final){

				combn_ParentVectorIndices <- cbind(combn_ParentVectorIndices,c(k,m))
			}
			init <- final+1
			final <- final +25
		}

		Parent_Combn_indices <- as.matrix(combn_ParentVectorIndices)

	}

	return(Parent_Combn_indices)



}

####################################################################################################
## Functions with meiosis code in Cpp ##############################################################
## Fn 20:
generateProgenyCpp <- function(nMarkers,nProgeny,RF,Parent1,Parent2){

		library(Rcpp)
		sourceCpp('mei.cpp')
		sourceCpp('mei2.cpp')

		m <- nMarkers
		n <- 2*nProgeny
		rr <- RF
		a1<- Parent1[,1]
		b1<- Parent1[,2]

		Pro <- mei(m,n,rr,a1,b1)

		a2<- Parent2[,1]
		b2<- Parent2[,2]

	    Pro2<- mei2(m,n,rr,a2,b2,Pro)

		return(Pro2)
	}

##################################################################################################################################################################################################### Functions  for Discrete selection - Fn 21-23 ###############################################################
#########################################################################################################################################################
## Fn 21:

runSimulations20X_DiscreteSelection <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap){

### Assign Variables #########################################################################

	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues
	  nFamilies <- 20
	  nSel_inFamily <- no_selected/20

	   if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }
	  nCrosses_inFamily <- nSel_inFamily
	  NAM_LinkMap_New <- NAM_LinkMap

###############################################################################################


	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny
	  nProgeny <- nProgeny/nSel_inFamily

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }


	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]

	  }

	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList

	  varE <- varEList

  	  nMarkers <- 4289
### Get predicted geno and pheno values for cycle1

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
      newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

#

## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()



	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }



### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)
	  # percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,nSel_inFamily*nCycles),nrow=nSel_inFamily,ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,nSel_inFamily)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,nSel_inFamily)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()


######## Functions on Cycle 1 Geno Data #################################################################

######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   } else if (selectionOnSimulated==TRUE) {

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	    }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  attainedGenoValues_List[[nFamily]][1] <- max(topNGenoValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


########### Functions on Cycle 1 Pheno Data#############################################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- max(topNPhenoValues)
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

###########################################################################

###########################################################################

# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################
   	for( i in 2:nCycles){


		selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }

	#ParentVectorIndices <- sample(c(1:(no_selected)),(no_selected),replace=FALSE)
	# combn_ParentVectorIndices <- combn(ParentVectorIndices,2)
	# Parent_Combn_indicies <- combn_ParentVectorIndices[,sample(c(1:ncol( combn_ParentVectorIndices)),nCrosses)]



	# Parent_Combn_indices <- getParentCombnIndices(no_selected)


	  Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nSel_inFamily,nMarkers,2,nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nSel_inFamily))

	  for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F2 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))


###################################################################################
		if(nSel_inFamily==2){
				nCrosses_inFamily <- 1
		}else{nCrosses_inFamily <- nSel_inFamily}

	    # nCrosses <- nSel_inFamily

		for( j in 1:(nCrosses_inFamily)){


		  if(j ==1 && nCrosses_inFamily==1){
			  Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=nCrosses_inFamily-1)) {
		       Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][1,,]
		  }

		Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)

		for(m in 1:nProgeny){
		  Parent1<- Cycle_Progeny_F1[j,,,1]
		  Parent2<- Cycle_Progeny_F1[j,,,1]
		  progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,m]  <- progeny1
		}

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

			for(nSel in 1:nSel_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format


		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nSel_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}




# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

###########################################################################################3

####################################################################################################################################################################################################

#nextGenGenoTable <- generate012GenoFormat(Cycle_Progeny_F5,no_selected)

### This section works with full data table

  		nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,cycle1_nProgeny,nFamilies,no_selected)

		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)



## Split Geno and Pheno values according to families  ########################

	    initIndex <- 1
	    finalIndex <- cycle1_nProgeny
	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()

	    for(nFamily in 1:20){

			genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
			phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
			genoValSimValues_List[[nFamily]] <- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues[initIndex:finalIndex]
			initIndex <- finalIndex+1
			finalIndex <- finalIndex+cycle1_nProgeny

	    }

### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
		}

       # if(i==2){
			# GSTable_afterMigration <- migrate2X(genoSelectionTable,1,2)
			# PSTable_afterMigration <- migrate2X(phenoSelectionTable,1,2)

			# genoSelectionTable_afterMigration  <- GSTable_afterMigration[[1]]
			# phenoSelectionTable_afterMigration <- PSTable_afterMigration[[1]]
			# emigrantGroups <- GSTable_afterMigration[[2]]
			# immigrantGroups<- GSTable_afterMigration[[3]]

		# } else if(i>2){

			# genoSelectionTable_afterMigration <- migrate(genoSelectionTable,1,emigrantGroups,immigrantGroups)
			# phenoSelectionTable_afterMigration <- migrate(phenoSelectionTable,1,emigrantGroups,immigrantGroups)
		# }

##############################################################################################

		
## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

 ####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}

	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List)

      return(simResults_List)


}

## Test Discrete selection simulation

 # i <-1
 # k <- 1

# simResultsTrial <- (runSimulations20X_DiscreteSelection(model_Sol_Pheno_List[[i]][[k]],model_Sol_Geno_List[[i]][[k]],F5_Progeny_ListReps[[k]][[i]],genotypicValuesSimListReps[[k]][[i]],phenotypeSimListReps[[k]][[i]],errorVarListReps[[k]][[i]],no_QTL[i],h2[i],numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues))

###################################################################################################

## Fn 22:

runSimulations20X_DiscreteSelection_GM <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap){

### Assign Variables #########################################################################

	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues
	  nFamilies <- 20
	  nSel_inFamily <- no_selected/20

	  nMarkers <-4289

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }
	  nCrosses_inFamily <- nSel_inFamily
	  NAM_LinkMap_New <- NAM_LinkMap

###############################################################################################

	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny
	  nProgeny <- nProgeny/nSel_inFamily

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }


	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]

	  }

	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList

	  varE <- varEList


### Get predicted geno and pheno values for cycle1

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
      newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

#############################################################################################

## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()



	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }



### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_NX_N_3c <- matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)
	  # percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()


######## Functions on Cycle 1 Geno Data #################################################################

######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   } else if (selectionOnSimulated==TRUE) {

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	    }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  attainedGenoValues_List[[nFamily]][1] <- max(topNGenoValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


########### Functions on Cycle 1 Pheno Data#############################################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- max(topNPhenoValues)
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

###########################################################################

###########################################################################

# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################
   	for(i in 2:nCycles){

	  	print(i)

		selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }

	#ParentVectorIndices <- sample(c(1:(no_selected)),(no_selected),replace=FALSE)
	# combn_ParentVectorIndices <- combn(ParentVectorIndices,2)
	# Parent_Combn_indicies <- combn_ParentVectorIndices[,sample(c(1:ncol( combn_ParentVectorIndices)),nCrosses)]



	# Parent_Combn_indices <- getParentCombnIndices(no_selected)


	  Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nSel_inFamily,nMarkers,2,nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nSel_inFamily))

	  for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F2 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))


###################################################################################
		if(nSel_inFamily==2){
				nCrosses_inFamily <- 1
		}else{nCrosses_inFamily <- nSel_inFamily}

	    # nCrosses <- nSel_inFamily

		for( j in 1:(nCrosses_inFamily)){


		  if(j ==1 && nCrosses_inFamily==1){
			  Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=nCrosses_inFamily-1)) {
		       Parent1<-  selectedGenoData_List[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][1,,]
		  }

##  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

##mod 12/16

		Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)


		 for(m in 1:nProgeny){

          Parent1<- Cycle_Progeny_F1[j,,,m]
		  Parent2<- Cycle_Progeny_F1[j,,,m]

		  progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,m]  <- progeny1

		 }

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

			for(nSel in 1:nSel_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format


		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nSel_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}




# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

###########################################################################################3

####################################################################################################################################################################################################

#nextGenGenoTable <- generate012GenoFormat(Cycle_Progeny_F5,no_selected)

### This section works with full data table

  		nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,cycle1_nProgeny,nFamilies,no_selected)

		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)



## Split Geno and Pheno values according to families  ########################

	    initIndex <- 1
	    finalIndex <- cycle1_nProgeny
	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()

	    for(nFamily in 1:20){

			genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
			phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
			genoValSimValues_List[[nFamily]] <- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues[initIndex:finalIndex]
			initIndex <- finalIndex+1
			finalIndex <- finalIndex+cycle1_nProgeny

	    }

### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
		}


##############################################################################################

		## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

 ####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}

	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List)

      return(simResults_List)


}


################################################################################################################
## Fn 22 A:

runSimulations20X_DiscreteSelection_GM_PG <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap){

### Assign Variables #########################################################################

	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues
	  nFamilies <- 20
	  nSel_inFamily <- no_selected/20

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }
	  nCrosses_inFamily <- nSel_inFamily
      nMarkers <- 4289
	  NAM_LinkMap <- NAM_LinkMap_New
###############################################################################################

	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny
	  nProgeny <- nProgeny/nSel_inFamily

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }


	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]

	  }

	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList

	  varE <- varEList


	  i<-1
	  nCyc <-1

      percentHeterozygous_List<- rep(list(rep(0,20)),nFamilies)
	  percentHeterozygous_QTL_List<- rep(list(rep(0,20)),nFamilies)
	  FavorableAlleleFreq <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableAllele <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableQTLAlleleFreq <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableQTLAllele <- rep(list(rep(list(list()),20)),nFamilies)
	  percentHet <- rep(0,20)

	  QTL_GenoTable_List <- list()


	  Diff_Stats_List_Global <- list()
	  Diff_Stats_List_Local <- list()




### Get predicted geno and pheno values for cycle1

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
      newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

#############################################################################################



######## Get percentHeterozygous & favorableAllele frequency
     if (selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[1])
	 } else if(selectionOnGeno==FALSE){ 
		MarkerEffects <- unlist(PredictionModel_Pheno[1])
	 }
	  
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:20){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[nFamily]],1,cycle1_nProgeny)
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

### get island PG parameters

		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")

		Diff_Stats <- diff_stats(nextGenGenoTable_GInd)
	    Diff_Stats_List_Local[[nCyc]] <- as.big.matrix(Diff_Stats[[1]])
	    Diff_Stats_List_Global[[nCyc]] <- Diff_Stats[[2]]

		rm(nextGenGenoTable_GInd)
		rm(Diff_Stats)
		gc()



## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()



	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }



### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_NX_N_3c <- matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)
	  # percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   } else if (selectionOnSimulated==TRUE) {

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	    }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  attainedGenoValues_List[[nFamily]][1] <- max(topNGenoValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


########### Functions on Cycle 1 Pheno Data#############################################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- max(topNPhenoValues)
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

###########################################################################

###########################################################################

# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################
   	for( i in 2:nCycles){

	  	print(i)
		nCyc <-i
		selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }

	#ParentVectorIndices <- sample(c(1:(no_selected)),(no_selected),replace=FALSE)
	# combn_ParentVectorIndices <- combn(ParentVectorIndices,2)
	# Parent_Combn_indicies <- combn_ParentVectorIndices[,sample(c(1:ncol( combn_ParentVectorIndices)),nCrosses)]



	# Parent_Combn_indices <- getParentCombnIndices(no_selected)


	  Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nSel_inFamily,nMarkers,2,nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nSel_inFamily))

	  for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F2 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))


###################################################################################
		if(nSel_inFamily==2){
				nCrosses_inFamily <- 1
		}else{nCrosses_inFamily <- nSel_inFamily}

	    # nCrosses <- nSel_inFamily

		for( j in 1:(nCrosses_inFamily)){


		  if(j ==1 && nCrosses_inFamily==1){
			  Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=nCrosses_inFamily-1)) {
		       Parent1<-  selectedGenoData_List[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][1,,]
		  }

##  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

##mod 12/16

		Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)


		 for(m in 1:nProgeny){

          Parent1<- Cycle_Progeny_F1[j,,,m]
		  Parent2<- Cycle_Progeny_F1[j,,,m]

		  progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,m]  <- progeny1

		 }

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

			for(nSel in 1:nSel_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format


		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nSel_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}




# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

###########################################################################################3

####################################################################################################################################################################################################

#nextGenGenoTable <- generate012GenoFormat(Cycle_Progeny_F5,no_selected)

### This section works with full data table

  		nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,cycle1_nProgeny,nFamilies,no_selected)

		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)



		for(nFamily in 1:20){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)
			QTL_GenoTable_List[[nFamily]] <- getQTLTable(nextGenGenoTable_List[[nFamily]],no_QTL)

			FavorableAlleleFreq[[nFamily]][[nCyc]] <- getFreq(nextGenGenoTable_List[[nFamily]],FavorableAllele[[nFamily]][[1]])
			percentHeterozygous_List[[nFamily]][nCyc] <- getPercentHeterozygous(nextGenGenoTable_List[[nFamily]])

			FavorableQTLAlleleFreq[[nFamily]][[nCyc]] <- getFreq(QTL_GenoTable_List[[nFamily]],FavorableQTLAllele[[nFamily]][[1]])
			percentHeterozygous_QTL_List[[nFamily]][nCyc] <- getPercentHeterozygous(QTL_GenoTable_List[[nFamily]])

		}



### get island PG parameters

		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")

		Diff_Stats <- diff_stats(nextGenGenoTable_GInd)
	    Diff_Stats_List_Local[[nCyc]] <- as.big.matrix(Diff_Stats[[1]])
	    Diff_Stats_List_Global[[nCyc]] <- Diff_Stats[[2]]

		rm(nextGenGenoTable_GInd)
		rm(Diff_Stats)
		gc()


#Split Geno and Pheno values according to families  ########################

	    initIndex <- 1
	    finalIndex <- cycle1_nProgeny
	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()

	    for(nFamily in 1:20){

			genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
			phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
			genoValSimValues_List[[nFamily]] <- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues[initIndex:finalIndex]
			initIndex <- finalIndex+1
			finalIndex <- finalIndex+cycle1_nProgeny

	    }


### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
		}


##############################################################################################
## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

 ####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}

	}


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List,QTL_GenoTable_List,FavorableAlleleFreq,percentHeterozygous_List,FavorableQTLAlleleFreq,percentHeterozygous_QTL_List,Diff_Stats_List_Local,Diff_Stats_List_Global)

      return(simResults_List)


}

## Fn 22 B


runSimulations20X_DiscreteSelection_GM_SelProgeny <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap){

### Assign Variables #########################################################################

	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues
	  nFamilies <- 20
	  nSel_inFamily <- no_selected/20

	  nMarkers <- 4289
	  
	  NAM_LinkMap_New <- NAM_LinkMap

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }
	  nCrosses_inFamily <- nSel_inFamily
      
###############################################################################################

	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny
	  nProgeny <- nProgeny/nSel_inFamily

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }


	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]

	  }

	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList

	  varE <- varEList

      NAM_LinkMap_New <- NAM_LinkMap

### Get predicted geno and pheno values for cycle1

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
      newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

#############################################################################################

## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()



	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }



### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_NX_N_3c <- matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)
	  # percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Sort genotype and phenotype data for Cycle 1 data #################################################################
   
    for(nFamily in 1:20){


        if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	    } else if (selectionOnSimulated==TRUE) {

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	    }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  attainedGenoValues_List[[nFamily]][1] <- max(topNGenoValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


########### Output variables for Cycle 1 Pheno Data #################################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- max(topNPhenoValues)
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

###########################################################################

###########################################################################

# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################
   	for(i in 2:nCycles){

	  	print(i)

		selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }

	
	  Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nSel_inFamily,nMarkers,2,nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nSel_inFamily))

	  for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F2 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))


###################################################################################
		if(nSel_inFamily==2){
				nCrosses_inFamily <- 1
		}else{nCrosses_inFamily <- nSel_inFamily}

	    # nCrosses <- nSel_inFamily

		for( j in 1:(nCrosses_inFamily)){


		  if(j ==1 && nCrosses_inFamily==1){
			  Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=nCrosses_inFamily-1)) {
		       Parent1<-  selectedGenoData_List[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][1,,]
		  }

##  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

##mod 12/16

		Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)


		 for(m in 1:nProgeny){

          Parent1<- Cycle_Progeny_F1[j,,,m]
		  Parent2<- Cycle_Progeny_F1[j,,,m]

		  progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,m]  <- progeny1

		 }

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

			for(nSel in 1:nSel_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format


		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nSel_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}




# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

###########################################################################################3

####################################################################################################################################################################################################

#nextGenGenoTable <- generate012GenoFormat(Cycle_Progeny_F5,no_selected)

### This section works with full data table

  		nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,cycle1_nProgeny,nFamilies,no_selected)

		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	    initIndex <- 1
	    finalIndex <- nProgeny*nCrosses_inFamily


	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()
		
		for(nFamily in 1:20){


		    genoValues.family <- genoValues[initIndex:finalIndex]
			phenoValues.family <- phenoValues[initIndex:finalIndex]
			genoValSimValues.family<- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues.family <- phenoValSimValues[initIndex:finalIndex]


			initIndex.Cross <-1
			finalIndex.Cross <- nProgeny

			genoValues.selected <- c()
			phenoValues.selected <- c()
			genoValSimValues.selected <- c()
			phenoValSimValues.selected <- c()



				genoValues.sorted <- sort(genoValues.family,decreasing=TRUE,index.return=TRUE)
				genoValues.selFamily <- genoValues.sorted[[1]][1:nProgeny]

				phenoValues.sorted <- sort(phenoValues.family,decreasing=TRUE,index.return=TRUE)
				phenoValues.selFamily <- phenoValues.sorted[[1]][1:nProgeny]

				genoValSimValues.selCross <- genoValSimValues.family[initIndex.Cross:finalIndex.Cross]
				genoValSimValues.sorted <- sort(genoValSimValues.family,decreasing=TRUE,index.return=TRUE)
				genoValSimValues.selFamily <- genoValSimValues.sorted[[1]][1:nProgeny]

				phenoValSimValues.selCross <- phenoValSimValues.family[initIndex.Cross:finalIndex.Cross]
				phenoValSimValues.sorted <- sort(phenoValSimValues.family,decreasing=TRUE,index.return=TRUE)
				phenoValSimValues.selFamily <- phenoValSimValues.sorted[[1]][1:nProgeny]


				genoValues_List[[nFamily]] <- genoValues.selFamily
				phenoValues_List[[nFamily]] <- phenoValues.selFamily

				genoValSimValues_List[[nFamily]] <- genoValSimValues.selFamily
				phenoValSimValues_List[[nFamily]] <- phenoValSimValues.selFamily

				initIndex <- finalIndex+1
				finalIndex <- finalIndex + (nProgeny*nCrosses_inFamily)
		}



### Apply selection for discrete families ######################################
### Selection based on simulated or predicted values

		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
		}


##############################################################################################
## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

 ####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}
					
				

	}
	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List)

      return(simResults_List)


}



## Fn 23:

runSimulations20X_DiscreteSelection_GM_SelProgeny_PG <- function(model_Sol_Pheno,model_Sol_Geno,F5_Progeny_List,genoValSimList,phenoValSimList,varEList,noQTLs,H2,numberSelected,noCrosses,no_Progeny,noCycles,selectionOnGenoValues,selectionOnSimulatedValues,NAM_LinkMap){

### Assign Variables #########################################################################

	  h2 <- H2
	  no_QTL <-noQTLs

	  no_selected <- numberSelected
	  nCrosses <- noCrosses
	  nProgeny <- no_Progeny
	  nCycles <- noCycles
	  selectionOnGeno <- selectionOnGenoValues
	  selectionOnSimulated <- selectionOnSimulatedValues
	  nFamilies <- 20
	  nSel_inFamily <- no_selected/20
	  nMarkers <- 4289
	  NAM_LinkMap_New <- NAM_LinkMap

	  if(nSel_inFamily ==1){
	     nSel_inFamily <- 2
	  }
	  nCrosses_inFamily <- nSel_inFamily

###############################################################################################

	  PredictionModel_Pheno <- model_Sol_Pheno
	  PredictionModel_Geno <- model_Sol_Geno

	  cycle1_Geno_Data <- F5_Progeny_List
	  cycle1_nProgeny <- nProgeny
	  nProgeny <- nProgeny/nSel_inFamily

	  F5_Progeny_List_Family<- F5_Progeny_List

	  F5_Progeny_List<- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		F5_Progeny_List[[nFamily]] <- F5_Progeny_List_Family[nFamily,,,]

	  }


	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)

	  for(nFamily in 1:20){

		cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]

	  }

	  genoSimValues <- genoValSimList
	  phenoSimValues <- phenoValSimList

	  varE <- varEList
	  NAM_LinkMap_New <- NAM_LinkMap

	  i<-1
	  nCyc <-1

	   

      percentHeterozygous_List<- rep(list(rep(0,20)),nFamilies)
	  percentHeterozygous_QTL_List<- rep(list(rep(0,20)),nFamilies)
	  FavorableAlleleFreq <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableAllele <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableQTLAlleleFreq <- rep(list(rep(list(list()),20)),nFamilies)
	  FavorableQTLAllele <- rep(list(rep(list(list()),20)),nFamilies)
	  percentHet <- rep(0,20)

	  QTL_GenoTable_List <- list()


	  Diff_Stats_List_Global <- list()
	  Diff_Stats_List_Local <- list()




### Get predicted geno and pheno values for cycle1

	  nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
      newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)
      genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
      phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

#############################################################################################



	######## Get percentHeterozygous & favorableAllele frequency
      if(selectionOnGeno==TRUE){
		MarkerEffects <- unlist(PredictionModel_Geno[1])
	  } else if(selectionOnGeno==FALSE){ 
		MarkerEffects <- unlist(PredictionModel_Pheno[1])
	  }
	  
	  cycle1_GenoData_List <- rep(list(array(0,c(nMarkers,2,cycle1_nProgeny))),nFamilies)
	  nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)

	  for(nFamily in 1:20){

			cycle1_GenoData_List[[nFamily]] <- cycle1_Geno_Data[nFamily,,,]
			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(cycle1_GenoData_List[[nFamily]],1,cycle1_nProgeny)
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

### get island PG parameters

		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")

		Diff_Stats <- diff_stats(nextGenGenoTable_GInd)
	    Diff_Stats_List_Local[[nCyc]] <- as.big.matrix(Diff_Stats[[1]])
	    Diff_Stats_List_Global[[nCyc]] <- Diff_Stats[[2]]

		rm(nextGenGenoTable_GInd)
		rm(Diff_Stats)
		gc()



## Split Geno and Pheno values according to families  ########################

	  initIndex <- 1
	  finalIndex <- cycle1_nProgeny
	  genoValues_List <- list()
	  phenoValues_List <- list()
	  genoSimValues_List <- list()
	  phenoSimValues_List <- list()



	  for(nFamily in 1:20){

		genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
		phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
		genoSimValues_List[[nFamily]] <- (genoSimValues)[initIndex:finalIndex]
		phenoSimValues_List[[nFamily]] <- (phenoSimValues)[initIndex:finalIndex]


		initIndex <- finalIndex+1
		finalIndex <- finalIndex+cycle1_nProgeny

	  }



### Initiate Arrays, Lists and Variables ######################################

	  GenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  GenoVal_NX_N_3c <- matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)
	  GenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)

	  PhenoVal_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)
	  PhenoVal_NX_N_3c <- matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)
	  PhenoVal_Sim_NX_2k_3c <- matrix(rep(0,(nCrosses*cycle1_nProgeny)*nCycles),nrow=(nCrosses*cycle1_nProgeny),ncol=nCycles)


	  selectedGenoIndividualIndices2 <- rep(0,no_selected)
	  attainedGenoValues <- rep(0,nCycles)
	  attainedPhenoValues <- rep(0,nCycles)
	  # percent_heterozygous_NX_3c <- rep(0,nCycles)

	  GenoVal_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  GenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)),nFamilies)
	  GenoVal_Sim_NX_2k_3c_List <- rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)
	  PhenoVal_NX_N_3c_List <- rep(list(matrix(rep(0,no_selected*nCycles),nrow=no_selected,ncol=nCycles)),nFamilies)
	  PhenoVal_Sim_NX_2k_3c_List <-  rep(list(matrix(rep(0,(cycle1_nProgeny)*nCycles),nrow=(cycle1_nProgeny),ncol=nCycles)),nFamilies)

	  selectedGenoIndividualIndices2_List <- list()
	  attainedPhenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  attainedGenoValues_List <- rep(list(rep(0,20)),nFamilies)
	  selectedGenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)
	  selectedPhenoIndividualIndices_List <- rep(list(rep(0,no_selected)),nFamilies)


	  sortedGenoValues_List <- list()
	  topNGenoValues_List <- list()
	  GenoTableIndices_topN_List <- list()

	  sortedPhenoValues_List <- list()
	  topNPhenoValues_List <- list()
	  PhenoTableIndices_topN_List <- list()



######## Functions on Cycle 1 Geno Data #################################################################
    for(nFamily in 1:20){


      if(selectionOnSimulated == FALSE){
		  sortedGenoValues <- sort(genoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	   } else if (selectionOnSimulated==TRUE) {

		  sortedGenoValues <- sort(genoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNGenoValues <- sortedGenoValues[[1]][1:nSel_inFamily]
		  GenoTableIndices_topN <- sortedGenoValues[[2]][1:nSel_inFamily]


		  sortedPhenoValues <- sort(phenoSimValues_List[[nFamily]],decreasing=TRUE,index.return=TRUE)
		  topNPhenoValues <- sortedPhenoValues[[1]][1:nSel_inFamily]
		  PhenoTableIndices_topN <- sortedPhenoValues[[2]][1:nSel_inFamily]

	    }

	  sortedGenoValues_List[[nFamily]] <- sortedGenoValues
	  topNGenoValues_List[[nFamily]] <- topNGenoValues
	  GenoTableIndices_topN_List[[nFamily]] <- GenoTableIndices_topN

	  GenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- genoSimValues_List[[nFamily]]
	  GenoVal_NX_2k_3c_List[[nFamily]][,1] <- genoValues_List[[nFamily]]
	  GenoVal_NX_N_3c_List[[nFamily]][,1] <- topNGenoValues


	  attainedGenoValues_List[[nFamily]][1] <- max(topNGenoValues)
	  selectedGenoIndividualIndices_List[[nFamily]] <- GenoTableIndices_topN


########### Functions on Cycle 1 Pheno Data#############################################################

      sortedPhenoValues_List[[nFamily]] <- sortedPhenoValues
	  topNPhenoValues_List[[nFamily]] <- topNPhenoValues
	  PhenoTableIndices_topN_List[[nFamily]] <- PhenoTableIndices_topN

	  PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,1] <- phenoSimValues_List[[nFamily]]
	  PhenoVal_NX_2k_3c_List[[nFamily]][,1] <- phenoValues_List[[nFamily]]
	  PhenoVal_NX_N_3c_List[[nFamily]][,1] <- topNPhenoValues

	  attainedPhenoValues_List[[nFamily]][1] <- max(topNPhenoValues)
	  selectedPhenoIndividualIndices_List[[nFamily]] <- PhenoTableIndices_topN

    }

###########################################################################

###########################################################################

# percent_heterozygous_NX_3c[1] <- percentHetList

############### Cycles 2*29 ###############################################
   	for( i in 2:nCycles){

	  	print(i)
		nCyc <-i
		selectedGenoData_List <- list()
		selectedGenoData_List <- list()

		for(nFamily in 1:20){

			if(i==2){

				if(selectionOnGeno == TRUE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],cycle1_GenoData_List[[nFamily]],nSel_inFamily)
				}

			}

			if(i>2){

				if(selectionOnGeno == TRUE){

					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedGenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}else if(selectionOnGeno == FALSE){
					selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
				}

		    }
	    }

	

	  Cycle_Progeny_F5_List_temp<- rep(list(array(0,c(nSel_inFamily,nMarkers,2,nProgeny))),nFamilies)

	  Cycle_Progeny_F5_List <- rep(list(array(0,c(nMarkers,2,nSel_inFamily*nProgeny))),nFamilies)

	  Cycle_Progeny_F5 <- array(0,c(nFamilies,nMarkers,2,nProgeny*nSel_inFamily))

	  for(nFamily in 1:20) {



		Cycle_Progeny_F1 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F2 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F3 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F4 <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))

		Cycle_Progeny_F5_temp <- array(0,c(nSel_inFamily,nMarkers,2,nProgeny))


###################################################################################
		if(nSel_inFamily==2){
				nCrosses_inFamily <- 1
		}else{nCrosses_inFamily <- nSel_inFamily}

	    # nCrosses <- nSel_inFamily

		for( j in 1:(nCrosses_inFamily)){


		  if(j ==1 && nCrosses_inFamily==1){
			  Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			  Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if ((nCrosses_inFamily>1)&&(j<=nCrosses_inFamily-1)) {
		       Parent1<-  selectedGenoData_List[[nFamily]][1,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][j+1,,]
		  }else if(j==nCrosses_inFamily){

		       Parent1<-  selectedGenoData_List[[nFamily]][j,,]
			   Parent2<-  selectedGenoData_List[[nFamily]][1,,]
		  }

##  Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)

##mod 12/16

		Cycle_Progeny_F1[j,,,] <- generateProgeny(Parent1,Parent2,nProgeny,NAM_LinkMap_New)


		 for(m in 1:nProgeny){

          Parent1<- Cycle_Progeny_F1[j,,,m]
		  Parent2<- Cycle_Progeny_F1[j,,,m]

		  progeny1<- generateProgeny(Parent1,Parent2,1,NAM_LinkMap_New)
		  Cycle_Progeny_F2[j,,,m]  <- progeny1

		 }

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

			for(nSel in 1:nSel_inFamily){

				for(nProg in 1:nProgeny){

					Cycle_Progeny_F5[nFamily,,,nSel*nProg] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSel,,,nProg]
				}
			}
		}


### Get Cycle_Progeny F5 list data in correct format


		for(nFamily in 1:nFamilies){

		    for(nSelFamily in 1:nSel_inFamily){

				for(nSel in 1:nProgeny){

					Cycle_Progeny_F5_List[[nFamily]][,,nSelFamily*nSel] <- Cycle_Progeny_F5_List_temp[[nFamily]][nSelFamily,,,nSel]


				}
			}
		}




# dim(F5_Progeny_List[[1]])
# [1]    1 4289    2  100

###########################################################################################3

####################################################################################################################################################################################################

#nextGenGenoTable <- generate012GenoFormat(Cycle_Progeny_F5,no_selected)

### This section works with full data table

  		nextGenGenoTable <- generateMinus101GenoFormat(Cycle_Progeny_F5,nSel_inFamily,cycle1_nProgeny,nFamilies,no_selected)

		genoValSimValues <- simulateGenotypicValues(nextGenGenoTable,no_QTL)
		phenoValSimValues <- simulatePhenotypicValues(nextGenGenoTable,varE,no_QTL)

		newNextGenGenoTable <- generateWeightedGenoTable(nextGenGenoTable,no_QTL)

		genoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Geno)
		phenoValues <- PredictRRBLUPPhenoValues(newNextGenGenoTable,PredictionModel_Pheno)

		nextGenGenoTable_List <- rep(list(matrix(rep(0,nrow=100,ncol=nMarkers))),nFamilies)


		for(nFamily in 1:20){

			nextGenGenoTable_List[[nFamily]] <- generateMinus101GenoFormat_V2(Cycle_Progeny_F5_List[[nFamily]],1,nCrosses_inFamily*nProgeny)
			QTL_GenoTable_List[[nFamily]] <- getQTLTable(nextGenGenoTable_List[[nFamily]],no_QTL)

			FavorableAlleleFreq[[nFamily]][[nCyc]] <- getFreq(nextGenGenoTable_List[[nFamily]],FavorableAllele[[nFamily]][[1]])
			percentHeterozygous_List[[nFamily]][nCyc] <- getPercentHeterozygous(nextGenGenoTable_List[[nFamily]])

			FavorableQTLAlleleFreq[[nFamily]][[nCyc]] <- getFreq(QTL_GenoTable_List[[nFamily]],FavorableQTLAllele[[nFamily]][[1]])
			percentHeterozygous_QTL_List[[nFamily]][nCyc] <- getPercentHeterozygous(QTL_GenoTable_List[[nFamily]])

		}



### get island PG parameters
		nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)
 		nextGenGenoTable_GInd <- df2genind(nextGenGenoTable_AlleleFormat,pop=populations,sep="")

		Diff_Stats <- diff_stats(nextGenGenoTable_GInd)
	    Diff_Stats_List_Local[[nCyc]] <- as.big.matrix(Diff_Stats[[1]])
	    Diff_Stats_List_Global[[nCyc]] <- Diff_Stats[[2]]

		rm(nextGenGenoTable_GInd)
		rm(Diff_Stats)
		gc()


### Split Geno and Pheno values according to families  ########################
### Sort and select top individuals per cross per family

	    initIndex <- 1
	    finalIndex <- nProgeny*nCrosses_inFamily


	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()
		genoValues_List.fam <- list()

		for(nFamily in 1:20){


		    genoValues.family <- genoValues[initIndex:finalIndex]
			phenoValues.family <- phenoValues[initIndex:finalIndex]
			genoValSimValues.family<- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues.family <- phenoValSimValues[initIndex:finalIndex]


			initIndex.Cross <-1
			finalIndex.Cross <- nProgeny

			genoValues.selected <- c()
			phenoValues.selected <- c()
			genoValSimValues.selected <- c()
			phenoValSimValues.selected <- c()



				genoValues.sorted <- sort(genoValues.family,decreasing=TRUE,index.return=TRUE)
				genoValues.selFamily <- genoValues.sorted[[1]][1:nProgeny]

				phenoValues.sorted <- sort(phenoValues.family,decreasing=TRUE,index.return=TRUE)
				phenoValues.selFamily <- phenoValues.sorted[[1]][1:nProgeny]

				genoValSimValues.selCross <- genoValSimValues.family[initIndex.Cross:finalIndex.Cross]
				genoValSimValues.sorted <- sort(genoValSimValues.family,decreasing=TRUE,index.return=TRUE)
				genoValSimValues.selFamily <- genoValSimValues.sorted[[1]][1:nProgeny]

				phenoValSimValues.selCross <- phenoValSimValues.family[initIndex.Cross:finalIndex.Cross]
				phenoValSimValues.sorted <- sort(phenoValSimValues.family,decreasing=TRUE,index.return=TRUE)
				phenoValSimValues.selFamily <- phenoValSimValues.sorted[[1]][1:nProgeny]


				genoValues_List[[nFamily]] <- genoValues.selFamily
				phenoValues_List[[nFamily]] <- phenoValues.selFamily

				genoValSimValues_List[[nFamily]] <- genoValSimValues.selFamily
				phenoValSimValues_List[[nFamily]] <- phenoValSimValues.selFamily

				initIndex <- finalIndex+1
				finalIndex <- finalIndex + (nProgeny*nCrosses_inFamily)
		}




### Initiate Arrays, Lists and Variables ######################################


		if(selectionOnSimulated ==TRUE){
			genoSelectionTable <- ApplySelection_Discrete(genoValSimValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValSimValues_List,no_selected)
		}else if (selectionOnSimulated==FALSE){
			genoSelectionTable <- ApplySelection_Discrete(genoValues_List,no_selected)
			phenoSelectionTable <- ApplySelection_Discrete(phenoValues_List,no_selected)
		}


##############################################################################################

## Assign output variables for cycle n

		for(nFamily in 1:20){

			GenoVal_Sim_NX_2k_3c_List[[nFamily]][,i] <- genoValSimValues_List[[nFamily]]
			GenoVal_NX_2k_3c_List[[nFamily]][,i] <- genoValues_List[[nFamily]]

			PhenoVal_Sim_NX_2k_3c_List[[nFamily]][,i]<- phenoValSimValues_List[[nFamily]]
			PhenoVal_NX_2k_3c_List[[nFamily]][,i] <- phenoValues_List[[nFamily]]

 ####   
		
			if(selectionOnGeno==TRUE){
				
				selectedGenoIndividualIndices <- as.vector(genoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedGenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedGenoIndividualIndices]
				
				selectedGenoIndividualIndices_List[[nFamily]]<- selectedGenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}else if(selectionOnGeno==FALSE){
			
				selectedPhenoIndividualIndices <- as.vector(phenoSelectionTable[[nFamily]][,2])
				genoSimSelectedValues <- genoSimValues[selectedPhenoIndividualIndices]
				phenoSimSelectedValues <- phenoSimValues[selectedPhenoIndividualIndices]
				
				selectedPhenoIndividualIndices_List[[nFamily]] <- selectedPhenoIndividualIndices
				attainedGenoValues_List[[nFamily]][i] <- max(genoSimSelectedValues)
				attainedPhenoValues_List[[nFamily]][i] <- max(phenoSimSelectedValues)
				
			}
			
			GenoVal_NX_N_3c_List[[nFamily]][,i] <- genoSelectionTable[[nFamily]][,1]
			PhenoVal_NX_N_3c_List[[nFamily]][,i] <- phenoSelectionTable[[nFamily]][,1]
		
		
		}

    }


	  simResults_List <- list(GenoVal_Sim_NX_2k_3c_List,GenoVal_NX_2k_3c_List,GenoVal_NX_N_3c_List, PhenoVal_Sim_NX_2k_3c_List,PhenoVal_NX_2k_3c_List,PhenoVal_NX_N_3c_List,attainedGenoValues_List, attainedPhenoValues_List,QTL_GenoTable_List,FavorableAlleleFreq,percentHeterozygous_List,FavorableQTLAlleleFreq,percentHeterozygous_QTL_List,Diff_Stats_List_Local,Diff_Stats_List_Global)

      return(simResults_List)

  }




#####################################################################################################################################
########################################### Functions for Island model selection Fn 24-30 ########################################
#####################################################################################################################################


## Fn 32:
extractGenoData_FamilyWise <- function(GenoTableIndices1,ProgenyGenoData,no_selected){

    GenoTableIndices <- GenoTableIndices1
	progenyGenoData <- ProgenyGenoData
    nMarkers <- 4289


	nSelectedIndividuals <- length(GenoTableIndices)

	genoData <- array(0,c(nSelectedIndividuals,nMarkers,2))


		for(j in 1:nSelectedIndividuals){

			individualIndex <- GenoTableIndices[j]

			genoData[j,,] <- progenyGenoData[,,individualIndex]

		}

	return(genoData)

}

#########

## test extractGenoData_FamilyWise()
## genoData_List <- list()
## for(nFamily in 1:20){

    # # GenoTableIndices1<- testSelection_Discrete[[i]][,2]
	# # ProgenyGenoData<-(cycle1_GenoData_List[[i]])
	# # nSel_inFamily <-1
	# # genoData_List[[nFamily]] <- extractGenoData_FamilyWise(GenoTableIndices1,ProgenyGenoData,nSel_inFamily)
	# # }



### Fn 36:
####### get Frequency and state of favorable allele

getFreq_BasePoln<- function(GenoTable,markerEffects){

   	nMarkers <- dim(GenoTable)[2]
	nIndividuals <- dim(GenoTable)[1]
	sign_markerEffects <- sign(markerEffects)
	n_alleles <- 2*nIndividuals

	Freq <- rep(0,nMarkers)
	FavAllele <- rep(0,nMarkers)
	count_1 <-rep(0,nMarkers)
	count_0 <-rep(0,nMarkers)

	for(j in 1:nMarkers){



		p_1_P <- length(which(GenoTable[,j]==1))
		p_0_H <- length(which(GenoTable[,j]==0))
		p_minus1_Q <- length(which(GenoTable[,j]==-1))


		# p_1_P <- length(which(GenoVector==1))/nIndividuals
		# p_0_H <- length(which(GenoVector==0))/nIndividuals
		# p_minus1_Q <- length(which(GenoVector==-1))/nIndividuals

		count_1[j] <- ((2*p_1_P)+p_0_H )
		count_0[j] <- ((2*p_minus1_Q)+p_0_H )

		if((count_1[j] >count_0[j]) && sign_markerEffects[j]==1){
			FavAllele[j] <- 1
			Freq[j] <- count_1[j]/n_alleles
		}else if((count_0[j] > count_1[j]) && sign_markerEffects[j]==1){
			FavAllele[j] <- 0
			Freq[j] <- count_0[j]/n_alleles

		}


	}

	return(list(Freq,FavAllele))

}

# GenoTable <- genoTable_Mod
# markerEffects <- TrainModel
#########################################################################################

## Fn 37: 
## get frequency of favorable allele given genotype table and state of favorable allele

getFreq <- function(GenoTable,FavAllele){


	nMarkers <- dim(GenoTable)[2]
	nIndividuals <- dim(GenoTable)[1]
	Freq <- rep(0,nMarkers)
	count_1 <-rep(0,nMarkers)
	count_0 <-rep(0,nMarkers)
	n_alleles <- 2*nIndividuals

	for(j in 1:nMarkers){


		p_1_P <- length(which(GenoTable[,j]==1))
		p_0_H <- length(which(GenoTable[,j]==0))
		p_minus1_Q <- length(which(GenoTable[,j]==-1))


		count_1[j] <- ((2*p_1_P)+p_0_H )
		count_0[j] <- ((2*p_minus1_Q)+p_0_H )

		if(FavAllele[j]==1){
			Freq[j] <- count_1[j]/n_alleles
		}else if(FavAllele[j]==0){
			Freq[j] <- count_0[j]/n_alleles
		}
	}

	return(Freq)

}


###########################################

## Fn 38:
##########################

 getAlleleFormat<- function(GenoTable,AlleleConversionTable_Combined){

	genoTable<- GenoTable
	nIndividuals <- nrow(GenoTable)
	IA3023_alleles <- as.character(AlleleConversionTable_Combined[,8])
	Alt_Homozygous_alleles<- as.character(AlleleConversionTable_Combined[,10])
	Het_alleles <- as.character(AlleleConversionTable_Combined[,22])

	for(i in 1:nIndividuals){

		oneIndices <- which(genoTable[i,] ==1)
		zeroIndices <- which(genoTable[i,] ==0)
		minusOneIndices <- which(genoTable[i,] ==-1)

		genoTable[i,oneIndices] <- IA3023_alleles[oneIndices]
		genoTable[i,zeroIndices] <- Het_alleles[zeroIndices]
		genoTable[i,minusOneIndices] <- Alt_Homozygous_alleles[minusOneIndices]
	}

	return(genoTable)


 }

 # genoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)

 ####################################################################################################################################

## Rewrite Migrate function based on policy 
## Crossing Function 
## BD1 
### Generate F5 RIL progeny for selected geno data
# Parent_Combn_indices <- getParentCombnIndices(no_selected)

# BD <- "BD1"

# F5RILs_Trial <- getF5RILs_BD (BD,selectedGenoData_List_afterExchange,nFamilies,nSel_inFamily,nProgeny,nMarkers,NAM_LinkMap)

 
# BD <- "GM" 
 # nFamilies <- 19
# F5RILs <- getF5RILs_GM(BD,selectedGenoData_List,selectedParents_NumProgeny_Geno_List,nFamilies,nMarkers,NAM_LinkMap)

# for(nFamily in 1:20){ 

	# selectedGenoData_List[[nFamily]] <- extractGenoData_FamilyWise(selectedPhenoIndividualIndices_List[[nFamily]],Cycle_Progeny_F5_List[[nFamily]],nSel_inFamily)
	
# }	
	
	# selectedPhenoIndividualIndices_List <- selectedParents_Pheno_List
	# nSel_inFamily <- length(selectedParents_Pheno_List)
	# selectedParents_NumProgeny_Geno_List <- selectedParents_NumProgeny_Pheno_List


################################################################################ 
###############
## 06/12/2018

# Crit_List <-getCriterionValueList(BD,genoValues,phenoValues,genoValSimValues,phenoValSimValues,nCrosses_inFamily,nProgeny,FALSE,FALSE)

## Fn 44:
 getCriterionValueList <- function(BD,genoValues,phenoValues,genoValSimValues,phenoValSimValues,nCrosses_inFamily,nProgeny,selectionOnGeno,selectionOnSimulated,nFamilies){

## Split Geno and Pheno values according to families  ########################
 
   if(BD =="BD1"){ 
## BD1 

        cycle_nProgeny <- nProgeny*nCrosses_inFamily
	    initIndex <- 1
	    finalIndex <- cycle_nProgeny
	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()

	    for(nFamily in 1:nFamilies){

			genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
			phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
			genoValSimValues_List[[nFamily]] <- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues[initIndex:finalIndex]
			initIndex <- finalIndex+1
			finalIndex <- finalIndex+cycle_nProgeny

	    }

   } else if(BD=="GM"){   
## GM 
### Sort and select top individuals per cross per family

	    initIndex <- 1
	    finalIndex <- nProgeny*nCrosses_inFamily


	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()

	    for(nFamily in 1:nFamilies){


		    genoValues.Family <- genoValues[initIndex:finalIndex]
			phenoValues.Family <- phenoValues[initIndex:finalIndex]
			genoValSimValues.Family<- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues.Family <- phenoValSimValues[initIndex:finalIndex]
			genoValues_List[[nFamily]] <- genoValues.Family
			phenoValues_List[[nFamily]] <- phenoValues.Family

			genoValSimValues_List[[nFamily]] <- genoValSimValues.Family
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues.Family



			initIndex <- finalIndex+1
			finalIndex <- finalIndex + (nProgeny*nCrosses_inFamily)

		}
	
	}else if(BD=="GM_SP"){
## GM_SP		
### Sort and select top individuals per cross per family

        cycle_nProgeny <- nProgeny*nCrosses_inFamily

	    initIndex <- 1
	    finalIndex <- nProgeny*nCrosses_inFamily


	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()
		
	    for(nFamily in 1:nFamilies){


		    genoValues.family <- genoValues[initIndex:finalIndex]
			phenoValues.family <- phenoValues[initIndex:finalIndex]
			genoSimValues.family<- genoValSimValues[initIndex:finalIndex]
			phenoSimValues.family <- phenoValSimValues[initIndex:finalIndex]

			
			genoValues.sorted <- sort(genoValues.family,decreasing=TRUE,index.return=TRUE)
		    genoSimValues.sorted <- sort(genoSimValues.family,decreasing=TRUE,index.return=TRUE)
            phenoValues.sorted <- sort(phenoValues.family,decreasing=TRUE,index.return=TRUE) 
			phenoSimValues.sorted <- sort(phenoSimValues.family,decreasing=TRUE,index.return=TRUE)
			
			if(selectionOnGeno==TRUE){ 
			
				if(selectionOnSimulated==FALSE){
					selectedIndices <- genoValues.sorted[[2]][1:nProgeny]
				}else if(selectionOnSimulated==TRUE){ 
				    selectedIndices <- genoSimValues.sorted[[2]][1:nProgeny]
				}
			}else if(selectionOnGeno==FALSE){ 
			
				if(selectionOnSimulated==FALSE){
					selectedIndices <- phenoValues.sorted[[2]][1:nProgeny]
					
				}else if(selectionOnSimulated==TRUE){ 
				   	selectedIndices <- phenoSimValues.sorted[[2]][1:nProgeny]
				}
			}

						
				genoValues.selFamily <- genoValues.sorted[[1]][selectedIndices]
				phenoValues.selFamily <- phenoValues.sorted[[1]][selectedIndices]
				genoSimValues.selFamily <- genoSimValues.sorted[[1]][selectedIndices]
				phenoSimValues.selFamily <- phenoSimValues.sorted[[1]][selectedIndices]

				genoValues_List[[nFamily]] <- genoValues.selFamily
				phenoValues_List[[nFamily]] <- phenoValues.selFamily
				genoValSimValues_List[[nFamily]] <- genoSimValues.selFamily
				phenoValSimValues_List[[nFamily]] <- phenoSimValues.selFamily

				initIndex <- finalIndex+1
				finalIndex <- finalIndex + (nProgeny*nCrosses_inFamily)
	    }
	
	}else if(BD =="BD2"){ 
## BD2 

        cycle_nProgeny <- nProgeny*nCrosses_inFamily
	    initIndex <- 1
	    finalIndex <- cycle_nProgeny
	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()

	    for(nFamily in 1:nFamilies){

			genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
			phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
			genoValSimValues_List[[nFamily]] <- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues[initIndex:finalIndex]
			initIndex <- finalIndex+1
			finalIndex <- finalIndex+cycle_nProgeny

	    }

    }else if(BD =="RM"){ 
## RM

        cycle_nProgeny <- nProgeny*nCrosses_inFamily
	    initIndex <- 1
	    finalIndex <- cycle_nProgeny
	    genoValues_List <- list()
	    phenoValues_List <- list()
	    genoValSimValues_List <- list()
	    phenoValSimValues_List <- list()

	    for(nFamily in 1:nFamilies){

			genoValues_List[[nFamily]] <- (genoValues)[initIndex:finalIndex]
			phenoValues_List[[nFamily]] <- (phenoValues)[initIndex:finalIndex]
			genoValSimValues_List[[nFamily]] <- genoValSimValues[initIndex:finalIndex]
			phenoValSimValues_List[[nFamily]] <- phenoValSimValues[initIndex:finalIndex]
			initIndex <- finalIndex+1
			finalIndex <- finalIndex+cycle_nProgeny

	    }
	}
    
	return(list(genoValues_List,phenoValues_List,genoValSimValues_List,phenoValSimValues_List))

}


####################################################################
######### Read in allele conversion table to get genotype table in Nt format
######### PG Table Format 

# snpAlleleFileName <- "SNP_allel_for_IA3023_and_other_genotypes_for_6k_chip.csv" 
# soyNAMParentsFileName <- "SoyNAM parents 4312 SNP genotypes Wm82.a2.csv"

# AlleleConversionTable_Combined <- getAlleleConversionTable(snpAlleleFileName,soyNAMParentsFileName)

getAlleleConversionTable <- function(snpAlleleFileName,soyNAMParentsFileName){


	AlleleConversionTable <- read.csv(snpAlleleFileName,header=TRUE)
	AlleleConversionTable_Form <- AlleleConversionTable[1:5305,]
	
	soyNAMParents <- read.csv(soyNAMParentsFileName,header=TRUE)
	soyNAM_ParentInfo <- soyNAMParents[,1:10]
	cnamesSoyNAM <- colnames(soyNAM_ParentInfo)
	cnamesSoyNAM[3] <- "Locus"
	colnames(soyNAM_ParentInfo) <- cnamesSoyNAM

########

	 AlleleConversionTable_Sorted <- c()

	 for(i in 1:20){
		if(i<=9){
		  lab <- paste("Gm0",i,sep="")
		}else if (i>9){ lab <- paste("Gm",i,sep="")}
		  a<- AlleleConversionTable_Form[as.character(AlleleConversionTable_Form[,2])==lab,3]
		  b<- sort(AlleleConversionTable_Form[as.character(AlleleConversionTable_Form[,2])==lab,3])
		  AlleleConversionTable_Form_Chr <-  AlleleConversionTable_Form[as.character(AlleleConversionTable_Form[,2])==lab,]
		  AlleleConversionTable_Form_ChrMod<- AlleleConversionTable_Form_Chr[match(b,a),]
		  AlleleConversionTable_Sorted<- rbind(AlleleConversionTable_Sorted,AlleleConversionTable_Form_ChrMod)
	 }

####
	  trial<- merge(AlleleConversionTable_Sorted,soyNAM_ParentInfo,by="Locus",sort=FALSE)
	  trial2<- merge(trial,NAM_LinkMap_New,by="Position",sort=FALSE)
	  dim(trial2)

	  HetColumn<-paste(trial2[,"Variant.A"],trial2[,"Variant.B"],sep="")

	  names(HetColumn)<- "Het.Variant"
	  AlleleConversionTable_Combined<- cbind(trial2,HetColumn)
	  
	  return(AlleleConversionTable_Combined)

  }
  
##############
 # populations <- rep(0,20*100)
 # init <- 1
 # for(i in 1:20){

	# final <- i*100
	# populations[init:final]<- rep(i,100)

	# init <- final+1

 # }

 # nextGenGenoTable <- generateMinus101GenoFormat_V1(cycle1_Geno_Data,nFamilies,cycle1_nProgeny)
 # nextGenGenoTable_AlleleFormat <- getAlleleFormat(nextGenGenoTable,AlleleConversionTable_Combined)

 
 
 	
### GM method for nextGen Geno data set
#for(nFamilyPairs in 1:10){
         
          # family1<- familyPairs[[nFamilyPairs]][1]
	      # family2<- familyPairs[[nFamilyPairs]][2]
		  # rNames_Family1 <- c() 
		  # rNames_Family2 <- c()
		
		  # for(famPairs in 1:nrow(nextGenGenoTable_GM_List[[family1]])){ 
		
		    # rNames_Family1 <-c(rNames_Family1,paste(family1," Ind",famPairs))
		    # rNames_Family2 <-c(rNames_Family2,paste(family2," Ind",famPairs))
		  # }
		
		 # rownames(nextGenGenoTable_GM_List[[family1]]) <- rNames_Family1
		 # rownames(nextGenGenoTable_GM_List[[family2]]) <- rNames_Family2
		 # M <- rbind(nextGenGenoTable_GM_List[[family1]],nextGenGenoTable_GM_List[[family2]]) 
	

#GM_SelGeno_List <-  getGM_selectedGenoList(nextGenGenoTable_GM_List,nextGenGenoTableMod_GM_List,PredictionModel_Geno,PredictionModel_Pheno,selectionOnGeno)



# ####
# Selected_GenoData_List <- selectedGenoData_List 
# Selected_PG_table_List <- genoSelectionTable
# MigrationSize <- 1
# EmigrantGroups <- emigrantGroups
# ImmigrantGroups <- immigrantGroups
# Direction <-1

# exchgTrial <- exchangeGenoData(Selected_GenoData_List,Selected_PG_table_List,MigrationSize,EmigrantGroups,ImmigrantGroups,Direction,Policy)
# Test for correct exchange
# for(i in 1:10){ 

  # a5<- (selected_Geno_table_List[[outFamily]][i,,])
  # a5GS <- (selectedGenoData_List[[outFamily]][i,,])
  # print(summary(a5-a5GS))

# }

# for(i in 1:10){ 

  # a5<- (selected_Geno_table_List[[inFamily]][i,,])
  # a5GS <- (selectedGenoData_List[[inFamily]][i,,])
  # print(summary(a5-a5GS))

# }
 
  
# selected_PG_table_List[[inFamily]][finalIndex, 1] :
####################
## Get percent heterozygous from genotype table (1,0,-1)
# FavorableAlleleFreq[[nFamily]][[nCyc]] <- getFreq(nextGenGenoTable_List[[nFamily]])
			# percentHeterozygous_List[[nFamily]][nCyc] <- getPercentHeterozygous(nextGenGenoTable_List[[nFamily]])

			# FavorableQTLAlleleFreq[[nFamily]][[nCyc]] <- getFreq(QTL_GenoTable_List[[nFamily]])
			# percentHeterozygous_QTL_List[[nFamily]][nCyc]

## Fn 47: 
################################################################
		
# ####  Exchange size list among family pairs	
		
		# parent1 <- rep(0,length(selectedParents_Geno)/2)
        # parent2 <- rep(0,length(selectedParents_Geno)/2)
	   
        # nParPairs <-1 
	   
	    # while(nParPairs<length(selectedParents_Geno)){ 
			# parent1[nParPairs] <-selectedParents_Geno[nParPairs]
			# parent2[nParPairs] <- selectedParents_Geno[nParPairs+1]
			# nParPairs <- nParPairs+2
		# }
	   
	    # parent1 <- parent1[seq(1,length(selectedParents_Geno),by=2)]
	    # parent2 <- parent2[seq(1,length(selectedParents_Geno),by=2)] 
	   
	    # exchSize <-0
	    # for(exchNum in 1:(length(parent1)-1)){if((parent1[exchNum]<=100 && parent2[exchNum] >=100)||(parent1[exchNum] >=100 && parent2[exchNum]<=100)){exchSize <- exchSize+1}}
	   
	    # exchSize_List[[nFamilyPairs]][[i]] <- exchSize
			
  
 
 
 
 
 
 
 
 
 
 
 
 
 
