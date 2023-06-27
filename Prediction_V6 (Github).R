
#########################################################
#install.packages("dtw")
library(dtw)

#########################################################
# The function to measure the distance between 2 patients.

PlotAlign<-function(Patient1, Patient2, steppattern){
align1<-dtw(Patient1,Patient2,step=steppattern,keep=TRUE)
align2<-dtw(Patient1,Patient2,step=steppattern,keep=TRUE,
		window.type = "sakoechiba", window.size=1)
A1<-align1$normalizedDistance
A2<-align2$normalizedDistance
return("AlignValues"=c(A1,A2))
} 

#########################################################
# Define the list of patients' Score TSs.

Patients<-list("P1"=SumP1adj, "P2"=SumP2adj, ..., "P14"=SumP14adj,
               "C1"=SumC1adj, "C2"=SumC2adj, ..., "C6"=SumC6adj)

length(Patients)	
names(Patients)	

#########################################################################
# Develop a matrix with the distances across the list of patients.
# The distances refer to the n time point at which the patients
# will be compared to each other.

DistanceMatrix<-function(steppattern, maindiagonal, n, PatientList){
NDistanceM<-matrix(NA,nrow=length(PatientList),ncol=length(PatientList))
# Complete the values of the lower triangular part of the matrix.
for(i in 2:length(PatientList) ){
for(j in 1: (i-1) ){
#Set the distance we prefer at the end (1 or 2)
NDistanceM[i,j]<-PlotAlign(PatientList[[i]][1:n],PatientList[[j]][1:n],steppattern)[2]
}
} ;NDistanceM

# Complete the values of the upper triangular part of the matrix.
for(i in 1:(length(PatientList)-1)){
for(j in (i+1): length(PatientList)){
#Set the distance we prefer at the end (1 or 2)
NDistanceM[i,j]<-PlotAlign(PatientList[[i]][1:n],PatientList[[j]][1:n],steppattern)[2]
}
};NDistanceM

# Main diagonal
for(i in 1:length(PatientList)){
NDistanceM[i,i]<-maindiagonal	#to easily exclude this distance below
};NDistanceM
rownames(NDistanceM)<-colnames(NDistanceM)<-names(PatientList)
return(NDistanceM)
}

# For example
NDistanceM<- DistanceMatrix(symmetric2,10,13,Patients);round(NDistanceM,2)

##################################################################################

Thresholds<-data.frame("PIDs"=names(Patients), 
                       "Thresholds"=(vector stating the binary threshold obtained 
                                     for each treated patient - for Controls
                                     use 100, in any case it is irrelevant for the analysis),
                       "IGHV"=(vector stating the IGHV status for each patient))

##################################################################################

# Define the vectors containing the treated patients IDs, and control IDs.
TreatedPNames<-names(Patients)[1:14];TreatedPNames
ControlPNames<-names(Patients)[15:20];ControlPNames


##################################################################################
# Develop the functions to find the group of similar patients.
# Three versions, corresponding to the three assessment criteria, 
# with two variations each, are presented below.

##############
# 1st version.
FindSimilarTreatedP<-function(steppattern, maindiagonal, n, PatientList, DistanceT, PatientID){
  PatientList<-PatientList[lengths(PatientList)>=n]
  NDistanceM<- DistanceMatrix(steppattern, maindiagonal, n, PatientList)
  print(dim(NDistanceM)[1])
  selectedRow<-NDistanceM[rownames(NDistanceM)==PatientID,];selectedRow
  ThresholdSort<-       sort(selectedRow[selectedRow<=DistanceT]);ThresholdSort
  GroupNamesSort<-names(sort(selectedRow[selectedRow<=DistanceT]));GroupNamesSort
  if ( length(GroupNamesSort)>=4 & sum(GroupNamesSort[1:4] %in% TreatedPNames)/4>=0.75 ) {
    GroupTNames<-GroupNamesSort[1:4][GroupNamesSort[1:4] %in% TreatedPNames]
  } else if ( length(GroupNamesSort)<4 & sum(GroupNamesSort %in% TreatedPNames)/length(GroupNamesSort)==1) {
    GroupTNames<-GroupNamesSort[GroupNamesSort %in% TreatedPNames]
  } else {
    GroupTNames<-c("Dont predict");print("Dont predict");print(ThresholdSort)
  }
  SimilarTreatedP<-ThresholdSort[GroupTNames]
  return(SimilarTreatedP)
}

#####################################
# 1st version (only for Mutated CLL).
Thresholds
MutatedPNames<-Thresholds$PIDs[c(2:6,10,11,14:19)];MutatedPNames

FindSimilarTreatedPM<-function(steppattern, maindiagonal, n, PatientList, DistanceT, PatientID){
  PatientList<-PatientList[lengths(PatientList)>=n]
  NDistanceM<- DistanceMatrix(steppattern, maindiagonal, n, PatientList)
  print(dim(NDistanceM)[1])
  selectedRow<-NDistanceM[rownames(NDistanceM)==PatientID,];selectedRow
  ThresholdSort<-       sort(selectedRow[selectedRow<=DistanceT]);ThresholdSort
  GroupNamesSort<-names(sort(selectedRow[selectedRow<=DistanceT]));GroupNamesSort
  GroupNamesSort<-GroupNamesSort[GroupNamesSort %in% MutatedPNames];GroupNamesSort
  print(GroupNamesSort)
  if ( length(GroupNamesSort)>=4 & sum(GroupNamesSort[1:4] %in% TreatedPNames)/4>=0.75 ) {
    GroupTNames<-GroupNamesSort[1:4][GroupNamesSort[1:4] %in% TreatedPNames]
  } else if ( length(GroupNamesSort)<4 & sum(GroupNamesSort %in% TreatedPNames)/length(GroupNamesSort)==1) {
    GroupTNames<-GroupNamesSort[GroupNamesSort %in% TreatedPNames]
  } else {
    GroupTNames<-c("Dont predict");print("Dont predict");print(ThresholdSort)
  }
  SimilarTreatedP<-ThresholdSort[GroupTNames]
  return(SimilarTreatedP)
}

##############
# 2nd version.
FindSimilarTreatedP<-function(steppattern, maindiagonal, n, PatientList, DistanceT, PatientID){
  PatientList<-PatientList[lengths(PatientList)>=n]
  NDistanceM<- DistanceMatrix(steppattern, maindiagonal, n, PatientList)
  print(dim(NDistanceM)[1])
  selectedRow<-NDistanceM[rownames(NDistanceM)==PatientID,];selectedRow
  ThresholdSort<-       sort(selectedRow[selectedRow<=DistanceT]);ThresholdSort
  GroupNamesSort<-names(sort(selectedRow[selectedRow<=DistanceT]));GroupNamesSort
  print(GroupNamesSort)
  if ( length(GroupNamesSort)>=6 & sum(GroupNamesSort[1:6] %in% TreatedPNames)/6>=0.65 ) {
    GroupTNames<-GroupNamesSort[1:6][GroupNamesSort[1:6] %in% TreatedPNames]
    GroupTNames<-GroupTNames[1:2] ; GroupTNames<-GroupTNames[!is.na(GroupTNames)]
  } else if ( length(GroupNamesSort)<6 & (sum(GroupNamesSort %in% TreatedPNames)/length(GroupNamesSort)>=0.65)) {
    GroupTNames<-GroupNamesSort[GroupNamesSort %in% TreatedPNames]
    GroupTNames<-GroupTNames[1:2] ; GroupTNames<-GroupTNames[!is.na(GroupTNames)]
  } else {
    GroupTNames<-c("Dont predict");print("Dont predict");print(ThresholdSort)
  }
  SimilarTreatedP<-ThresholdSort[GroupTNames]
  return(SimilarTreatedP)
}

#####################################
# 2nd version (only for Mutated CLL).
FindSimilarTreatedPM<-function(steppattern, maindiagonal, n, PatientList, DistanceT, PatientID){
  PatientList<-PatientList[lengths(PatientList)>=n]
  NDistanceM<- DistanceMatrix(steppattern, maindiagonal, n, PatientList)
  print(dim(NDistanceM)[1])
  selectedRow<-NDistanceM[rownames(NDistanceM)==PatientID,];selectedRow
  ThresholdSort<-       sort(selectedRow[selectedRow<=DistanceT]);ThresholdSort
  GroupNamesSort<-names(sort(selectedRow[selectedRow<=DistanceT]));GroupNamesSort
  GroupNamesSort<-GroupNamesSort[GroupNamesSort %in% MutatedPNames];GroupNamesSort
  print(GroupNamesSort)
  if ( length(GroupNamesSort)>=6 & sum(GroupNamesSort[1:6] %in% TreatedPNames)/6>=0.65 ) {
    GroupTNames<-GroupNamesSort[1:6][GroupNamesSort[1:6] %in% TreatedPNames]
    GroupTNames<-GroupTNames[1:2] ; GroupTNames<-GroupTNames[!is.na(GroupTNames)]
  } else if ( length(GroupNamesSort)<6 & (sum(GroupNamesSort %in% TreatedPNames)/length(GroupNamesSort)>=0.65)) {
    GroupTNames<-GroupNamesSort[GroupNamesSort %in% TreatedPNames]
    GroupTNames<-GroupTNames[1:2] ; GroupTNames<-GroupTNames[!is.na(GroupTNames)]
  } else {
    GroupTNames<-c("Dont predict");print("Dont predict");print(ThresholdSort)
  }
  SimilarTreatedP<-ThresholdSort[GroupTNames]
  return(SimilarTreatedP)
}

##############
# 3rd version.
FindSimilarTreatedP<-function(steppattern, maindiagonal, n, PatientList, DistanceT, PatientID){
  PatientList<-PatientList[lengths(PatientList)>=n]
  NDistanceM<- DistanceMatrix(steppattern, maindiagonal, n, PatientList)
  print(dim(NDistanceM)[1])
  selectedRow<-NDistanceM[rownames(NDistanceM)==PatientID,];selectedRow
  ThresholdSort<-       sort(selectedRow[selectedRow<=DistanceT]);ThresholdSort
  GroupNamesSort<-names(sort(selectedRow[selectedRow<=DistanceT]));GroupNamesSort
  print(GroupNamesSort)
  if ( length(GroupNamesSort)>=6 & sum(GroupNamesSort[1:6] %in% TreatedPNames)/6>=0.5 ) {
    GroupTNames<-GroupNamesSort[1:6][GroupNamesSort[1:6] %in% TreatedPNames]
    GroupTNames<-GroupTNames[1:3] ; GroupTNames<-GroupTNames[!is.na(GroupTNames)]
  } else if ( length(GroupNamesSort)<6 & (sum(GroupNamesSort %in% TreatedPNames)/length(GroupNamesSort)>=0.5)) {
    GroupTNames<-GroupNamesSort[GroupNamesSort %in% TreatedPNames]
    GroupTNames<-GroupTNames[1:3] ; GroupTNames<-GroupTNames[!is.na(GroupTNames)]
  } else {
    GroupTNames<-c("Dont predict");print("Dont predict");print(ThresholdSort)
  }
  SimilarTreatedP<-ThresholdSort[GroupTNames]
  return(SimilarTreatedP)
}

#####################################
# 3rd version (only for Mutated CLL).
FindSimilarTreatedPM<-function(steppattern, maindiagonal, n, PatientList, DistanceT, PatientID){
  PatientList<-PatientList[lengths(PatientList)>=n]
  NDistanceM<- DistanceMatrix(steppattern, maindiagonal, n, PatientList)
  print(dim(NDistanceM)[1])
  selectedRow<-NDistanceM[rownames(NDistanceM)==PatientID,];selectedRow
  ThresholdSort<-       sort(selectedRow[selectedRow<=DistanceT]);ThresholdSort
  GroupNamesSort<-names(sort(selectedRow[selectedRow<=DistanceT]));GroupNamesSort
  GroupNamesSort<-GroupNamesSort[GroupNamesSort %in% MutatedPNames];GroupNamesSort
  print(GroupNamesSort)
  if ( length(GroupNamesSort)>=6 & sum(GroupNamesSort[1:6] %in% TreatedPNames)/6>=0.5 ) {
    GroupTNames<-GroupNamesSort[1:6][GroupNamesSort[1:6] %in% TreatedPNames]
    GroupTNames<-GroupTNames[1:3] ; GroupTNames<-GroupTNames[!is.na(GroupTNames)]
  } else if ( length(GroupNamesSort)<6 & (sum(GroupNamesSort %in% TreatedPNames)/length(GroupNamesSort)>=0.5)) {
    GroupTNames<-GroupNamesSort[GroupNamesSort %in% TreatedPNames]
    GroupTNames<-GroupTNames[1:3] ; GroupTNames<-GroupTNames[!is.na(GroupTNames)]
  } else {
    GroupTNames<-c("Dont predict");print("Dont predict");print(ThresholdSort)
  }
  SimilarTreatedP<-ThresholdSort[GroupTNames]
  return(SimilarTreatedP)
}

##################################################################################
# For each “new” patient, separately, the following procedure should be applied.
# Find the treated patients within the group of patients to which the “new” patient 
# is mostly similar to, based on the measured distances. The procedure begins at
# e.g., nstart=13, and continues, if necessary, with nstart=14, 15, etc.
# The "DistanceT" should be selected to be 2.5, 3, or 3.5, depending on the 
# assessment criterion used (corresponding to each of the three function versions).

# Employing the whole reference pool of patients.
PatientsDist<-FindSimilarTreatedP(symmetric2,10,nstart,Patients,DistanceT, "PatientID")
# Employing only mutated patients within the whole reference pool of patients.
PatientsDist<-FindSimilarTreatedPM(symmetric2,10,nstart,Patients,DistanceT, "PatientID")


# Find the mean threshold of the group. Use it for each patient (within the group) 
# to find the time point (after n) at which the predictions/steps for the 
# Binary Score are consecutively 1.

GroupNames<-names(PatientsDist);GroupNames
Thresholds[Thresholds$PIDs %in% GroupNames,]
MeanThreshold<-mean(Thresholds$Thresholds[Thresholds$PIDs %in% GroupNames]);MeanThreshold
ProbThreshold<-0.5

# Find the steps with 1s until 1st treatment
Steps<-c(NA)
for(i in 1:length(GroupNames)){
  print(i)
  PatientRef <-ifelse(Patients[[which(names(Patients) %in% GroupNames)[i]]]>=MeanThreshold,1,0);PatientRef 
  n<-nstart
  print(n)
  while(n <length(PatientRef)){
    trainB<-PatientRef[1:n]  ;trainB
    #Nofeedback models with lag=m.
    #Results via glm().
    size<-length(trainB);size
    M3 <- glm(trainB[4:size]~trainB[3:(size-1)]+trainB[2:(size-2)]+trainB[1:(size-3)], family=binomial);M3
    #M3 coefficients.
    M3d <- M3$coef[1]	;M3d
    M3b1<- M3$coef[2]	;M3b1
    M3b2<- M3$coef[3]	;M3b2
    M3b3<- M3$coef[4]	;M3b3
    
    if ( sum(is.na(M3$coef))>0 ) {
      print("glm did not properly fit")
      n<-n+1
      print(n)
    } else {
      #Prediction.
      M3predprob<-c(NA)
      M3testhat <-c(NA)
      if ( (length(PatientRef)-n)>=1 ) {
        M3predprob[1] <- plogis(M3d+M3b1*trainB[length(trainB)]+M3b2*trainB[length(trainB)-1]+M3b3*trainB[length(trainB)-2])
        M3testhat[1]<-ifelse(M3predprob[1]>=ProbThreshold,1,0)
      } else {
        M3testhat[1]<-NA
      }
      if ( (length(PatientRef)-n)>=2 ) {
        M3predprob[2] <- plogis(M3d+M3b1*M3testhat[1]+M3b2*trainB[length(trainB)]+M3b3*trainB[length(trainB)-1])
        M3testhat[2]<-ifelse(M3predprob[2]>=ProbThreshold,1,0)
      } else {
        M3testhat[2]<-NA
      }
      if ( (length(PatientRef)-n)>=3 ) {
        M3predprob[3] <- plogis(M3d+M3b1*M3testhat[2]+M3b2*M3testhat[1]+M3b3*trainB[length(trainB)])
        M3testhat[3]<-ifelse(M3predprob[3]>=ProbThreshold,1,0)
      } else {
        M3testhat[3]<-NA
      }
      if ( (length(PatientRef)-n)>=4 ) {
        M3predprob[4] <- plogis(M3d+M3b1*M3testhat[3]+M3b2*M3testhat[2]+M3b3*M3testhat[1])
        M3testhat[4]<-ifelse(M3predprob[4]>=ProbThreshold,1,0)
      } else {
        M3testhat[4]<-NA
      }
      if ( (length(PatientRef)-n)>=5 ) {
        M3predprob[5] <- plogis(M3d+M3b1*M3testhat[4]+M3b2*M3testhat[3]+M3b3*M3testhat[2])
        M3testhat[5]<-ifelse(M3predprob[5]>=ProbThreshold,1,0)
      } else {
        M3testhat[5]<-NA
      }
      M3predprob;M3testhat
      print(M3predprob)
      print(M3testhat)
      M3testhat<-M3testhat[!is.na(M3testhat)]
      if ( sum(M3testhat)<length(M3testhat) ) {
        print("Not all 1's")
        n<-n+1
        print(n)
      } else {
        Steps[i]<- length(PatientRef)-n
        break
      }
    }
  }
}
Steps
#  Their mean value.
mean(Steps, na.rm=T)


# Set the specific "PatientNumber" (of the “new” patient) in the list of Patients 
# in order to find, if appropriate, the estimated number of steps until 1st treatment
PatientRef <-ifelse(Patients[[PatientNumber]] >=MeanThreshold,1,0);PatientRef 
n<-nstart;n
trainB<-PatientRef[1:n]  ;trainB  ;length(trainB);size<-length(trainB);size
#Nofeedback models with lag=m.
#Results via glm().
M3 <- glm(trainB[4:size]~trainB[3:(size-1)]+trainB[2:(size-2)]+trainB[1:(size-3)], family=binomial);M3
#M3 coefficients.
M3d <- M3$coef[1]	;M3d;M3b1<- M3$coef[2]	;M3b1;M3b2<- M3$coef[3]	;M3b2;M3b3<- M3$coef[4]	;M3b3

if ( sum(is.na(M3$coef))>0 ) {
print("glm did not properly fit")
} else {
#Prediction.
M3predprob<-c(NA)
M3testhat<-c(NA)
if ( (length(PatientRef)-n)>=1 ) {
M3predprob[1] <- plogis(M3d+M3b1*trainB[length(trainB)]+M3b2*trainB[length(trainB)-1]+M3b3*trainB[length(trainB)-2])
M3testhat[1]<-ifelse(M3predprob[1]>=ProbThreshold,1,0)
} else {
M3testhat[1]<-NA
}
if ( (length(PatientRef)-n)>=2 ) {
M3predprob[2] <- plogis(M3d+M3b1*M3testhat[1]+M3b2*trainB[length(trainB)]+M3b3*trainB[length(trainB)-1])
M3testhat[2]<-ifelse(M3predprob[2]>=ProbThreshold,1,0)
} else {
M3testhat[2]<-NA
}
if ( (length(PatientRef)-n)>=3 ) {
M3predprob[3] <- plogis(M3d+M3b1*M3testhat[2]+M3b2*M3testhat[1]+M3b3*trainB[length(trainB)])
M3testhat[3]<-ifelse(M3predprob[3]>=ProbThreshold,1,0)
} else {
M3testhat[3]<-NA
}
if ( (length(PatientRef)-n)>=4 ) {
M3predprob[4] <- plogis(M3d+M3b1*M3testhat[3]+M3b2*M3testhat[2]+M3b3*M3testhat[1])
M3testhat[4]<-ifelse(M3predprob[4]>=ProbThreshold,1,0)
} else {
M3testhat[4]<-NA
}
if ( (length(PatientRef)-n)>=5 ) {
M3predprob[5] <- plogis(M3d+M3b1*M3testhat[4]+M3b2*M3testhat[3]+M3b3*M3testhat[2])
M3testhat[5]<-ifelse(M3predprob[5]>=ProbThreshold,1,0)
} else {
M3testhat[5]<-NA
}
round(M3predprob,3);round(M3testhat,1)
M3testhat<-M3testhat[!is.na(M3testhat)]
#Steps with 1s until 1st treatment
if ( sum(M3testhat)<length(M3testhat) ) {
print("Not all 1's");print(round(M3testhat,1))
} else {
StepsPN<- length(PatientRef)-n;print(round(M3testhat,1))
}
}

StepsPN    #Observed
mean(Steps)#Expected
mean(Steps, na.rm=T)
n

# In case the procedure does not lead to obtain the "StepsPN", the next value
# for "nstart" should be considered.











