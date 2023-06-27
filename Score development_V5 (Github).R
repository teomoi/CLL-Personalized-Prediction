
################
# Load the data.

getwd()

###################
# Treated Patients.

dataP1<-read.csv("P1.csv",row.names=1)
dim(dataP1)#29 19 #Just to check the dimensions
#
dataP2<-read.csv("P2.csv",row.names=1)
dim(dataP2)#21 14
#
dataP3<-read.csv("P3.csv",row.names=1)
dim(dataP3)#19 20
#
dataP4<-read.csv("P4.csv",row.names=1)
dim(dataP4)#24 16
#
dataP5<-read.csv("P5.csv",row.names=1)
dim(dataP5)#23 20
#
dataP6<-read.csv("P6.csv",row.names=1)
dim(dataP6)#21 21
#
dataP7<-read.csv("P7.csv",row.names=1)
dim(dataP7)#22 16
#
dataP8<-read.csv("P8.csv",row.names=1)
dim(dataP8)#36 20
#
dataP9<-read.csv("P9.csv",row.names=1)
dim(dataP9)#17 20
#
dataP10<-read.csv("P10.csv",row.names=1)
dim(dataP10)#18 15
#
dataP11<-read.csv("P11.csv",row.names=1)
dim(dataP11)#35 20
#
dataP12<-read.csv("P12.csv",row.names=1)
dim(dataP12)#17 19
#
dataP13<-read.csv("P13.csv",row.names=1)
dim(dataP13)#25 19
#
dataP14<-read.csv("P14.csv",row.names=1)
dim(dataP14)#18 21
########################


###################
# Control Patients.

dataC1<-read.csv("C1.csv",row.names=1)
dim(dataC1)#62 20
#
dataC2<-read.csv("C2.csv",row.names=1)
dim(dataC2)#48 20
#
dataC3<-read.csv("C3.csv",row.names=1)
dim(dataC3)#74 18
#
dataC4<-read.csv("C4.csv",row.names=1)
dim(dataC4)#50 20
#
dataC5<-read.csv("C5.csv",row.names=1)
dim(dataC5)#33 17
#
dataC6<-read.csv("C6.csv",row.names=1)
dim(dataC6)#32 18
########################


########################
# Assign to "datadf" the dataframe of each Patient, in order
# to compute the corresponding Score TS.
# This should be separately performed for each patient.

# E.g. P5.
datadf<-dataP5
dim(datadf)
names(datadf)
########################


#######################
# The function C_{j,F_i}(t), when the positive change is of negative nature. 

previous3posthresholdV4<-function(variable, L, U, percentage, percentile, percentileweight){
range<- (U-L)
threshold<-range*percentage
Q<-quantile(seq(L,U,by=0.1), probs = percentile)
Vprevious3<-c(NA)
#1st value
if( variable[2]>=variable[1] & variable[2]>=Q ) {
Vprevious3[1]<- min( (variable[2]-variable[1])/threshold, 1) * percentileweight
} else if ( variable[2]>=variable[1] & variable[2]<Q ) {
Vprevious3[1]<- min( (variable[2]-variable[1])/threshold, 1)
} else {
Vprevious3[1]<- 0
}
#2nd value
if( variable[3]>=mean(variable[1:2]) & variable[3]>=Q ) {
Vprevious3[2]<- min( (variable[3]-mean(variable[1:2]))/threshold, 1) * percentileweight
} else if ( variable[3]>=mean(variable[1:2]) & variable[3]<Q ) {
Vprevious3[2]<- min( (variable[3]-mean(variable[1:2]))/threshold, 1)
} else {
Vprevious3[2]<- 0
}
#i-th value
for(i in 3:(dim(datadf)[1]-1)){
if( variable[i+1]>=mean(variable[(i-2):(i)]) & variable[i+1]>=Q ) {
Vprevious3[i]<- min( (variable[i+1]-mean(variable[(i-2):(i)]))/threshold, 1) * percentileweight
} else if ( variable[i+1]>=mean(variable[(i-2):(i)]) & variable[i+1]<Q ) {
Vprevious3[i]<- min( (variable[i+1]-mean(variable[(i-2):(i)]))/threshold, 1)
} else {
Vprevious3[i]<- 0
}
}
return(Vprevious3)
}
#######################


#######################
# The function C_{j,F_i}(t), when the negative change is of negative nature.

previous3negthresholdV4<-function(variable, L, U, percentage, percentile, percentileweight){
range<- (U-L)
threshold<-range*percentage
Q<-quantile(seq(L,U,by=0.1), probs = percentile)
Vprevious3<-c(NA)
#1st value
if( variable[2]<=variable[1] & variable[2]<=Q ) {
Vprevious3[1]<- min( abs(variable[2]-variable[1])/threshold, 1) * percentileweight
} else if ( variable[2]<=variable[1] & variable[2]>Q ) {
Vprevious3[1]<- min( abs(variable[2]-variable[1])/threshold, 1)
} else {
Vprevious3[1]<- 0
}
#2nd value
if( variable[3]<=mean(variable[1:2]) & variable[3]<=Q ) {
Vprevious3[2]<- min( abs(variable[3]-mean(variable[1:2]))/threshold, 1) * percentileweight
} else if ( variable[3]<=mean(variable[1:2]) & variable[3]>Q ) {
Vprevious3[2]<- min( abs(variable[3]-mean(variable[1:2]))/threshold, 1)
} else {
Vprevious3[2]<- 0
}
#i-th value
for(i in 3:(dim(datadf)[1]-1)){
if( variable[i+1]<=mean(variable[(i-2):(i)]) & variable[i+1]<=Q ) {
Vprevious3[i]<- min( abs(variable[i+1]-mean(variable[(i-2):(i)]))/threshold, 1) * percentileweight
} else if ( variable[i+1]<=mean(variable[(i-2):(i)]) & variable[i+1]>Q ) {
Vprevious3[i]<- min( abs(variable[i+1]-mean(variable[(i-2):(i)]))/threshold, 1)
} else {
Vprevious3[i]<- 0
}
}
return(Vprevious3)
}
#######################



#######################################################################
# The function that will produce the Score for each patient.
# In case one, or more of the features are not included in the dataframe
# of a patient, the corresponding lines in the Score function below
# should be excluded using "#".

Score<-function(percentage, Lpercentile, Upercentile, percentileweight){
# WBC (limits 4000-11000).
WBCMilestone<-c(NA)
for(i in 1:dim(datadf)[1]){
if ( datadf$WBC[i] < 50000) {
WBCMilestone[i]=0
} else if ( datadf$WBC[i] >= 50000 & datadf$WBC[i] < 75000) {
WBCMilestone[i]=1
} else if ( datadf$WBC[i] >= 75000 & datadf$WBC[i] < 100000) {
WBCMilestone[i]=2
} else {
WBCMilestone[i]=3
}
}
#plot(WBCMilestone)
# Lymph (1100-4000).
LymphRatioToDiagnosis<-c(NA)
for(i in 2:dim(datadf)[1]){
if (datadf$Lymph[i]<2*datadf$Lymph[1]) {
LymphRatioToDiagnosis[i-1]=0
} else if ( datadf$Lymph[i]>=2*datadf$Lymph[1] & datadf$Lymph[i]<3*datadf$Lymph[1] ) {
LymphRatioToDiagnosis[i-1]=1
} else if ( datadf$Lymph[i]>=3*datadf$Lymph[1] & datadf$Lymph[i]<4*datadf$Lymph[1] ) {
LymphRatioToDiagnosis[i-1]=2
} else if ( datadf$Lymph[i]>=4*datadf$Lymph[1] & datadf$Lymph[i]<5*datadf$Lymph[1] ) {
LymphRatioToDiagnosis[i-1]=3
} else if ( datadf$Lymph[i]>=5*datadf$Lymph[1] & datadf$Lymph[i]<6*datadf$Lymph[1] ) {
LymphRatioToDiagnosis[i-1]=4
} else {
LymphRatioToDiagnosis[i-1]=5
}
}
#plot(LymphRatioToDiagnosis)
Lymphprevious3<-previous3posthresholdV4(datadf$Lymph,1100,4000,percentage+5, Upercentile, percentileweight)
#plot(Lymphprevious3)
LymphPerc<-round(datadf$Lymph/datadf$WBC ,2)
LymphPerc<-LymphPerc*2
#plot(LymphPerc)
# RBC (4500-5900).
# RBCprevious3<-previous3negthresholdV4(datadf$RBC,4500,5900,percentage, Lpercentile, percentileweight)
#plot(RBCprevious3)
# Hb (limits 12-17).
HbMilestone<-c(NA)
for(i in 1:dim(datadf)[1]){
if ( datadf$Hb[i] > 12) {
HbMilestone[i]=0
} else if ( datadf$Hb[i] <= 12 & datadf$Hb[i] > 11) {
HbMilestone[i]=1
} else if ( datadf$Hb[i] <= 11 & datadf$Hb[i] > 10.5) {
HbMilestone[i]=2
} else if ( datadf$Hb[i] <= 10.5 & datadf$Hb[i] > 10) {
HbMilestone[i]=3
} else {
HbMilestone[i]=4
}
}
#plot(HbMilestone)
Hbprevious3<-previous3negthresholdV4(datadf$Hb,12,17,percentage, Lpercentile, percentileweight)
#plot(Hbprevious3)
# PLT (limits 150000-400000).
PLTRatioToDiagnosis<-c(NA)
for(i in 2:dim(datadf)[1]){
if (datadf$PLT[i]>0.667*datadf$PLT[1]) {
PLTRatioToDiagnosis[i-1]=0
} else if ( datadf$PLT[i]<0.667*datadf$PLT[1] & datadf$PLT[i]>0.5*datadf$PLT[1] ) {
PLTRatioToDiagnosis[i-1]=1
} else {
PLTRatioToDiagnosis[i-1]=2
}
}
#plot(PLTRatioToDiagnosis)
PLTMilestone<-c(NA)
for(i in 1:dim(datadf)[1]){
if ( datadf$PLT[i] > 200000) {
PLTMilestone[i]=0
} else if ( datadf$PLT[i] <= 200000 & datadf$PLT[i] > 150000) {
PLTMilestone[i]=1
} else if ( datadf$PLT[i] <= 150000 & datadf$PLT[i] > 100000) {
PLTMilestone[i]=2
} else {
PLTMilestone[i]=3
}
}
#plot(PLTMilestone)
PLTprevious3<-previous3negthresholdV4(datadf$PLT,150000,400000,percentage, Lpercentile, percentileweight)
#plot(PLTprevious3)
# ESR (limits 0-20).
ESRMilestone	<-ifelse(datadf$ESR>=100,1,0)
#plot(ESRMilestone)
ESRprevious3<-previous3posthresholdV4(datadf$ESR, 0, 20,percentage, Upercentile, percentileweight)
#plot(ESRprevious3)
# Ur (limits 20-50).
Urprevious3<-previous3posthresholdV4(datadf$Ur,20,50,percentage, Upercentile, percentileweight)
#plot(Urprevious3)
# Cr (0.8-1.4).
CrMilestone	<-ifelse(datadf$Cr>=2,1,0)
#plot(CrMilestone)
Crprevious3<-previous3posthresholdV4(datadf$Cr,0.8,1.4,percentage, Upercentile, percentileweight)
#plot(Crprevious3)
# Glu (70-105).
Gluprevious3<-previous3posthresholdV4(datadf$Glu,70,105,percentage, Upercentile, percentileweight)
#plot(Gluprevious3)
# ALP (25-125).
ALPprevious3<-previous3posthresholdV4(datadf$ALP,25,125,percentage, Upercentile, percentileweight)
#plot(ALPprevious3)
# SGOT (limits 10-40).
SGOTprevious3<-previous3posthresholdV4(datadf$SGOT,10,40,percentage, Upercentile, percentileweight)
#plot(SGOTprevious3)
# SGPT (limits 10-40).
SGPTprevious3<-previous3posthresholdV4(datadf$SGPT,10,40,percentage, Upercentile, percentileweight)
#plot(SGPTprevious3)
# CRP (0-1).
CRPMilestone	<-ifelse(datadf$CRP>=6,1,0)
#plot(CRPMilestone)
CRPprevious3<-previous3posthresholdV4(datadf$CRP,0,1,percentage, Upercentile, percentileweight)
#plot(CRPprevious3)
# gGT (limits 10-49).
gGTprevious3<-previous3posthresholdV4(datadf$gGT,10,49,percentage, Upercentile, percentileweight)
#plot(gGTprevious3)
# UA (limits 3.6-8.4).
UAMilestone	<-ifelse(datadf$UA>=10,1,0)
#plot(UAMilestone)
UAprevious3<-previous3posthresholdV4(datadf$UA,3.6,8.4,percentage, Upercentile, percentileweight)
#plot(UAprevious3)
# LDH (limits 140-280).
LDHMilestone	<-ifelse(datadf$LDH>=246,1,0)
#plot(LDHMilestone)
LDHprevious3<-previous3posthresholdV4(datadf$LDH,140,280,percentage, Upercentile, percentileweight)
#plot(LDHprevious3)
# Ca (limits 8.5-10.5).
Caprevious3<-previous3posthresholdV4(datadf$Ca,8.5,10.5,percentage, Upercentile, percentileweight)
#plot(Caprevious3)
# TBil (0.3-1.2).
TBilprevious3<-previous3posthresholdV4(datadf$TBil,0.3,1.2,percentage, Upercentile, percentileweight)
#plot(TBilprevious3)
# alb (limits 3.5-5.2).
albMilestone	<-ifelse(datadf$alb<=3.5,1,0)
#plot(albMilestone)
albprevious3<-previous3negthresholdV4(datadf$alb,3.5,5.2,percentage, Lpercentile, percentileweight)
#plot(albprevious3)
# Glob (limits 2.7-3.3).
Globprevious3<-previous3negthresholdV4(datadf$Glob,2.7,3.3,percentage, Lpercentile, percentileweight)
#plot(Globprevious3)
SpleenSize<-datadf$SpleenSize/5
#plot(SpleenSize)

Sum <-WBCMilestone[-1]+
	LymphRatioToDiagnosis+Lymphprevious3+LymphPerc[-1]+
	# RBCprevious3 +
	HbMilestone[-1]+Hbprevious3+
	PLTRatioToDiagnosis+PLTMilestone[-1]+PLTprevious3+
  ESRprevious3 + ESRMilestone[-1] + Urprevious3 + 
  SpleenSize[-1] +
  CrMilestone[-1] + Crprevious3+
  Gluprevious3+
	ALPprevious3 + SGOTprevious3 + SGPTprevious3+
  CRPMilestone[-1] + CRPprevious3+
	gGTprevious3+ 
	UAMilestone[-1] + UAprevious3+
	LDHMilestone[-1]+ LDHprevious3+
	Caprevious3 + TBilprevious3 +
	albMilestone[-1] + albprevious3 + 
	Globprevious3 

plot(Sum, type="o", main="Total Score", ylab="Score",
	xlab="Follow-up visits until FT", ylim=c(0, max(Sum)+1),lwd=1.2)
lines(lowess(Sum, f=0.3),col=11,lwd=1.2)
return(Sum)
}

####################################################################################
# By setting percentage=0.2; Lpercentile=0.3; Upercentile=0.7; percentileweight=2,
# then we have the Score TS for the individual patient, represented by "Sum"

P0.2LP0.3UP0.7PW2 <-Score(percentage=0.2, Lpercentile=0.3, Upercentile=0.7, percentileweight=2)
Sum <-P0.2LP0.3UP0.7PW2

plot(Sum, type="o", main="P5 with a Smooth curve", ylab="Score",
	xlab="Follow-up visits until FT", ylim=c(0, max(Sum)+1))
lines(lowess(Sum, f=0.3),lwd=2,col=11)

# For each patient, the corresponding "Sum" should be computed, e.g. SumP1, SumP2,
# SumP3, ..., SumP14, SumC1, ..., SumC6.


############################################################################
# After computing the Sum, based on the Score function, for each patient, 
# the Score values for each patient should be adjusted for the number of 
# components the corresponding Score TS included. E.g. P4, did not have
# Glucose values (maximum contribution to the Score=2), thus, the Score 
# should be multiplied by (64/62), where 64 is the maximum value of the Score, 
# assumed on the maximum values of its components, and with the weight set to 2,
# and 62 is the actual maximum value of the Score that could have been
# received in the example of this patient.
# Then, for each individual treated patient, the adjusted Score, represented
# by, lets say, e.g., SumP1adj, SumP2adj, SumP3adj, ..., SumP14adj, should be 
# used in the function "AssessPrediction5", to define the corresponding threshold.
# For control patients, no threshold is computed, however, the adjustment is
# performed similarly as in the case of treated patients, resulting in the 
# adjusted Score, lets say, e.g., SumC1adj, ..., SumC6adj.
# Please note, that the adjustment will not necessarily be performed in all
# patients, however, for uniformity, we have used the addition of "adj" in the
# notation corresponding to all patients.








