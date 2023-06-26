
###########################################################################
# The function used to estimate a threshold, specific for each treated patient, 
# to transform the continuous Score TS into a binary one, with values 0/1, 
# representing lower/higher Score, respectively.

AssessPrediction5<-function(timeseries,ProbThreshold){
k<-5
train<-timeseries[1:(length(timeseries)-k)]
test <-timeseries[(length(timeseries)-k+1):length(timeseries)]
#Create the arrays to store the results.
M1results <- array(c(NA),dim = c(5,5,5),dimnames = list(c("Step 1","Step 2","Step 3","Step 4","Step 5"),
c("PredProb","Predictions","Observations","Criterion1","Criterion2"), 
c("Model1-Threshold1","Model1-Threshold2","Model1-Threshold3","Model1-Threshold4","Model1-Threshold5") ))
M2results <- array(c(NA),dim = c(5,5,5),dimnames = list(c("Step 1","Step 2","Step 3","Step 4","Step 5"),
c("PredProb","Predictions","Observations","Criterion1","Criterion2"), 
c("Model2-Threshold1","Model2-Threshold2","Model2-Threshold3","Model2-Threshold4","Model2-Threshold5") ))
M3results <- array(c(NA),dim = c(5,5,5),dimnames = list(c("Step 1","Step 2","Step 3","Step 4","Step 5"),
c("PredProb","Predictions","Observations","Criterion1","Criterion2"), 
c("Model3-Threshold1","Model3-Threshold2","Model3-Threshold3","Model3-Threshold4","Model3-Threshold5") ))
Thresholds<-quantile(train, probs = c(0.6,0.65,0.7,0.75,0.8))
for(i in 1:5){
trainB<-ifelse(train>=Thresholds[i],1,0)	;trainB
testB <-ifelse(test >=Thresholds[i],1,0)	;testB 

#No feedback models with lag=m.
#Results via glm().
size<-length(trainB)
M1 <- glm(trainB[2:size]~trainB[1:(size-1)], family=binomial) #m=1.
M2 <- glm(trainB[3:size]~trainB[2:(size-1)]+trainB[1:(size-2)], family=binomial)#m=2.
M3 <- glm(trainB[4:size]~trainB[3:(size-1)]+trainB[2:(size-2)]+trainB[1:(size-3)], family=binomial)#m=3.

#M1 coefficients.
M1d <- M1$coef[1]	;M1d
M1b1<- M1$coef[2]	;M1b1
#M2 coefficients.
M2d <- M2$coef[1]	;M2d
M2b1<- M2$coef[2]	;M2b1
M2b2<- M2$coef[3]	;M2b2
#M3 coefficients.
M3d <- M3$coef[1]	;M3d
M3b1<- M3$coef[2]	;M3b1
M3b2<- M3$coef[3]	;M3b2
M3b3<- M3$coef[4]	;M3b3

#Prediction.
M1predprob<-rep(NA,k)
M2predprob<-rep(NA,k)
M3predprob<-rep(NA,k)
M1testhat<-rep(NA,k)
M2testhat<-rep(NA,k)
M3testhat<-rep(NA,k)
M1Criterion1<-rep(NA,k)
M2Criterion1<-rep(NA,k)
M3Criterion1<-rep(NA,k)
M1Criterion2<-rep(NA,k)
M2Criterion2<-rep(NA,k)
M3Criterion2<-rep(NA,k)

#M1.
M1predprob[1] <- plogis(M1d+M1b1*trainB[length(trainB)])
M1testhat[1]<-ifelse(M1predprob[1]>=ProbThreshold,1,0)
M1predprob[2] <- plogis(M1d+M1b1*M1testhat[1])
M1testhat[2]<-ifelse(M1predprob[2]>=ProbThreshold,1,0)
M1predprob[3] <- plogis(M1d+M1b1*M1testhat[2])
M1testhat[3]<-ifelse(M1predprob[3]>=ProbThreshold,1,0)
M1predprob[4] <- plogis(M1d+M1b1*M1testhat[3])
M1testhat[4]<-ifelse(M1predprob[4]>=ProbThreshold,1,0)
M1predprob[5] <- plogis(M1d+M1b1*M1testhat[4])
M1testhat[5]<-ifelse(M1predprob[5]>=ProbThreshold,1,0)

M1Criterion1[1]<- abs(M1predprob[1]-testB[1])
M1Criterion1[2]<- abs(M1predprob[2]-testB[2])
M1Criterion1[3]<- abs(M1predprob[3]-testB[3])
M1Criterion1[4]<- abs(M1predprob[4]-testB[4])
M1Criterion1[5]<- abs(M1predprob[5]-testB[5])
M1Criterion2[1]<- (M1predprob[1]-testB[1])^2/(M1predprob[1]*(1-M1predprob[1]))
M1Criterion2[2]<- (M1predprob[2]-testB[2])^2/(M1predprob[2]*(1-M1predprob[2]))
M1Criterion2[3]<- (M1predprob[3]-testB[3])^2/(M1predprob[3]*(1-M1predprob[3]))
M1Criterion2[4]<- (M1predprob[4]-testB[4])^2/(M1predprob[4]*(1-M1predprob[4]))
M1Criterion2[5]<- (M1predprob[5]-testB[5])^2/(M1predprob[5]*(1-M1predprob[5]))

#M2.
M2predprob[1] <- plogis(M2d+M2b1*trainB[length(trainB)]+M2b2*trainB[length(trainB)-1])
M2testhat[1]<-ifelse(M2predprob[1]>=ProbThreshold,1,0)
M2predprob[2] <- plogis(M2d+M2b1*M2testhat[1]+M2b2*trainB[length(trainB)])
M2testhat[2]<-ifelse(M2predprob[2]>=ProbThreshold,1,0)
M2predprob[3] <- plogis(M2d+M2b1*M2testhat[2]+M2b2*M2testhat[1])
M2testhat[3]<-ifelse(M2predprob[3]>=ProbThreshold,1,0)
M2predprob[4] <- plogis(M2d+M2b1*M2testhat[3]+M2b2*M2testhat[2])
M2testhat[4]<-ifelse(M2predprob[4]>=ProbThreshold,1,0)
M2predprob[5] <- plogis(M2d+M2b1*M2testhat[4]+M2b2*M2testhat[3])
M2testhat[5]<-ifelse(M2predprob[5]>=ProbThreshold,1,0)


M2Criterion1[1]<- abs(M2predprob[1]-testB[1])
M2Criterion1[2]<- abs(M2predprob[2]-testB[2])
M2Criterion1[3]<- abs(M2predprob[3]-testB[3])
M2Criterion1[4]<- abs(M2predprob[4]-testB[4])
M2Criterion1[5]<- abs(M2predprob[5]-testB[5])
M2Criterion2[1]<- (M2predprob[1]-testB[1])^2/(M2predprob[1]*(1-M2predprob[1]))
M2Criterion2[2]<- (M2predprob[2]-testB[2])^2/(M2predprob[2]*(1-M2predprob[2]))
M2Criterion2[3]<- (M2predprob[3]-testB[3])^2/(M2predprob[3]*(1-M2predprob[3]))
M2Criterion2[4]<- (M2predprob[4]-testB[4])^2/(M2predprob[4]*(1-M2predprob[4]))
M2Criterion2[5]<- (M2predprob[5]-testB[5])^2/(M2predprob[5]*(1-M2predprob[5]))

#M3.
M3predprob[1] <- plogis(M3d+M3b1*trainB[length(trainB)]+M3b2*trainB[length(trainB)-1]+M3b3*trainB[length(trainB)-2])
M3testhat[1]<-ifelse(M3predprob[1]>=ProbThreshold,1,0)
M3predprob[2] <- plogis(M3d+M3b1*M3testhat[1]+M3b2*trainB[length(trainB)]+M3b3*trainB[length(trainB)-1])
M3testhat[2]<-ifelse(M3predprob[2]>=ProbThreshold,1,0)
M3predprob[3] <- plogis(M3d+M3b1*M3testhat[2]+M3b2*M3testhat[1]+M3b3*trainB[length(trainB)])
M3testhat[3]<-ifelse(M3predprob[3]>=ProbThreshold,1,0)
M3predprob[4] <- plogis(M3d+M3b1*M3testhat[3]+M3b2*M3testhat[2]+M3b3*M3testhat[1])
M3testhat[4]<-ifelse(M3predprob[4]>=ProbThreshold,1,0)
M3predprob[5] <- plogis(M3d+M3b1*M3testhat[4]+M3b2*M3testhat[3]+M3b3*M3testhat[2])
M3testhat[5]<-ifelse(M3predprob[5]>=ProbThreshold,1,0)

M3Criterion1[1]<- abs(M3predprob[1]-testB[1])
M3Criterion1[2]<- abs(M3predprob[2]-testB[2])
M3Criterion1[3]<- abs(M3predprob[3]-testB[3])
M3Criterion1[4]<- abs(M3predprob[4]-testB[4])
M3Criterion1[5]<- abs(M3predprob[5]-testB[5])
M3Criterion2[1]<- (M3predprob[1]-testB[1])^2/(M3predprob[1]*(1-M3predprob[1]))
M3Criterion2[2]<- (M3predprob[2]-testB[2])^2/(M3predprob[2]*(1-M3predprob[2]))
M3Criterion2[3]<- (M3predprob[3]-testB[3])^2/(M3predprob[3]*(1-M3predprob[3]))
M3Criterion2[4]<- (M3predprob[4]-testB[4])^2/(M3predprob[4]*(1-M3predprob[4]))
M3Criterion2[5]<- (M3predprob[5]-testB[5])^2/(M3predprob[5]*(1-M3predprob[5]))

M1results[,,i] <- c(M1predprob,M1testhat,testB,M1Criterion1,M1Criterion2)
M2results[,,i] <- c(M2predprob,M2testhat,testB,M2Criterion1,M2Criterion2)
M3results[,,i] <- c(M3predprob,M3testhat,testB,M3Criterion1,M3Criterion2)
}

#Assemble results
#M1.
S11C1<-apply(M1results[,,1], 2, sum)[4]
S12C1<-apply(M1results[,,2], 2, sum)[4]
S13C1<-apply(M1results[,,3], 2, sum)[4]
S14C1<-apply(M1results[,,4], 2, sum)[4]
S15C1<-apply(M1results[,,5], 2, sum)[4]
S11C2<-apply(M1results[,,1], 2, sum)[5]
S12C2<-apply(M1results[,,2], 2, sum)[5]
S13C2<-apply(M1results[,,3], 2, sum)[5]
S14C2<-apply(M1results[,,4], 2, sum)[5]
S15C2<-apply(M1results[,,5], 2, sum)[5]
M1resultsA<-data.frame("Thresholds"=c("Perc 60%","Perc 65%","Perc 70%","Perc 75%","Perc 80%"),
	"ThresholdValues"=round(as.numeric(Thresholds),2), 
	"Criterion1"=round(c(S11C1,S12C1,S13C1,S14C1,S15C1),2),
	"MinCriterion1"=c(S11C1,S12C1,S13C1,S14C1,S15C1)==min(c(S11C1,S12C1,S13C1,S14C1,S15C1)),
	"Criterion2"=round(c(S11C2,S12C2,S13C2,S14C2,S15C2),2),
	"MinCriterion2"=c(S11C2,S12C2,S13C2,S14C2,S15C2)==min(c(S11C2,S12C2,S13C2,S14C2,S15C2)))

#M2.
S21C1<-apply(M2results[,,1], 2, sum)[4]
S22C1<-apply(M2results[,,2], 2, sum)[4]
S23C1<-apply(M2results[,,3], 2, sum)[4]
S24C1<-apply(M2results[,,4], 2, sum)[4]
S25C1<-apply(M2results[,,5], 2, sum)[4]
S21C2<-apply(M2results[,,1], 2, sum)[5]
S22C2<-apply(M2results[,,2], 2, sum)[5]
S23C2<-apply(M2results[,,3], 2, sum)[5]
S24C2<-apply(M2results[,,4], 2, sum)[5]
S25C2<-apply(M2results[,,5], 2, sum)[5]
M2resultsA<-data.frame("Thresholds"=c("Perc 60%","Perc 65%","Perc 70%","Perc 75%","Perc 80%"),
	"ThresholdValues"=round(as.numeric(Thresholds),2), 
	"Criterion1"=round(c(S21C1,S22C1,S23C1,S24C1,S25C1),2),
	"MinCriterion1"=c(S21C1,S22C1,S23C1,S24C1,S25C1)==min(c(S21C1,S22C1,S23C1,S24C1,S25C1)),
	"Criterion2"=round(c(S21C2,S22C2,S23C2,S24C2,S25C2),2),
	"MinCriterion2"=c(S21C2,S22C2,S23C2,S24C2,S25C2)==min(c(S21C2,S22C2,S23C2,S24C2,S25C2)))

#M3.
S31C1<-apply(M3results[,,1], 2, sum)[4]
S32C1<-apply(M3results[,,2], 2, sum)[4]
S33C1<-apply(M3results[,,3], 2, sum)[4]
S34C1<-apply(M3results[,,4], 2, sum)[4]
S35C1<-apply(M3results[,,5], 2, sum)[4]
S31C2<-apply(M3results[,,1], 2, sum)[5]
S32C2<-apply(M3results[,,2], 2, sum)[5]
S33C2<-apply(M3results[,,3], 2, sum)[5]
S34C2<-apply(M3results[,,4], 2, sum)[5]
S35C2<-apply(M3results[,,5], 2, sum)[5]
M3resultsA<-data.frame("Thresholds"=c("Perc 60%","Perc 65%","Perc 70%","Perc 75%","Perc 80%"),
	"ThresholdValues"=round(as.numeric(Thresholds),2), 
	"Criterion1"=round(c(S31C1,S32C1,S33C1,S34C1,S35C1),2),
	"MinCriterion1"=c(S31C1,S32C1,S33C1,S34C1,S35C1)==min(c(S31C1,S32C1,S33C1,S34C1,S35C1)),
	"Criterion2"=round(c(S31C2,S32C2,S33C2,S34C2,S35C2),2),
	"MinCriterion2"=c(S31C2,S32C2,S33C2,S34C2,S35C2)==min(c(S31C2,S32C2,S33C2,S34C2,S35C2)))

return(list("M1"=round(M1results,2),"M2"=round(M2results,2),
            "M3"=round(M3results,2),
	"M1A"=M1resultsA,"M2A"=M2resultsA,"M3A"=M3resultsA,
	"Thresholds"=Thresholds))
}

####################################################################################################

# Sum<-Sum*(64/62)

# In this study, the probability threshold was set to 0.5.
AssessPrediction5(Sum,0.5)
# However, other probability thresholds may be used as well, such as 0.4.
AssessPrediction5(Sum,0.4)

# At the end, a specific threshold should be computed for each treated patient.







