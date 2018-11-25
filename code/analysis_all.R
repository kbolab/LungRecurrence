library(OptimalCutpoints)

load("./LUNG DATA/GTV_NOFILTER.RData")
data.gtv <- as.data.frame(matrice,row.names = F,stringsAsFactors = F)
data.gtv$V1 <- sub('.*\\/', '', data.gtv$V1)
data.gtv[, c(1:length(data.gtv))] <- sapply(data.gtv[, c(1:length(data.gtv))], as.numeric)
names(data.gtv)[1] <- "CS"
load("./LUNG DATA/PERITUMOR_NOFILTER.RData")
data.peritumor <- as.data.frame(matrice,row.names = F,stringsAsFactors = F)
data.peritumor$V1 <- sub('.*\\/', '', data.peritumor$V1)
data.peritumor[, c(1:length(data.peritumor))] <- sapply(data.peritumor[, c(1:length(data.peritumor))], as.numeric)
names(data.peritumor)[1] <- "CS"
load("./LUNG DATA/LOBE_NOFILTER.RData")
data.lobe <- as.data.frame(matrice,row.names = F,stringsAsFactors = F)
data.lobe$V1 <- sub('.*\\/', '', data.lobe$V1)
data.lobe[, c(1:length(data.lobe))] <- sapply(data.lobe[, c(1:length(data.lobe))], as.numeric)
names(data.lobe)[1] <- "CS"
rm(matrice)


bootstrap <- function(n, data, model, outcome, outcomeTime,model_roi) {

  names(model$assign)[1]
  x.boot <- list()
  allAuc <- numeric()
  COX_ROC <- list()
  covariates <- cbind(outcome, outcomeTime, names(model$assign))

  for (i in 1:n){
    ind <- c()
    COX_AUC <- c()
    predittore_cox_boot <- c()
    ind <- sample(nrow(data), size=nrow(data), replace = T)
    x.boot[[i]] <- data[ind,covariates]
    x.boot[[i]] <-  x.boot[[i]][which(complete.cases(x.boot[[i]])),]
    predittore_cox_boot <- predict(object = model, newdata = x.boot[[i]])
    COX_ROC[[i]] <- (roc(x.boot[[i]][,outcome], predittore_cox_boot))
    allAuc[i] <- COX_ROC[[i]]$auc[1]
  }

  medianAuc <- median(allAuc)
  diffAuc<- numeric()

  for(i in 1:length(allAuc)){
    diffAuc[i] <- abs(COX_ROC[[i]]$auc[1] - medianAuc)
  }
  medianROC <- COX_ROC[[which(diffAuc==min(diffAuc))[1]]]
  plot(medianROC, main=paste(model_roi,"Validation ROC, (AUC=", round(x = medianROC$auc[1], digits = 3),")"))
  return(medianROC)
}


#####################################################################
library(readxl)
X2_5MDC_final_version <- read_excel("./2_5MDC final version.xlsx")

datasetLung <- X2_5MDC_final_version
datasetLung$CS <- as.character(datasetLung$CS)
datasetLung$`rec date` <- datasetLung$`rec date (if any, if not last fup)`
datasetLung$`rec date` <- as.Date(as.character(datasetLung$`rec date`), format = "%Y-%m-%d")
datasetLung$`Birth date` <- as.Date(as.character(datasetLung$`Birth date`), format = "%Y-%m-%d")
datasetLung$`surgery date` <- as.Date(as.character(datasetLung$`surgery date`), format = "%Y-%m-%d")
datasetLung$`last fup date (all pts)` <- as.Date(as.character(datasetLung$`last fup date (all pts)`), format = "%Y-%m-%d")
datasetLung$`Date of death` <- as.Date(as.character(datasetLung$`Date of death`), format = "%Y-%m-%d")
datasetLung$`TC (pre-surgery)` <- as.Date(as.character(datasetLung$`TC (pre-surgery)`), format = "%Y-%m-%d")

recurrence_days <- difftime(datasetLung$`rec date`, datasetLung$`surgery date`, units = "days")
recurrence_months <- floor(recurrence_days * 0.0328767)
datasetLung$recurrence_months <- as.numeric(recurrence_months)

age_days <- difftime(datasetLung$`surgery date`, datasetLung$`Birth date`, units = "days")
age_years <- floor(age_days * 0.00273973)
datasetLung$age_years <- as.numeric(age_years)

fup_days <- difftime(datasetLung$`last fup date (all pts)`, datasetLung$`surgery date`, units = "days")
fup_months <- floor(fup_days * 0.0328767)
datasetLung$fup_months <- as.numeric(fup_months)

death_days <- difftime(datasetLung$`Date of death`, datasetLung$`surgery date`, units = "days")
death_months <- floor(death_days * 0.0328767)
datasetLung$death_months <- as.numeric(death_months)

tc_days <- difftime(datasetLung$`surgery date`, datasetLung$`TC (pre-surgery)`, units = "days")
tc_months <- floor(tc_days * 0.0328767)
datasetLung$tc_months <- as.numeric(tc_months)



datasetLung$OverallStage <- character(length = dim(datasetLung)[1])
datasetLung$OverallStage[which(datasetLung$pT=="1a" & datasetLung$pN=="0" & datasetLung$cM=="0")] <- "IA"
datasetLung$OverallStage[which(datasetLung$pT=="1b" & datasetLung$pN=="0" & datasetLung$cM=="0")] <- "IA"
datasetLung$OverallStage[which(datasetLung$pT=="1" & datasetLung$pN=="0" & datasetLung$cM=="0")] <- "IA"
datasetLung$OverallStage[which(datasetLung$pT=="2a" & datasetLung$pN=="0" & datasetLung$cM=="0")] <- "IB"
datasetLung$OverallStage[which(datasetLung$pT=="2" & datasetLung$pN=="0" & datasetLung$cM=="0")] <- "IB"
datasetLung$OverallStage[which(datasetLung$pT=="2b" & datasetLung$pN=="0" & datasetLung$cM=="0")] <- "IIA"
datasetLung$OverallStage[which(datasetLung$pT=="2a" & datasetLung$pN=="1" & datasetLung$cM=="0")] <- "IIA"
datasetLung$OverallStage[which(datasetLung$pT=="1a" & datasetLung$pN=="1" & datasetLung$cM=="0")] <- "IIA"
datasetLung$OverallStage[which(datasetLung$pT=="1a " & datasetLung$pN=="1" & datasetLung$cM=="0")] <- "IIA"
datasetLung$OverallStage[which(datasetLung$pT=="1b" & datasetLung$pN=="1" & datasetLung$cM=="0")] <- "IIA"
datasetLung$OverallStage[which(datasetLung$pT=="2b" & datasetLung$pN=="1" & datasetLung$cM=="0")] <- "IIB"
datasetLung$OverallStage[which(datasetLung$pT=="3" & datasetLung$pN=="0" & datasetLung$cM=="0")] <- "IIB"
datasetLung$OverallStage[which(datasetLung$pT=="1a" & datasetLung$pN=="2" & datasetLung$cM=="0")] <- "IIIA"
datasetLung$OverallStage[which(datasetLung$pT=="1b" & datasetLung$pN=="2" & datasetLung$cM=="0")] <- "IIIA"
datasetLung$OverallStage[which(datasetLung$pT=="2a" & datasetLung$pN=="2" & datasetLung$cM=="0")] <- "IIIA"
datasetLung$OverallStage[which(datasetLung$pT=="2b" & datasetLung$pN=="2" & datasetLung$cM=="0")] <- "IIIA"
datasetLung$OverallStage[which(datasetLung$pT=="2" & datasetLung$pN=="2" & datasetLung$cM=="0")] <- "IIIA"
datasetLung$OverallStage[which(datasetLung$pT=="3" & datasetLung$pN=="1" & datasetLung$cM=="0")] <- "IIIA"
datasetLung$OverallStage[which(datasetLung$pT=="3" & datasetLung$pN=="2" & datasetLung$cM=="0")] <- "IIIA"
datasetLung$OverallStage[which(datasetLung$pT=="4" & datasetLung$pN=="0" & datasetLung$cM=="0")] <- "IIIA"
datasetLung$OverallStage[which(datasetLung$pT=="4" & datasetLung$pN=="1" & datasetLung$cM=="0")] <- "IIIA"

datasetLung$OverallStage.factor <- as.factor(datasetLung$OverallStage)

datasetLung$histopathology <- character(length = dim(datasetLung)[1])
datasetLung$histopathology[datasetLung$`Histo-pathology`== 1] <- "adenoca"
datasetLung$histopathology[datasetLung$`Histo-pathology`== 2] <- "non adenoca"
datasetLung$histopathology[datasetLung$`Histo-pathology`== 3] <- "non adenoca"
datasetLung$histopathology <- as.factor(datasetLung$histopathology)

datasetLung$smoker <- character(length = dim(datasetLung)[1])
datasetLung$smoker[datasetLung$`smoking history pack/year`!= 0] <- "smoker"
datasetLung$smoker[datasetLung$`smoking history pack/year`== 0] <- "non smoker"
datasetLung$smoker <- as.factor(datasetLung$smoker)

datasetLung$surgerytype <- character(length = dim(datasetLung)[1])
datasetLung$surgerytype[datasetLung$`Type of Surgery`== 1] <- "lobectomy"
datasetLung$surgerytype[datasetLung$`Type of Surgery`== 2] <- "lobectomyy"
datasetLung$surgerytype[datasetLung$`Type of Surgery`== 3] <- "lobectomy"
datasetLung$surgerytype[datasetLung$`Type of Surgery`== 4] <- "lobectomy"
datasetLung$surgerytype[datasetLung$`Type of Surgery`== 5] <- "lobectomy"
datasetLung$surgerytype[datasetLung$`Type of Surgery`== 6] <- "pneumonectomy"
datasetLung$surgerytype[datasetLung$`Type of Surgery`== 7] <- "pneumonectomy"
datasetLung$surgerytype[datasetLung$`Type of Surgery`== 8] <- "atypic resection"
datasetLung$surgerytype[datasetLung$`Type of Surgery`== 9] <- "lobectomy"
datasetLung$surgerytype <- as.factor(datasetLung$surgerytype)


datasetLung$Gender[datasetLung$Gender== 1] <- "male"
datasetLung$Gender[datasetLung$Gender== 2] <- "female"
datasetLung$Gender <- as.factor(datasetLung$Gender)

datasetLung$GTV.factor <- character(length = dim(datasetLung)[1])
datasetLung$GTV.factor[datasetLung$GTV== 1] <- "gg"
datasetLung$GTV.factor[datasetLung$GTV== 2] <- "subsolid"
datasetLung$GTV.factor[datasetLung$GTV== 3] <- "solid"
datasetLung$GTV.factor <- as.factor(datasetLung$GTV.factor)

datasetLung$age_factor <- 0
datasetLung$age_factor[datasetLung$age_years<=60] <- "<=60"
datasetLung$age_factor[datasetLung$age_years>60] <- ">60"
datasetLung$age_factor <- as.factor(datasetLung$age_factor)
######################CREATE OUTCOME COLS############
datasetLung$LR <- 0
datasetLung$LR[which(datasetLung$`recurrence type`==1)] <- 1
datasetLung$M <- 0
datasetLung$M[which(datasetLung$`recurrence type`==2)] <- 1

########
datasetLung$pT <- as.numeric(as.factor(datasetLung$pT))
datasetLung$pN <- as.numeric(as.factor(datasetLung$pN))
subsetLung<- datasetLung[,c("CS","LR","M","recurrence_months","recurrence","pT","OverallStage")]
subsetLung$CS <- as.numeric(subsetLung$CS)


#####################################################################

df.gtv.complete <- merge(data.gtv,subsetLung,by = "CS")
df.peritumor.complete <- merge(data.peritumor,subsetLung,by = "CS")
df.lobe.complete <- merge(data.lobe,subsetLung,by = "CS")

df.gtv.complete <- df.gtv.complete[,-which(names(df.gtv.complete) %in% c("F_szm.lgze","F_szm.szlge","F_szm.lzlge","meanFD","medianFD","minFD","maxFD","sdFD"))]
df.peritumor.complete <- df.peritumor.complete[,-which(names(df.peritumor.complete) %in% c("F_szm.lgze","F_szm.szlge","F_szm.lzlge"))]
df.lobe.complete <- df.lobe.complete[,-which(names(df.lobe.complete) %in% c("F_szm.lgze","F_szm.szlge","F_szm.lzlge"))]

################## GTV

# #########FEATURE SELECTION AND MODEL GTV ##############################################

library(caret)
library(fscaret)
library(rms)

set.seed(1234)
df.gtv.complete <- df.gtv.complete[,c(which(colnames(df.gtv.complete)!="recurrence"),which(colnames(df.gtv.complete)=="recurrence"))]
splitIndex <- createDataPartition(df.gtv.complete$recurrence, p = .75, list = FALSE, times = 1)
trainDF.gtv <- df.gtv.complete[ splitIndex,-which(names(df.gtv.complete) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]
testDF.gtv  <- df.gtv.complete[-splitIndex,-which(names(df.gtv.complete) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]

#fsModels <- c("glm", "gbm", "treebag", "ridge", "lasso")
myFS<-fscaret(trainDF.gtv, testDF.gtv, preprocessData=TRUE,
              Used.funcRegPred = "all", with.labels=TRUE,
              supress.output=FALSE, no.cores=2)

results.gtv <- myFS$VarImp$matrixVarImp.MSE
results.gtv$Input_no <- as.numeric(results.gtv$Input_no)
results.gtv <- results.gtv[c("SUM","SUM%","ImpGrad","Input_no")]
myFS$PPlabels$Input_no <-  as.numeric(rownames(myFS$PPlabels))
results.gtv <- merge(x=results.gtv, y=myFS$PPlabels, by="Input_no", all.x=T)
results.gtv <- results.gtv[c('Labels', 'SUM')]
results.gtv <- subset(results.gtv,results.gtv$SUM !=0)
results.gtv <- results.gtv[order(-results.gtv$SUM),]
print(results.gtv)

################################

df.gtv.complete$recurrence_months <- as.numeric(df.gtv.complete$recurrence_months)
f.gtv.cox <- as.formula(paste("Surv(recurrence_months, recurrence) ~", paste(results.gtv$Labels[1:10],collapse="+")))
cox_multi.gtv <- coxph(formula = f.gtv.cox, data=df.gtv.complete)
cox_multi_step.gtv <- step(cox_multi.gtv)
cox.zph(cox_multi_step.gtv, transform="km", global=TRUE)


library("pROC")
predittore_cox <- predict(object = cox_multi_step.gtv)
COX_ROC.gtv <- roc(df.gtv.complete$recurrence, predittore_cox,ci = T)
plot(COX_ROC.gtv,main = paste("COX GTV", " AUC:",paste(round(pROC::auc(COX_ROC.gtv),digits = 3)),sep = " "))
#bootstrap
bootResults.cox.gtv <- bootstrap(n = 100,data = df.gtv.complete,model = cox_multi_step.gtv,outcome = "recurrence",outcomeTime = "recurrence_months",model_roi = "COX GTV ")

cox.for.nomo <- rad.sign(cox_multi_step.gtv,status = "recurrence",time = "recurrence_months",dataset = df.gtv.complete)
plotNomogramCoxFromGLM(x = 10,coxResult = cox.for.nomo,outcomeTime = "recurrence_months",outcome = "recurrence")



################## PERITUMOR

# #########FEATURE SELECTION AND MODEL PERITUMOR ##############################################

library(caret)
library(fscaret)

set.seed(1234)
df.peritumor.complete <- df.peritumor.complete[,c(which(colnames(df.peritumor.complete)!="recurrence"),which(colnames(df.peritumor.complete)=="recurrence"))]
splitIndex <- createDataPartition(df.peritumor.complete$recurrence, p = .75, list = FALSE, times = 1)
trainDF.peritumor <- df.peritumor.complete[ splitIndex,-which(names(df.peritumor.complete) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]
testDF.peritumor <- df.peritumor.complete[-splitIndex,-which(names(df.peritumor.complete) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]

#fsModels <- c("glm", "gbm", "treebag", "ridge", "lasso")
myFS<-fscaret(trainDF.peritumor, testDF.peritumor, preprocessData=TRUE,
              Used.funcRegPred = "all", with.labels=TRUE,
              supress.output=FALSE, no.cores=2)

results.peritumor <- myFS$VarImp$matrixVarImp.MSE
results.peritumor$Input_no <- as.numeric(results.peritumor$Input_no)
results.peritumor <- results.peritumor[c("SUM","SUM%","ImpGrad","Input_no")]
myFS$PPlabels$Input_no <-  as.numeric(rownames(myFS$PPlabels))
results.peritumor <- merge(x=results.peritumor, y=myFS$PPlabels, by="Input_no", all.x=T)
results.peritumor <- results.peritumor[c('Labels', 'SUM')]
results.peritumor <- subset(results.peritumor,results.peritumor$SUM !=0)
results.peritumor <- results.peritumor[order(-results.peritumor$SUM),]
print(results.peritumor)

################################

##COX

df.peritumor.complete$recurrence_months <- as.numeric(df.peritumor.complete$recurrence_months)
f.peritumor.cox <- as.formula(paste("Surv(recurrence_months, recurrence) ~", paste(results.peritumor$Labels[1:10],collapse="+")))
cox_multi.peritumor <- coxph(formula = f.peritumor.cox, data=df.peritumor.complete)
cox_multi_step.peritumor <- step(cox_multi.peritumor)

library("pROC")
predittore_cox <- predict(object = cox_multi_step.peritumor)
COX_ROC.peritumor <- roc(df.peritumor.complete$recurrence, predittore_cox)
plot(COX_ROC.peritumor,main = paste("COX Peritumor", " AUC:",paste(round(pROC::auc(COX_ROC.peritumor),digits = 3)),sep = " "))
bootResults.cox.peritumor <- bootstrap(n = 100,data = df.peritumor.complete,model = cox_multi_step.peritumor,outcome = "recurrence",outcomeTime = "recurrence_months",model_roi = "COX PERITUMOR ")


################## GTV + PERITUMOR

###########FEATURE SELECTION AND MODEL PERITUMOR ##############################################

df.all.complete <- merge(df.gtv.complete[,-c(88:92)],df.peritumor.complete,by = "CS",suffixes = c(".gtv",".peritumor"))
names(df.all.complete)[174:178] <- c("meanFD.peritumor" , "medianFD.peritumor" ,"minFD.peritumor"  ,  "maxFD.peritumor"  ,  "sdFD.peritumor" )

library(caret)
library(fscaret)

set.seed(1234)
df.all.complete <- df.all.complete[,c(which(colnames(df.all.complete)!="recurrence"),which(colnames(df.all.complete)=="recurrence"))]
splitIndex <- createDataPartition(df.all.complete$recurrence, p = .75, list = FALSE, times = 1)
trainDF<- df.all.complete[ splitIndex,-which(names(df.all.complete) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]
testDF <- df.all.complete[-splitIndex,-which(names(df.all.complete) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]


fsModels <- c("glm", "gbm", "treebag", "ridge", "lasso")
myFS<-fscaret(trainDF, testDF, preprocessData=TRUE, with.labels=TRUE,
              supress.output=FALSE, no.cores=2)

results <- myFS$VarImp$matrixVarImp.MSE
results$Input_no <- as.numeric(results$Input_no)
results <- results[c("SUM","SUM%","ImpGrad","Input_no")]
myFS$PPlabels$Input_no <-  as.numeric(rownames(myFS$PPlabels))
results <- merge(x=results, y=myFS$PPlabels, by="Input_no", all.x=T)
results <- results[c('Labels', 'SUM')]
results <- subset(results,results$SUM !=0)
results <- results[order(-results$SUM),]
print(results)

################################

##COX

df.all.complete$recurrence_months <- as.numeric(df.all.complete$recurrence_months)
f.cox <- as.formula(paste("Surv(recurrence_months, recurrence) ~", paste(results$Labels[1:10],collapse="+")))
cox_multi <- coxph(formula = f.cox, data=df.all.complete)
cox_multi_step <- step(cox_multi)

library("pROC")
predittore_cox <- predict(object = cox_multi_step)
COX_ROC <- roc(df.all.complete$recurrence, predittore_cox)
plot(COX_ROC,main = paste("COX Peritumor + GTV", " AUC:",paste(round(pROC::auc(COX_ROC),digits = 3)),sep = " "))
bootResults.cox.gtvperitumor <- bootstrap(n = 100,data = df.all.complete,model = cox_multi_step,outcome = "recurrence",outcomeTime = "recurrence_months",model_roi = "COX PERITUMOR + GTV ")
plot(COX_ROC,main = paste("COX Peritumor + GTV", " AUC:",paste(round(pROC::auc(bootResults.cox.gtvperitumor),digits = 3)),sep = " "))
######################################

################## LOBE

# #########FEATURE SELECTION AND MODEL LOBE ##############################################

library(caret)
library(fscaret)

set.seed(1234)
df.lobe.complete <- df.lobe.complete[,c(which(colnames(df.lobe.complete)!="recurrence"),which(colnames(df.lobe.complete)=="recurrence"))]
splitIndex <- createDataPartition(df.lobe.complete$recurrence, p = .75, list = FALSE, times = 1)
trainDF.lobe <- df.peritumor.complete[ splitIndex,-which(names(df.lobe.complete) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]
testDF.lobe <- df.peritumor.complete[-splitIndex,-which(names(df.lobe.complete) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]

#fsModels <- c("glm", "gbm", "treebag", "ridge", "lasso")
myFS<-fscaret(trainDF.lobe, testDF.lobe, preprocessData=TRUE,
              Used.funcRegPred = "all", with.labels=TRUE,
              supress.output=FALSE, no.cores=2)

results.lobe <- myFS$VarImp$matrixVarImp.MSE
results.lobe$Input_no <- as.numeric(results.lobe$Input_no)
results.lobe <- results.lobe[c("SUM","SUM%","ImpGrad","Input_no")]
myFS$PPlabels$Input_no <-  as.numeric(rownames(myFS$PPlabels))
results.lobe <- merge(x=results.lobe, y=myFS$PPlabels, by="Input_no", all.x=T)
results.lobe <- results.lobe[c('Labels', 'SUM')]
results.lobe <- subset(results.lobe,results.lobe$SUM !=0)
results.lobe <- results.lobe[order(-results.lobe$SUM),]
print(results.lobe)

################################

##COX

df.lobe.complete$recurrence_months <- as.numeric(df.lobe.complete$recurrence_months)
f.lobe.cox <- as.formula(paste("Surv(recurrence_months, recurrence) ~", paste(results.lobe$Labels[1:10],collapse="+")))
cox_multi.lobe <- coxph(formula = f.lobe.cox, data=df.lobe.complete)
cox_multi_step.lobe <- step(cox_multi.lobe)

library("pROC")
predittore_cox <- predict(object = cox_multi_step.lobe)
COX_ROC.lobe <- roc(df.lobe.complete$recurrence, predittore_cox)
plot(COX_ROC.lobe,main = paste("COX Lobe", " AUC:",paste(round(pROC::auc(COX_ROC.lobe),digits = 3)),sep = " "))
bootResults.cox.lobe <- bootstrap(n = 100,data = df.lobe.complete,model = cox_multi_step.lobe,outcome = "recurrence",outcomeTime = "recurrence_months",model_roi = "COX LOBE ")


################## GTV + LOBE

# #########FEATURE SELECTION AND MODEL GTV + LOBE ##############################################

df.all.gtvl <- merge(df.gtv.complete[,-c(88:92)],df.lobe.complete,by = "CS",suffixes = c(".gtv",".lobe"))
names(df.all.gtvl)[174:178] <- c("meanFD.lobe"  , "medianFD.lobe" ,"minFD.lobe"  ,  "maxFD.lobe"  ,  "sdFD.lobe" )

library(caret)
library(fscaret)

set.seed(1234)
df.all.gtvl <- df.all.gtvl[,c(which(colnames(df.all.gtvl)!="recurrence"),which(colnames(df.all.gtvl)=="recurrence"))]
splitIndex <- createDataPartition(df.all.gtvl$recurrence, p = .75, list = FALSE, times = 1)
trainDF<- df.all.gtvl[ splitIndex,-which(names(df.all.gtvl) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]
testDF <- df.all.gtvl[-splitIndex,-which(names(df.all.gtvl) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]


fsModels <- c("glm", "gbm", "treebag", "ridge", "lasso")
myFS<-fscaret(trainDF, testDF, preprocessData=TRUE,
              Used.funcRegPred = "all", with.labels=TRUE,
              supress.output=FALSE, no.cores=2)

results.gtvl <- myFS$VarImp$matrixVarImp.MSE
results.gtvl$Input_no <- as.numeric(results.gtvl$Input_no)
results.gtvl <- results.gtvl[c("SUM","SUM%","ImpGrad","Input_no")]
myFS$PPlabels$Input_no <-  as.numeric(rownames(myFS$PPlabels))
results.gtvl <- merge(x=results.gtvl, y=myFS$PPlabels, by="Input_no", all.x=T)
results.gtvl <- results.gtvl[c('Labels', 'SUM')]
results.gtvl <- subset(results.gtvl,results.gtvl$SUM !=0)
results.gtvl <- results.gtvl[order(-results.gtvl$SUM),]
print(results.gtvl)

################################

##COX

df.all.gtvl$recurrence_months <- as.numeric(df.all.gtvl$recurrence_months)
f.cox.gtvl <- as.formula(paste("Surv(recurrence_months, recurrence) ~", paste(results.gtvl$Labels[1:10],collapse="+")))
cox_multi.gtvl <- coxph(formula = f.cox.gtvl, data=df.all.gtvl)
cox_multi_step.gtvl <- step(cox_multi.gtvl)

library("pROC")
predittore_cox <- predict(object = cox_multi_step.gtvl)
COX_ROC.gtvl <- roc(df.all.gtvl$recurrence, predittore_cox)
plot(COX_ROC.gtvl,main = paste("COX GTV + LOBE", " AUC:",paste(round(pROC::auc(COX_ROC.gtvl),digits = 3)),sep = " "))
bootResults.cox.gtvl <- bootstrap(n = 100,data = df.all.gtvl,model = cox_multi_step.gtvl,outcome = "recurrence",outcomeTime = "recurrence_months",model_roi = "COX GTV + LOBE ")


###############

################## GTV + PERITUMOR + LOBE

# #########FEATURE SELECTION AND MODEL PERITUMOR ##############################################

df.allall.complete <- merge(df.all.complete[,-c(179:183)],df.lobe.complete,by = "CS",suffixes = c("",".lobe"))
names(df.allall.complete)[179:269] <- paste(names(df.allall.complete)[179:269] , "lobe", sep = ".")

library(caret)
library(fscaret)

set.seed(1234)
df.allall.complete <- df.allall.complete[,c(which(colnames(df.allall.complete)!="recurrence"),which(colnames(df.allall.complete)=="recurrence"))]
splitIndex <- createDataPartition(df.allall.complete$recurrence, p = .75, list = FALSE, times = 1)
trainDF.allall<- df.allall.complete[ splitIndex,-which(names(df.allall.complete) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]
testDF.allall <- df.allall.complete[-splitIndex,-which(names(df.allall.complete) %in% c("CS","recurrence_months","recurrence","LR","pT","OverallStage"))]


fsModels <- c("glm", "gbm", "treebag", "ridge", "lasso")
myFS<-fscaret(trainDF.allall, testDF.allall, preprocessData=TRUE,
              Used.funcRegPred = "all", with.labels=TRUE,
              supress.output=FALSE, no.cores=2)

results.allall <- myFS$VarImp$matrixVarImp.MSE
results.allall$Input_no <- as.numeric(results.allall$Input_no)
results.allall <- results.allall[c("SUM","SUM%","ImpGrad","Input_no")]
myFS$PPlabels$Input_no <-  as.numeric(rownames(myFS$PPlabels))
results.allall <- merge(x=results.allall, y=myFS$PPlabels, by="Input_no", all.x=T)
results.allall <- results.allall[c('Labels', 'SUM')]
results.allall <- subset(results.allall,results.allall$SUM !=0)
results.allall <- results.allall[order(-results.allall$SUM),]
print(results.allall)

################################

##COX

df.allall.complete$recurrence_months <- as.numeric(df.allall.complete$recurrence_months)
f.cox.allall <- as.formula(paste("Surv(recurrence_months, recurrence) ~", paste(results.allall$Labels[1:10],collapse="+")))
cox_multi.allall <- coxph(formula = f.cox.allall, data=df.allall.complete)
cox_multi_step.allall <- step(cox_multi.allall)

library("pROC")
predittore_cox <- predict(object = cox_multi_step.allall)
COX_ROC.allall <- roc(df.allall.complete$recurrence, predittore_cox)
plot(COX_ROC.allall,main = paste("COX Peritumor + GTV + Lobe", " AUC:",paste(round(pROC::auc(COX_ROC.allall),digits = 3)),sep = " "))
bootResults.cox.gtvlobeper <- bootstrap(n = 100,data = df.allall.complete,model = cox_multi_step.allall,outcome = "recurrence",outcomeTime = "recurrence_months",model_roi = "COX PERITUMOR + GTV + LOBE ")


##TABELLA CON TUUE LE ROI E IMODELLI MILGIORI

tableROIs_names <- c("GTV","Peritumor","Lobe","GTV+Peritumor","GTV+Lobe","GTV+Peritumor+Lobe")
tableROIs_auc <- c(COX_ROC.gtv$auc,COX_ROC.peritumor$auc,COX_ROC.lobe$auc,COX_ROC$auc,COX_ROC.gtvl$auc,COX_ROC.allall$auc)
tableROIs_auc <- round(as.numeric(tableROIs_auc),3)
tableROIs <- cbind("ROI"=tableROIs_names,"AUC"=tableROIs_auc)






####PRIMA VERIFICARE QUALE E' LA COMBINAZIONE MIGLIORE DI GTV PEITUM E LOBE!

cox.for.nomo <- rad.sign(cox_multi_step,status = "recurrence",time = "recurrence_months",dataset = df.all.complete)
plotNomogramCoxFromGLM(x = 12,coxResult = cox.for.nomo,outcomeTime = "recurrence_months",outcome = "recurrence")

cox.for.nomo <- rad.sign(cox_multi_step,status = "recurrence",time = "recurrence_months",dataset = df.all.complete)
plotNomogramCoxFromGLM(x = 24,coxResult = cox.for.nomo,outcomeTime = "recurrence_months",outcome = "recurrence")

cox.for.nomo <- rad.sign(cox_multi_step,status = "recurrence",time = "recurrence_months",dataset = df.all.complete)
plotNomogramCoxFromGLM(x = 36,coxResult = cox.for.nomo,outcomeTime = "recurrence_months",outcome = "recurrence")

optimal.cutpoint.Youden <- optimal.cutpoints(X = "rad.signature", status = "recurrence", tag.healthy = 0,
                                             methods = "Youden", data = cox.for.nomo[[2]], pop.prev = NULL,
                                             control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95, trace = FALSE)

df.all.complete$rad.signature <- cox.for.nomo[[2]]$rad.signature
df.all.complete$rad.signature.factor <- 0
df.all.complete$rad.signature.factor[which(df.all.complete$rad.signature >= optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$cutoff)] <- 1
df.all.complete$rad.signature.factor <- as.factor(df.all.complete$rad.signature.factor)
levels(df.all.complete$rad.signature.factor)[[1]] <- "Low risk"
levels(df.all.complete$rad.signature.factor)[[2]] <- "High risk"

############################

KM <- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=df.all.complete)
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=df.all.complete)
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM, main = paste("Radiomics signature", 'p value=',round(p.val, digits = 2)), col = c('red', 'blue','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrenceat all")
legend(x = 5, y = 0.28, legend = levels(df.all.complete$rad.signature.factor), lty = 1, col = c('red', 'blue','black', 'green', 'yellow', 'grey'), cex=0.55)

######################### CONFRONTO TRA RADIOMICS * OVERALLSTAGE E OVERALL STAGE DA SOLO

datasetLung$OverallStage <- as.numeric(datasetLung$OverallStage.factor)
datasetLung$histopathology.factor <- datasetLung$histopathology
datasetLung$histopathology <- as.numeric(datasetLung$histopathology)
datasetLung$surgerytype.factor <- datasetLung$surgerytype
datasetLung$surgerytype <- as.numeric(datasetLung$surgerytype)
datasetLung$gender.factor <- datasetLung$Gender
datasetLung$gender <- as.numeric(datasetLung$Gender)


subsetLung<- datasetLung[,c("CS","LR","recurrence_months","recurrence","pT","OverallStage","histopathology","age_years","GTV","surgerytype","gender")]
subsetLung$CS <- as.numeric(subsetLung$CS)
subsetRad <- df.all.complete[,c("CS","rad.signature","rad.signature.factor")]
finaldata <- merge(subsetLung,subsetRad,by = "CS")

#RAD only
cox_rad <- coxph(formula = Surv(recurrence_months, recurrence) ~ rad.signature, data=finaldata)

#RAD + CLINICAL
cox_rad_clin <- coxph(formula = Surv(recurrence_months, recurrence) ~ rad.signature + OverallStage + histopathology + age_years + surgerytype, data=finaldata)
predittore_cox_rad_clin <- predict(object = cox_rad_clin )
COX_ROC.cox_rad_clin <- roc(finaldata$recurrence, predittore_cox_rad_clin)
plot(COX_ROC.cox_rad_clin)


## RAD + CLIN STEP
cox_rad_clin.step <- step(cox_rad_clin)
predittore_cox_rad_clin.step <- predict(object = cox_rad_clin.step )
COX_ROC.cox_rad_clin.step <- roc(finaldata$recurrence, predittore_cox_rad_clin.step)
plot(COX_ROC.cox_rad_clin.step,main = paste("COX RAD + Clinical for total recurrence", " AUC:",paste(round(pROC::auc(COX_ROC.cox_rad_clin.step),digits = 3)),sep = " "))


#cox.for.nomo.rad_clin <- rad.sign(cox_rad_clin.step,status = "recurrence",time = "recurrence_months",dataset = finaldata)
cox.for.nomo.rad_clin <- list(cox_rad_clin.step,finaldata,finaldata[,c("rad.signature","OverallStage","recurrence_months","recurrence")])
plotNomogramCoxFromGLM(x = 12,coxResult = cox.for.nomo.rad_clin,outcomeTime = "recurrence_months",outcome = "recurrence")


# JUST OVERALL STAGE
cox_juct.overalstage<- coxph(formula = Surv(recurrence_months, recurrence) ~ OverallStage, data=finaldata)
predittore_cox_juct.overalstage <- predict(object = cox_juct.overalstage )
COX_ROC.juct.overalstage <- roc(finaldata$recurrence, predittore_cox_juct.overalstage)
COX_ROC.juct.overalstage
#cox.for.nomo.justoverall <- rad.sign(cox_juct.overalstage,status = "recurrence",time = "recurrence_months",dataset = finaldata)
#plotNomogramCoxFromGLM(x = 12,coxResult = cox.for.nomo.justoverall,outcomeTime = "recurrence_months",outcome = "recurrence")

typemodel <- c("RAD only","RAD + clinical (overall staging)","Just overall staging")
aucs <- c(COX_ROC$auc,COX_ROC.cox_rad_clin.step$auc,COX_ROC.juct.overalstage$auc)
aics <- c(AIC(coxph(formula = Surv(recurrence_months, recurrence) ~ rad.signature, data=finaldata)),AIC(cox_rad_clin.step),AIC(cox_juct.overalstage))
models_aucs <- cbind("model"=typemodel,"AUC"=round(as.numeric(aucs),3),"AIC"=round(as.numeric(aics),0))

#
# library(DynNom)
# library(PASWR)
# library("rms")
# DynNom(cox_rad_clin.step.numeric, finaldata)
# plot(nomogram(cox_rad_clin.step.numeric)

##### calibration plot
form <- as.formula(Surv(recurrence_months, recurrence) ~ rad.signature + OverallStage)
cph_multi.12<- cph(formula = form,time.inc = 12,
                 data=finaldata, model = TRUE, surv = TRUE, x = TRUE, y = TRUE)
cph_multi.24<- cph(formula = form,time.inc = 24,
                   data=finaldata, model = TRUE, surv = TRUE, x = TRUE, y = TRUE)
cph_multi.36<- cph(formula = form,time.inc = 36,
                   data=finaldata, model = TRUE, surv = TRUE, x = TRUE, y = TRUE)


# Creo l'oggetto calibrate utilizzando 20 pazienti per ogni subset per avere 10 quartili
cal_harrel.12 <- calibrate(fit = cph_multi.12, cmethod = 'KM', u = 12, m = 20)
cal_harrel.24 <- calibrate(fit = cph_multi.24, cmethod = 'KM', u = 24, m = 20)
cal_harrel.36 <- calibrate(fit = cph_multi.36, cmethod = 'KM', u = 36, m = 20)

# Calibration plot
plot(cal_harrel.12, xlim = c(0,1), ylim = c(0,1), xlab="Predicted outcome", ylab="Actual outcome")
abline(a = 0, b = 1, col = 'red', lwd = 2, lty = 2)






###KAPLANS
datasetLung$pT <- as.factor(datasetLung$pT)
subsetLung.factor<- datasetLung[,c("CS","LR","recurrence_months","recurrence","M","OverallStage.factor","histopathology.factor","age_factor","GTV.factor","surgerytype","gender.factor","smoker","pT")]
subsetRad <- df.all.complete[,c("CS","rad.signature","rad.signature.factor")]
#subsetLung.factor$CS <- as.numeric(subsetLung.factor$CS)
finaldata.factor <- merge(subsetLung.factor,subsetRad,by = "CS")
finaldata.factor$surgerytype[which(finaldata.factor$surgerytype=="1"
                                   | finaldata.factor$surgerytype=="2"
                                   | finaldata.factor$surgerytype=="3"
                                   | finaldata.factor$surgerytype=="4"
                                   | finaldata.factor$surgerytype=="5"
                                   | finaldata.factor$surgerytype=="9")] <- "lobectomy"


KM.male.rec <- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$gender.factor=="male"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$gender.factor=="male"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.male.rec, main = paste("Gender: Male", 'p value=',round(p.val, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$gender.factor=="male"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.7)

KM.female <- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$gender.factor=="female"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$gender.factor=="female"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.female, main = paste("Gender: female", 'p value=',round(p.val, digits = 4)), col = c('red', 'blue','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$gender.factor=="female"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.55)

KM.adenoca <- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$histopathology.factor=="adenoca"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$histopathology.factor=="adenoca"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.adenoca, main = paste("Histology: adenoca", 'p value=',round(0.0002, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$histopathology.factor=="adenoca"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.7)

KM.nonadenoca <- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$histopathology.factor=="adenoca"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$histopathology.factor=="non adenoca"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.nonadenoca, main = paste("Histology: non adenoca", 'p value=',round(p.val, digits = 4)), col = c('red', 'blue','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$histopathology.factor=="non adenoca"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.55)

KM.lobectomy <- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$surgerytype=="lobectomy"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$surgerytype=="lobectomy"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.lobectomy, main = paste("Surgery: lobectomy", 'p value=',"0.003"), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$surgerytype=="lobectomy"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.7)

KM.atypicresection <- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$surgerytype=="atypic resection"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$surgerytype=="atypic resection"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.atypicresection, main = paste("Surgery: atypic resection", 'p value=',round(p.val, digits = 4)), col = c('red', 'blue','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$surgerytype=="atypic resection"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.55)

KM.ageold <- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$age_factor==">60"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$age_factor==">60"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.ageold, main = paste("Age: >60", 'p value=',round(p.val, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$age_factor==">60"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.7)

KM.ageyoung <- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$age_factor=="<=60"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$age_factor=="<=60"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.ageyoung, main = paste("Age: <=60", 'p value=',round(p.val, digits = 4)), col = c('red', 'blue','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$age_factor=="<=60"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.55)

### last KMs

KM.GTVsolid<- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$GTV.factor=="solid"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$GTV.factor=="solid"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.GTVsolid, main = paste("GTV: solid", 'p value=',round(p.val, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$GTV.factor=="solid"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.7)

KM.GTVnonsolid<- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$GTV.factor=="gg" | finaldata.factor$GTV.factor=="subsolid"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$GTV.factor=="gg" | finaldata.factor$GTV.factor=="subsolid"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.GTVnonsolid, main = paste("GTV: non solid", 'p value=',round(p.val, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$GTV.factor=="gg" | finaldata.factor$GTV.factor=="subsolid"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.55)

KM.t1<- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$pT=="1a" | finaldata.factor$pT=="1a " | finaldata.factor$pT=="1b"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$pT=="1a" | finaldata.factor$pT=="1a " | finaldata.factor$pT=="1b"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.t1, main = paste("T stage: 1", 'p value=',round(p.val, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$pT=="1a" | finaldata.factor$pT=="1a " | finaldata.factor$pT=="1b"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.55)

KM.t2<- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$pT=="2a" | finaldata.factor$pT=="2b"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$pT=="2a" | finaldata.factor$pT=="2b"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.t2, main = paste("T stage: 2", 'p value=',round(p.val, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$pT=="2a" | finaldata.factor$pT=="2b"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.55)

KM.t3<- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$pT=="3"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$pT=="2a" | finaldata.factor$pT=="3"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.t3, main = paste("T stage: 3", 'p value=',round(p.val, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$pT=="3"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.55)

KM.t4<- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$pT=="4"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$pT=="2a" | finaldata.factor$pT=="4"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.t4, main = paste("T stage: 4", 'p value=',round(p.val, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$pT=="4"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.55)

KM.smoker<- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$smoker=="smoker"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$pT=="2a" | finaldata.factor$smoker=="smoker"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.smoker, main = paste("Smoker", 'p value=',round(0.0004, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$smoker=="smoker"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.7)

KM.nonsmoker<- survfit(Surv(recurrence_months, recurrence) ~ rad.signature.factor,  type="kaplan-meier", conf.type="log", data=finaldata.factor[which(finaldata.factor$smoker=="non smoker"),])
KAPLANcovSign <- survdiff(Surv(recurrence_months, recurrence) ~ rad.signature.factor, data=finaldata.factor[which(finaldata.factor$pT=="2a" | finaldata.factor$smoker=="non smoker"),])
p.val<- 1 - pchisq(KAPLANcovSign$chisq, length(KAPLANcovSign$n) - 1)
plot(KM.nonsmoker, main = paste("Non smoker", 'p value=',round(p.val, digits = 4)), col = c('blue', 'red','black', 'green', 'yellow', 'grey'),xlab = "t (months)", ylab="Recurrence at all")
legend(x = 5, y = 0.28, legend = levels(finaldata.factor[which(finaldata.factor$smoker=="non smoker"),]$rad.signature.factor), lty = 1, col = c('blue', 'red','black', 'green', 'yellow', 'grey'), cex=0.55)


##### WRITE REPORT PDF #########################
fileName.PDF<-'./totalR.pdf'
render("./total_recurrence_final.Rmd", "pdf_document", fileName.PDF, output_dir="report")
