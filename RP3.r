library(glmnet)
library(Matrix)
library(foreach)
library(pROC)
library(PRROC)
library(caret)
library(lattice)

#Input data
data <- read.csv("./data/example_data.csv",header = T,row.names=1)
ppi <- read.csv("./data/PPI.csv",header = F)
cross_tissue_expression <- read.csv("./data/cross_tissue_gene_expression.csv",header = T)

#Set cutoff
cutoff <- 0.75
TF <- colnames(cross_tissue_expression)
TF_number <- ncol(cross_tissue_expression)

ppi <- as.matrix(ppi)
ppi_matrix <- matrix(0,TF_number,TF_number)
rownames(ppi_matrix) <- TF
colnames(ppi_matrix) <- TF

for(i in 1:TF_number){
  index1 <- which(TF[i]==ppi[,1])
  index_vector1 <- ppi[index1,2]
  index2 <- which(TF[i]==ppi[,2])
  index_vector2 <- ppi[index2,1]
  index_vector <- c(index_vector1,index_vector2)
  index_final <- c()
  for(j in index_vector){
    index3 <- which(j==TF)
    index_final <- c(index_final,index3)
  }
  ppi_matrix[index_final,i] <- 1
}

co_expression <- cor(cross_tissue_expression)
co_expression[is.na(co_expression)] <- 0
PPIandCoExPression <- ppi_matrix* co_expression
PPIandCoExPression[which(PPIandCoExPression <= cutoff)]<-0
TFweight <- colSums(PPIandCoExPression)
candidateTF_index <- which(TFweight > 0)
candidateTF <- TF[candidateTF_index]
candidateTF_weight <- TFweight[candidateTF_index]

TrainingData <- data.frame(data[,candidateTF_index],data[,(ncol(data)-11):ncol(data)])
other_names <- colnames(data[,(ncol(data)-11):ncol(data)])
name <- c(candidateTF,other_names)
colnames(TrainingData) <- name
other_w <- rep(0,11)
w <- c(candidateTF_weight,other_w)

set.seed(2)
require(caret)
folds <- createFolds(y=TrainingData$Label,k=10)
fold_predict <- c()
true_value <- c()
a <- ncol(TrainingData) -1
auc_record<-c()
model_recond <- list()
coefficients_recond <- list()

for(i in 1:10){
  fold_test <- TrainingData[folds[[i]],]   
  fold_train <- TrainingData[-folds[[i]],] 
  y1<-fold_test$Label
  y2<-fold_train$Label
  x1<-fold_test[,1:a]
  x2<-fold_train[,1:a]
  x1<-as.matrix(x1)
  x2<-as.matrix(x2)
  fold_pre <-cv.glmnet(x2,y2,family = "binomial",penalty.factor=w)
  fold_predict1 <- predict(fold_pre,type='response',x1)
  fold_predict1<- as.numeric(fold_predict1)
  true_value1 <- y1
  fold_predict <- append(fold_predict,fold_predict1)
  true_value <- append(true_value,true_value1)
  fit<-fold_pre
  coefficients <- coef(fit,s=fit$lambda.min)
  coefficients_recond[[i]] <- coefficients
  TFcoefficents <- coefficients[,1]
  Active.Index<-which(TFcoefficents!=0)
  GuidingTF <- names(Active.Index)
  model_recond[[i]] <- GuidingTF
  modelroc_fold <- roc(true_value,fold_predict)
  auc_fold<-auc(modelroc_fold)
  auc_record[i] <- auc_fold
}
result_rp3 <- data.frame(true_value,fold_predict)

modelroc <- roc(result_rp3$true_value,result_rp3$fold_predict)
auroc_value<-auc(modelroc)
auroc_value<-round(auroc_value,4)
#set pdf output path 
pdf("./output_path/roc_pr.pdf",width = 9,height = 9)
par(pty="s")
plot(modelroc, col='#377eb8',main = "ROC curve",auc.polygon=TRUE,auc.polygon.col="white",legacy.axes=TRUE,
     grid.col=c("black", "black"),xaxs="i",yaxs="i",font.lab=2,cex.main=2,cex.lab=2,lwd=3)
feature_name<-c("RP3, AUC =")
auc_legend <-paste(feature_name,auroc_value)
legend("bottomright",legend = auc_legend,col='#377eb8',lwd = 4,cex=1.5)

fg = result_rp3$fold_predict[which(result_rp3$true_value==1)]
bg = result_rp3$fold_predict[which(result_rp3$true_value==0)]
pr <-pr.curve(scores.class0 = fg, scores.class1 = bg,curve = TRUE)
aupr_value <- round(pr$auc.integral,4)
plot(pr, col='#4daf4a',auc.polygon=TRUE,auc.polygon.col="white",legacy.axes=TRUE,auc.main = F,
     grid.col=c("black", "black"),xaxs="i",yaxs="i",font.lab=2,cex.main=2,cex.lab=2,lwd=3)
aupr_legend<-paste(feature_name,aupr_value)
legend("bottomleft",legend = aupr_legend,col='#4daf4a',lwd = 4,cex=1.5)
dev.off()

best_model_index<- which.max(auc_record)
best_fetures <- model_recond[[best_model_index]]
fetures_coefficients <- coefficients_recond[[best_model_index]]
fetures_coefficients <- as.matrix(fetures_coefficients)
fetures_coefficients <- as.data.frame(fetures_coefficients)

write.csv(best_fetures,"output_path/output_tf.csv")
