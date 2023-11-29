library(glmnet)
library(Matrix)
library(foreach)
library(pROC)
library(caret)
library(lattice)

# Input data
data <- data <- read.csv("./data/example_histone.csv",header=T)
histone_name <- colnames(data)

# Generate the nonlinear feature library
recond_2 <- c()
x <- 1:5
for (i in 1:5) {
  for (j in i:5) {
    a <- c(x[i],x[j])
    recond_2 <-rbind(recond_2,a)
  }
}
recond_3 <- c()
x <- 1:5
for (i in 1:5) {
  for (j in i:5) {
    for (k in j:5) {
      
      a <- c(x[i],x[j],x[k])
      recond_3 <-rbind(recond_3,a)}
    
  }
}
recond_4 <- c()
x<-1:5
for (i in 1:5) {
  for (j in i:5) {
    for (k in j:5) {
      for (m in k:5) {
        a <- c(x[i],x[j],x[k],x[m])
        recond_4 <- rbind(recond_4,a)
      }
      
    }
  }
}
recond_5 <- c()
x<-1:5
for (i in 1:5) {
  for (j in i:5) {
    for (k in j:5) {
      for (m in k:5) {
        for (n in m:5) {
          
          
          a <- c(x[i],x[j],x[k],x[m],x[n])
          recond_5 <- rbind(recond_5,a)
        }
      }
    }
  }
}
# Combination of two histone features
combin_matrix <- t(recond_2)
recond <- matrix(NA,nrow(data),ncol(combin_matrix))
recond_name <- c()
for(i in 1:ncol(combin_matrix)){
  feature_vec <- combin_matrix[,i]
  data_sub = data[,feature_vec[1]]*data[,feature_vec[2]]
  sub_name = histone_name[feature_vec]
  merge_name<- paste(sub_name[1],"*",sub_name[2],sep = "")
  recond[,i] <- data_sub
  recond_name[i] <- merge_name
}
colnames(recond)<-recond_name
histone_2 <- recond

# Combination of three histone features
combin_matrix <- t(recond_3)
recond <- matrix(NA,nrow(data),ncol(combin_matrix))
recond_name <- c()
for(i in 1:ncol(combin_matrix)){
  feature_vec <- combin_matrix[,i]
  data_sub = data[,feature_vec[1]]*data[,feature_vec[2]]*data[,feature_vec[3]]
  sub_name = histone_name[feature_vec]
  merge_name<- paste(sub_name[1],"*",sub_name[2],"*",sub_name[3],sep = "")
  recond[,i] <- data_sub
  recond_name[i] <- merge_name
}
colnames(recond)<-recond_name
histone_3 <- recond

# Combination of four histone features
combin_matrix <- t(recond_4)
recond <- matrix(NA,nrow(data),ncol(combin_matrix))
recond_name <- c()
for(i in 1:ncol(combin_matrix)){
  feature_vec <- combin_matrix[,i]
  data_sub = data[,feature_vec[1]]*data[,feature_vec[2]]*data[,feature_vec[3]]*data[,feature_vec[4]]
  sub_name = histone_name[feature_vec]
  merge_name<- paste(sub_name[1],"*",sub_name[2],"*",sub_name[3],"*",sub_name[4],sep = "")
  recond[,i] <- data_sub
  recond_name[i] <- merge_name
}
colnames(recond)<-recond_name
histone_4 <- recond

# Combination of five histone features
combin_matrix <- t(recond_5)
recond <- matrix(NA,nrow(data),ncol(combin_matrix))
recond_name <- c()
for(i in 1:ncol(combin_matrix)){
  feature_vec <- combin_matrix[,i]
  data_sub = data[,feature_vec[1]]*data[,feature_vec[2]]*data[,feature_vec[3]]*data[,feature_vec[4]]*data[,feature_vec[5]]
  sub_name = histone_name[feature_vec]
  merge_name<- paste(sub_name[1],"*",sub_name[2],"*",sub_name[3],"*",sub_name[4],sub_name[5],sep = "")
  recond[,i] <- data_sub
  recond_name[i] <- merge_name
}
colnames(recond)<-recond_name
histone_5 <- recond

data <- cbind(histone_2,histone_3,histone_4,histone_5,data)


# Ordinary logistic regression
OLR <- function(data, folds = 5) {
  set.seed(2)
  require(caret)
  all_fold_predict <- c()
  all_true_value <- c()
  for (i in 1:folds) {
    fold_test <- data[createFolds(y = data$Label, k = folds)[[i]],]
    fold_train <- data[-createFolds(y = data$Label, k = folds)[[i]],]
    fold_pre <- glm(Label ~ ., family = binomial(link = 'logit'), data = fold_train)
    fold_predict <- predict(fold_pre, type = 'response', newdata = fold_test)
    true_value <- fold_test$Label
    all_fold_predict <- c(all_fold_predict, fold_predict)
    all_true_value <- c(all_true_value, true_value)
  }
  
  modelroc <- roc(all_true_value, all_fold_predict)
  auc_value <- auc(modelroc)
  auc_value<-round(auc_value,4)
  return(auc_value)
}

# Evaluate individual histone prediction performance
feature_recond <- c()
auc_recond<-c()
feature <- data[,-ncol(data)]
feature_name <- colnames(feature)
feature_length <- length(feature_name)
recond <- matrix(NA,feature_length,2)
for (i in 1:feature_length) {
  print(i)
  data_sub <- data.frame(feature[,i],data$Label)
  colnames(data_sub) <- c(feature_name[i],"Label")
  auc_value <- OLR(data_sub)
  recond[i,1] <- feature_name[i]
  recond[i,2] <- auc_value
}
value <- as.numeric(recond[,2])
index <- which(value==max(value))
feature_base <- feature[,index]
feature_all <- feature[,-(index)]
feature_recond[1]<- feature_name[index]
auc_recond[1] <- max(value)
# Forward feature search
steps = 30
for (k in 2:steps) {
  feature <- feature_all
  feature_name <- colnames(feature)
  feature_length <- length(feature_name)
  recond <- matrix(NA,feature_length,2)
  for (i in 1:feature_length) {
    data_sub <- data.frame(feature_base,feature[,i],data$Label)
    colnames(data_sub) <- c(feature_recond,feature_name[i],"Label")
    auc_value <- OLR(data_sub)
    recond[i,1] <- feature_name[i]
    recond[i,2] <- auc_value
  }
  value <- as.numeric(recond[,2])
  index <- which(value==max(value))
  feature_base <- cbind(feature_base,feature[,index])
  feature_all <- feature[,-(index)]
  feature_recond[k]<- feature_name[index]
  auc_recond[k] <- max(value)
  if(auc_recond[k]<auc_recond[k-1]){
    break
  }
  
}
result <- cbind(feature_recond,auc_recond)
write.csv(feature_recond,"output_path/result.csv")
