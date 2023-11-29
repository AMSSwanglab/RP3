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
record_2 <- c()
x <- 1:5
for (i in 1:5) {
  for (j in i:5) {
    a <- c(x[i],x[j])
    record_2 <-rbind(record_2,a)
  }
}
record_3 <- c()
x <- 1:5
for (i in 1:5) {
  for (j in i:5) {
    for (k in j:5) {
      
      a <- c(x[i],x[j],x[k])
      record_3 <-rbind(record_3,a)}
    
  }
}
record_4 <- c()
x<-1:5
for (i in 1:5) {
  for (j in i:5) {
    for (k in j:5) {
      for (m in k:5) {
        a <- c(x[i],x[j],x[k],x[m])
        record_4 <- rbind(record_4,a)
      }
      
    }
  }
}
record_5 <- c()
x<-1:5
for (i in 1:5) {
  for (j in i:5) {
    for (k in j:5) {
      for (m in k:5) {
        for (n in m:5) {
          
          
          a <- c(x[i],x[j],x[k],x[m],x[n])
          record_5 <- rbind(record_5,a)
        }
      }
    }
  }
}
# Combination of two histone features
combin_matrix <- t(record_2)
record <- matrix(NA,nrow(data),ncol(combin_matrix))
record_name <- c()
for(i in 1:ncol(combin_matrix)){
  feature_vec <- combin_matrix[,i]
  data_sub = data[,feature_vec[1]]*data[,feature_vec[2]]
  sub_name = histone_name[feature_vec]
  merge_name<- paste(sub_name[1],"*",sub_name[2],sep = "")
  record[,i] <- data_sub
  record_name[i] <- merge_name
}
colnames(record)<-record_name
histone_2 <- record

# Combination of three histone features
combin_matrix <- t(record_3)
record <- matrix(NA,nrow(data),ncol(combin_matrix))
record_name <- c()
for(i in 1:ncol(combin_matrix)){
  feature_vec <- combin_matrix[,i]
  data_sub = data[,feature_vec[1]]*data[,feature_vec[2]]*data[,feature_vec[3]]
  sub_name = histone_name[feature_vec]
  merge_name<- paste(sub_name[1],"*",sub_name[2],"*",sub_name[3],sep = "")
  record[,i] <- data_sub
  record_name[i] <- merge_name
}
colnames(record)<-record_name
histone_3 <- record

# Combination of four histone features
combin_matrix <- t(record_4)
record <- matrix(NA,nrow(data),ncol(combin_matrix))
record_name <- c()
for(i in 1:ncol(combin_matrix)){
  feature_vec <- combin_matrix[,i]
  data_sub = data[,feature_vec[1]]*data[,feature_vec[2]]*data[,feature_vec[3]]*data[,feature_vec[4]]
  sub_name = histone_name[feature_vec]
  merge_name<- paste(sub_name[1],"*",sub_name[2],"*",sub_name[3],"*",sub_name[4],sep = "")
  record[,i] <- data_sub
  record_name[i] <- merge_name
}
colnames(record)<-record_name
histone_4 <- record

# Combination of five histone features
combin_matrix <- t(record_5)
record <- matrix(NA,nrow(data),ncol(combin_matrix))
record_name <- c()
for(i in 1:ncol(combin_matrix)){
  feature_vec <- combin_matrix[,i]
  data_sub = data[,feature_vec[1]]*data[,feature_vec[2]]*data[,feature_vec[3]]*data[,feature_vec[4]]*data[,feature_vec[5]]
  sub_name = histone_name[feature_vec]
  merge_name<- paste(sub_name[1],"*",sub_name[2],"*",sub_name[3],"*",sub_name[4],sub_name[5],sep = "")
  record[,i] <- data_sub
  record_name[i] <- merge_name
}
colnames(record)<-record_name
histone_5 <- record

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
feature_record <- c()
auc_record<-c()
feature <- data[,-ncol(data)]
feature_name <- colnames(feature)
feature_length <- length(feature_name)
record <- matrix(NA,feature_length,2)
for (i in 1:feature_length) {
  print(i)
  data_sub <- data.frame(feature[,i],data$Label)
  colnames(data_sub) <- c(feature_name[i],"Label")
  auc_value <- OLR(data_sub)
  record[i,1] <- feature_name[i]
  record[i,2] <- auc_value
}
value <- as.numeric(record[,2])
index <- which(value==max(value))
feature_base <- feature[,index]
feature_all <- feature[,-(index)]
feature_record[1]<- feature_name[index]
auc_record[1] <- max(value)
# Forward feature search
steps = 30
for (k in 2:steps) {
  feature <- feature_all
  feature_name <- colnames(feature)
  feature_length <- length(feature_name)
  record <- matrix(NA,feature_length,2)
  for (i in 1:feature_length) {
    data_sub <- data.frame(feature_base,feature[,i],data$Label)
    colnames(data_sub) <- c(feature_record,feature_name[i],"Label")
    auc_value <- OLR(data_sub)
    record[i,1] <- feature_name[i]
    record[i,2] <- auc_value
  }
  value <- as.numeric(record[,2])
  index <- which(value==max(value))
  feature_base <- cbind(feature_base,feature[,index])
  feature_all <- feature[,-(index)]
  feature_record[k]<- feature_name[index]
  auc_record[k] <- max(value)
  if(auc_record[k]<auc_record[k-1]){
    break
  }
  
}
result <- cbind(feature_record,auc_record)
write.csv(feature_record,"output_path/result.csv")
