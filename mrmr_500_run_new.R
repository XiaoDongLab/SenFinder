#clear the R environment
rm(list=ls())

#call library
library(dplyr)
library(tidyr)
library(e1071)
library(caret)
library(mRMRe)
library(randomForest)
library(rpart)
library(rrcovHD)


# read the dataset
# set seed
set.seed(12)

#read files
dat1=read.csv("./final_dat.csv")
dat1=dat1[,-1]
colnames(dat1)[1]<-"senescence"

# remove the columns which has zero values
dat2=dat1[,colSums(dat1 != 0)>0]


# create a list for adding dataset
xnum<-c()
for(m in 1:9)
  xnum<-c(xnum,rep(m,11))

xnum<-c(xnum,rep(10,15))
df<-as.data.frame(xnum)

#Check the number
table(xnum)

#data modification
dat2$sample_no=factor(dat2$sample_no)

# create mrmr object
dat2[,1:45394] <- lapply(dat2[,1:45394],as.numeric)
dd <- mRMR.data(data=dat2)

#Call value storage
dst<-c()
rft<-c()
svt<-c()
rst<-c()

# Iterative feature selection with 10 fold cross validation
# Loop for 500 iteration
for(i in 1:500){
  
  # 10 times value store
  accu<-c()
  bccu<-c()
  cccu<-c()
  dccu<-c()
  
  # Create data set for MRMR
  # Classify the features
  xx<-mRMR.classic(data=dd,target_indices=c(1),feature_count=i)
  
  # Identify the index
  imp_index <-solutions(xx)[[1]]
  
  # Create a new dataset
  newnames<-c(1,imp_index)
  newdat1<-dat2[,newnames]
  newdat1<-cbind(newdat1,df)
  
  # Loop for 10 fold cross validation
  for (j in 1:10){
    #create the training dataset
    intest<-newdat1[xnum==j,]
    intest<-subset(intest,select=-c(xnum))
   
    #create the test dataset
    intrain<-newdat1[!(xnum==j),]
    intrain<-subset(intrain,select=-c(xnum))
    
    
    # run for ds
    # train the data with decision tree
    model_dt <-rpart(senescence ~ . , intrain, method="class")
    # test the data
    pred1 <- predict(model_dt,intest[-1],type="class")
    # store results
    conf1=confusionMatrix(as.factor(intest$senescence),as.factor(pred1))
    accu[j] <- conf1$overall[1]
    
    #run for rf
    # train the data with RF
    model_rf <-randomForest(as.factor(senescence) ~ . , intrain, ntree=50)
    # test the data
    pred2 <- predict(model_rf,intest[-1],type="class")
    # store results
    conf2=confusionMatrix(as.factor(intest$senescence),as.factor(pred2))
    bccu[j] <- conf2$overall[1]
    
    # run for SVM
    # train the data with svm
    model_svm <- svm(as.factor(senescence) ~ .,data=intrain, kernel="linear")
    # test the data
    pred3 <- predict(model_svm,intest[-1],type="class")
    # store results
    conf3=confusionMatrix(as.factor(intest$senescence),as.factor(pred3))
    cccu[j] <- conf3$overall[1]
    
    # run for rs
    # train the model with rsimca
    model_rs <-RSimca(as.factor(senescence) ~ . , intrain)
    # test the data
    pred4 <- predict(model_rs,intest[-1],type="class")
    
    # store results
    conf=confusionMatrix(as.factor(intest$senescence),as.factor(pred4@classification))
    dccu[j] <- conf$overall[1]
    
  }
  dst[i]=mean(accu)
  rft[i]=mean(bccu)
  svt[i]=mean(cccu)
  rst[i]=mean(dccu)
  
}

# use all save variable to make dataframe
simout<-cbind(dst,rft,svt,rst)

# write the output
write.csv(simout,"./simout.csv")


