source("energytree.R")
require(igraph)
require(reshape)
require(plyr)
require(Metrics)
require(pracma)

load("../data/data_carcrash.Rdata")
data2 <- vector("list", 2)
names(data2) <- c("Y", "V1")
data2[["Y"]] <- data$y
data2[["V1"]] <- data$X
data <- data2
data2 <- NULL

set.seed(123)
n = length(data[[1]])
rnd.ind = sample(x = seq(n), size = 0.8*n)
train = test = data

# split the data
for(i in 1:length(data)) {
  if(class(train[[i]])=="numeric" | class(train[[i]])=="integer" | class(train[[i]])=="double" | class(train[[i]])=="list") {
    train[[i]] = train[[i]][c(rnd.ind)]
    test[[i]] = test[[i]][-c(rnd.ind)] 
  }
  else if(class(train[[i]])=="data.frame") {
    train[[i]] = train[[i]][c(rnd.ind),]
    test[[i]] = test[[i]][-c(rnd.ind),] 
  }
}

# this line is used to generate the data to me used in the comparison with other algorithms
# it's not used by ours
m.data <- list2matrix(lst = data[-1], struct2copy = NULL, Y = list("Y" = data$Y), nb = NULL)

Y.test <- test[[which(names(test) == "Y")]]
test <- test[which(names(test) != "Y")]
for(a in seq(0.01, 0.4, 0.01)) {
  set.seed(123)
  myS  <-  mytree     ("Y", data = train, weights = NULL,
                       minbucket = 3,
                       alpha = a)
  ctree.model <- ctree(formula = Y~., data = train, control = ctree_control(alpha = a))
  
  #plot(myS)
  
  ###### PREDICTION
  my.pred <- my.predict(model = myS, newdata = test)
  ctree.pred <- predict(object = ctree.model, newdata = test)
  
  print(paste("alpha:" ,a, " - RMSE Energy Tree:", rmse(Y.test, my.pred)))
  print(paste("alpha:" ,a, " - RMSE Conditional Inference Tree:", rmse(Y.test, ctree.pred)))
}

###### GET THE TABULAR DATA TO BE USED BY OTHER ALGORITHMS
#change name of variables (it was not accepting numbers as header)
header <- names(m.data)
header <- unlist(lapply(seq(1,length(header)), function(i) {
  header[i] <- paste("V1.",i,sep="")
}))
header[length(header)] <- 'Y'
names(m.data) <- header

train <- m.data[c(rnd.ind),]
test <- m.data[-c(rnd.ind),]
test <- test[which(names(test) != "Y")]

###### COMPARISON WITH DECISION TREES (CART)
library(rpart)
dtree.model <- rpart(Y~., data=train)
dtree.pred <- predict(object = dtree.model, newdata = test)

###### COMPARISON WITH CONDITIONAL TREES
ctree.model <- ctree(formula = Y~., data = train, control = ctree_control(alpha = 0.2))
ctree.pred <- predict(object = ctree.model, newdata = test)

###### COMPARISON
print(paste("RMSE Energy Tree:", rmse(Y.test, my.pred)))
print(paste("RMSE Decision Tree:", rmse(Y.test, dtree.pred)))
print(paste("RMSE Conditional Inference Tree:", rmse(Y.test, ctree.pred)))

