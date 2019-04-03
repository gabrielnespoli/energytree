source("energytree.R")
require(igraph)
require(reshape)
require(plyr)
require(Metrics)
require(pracma)

load("../data/data.Rdata")

set.seed(123)
n = length(data[[1]])
rnd.ind = sample(x = seq(n), size = 0.8*n)
train = test = data

# split the data
for(i in 1:length(data)) {
  if(class(train[[i]])=="numeric" | class(train[[i]])=="double" | class(train[[i]])=="list") {
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
m.data <- list2matrix(data[-1], NULL, list("Y" = data$Y))

train = train[-3] # remove the functional connectivity variable

myS  <-  mytree     ("Y", data = train, weights = NULL,
                     minbucket = 1,
                     alpha = 0.85)
plot(myS)

###### PREDICTION
Y.test <- test[[which(names(test) == "Y")]]
test <- test[which(names(test) != "Y")]
my.pred <- my.predict(model = myS, newdata = test)

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

###### PREDICTION WITH LINEAR REGRESSION
lm.model <- lm(formula = Y~., data = train)
lm.pred <- predict(lm.model, newdata = test)

###### COMPARISON WITH CONDITIONAL TREES
ctree.model <- ctree(formula = Y~., data = train)
ctree.pred <- predict(object = ctree.model, newdata = test)

###### COMPARISON
print(paste("RMSE Energy Tree:", rmse(Y.test, my.pred)))
print(paste("RMSE Decision Tree:", rmse(Y.test, dtree.pred)))
print(paste("RMSE Linear Regression:", rmse(Y.test, lm.pred)))
print(paste("RMSE Conditional Inference Tree:", rmse(Y.test, ctree.pred)))
