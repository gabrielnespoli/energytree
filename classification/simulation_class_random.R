source("energytree.R")
require(igraph)
require(reshape)
require(plyr)
require(Metrics)
require(pracma)

set.seed(1)
n_simulations <- 50
acc_v <- list()
acc_v[["etree"]] <- rep(NA, n_simulations)
acc_v[["dtree"]] <- rep(NA, n_simulations)
acc_v[["ctree"]] <- rep(NA, n_simulations)

n.graphs <- 400
all.graphs <- list()
Y.s <- rep(NA, n.graphs)
for(j in 1:n_simulations) {
  print(paste("SIMULATION ", j))
  
  # creates 200 random graphs
  for(i in 1:134) {
    g <- sample_gnp(directed = F, loops = F, p = 0.01, n = 200)
    all.graphs[[i]] <- g
    Y.s[[i]] <- 1
  }
  for(i in 135:267) {
    g <- sample_gnp(directed = F, loops = F, p = 0.015, n = 200)
    all.graphs[[i]] <- g
    Y.s[[i]] <- 2
  }
  for(i in 268:400) {
    g <- sample_gnp(directed = F, loops = F, p = 0.02, n = 200)
    all.graphs[[i]] <- g
    Y.s[[i]] <- 3
  }
  
  # for(i in 1:100) {
  #   g <- sample_gnp(directed = F, loops = F, p = 0.01, n = 200)
  #   all.graphs[[i]] <- g
  #   Y.s[[i]] <- 1
  # }
  # for(i in 101:200) {
  #   g <- sample_gnp(directed = F, loops = F, p = 0.05, n = 200)
  #   all.graphs[[i]] <- g
  #   Y.s[[i]] <- 2
  # }
  
  # shuffle the data with the same seed
  set.seed(123)
  Y.s <- sample(Y.s)
  set.seed(123)
  all.graphs <- sample(all.graphs)
  
  data <- vector("list", 2)
  names(data) <- c("Y", "V1")
  data$Y <- Y.s
  data$V1 <- all.graphs
  
  ###### split the data into train and test sets
  n = length(data[[1]])
  rnd.ind = sample(x = seq(n), size = 0.8*n)
  train = test = data
  n.var <- which(names(data) != "Y")
  
  # this line is used to generate the data to me used in the comparison with other algorithms
  # it's not used by ours
  m.data <- list2matrix(data[-1], NULL, list("Y" = data$Y))
  
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
  train$Y <- as.factor(train$Y)
  myS  <-  mytree     ("Y", data = train, weights = NULL,
                       minbucket = 17,
                       alpha = 0.05)
  plot(myS)
  
  ###### PREDICTION
  Y.test <- test[[which(names(test) == "Y")]]
  test <- test[which(names(test) != "Y")]
  my.pred <- my.predict(model = myS, newdata = test)
  # my.pred <- round(my.pred)
  
  ###### GET THE TABULAR DATA TO BE USED BY OTHER ALGORITHMS
  train <- m.data[c(rnd.ind),]
  test <- m.data[-c(rnd.ind),]
  test <- test[which(names(test) != "Y")]
  train$Y <- as.factor(train$Y)
  
  ###### COMPARISON WITH DECISION TREES (CART)
  library(rpart)
  dtree.model <- rpart(Y~., data=train)
  dtree.pred <- predict(object = dtree.model, newdata = test)
  dtree.pred <- apply(data.frame(dtree.pred), 1, which.max) #convert the probability matrix to one single vector
  
  
  ###### COMPARISON WITH CONDITIONAL TREES
  # for an unknown reason, ctree is not accepting numerical names for columns, so we 
  # converted it to strings
  names(train) <- lapply(names(train), function(i) {
    if(i != "Y")
      paste("V1.",i,sep="")
    else
      "Y"
  })
  names(test) <- lapply(names(test), function(i) {
    if(i != "Y")
      paste("V1.",i,sep="")
    else
      "Y"
  })
  ctree.model <- ctree(Y~., data = train, control = ctree_control(splitstat = "maximum"))
  ctree.pred <- predict(object = ctree.model, newdata = test)
  # ctree.pred <- apply(data.frame(ctree.pred), 1, which.max) #convert the probability matrix to one single vector
  
  acc_v[["etree"]][j] <- accuracy(Y.test, my.pred)
  acc_v[["dtree"]][j] <- accuracy(Y.test, dtree.pred)
  acc_v[["ctree"]][j] <- accuracy(Y.test, ctree.pred)
}

###### COMPARISON
print(paste("Acc Energy Tree:", round(mean(acc_v[["etree"]]),4)))
print(paste("Acc Decision Tree:", round(mean(acc_v[["dtree"]]),4)))
print(paste("Acc Conditional Tree:", round(mean(acc_v[["ctree"]]),4)))
print("------------------------")
print(paste("SD Energy Tree:", round(sd(acc_v[["etree"]]),4)))
print(paste("SD Decision Tree:", round(sd(acc_v[["dtree"]]),4)))
print(paste("SD Conditional Tree:", round(sd(acc_v[["ctree"]]),4)))

require(caret)
print(confusionMatrix(as.factor(Y.test),as.factor(my.pred)))

