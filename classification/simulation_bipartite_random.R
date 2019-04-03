source("energytree.R")
require(igraph)
require(reshape)
require(plyr)
require(Metrics)

set.seed(123)
n_simulations <- 50
acc_v <- list()
acc_v[["etree"]] <- rep(NA, n_simulations)
acc_v[["dtree"]] <- rep(NA, n_simulations)
acc_v[["ctree"]] <- rep(NA, n_simulations)

n.graphs <- 200
for(j in 1:n_simulations) {
  print(paste("SIMULATION ", j))
  
  all.graphs <- list()
  Y.s <- rep(NA, n.graphs)
  no.of.nodes.s <- seq(20, 100)
  
  # creates 100 bipartite graphs
  for(i in 1:(n.graphs/2)) {
    g <- sample_bipartite(75, 125, p=.07)
    all.graphs[[i]] <- g
    Y.s[[i]] <- 0
  }
  
  # creates 100 random graphs
  for(i in (n.graphs/2 + 1):n.graphs) {
    g <- sample_gnp(directed = F, loops = F, p = 0.03, n = 200)
    all.graphs[[i]] <- g
    Y.s[[i]] <- 1
  }
  
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
                       minbucket = 5,
                       alpha = 0.05)
  plot(myS)
  
  ###### PREDICTION
  Y.test <- test[[which(names(test) == "Y")]]
  test <- test[which(names(test) != "Y")]
  
  my.pred <- my.predict(model = myS, newdata = test)
  my.pred <- round(my.pred)
  
  Y.train <- list("Y" = train$Y)
  train <- cbind(myS$data, Y.train)
  test <- list2matrix(test, myS$data, NULL)
  
  ###### COMPARISON WITH DECISION TREES (CART)
  library(rpart)
  dtree.model <- rpart(as.factor(Y)~., data=train, method="class")
  dtree.pred <- predict(object = dtree.model, newdata = test, type="class")
  
  ###### COMPARISON WITH CONDITIONAL TREES
  ctree.model <- ctree(as.factor(Y)~., data=train)
  ctree.pred <- predict(object = ctree.model, newdata = test)
  
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