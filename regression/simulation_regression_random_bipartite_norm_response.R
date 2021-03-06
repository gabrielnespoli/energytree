source("energytree.R")
require(igraph)
require(reshape)
require(plyr)
require(Metrics)
require(pracma)

set.seed(123)

n_simulations <- 50
rmse_v <- list()
rmse_v[["etree"]] <- rep(NA, n_simulations)
rmse_v[["dtree"]] <- rep(NA, n_simulations)
rmse_v[["ctree"]] <- rep(NA, n_simulations)

n.graphs <- 200
all.graphs <- list()
Y.s <- rep(NA, n.graphs)
for(j in 1:n_simulations) {
  print(paste("SIMULATION ", j))
  
  # creates 100 bipartite graphs
  for(i in 1:(n.graphs/2)) {
    g <- sample_bipartite(75, 125, p=.07)
    all.graphs[[i]] <- g
    Y.s[[i]] <- rnorm(n = 1, mean = 10,sd = 2)
  }
  # creates 100 random graphs
  for(i in (n.graphs/2 + 1):n.graphs) {
    g <- sample_gnp(directed = F, loops = F, p = 0.03, n = 200)
    all.graphs[[i]] <- g
    Y.s[[i]] <- rnorm(n = 1, mean = 30,sd = 3)
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
  
  myS  <-  mytree     ("Y", data = train, weights = NULL,
                       minbucket = 10,
                       alpha = 0.05)
  plot(myS)
  
  ###### PREDICTION
  Y.test <- test[[which(names(test) == "Y")]]
  test <- test[which(names(test) != "Y")]
  my.pred <- my.predict(model = myS, newdata = test)
  
  ###### GET THE TABULAR DATA TO BE USED BY OTHER ALGORITHMS
  train <- m.data[c(rnd.ind),]
  test <- m.data[-c(rnd.ind),]
  test <- test[which(names(test) != "Y")]
  
  ###### COMPARISON WITH DECISION TREES (CART)
  library(rpart)
  dtree.model <- rpart(Y~., data=train)
  dtree.pred <- predict(object = dtree.model, newdata = test)
  
  ###### COMPARISON WITH CONDITIONAL TREES
  ctree.model <- ctree(formula = Y~., data = train)
  ctree.pred <- predict(object = ctree.model, newdata = test)
  
  rmse_v[["etree"]][j] <-  rmse(Y.test, my.pred)
  rmse_v[["dtree"]][j] <- rmse(Y.test, dtree.pred)
  rmse_v[["ctree"]][j] <- rmse(Y.test, ctree.pred)
}

###### COMPARISON
print(paste("MRMSE Energy Tree:", round(mean(rmse_v[["etree"]]),4)))
print(paste("MRMSE Decision Tree:", round(mean(rmse_v[["dtree"]]),4)))
print(paste("MRMSE Conditional Inference Tree:", round(mean(rmse_v[["ctree"]]),4)))

print(paste("SD RMSE Energy Tree:", round(sd(rmse_v[["etree"]]),4)))
print(paste("SD RMSE Decision Tree:", round(sd(rmse_v[["dtree"]]),4)))
print(paste("SD RMSE Conditional Inference Tree:", round(sd(rmse_v[["ctree"]]),4)))
