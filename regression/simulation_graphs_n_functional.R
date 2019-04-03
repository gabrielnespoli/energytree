source("energytree.R")
source("gen_data_REG.R")
require(igraph)
require(reshape)
require(plyr)
require(roahd)
library(fda.usc)
library(energy)
library(entropy)
library(partykit)
require(Metrics)
require(pracma)

set.seed(100)
n.graphs <- 200
Y.s <- rep(NA, n.graphs)
all.graphs <- list()

n_simulations <- 50
rmse_v <- list()
rmse_v[["etree"]] <- rep(NA, n_simulations)
rmse_v[["dtree"]] <- rep(NA, n_simulations)
rmse_v[["ctree"]] <- rep(NA, n_simulations)
for(j in 1:n_simulations) {
  print(paste("SIMULATION ", j)) 
  
  # creates 200 random graphs
  for(i in 1:n.graphs) {
    g <- sample_gnp(directed = F, loops = F, p = 0.03, n = 200)
    all.graphs[[i]] <- g
    
    coreness.distr = count(coreness(g)) # aggr. by count
    rownames(coreness.distr) <- coreness.distr$x # re-index the df by the shellness number
    coreness.distr = coreness.distr[c('freq')] # keep just the frequency column
    coreness.distr = t(coreness.distr) # transpose the df. Convert the column-df into row-df. This will ease the join with df.shellness.distr
    
    # sum of multiplication of the row number by the value in the row
    Y.s[[i]] <- dot(coreness.distr[1,], seq(0,ncol(coreness.distr)-1))
  }
  
  # creates functional dataset
  size <- c(66,67,67) # 200 datapoints
  P <- 1e2
  n.var <- 1
  data=list()
  y_mean=c(10,30,50)
  y_sd=c(2,3,2.5)
  alpha <- matrix( round(runif(3,0.1,1),2),
                   nrow=3, ncol=1)
  beta <- matrix( round(runif(3,0.1,1),2),
                  nrow=3, ncol=1)
  b <- matrix( sample(c(1,1.5,2,2.5,3,3.5,4),size=3, replace=T),
               nrow=3, ncol=1)
  a <- matrix( sample(c(0,1,2,3),size=3, replace=T),
               nrow=3, ncol=1)
  func.data = gen_data_reg(y_mean, y_sd ,size=size,P=P,n.var=n.var,alpha=alpha,beta=beta,
                           a=a,b=b)
  
  all.graphs <- sample(all.graphs)
  
  data <- vector("list", 3)
  names(data) <- c("Y", "V1", "V2")
  
  # data$Y[1:100] <- func.data$Y[1:100]
  # data$Y[101:200] <- runif(100,min = 8, max = 35)
  data$Y[1:200] <- func.data$Y[1:200] + Y.s
  data$V1 <- all.graphs
  data$V2 <- func.data$V1 # selects just the covariates of the functional data (not Y)
  
  # this line is used to generate the data to me used in the comparison with other algorithms
  # it's not used by ours
  nb <- 15
  m.data <- list2matrix(data[-1], NULL, list("Y" = data$Y), nb=nb)
  
  ###### split the data into train and test sets
  n = length(data[[1]])
  rnd.ind = sample(x = seq(n), size = 0.8*n)
  train = test = data
  
  # convert the graph to its respective shellness distribution dataframe
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
    else if(class(train[[i]])=="fdata") {
      train[[i]] = train[[i]][c(rnd.ind),]
      test[[i]] = test[[i]][-c(rnd.ind),] 
    }
  }
  
  myS  <-  mytree     ("Y", data = train, weights = NULL,
                       minbucket = 3,
                       alpha = 0.2,
                       R = 1000, 
                       rnd.sel=T,
                       rnd.splt=TRUE,
                       nb=nb)
  
  plot(myS)
  
  ###### PREDICTION
  Y.test <- test[[which(names(test) == "Y")]]
  test <- test[which(names(test) != "Y")]
  my.pred <- my.predict(model = myS, newdata = test, nb = nb)
  
  ###### GET THE TABULAR DATA TO BE USED BY OTHER ALGORITHMS
  train <- m.data[c(rnd.ind),]
  test <- m.data[-c(rnd.ind),]
  test <- test[which(names(test) != "Y")]
  
  ###### COMPARISON WITH DECISION TREES (CART)
  library(rpart)
  dtree.model <- rpart(Y~., data=train)
  dtree.pred <- predict(object = dtree.model, newdata = test)
  
  ###### COMPARISON WITH CONDITIONAL TREES
  ctree.model <- ctree(Y~., data=train)
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