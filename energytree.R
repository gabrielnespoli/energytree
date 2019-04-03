#### ver. "2019-03-17"

variable.importance <- function(response_name = "Y",
                                data,
                                deep = TRUE,
                                nb = 15,
                                R = 1000) {
  response <- data[[response_name]]
  datanew = data
  # extract the features from each variable of the list, otherwise
  # just check importance of the variables of the list
  if(deep) {
    if(class(data) == "list") {
      datanew <- list2matrix(
        lst = data[which(names(data) != response_name)],
        struct2copy = NULL,
        Y = NULL,
        nb = nb)
    }
  }
  else {
    n.var <- which(names(data) != response_name)
    if (class(data) == "list") {
      datanew <- list("response" = response)
      for (j in n.var) {
        if (class(data[[j]]) == "fdata") {
          foo <- min.basis(data[[j]], numbasis = nb)
          fd3 <-
            fdata2fd(foo$fdata.est,
                     type.basis = "bspline",
                     nbasis = foo$numbasis.opt)
          foo$coef <- t(fd3$coefs)
          datanew[[j]] <- foo
        }
        else if (class(data[[j]]) == "list" &
                 class(data[[j]][[1]]) == "igraph") {
          datanew[[j]] <- graph.to.shellness.distr.df(data[[j]])
        }
        else if (class(data[[j]]) == "data.frame") {
          datanew[[j]] = data[[j]]
        }
      }
      names(datanew)[-1] <- names(data)[-1]
    }
    else if (class(data) == "data.frame") {
      colnames(datanew)[
        colnames(datanew) == response_name] <- "response"
    }
    
  }
  data <- NULL
  
  n.var <- which(names(datanew) != "response")
  p = sapply(n.var, function(i) {
    if (class(datanew[[i]]) == "list") {
      mytestREG(x = datanew[[i]]$fdata.est,
                y = response,
                R = R)
    }
    
    else if (class(datanew[[i]]) == "data.frame") {
      mytestREG(x = datanew[[i]], y = response, R = R)
    }
    
    else if (class(datanew[[i]]) == "numeric" | class(datanew[[i]]) == "integer") {
      mytestREG(x = datanew[[i]], y = response, R = R)
    }
  })
  feature <- names(datanew)[which(names(datanew) != "response")]
  p_value <- p[2,]
  
  # bonferroni correction
  p_value <- 1 - (1 - p_value) ^ sum(!is.na(p_value))
  dist.corr <- p[1,]
  result <- data.frame(feature, p_value, dist.corr)
  return(result[order(result$p_value, -rank(result$dist.corr)), ])
}

#Testing equal distributions
mytestREG <- function(x, y, R = 1000) {
  library(cluster)
  library(fda.usc)
  
  if (is.factor(x)) {
    d1 = daisy(as.data.frame(x))
  }
  if (is.numeric(x)) {
    d1 = metric.dist(as.data.frame(x))
  }
  if (is.data.frame(x)) {
    d1 = metric.dist(x)
  }
  if (is.fdata(x)) {
    d1 = metric.lp(x)
  }
  y = data.frame(y)
  d2 <- metric.dist(y)
  ct <- energy::dcor.test(d1, d2, R = R)
  if (!is.na(ct$statistic)) {
    return(c(ct$statistic, ct$p.value))
  } else{
    c(NA, NA)
  }
}

##test for split nominal variables
##perform chi-squared test of yvs.x
mychisqtest <- function(x, y) {
  x <- factor(x)
  if (length(levels(x)) < 2)
    return(NA)
  ct <- suppressWarnings(chisq.test(table(y, x), correct = FALSE))
  pchisq(ct$statistic,
         ct$parameter,
         log = TRUE,
         lower.tail = FALSE)
}

# See sample()'s surprise -- example in help file
resample <- function(x, ...)
  x[sample.int(length(x), ...)]

graph.to.shellness.distr.df <- function(data, shell.limit = NULL) {
  tot.graphs = length(data)
  
  list.df <- list()
  max.shellness = 0
  
  for (i in 1:tot.graphs) {
    g = data[[i]]
    coreness.distr = count(coreness(g)) # aggr. by count
    rownames(coreness.distr) <-
      coreness.distr$x # re-index the df by the shellness number
    
    # keep just the frequency column
    coreness.distr = coreness.distr[c('freq')]
    
    # transpose the df. Convert the column-df into row-df. 
    #This will ease the join with df.shellness.distr
    coreness.distr = t(coreness.distr)
    list.df[[i]] <- coreness.distr
    
    this.max.shellness = colnames(coreness.distr)[
      ncol(coreness.distr)]
    
    # update the maximum shellness found so far (used to build 
    #the df of shellness distr)
    if (this.max.shellness > max.shellness) {
      max.shellness = this.max.shellness
    }
  }
  
  if (!is.null(shell.limit)) {
    # calculates the max shellness between the number of 
    # predictors used in train set and the one calculated in 
    # test set
    max.shellness <-
      if (as.numeric(max.shellness) < shell.limit - 1)
        shell.limit - 1
    else
      max.shellness
  }
  
  col.names = seq(0, max.shellness, 1)
  col.names = lapply(col.names, function(x)
    as.character(x))  # convert to char
  df.shellness.distr = data.frame(matrix(
    data = NA_integer_,
    nrow = tot.graphs,
    ncol = length(col.names)
  )) #df with all graphs shellness distribution
  colnames(df.shellness.distr) <- col.names
  
  # fill in the df with the shellness distribution of each graph
  for (i in 1:tot.graphs) {
    updated.cols = colnames(list.df[[i]])
    
    for (x in updated.cols) {
      df.shellness.distr[i, x] = list.df[[i]][, x]
    }
  }
  
  df.shellness.distr[is.na(df.shellness.distr)] <-
    0 # replace NA by 0
  
  # converted the df columns to integer
  df.shellness.distr[, seq(1, ncol(df.shellness.distr))] <-
    sapply(df.shellness.distr[, seq(1, ncol(df.shellness.distr))], 
           as.integer)
  
  return(df.shellness.distr)
}

mytree <- function(Y,
                   data,
                   weights = NULL,
                   minbucket = 1,
                   alpha = 0.05,
                   R = 1000,
                   rnd.sel = T,
                   rnd.splt = TRUE,
                   nb = 5) {
  # name of the response variable
  response <- data[[which(names(data) == Y)]]
  
  if (is.null(weights))
    weights <- rep(1L, length(response))
  
  n.var <- which(names(data) != Y)
  if (class(data) == "list") {
    datanew <- list("response" = response)
    for (j in n.var) {
      if (class(data[[j]]) == "fdata") {
        foo <- min.basis(data[[j]], numbasis = nb)
        fd3 <-
          fdata2fd(foo$fdata.est,
                   type.basis = "bspline",
                   nbasis = foo$numbasis.opt)
        foo$coef <- t(fd3$coefs)
        datanew[[j]] <- foo
      }
      else if (class(data[[j]]) == "list" &
               class(data[[j]][[1]]) == "igraph") {
        datanew[[j]] <- graph.to.shellness.distr.df(data[[j]])
      }
      else if (class(data[[j]]) == "data.frame") {
        datanew[[j]] = data[[j]]
      }
    }
    names(datanew)[-1] <- names(data)[-1]
  }
  else if (class(data) == "data.frame") {
    datanew = data
    colnames(datanew)[colnames(datanew) == "Y"] <- "response"
  }
  
  nodes <-
    growtree(
      id = 1L,
      response = datanew$response,
      data = datanew,
      weights,
      minbucket = minbucket,
      alpha = alpha,
      R = R,
      rnd.sel = rnd.sel,
      rnd.splt = rnd.splt,
      n.var = n.var
    )
  
  # compute terminal node number for each observation
  response <- response
  response <- data.frame(response)
  y = response
  m.data <- c()
  
  for (j in n.var) {
    if (class(data[[j]]) == "fdata") {
      foo <- datanew[[j]]$coef
      colnames(foo) <-
        paste(names(data)[j], colnames(datanew[[j]]$coef), 
              sep = ".")
    }
    else if (class(data[[j]]) == "data.frame" |
             (class(data[[j]]) == "list" &
              class(data[[j]][[1]]) == "igraph")) {
      foo <- datanew[[j]]
      colnames(foo) <-
        paste(names(data)[j], colnames(datanew[[j]]), sep = ".")
    }
    
    # fill in m.data or initialize it
    if (!is.null(m.data)) {
      m.data <- cbind(m.data, foo)
    }
    else {
      m.data <- foo
    }
  }
  
  data1 = cbind(response, m.data)
  m.data = m.data
  data1 = data1
  
  fitted <- fitted_node(nodes, data = data.frame(m.data))
  formula = response ~ .
  
  # return rich constparty object
  ret <- party(
    nodes,
    data = data.frame(m.data),
    fitted = data.frame(
      "(fitted)" = fitted,
      "(response)" = data1$response,
      "(weights)" = weights,
      check.names = FALSE
    ),
    terms = terms(formula, data = data1)
  )
  
  as.constparty(ret)
}

####

growtree <- function(id = 1L,
                     response,
                     data,
                     weights,
                     minbucket,
                     alpha,
                     R,
                     rnd.sel,
                     rnd.splt,
                     n.var) {
  # for less than <minbucket> observations stop here (for ctree()
  # is 7 in ?ctree_control)
  if (sum(weights) <= minbucket) {
    return(partynode(id = id))
  }
  
  # find best split
  res <- findsplit(
    response,
    data,
    weights,
    alpha = alpha,
    R = R,
    rnd.sel = rnd.sel,
    rnd.splt = rnd.splt,
    n.var = n.var
  )
  
  sp <- res$sp
  varselect <- res$varselect
  
  # no split found, stop here
  if (is.null(sp)) {
    return(partynode(id = id))
  }
  
  kidids <- c()
  
  if (class(data[[varselect]]) == "list") {
    kidids[which(data[[
      varselect]]$coef[, sp$varid] <= sp$breaks)] <- 1
    kidids[which(data[[
      varselect]]$coef[, sp$varid] > sp$breaks)] <- 2
    
    sum1 <- length(which(data[[varselect]]$coef[
      which(weights == 1), sp$varid] <= sp$breaks))
    sum2 <-
      length(which(data[[varselect]]$coef[
        which(weights == 1), sp$varid] > sp$breaks))
    
    nb = 0
    for (i in n.var) {
      k = data[[i]]$numbasis.opt
      nb = c(nb, k)
    }
    
    # shift the varid of the tree based on the quantity of the 
    # previous features/basis
    # Ex: if variable 3 is selected for splitting (variable 1 
    # is the response, it's ignored), then shift varid by the 
    # number of basis of variable 2 (if it's functional) or the 
    # maximum k_core found in the graphs (if it's a graph)
    if (varselect != min(n.var)) {
      total_features <- c()
      lapply(n.var, function(v) {
        if (class(data[[v]]) == 'list')
          total_features[[v]] <<- data[[v]]$numbasis.opt
        if (class(data[[v]]) == 'data.frame')
          total_features[[v]] <<- ncol(data[[v]])
        if (class(data[[v]]) == 'numeric')
          total_features[[v]] <<- 1
      })
      step <-
        sum(total_features[n.var[
          which(n.var < varselect)]], na.rm = T)
      sp$varid = sp$varid + as.integer(step)
    }
  }
  else if (class(data[[varselect]]) == "data.frame") {
    kidids[(which(data[[
      varselect]][, sp$varid] <= sp$breaks))] <- 1
    kidids[(which(data[[
      varselect]][, sp$varid] > sp$breaks))] <- 2
    
    sum1 <-
      length(which(data[[
        varselect]][, sp$varid][which(weights == 1)] <= sp$breaks))
    sum2 <-
      length(which(data[[
        varselect]][, sp$varid][which(weights == 1)] > sp$breaks))
  }
  else if (class(data[[varselect]]) == "numeric") {
    kidids[(which(data[[varselect]] <= sp$breaks))] <- 1
    kidids[(which(data[[varselect]] > sp$breaks))] <- 2
    
    sum1 <-
      length(which(data[[
        varselect]][which(weights == 1)] <= sp$breaks))
    sum2 <-
      length(which(data[[
        varselect]][which(weights == 1)] > sp$breaks))
  }
  
  if (all(kidids == 1) | all(kidids == 2))
    return(partynode(id = id))
  
  if ((sum1 == 0 | sum2 == 0)) {
    return(partynode(id = id))
  }
  
  # setup all daugther nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))
  
  for (kidid in 1:length(kids)) {
    # select observations for current node
    w <- weights
    w[kidids != kidid] <- 0
    
    # get next node id
    if (kidid > 1) {
      myid <- max(nodeids(kids[[kidid - 1]]))
    } else{
      myid <- id
    }
    
    # Start recursion on this child node
    kids[[kidid]] <-
      growtree(
        id = as.integer(myid + 1),
        response,
        data,
        w,
        minbucket,
        alpha,
        R,
        rnd.sel,
        rnd.splt ,
        n.var = n.var
      )
  }
  
  # return nodes
  return(partynode(
    id = as.integer(id),
    split = sp,
    kids = kids,
    info = list(p.value = 
                  min(info_split(sp)$p.value, na.rm = TRUE))
  ))
}

findsplit <- function(response,
                      data,
                      weights,
                      alpha,
                      R,
                      rnd.sel,
                      rnd.splt,
                      n.var) {
  require(partykit)
  
  # extract response values from data
  y <- response[which(weights == 1)]
  
  p <- matrix(NA, nrow = length(n.var), ncol = 2)
  colnames(p) <- c("statistic", "p-value")
  p = sapply(n.var, function(i) {
    if (class(data[[i]]) == "list") {
      mytestREG(x = data[[i]]$fdata.est[which(weights == 1)],
                y = y,
                R = R)
    }
    
    else if (class(data[[i]]) == "data.frame") {
      mytestREG(x = data[[i]][which(weights == 1),], y = y, R = R)
    }
    
    else if (class(data[[i]]) == "numeric") {
      mytestREG(x = data[[i]][which(weights == 1)], y = y, R = R)
    }
  })
  
  # Bonferroni-adjusted p-value small enough?
  if (all(is.na(p[2,])))
    return(NULL)
  
  minp <- min(p[2,], na.rm = TRUE)
  minp <- 1 - (1 - minp) ^ sum(!is.na(p[2,]))
  if (minp > alpha)
    return(NULL)
  
  xselect <- n.var
  if (length(which(p[2,] == min(p[2,], na.rm = T))) > 1) {
    xselect <- which.max(p[1,]) + 1
  } else{
    xselect <- which.min(p[2,]) + 1
  }
  
  x <-  data[[xselect]]
  if (is.list(x)) {
    if (is.fdata(x$fdata.est))
      x1 = x$coef[which(weights == 1),]
  }
  
  # split into two groups minimizing entropy
  if (is.factor(x)) {
    # setup all possible splits in two kid nodes
    lev <- levels(x[drop = TRUE])
    if (length(lev) == 2) {
      splitpoint <- lev[1]
    } else{
      comb <- do.call("c", lapply(1:(length(lev) - 2),
                                  function(x)
                                    combn(lev, 
                                          x, 
                                          simplify = FALSE)))
      xlogp <- sapply(comb, function(q)
        mychisqtest(x %in% q, y))
      splitpoint <- comb[[which.min(xlogp)]]
    }
    
    # split into two groups (setting groups that do not occur 
    # to NA)
    splitindex <- !(levels(data[[xselect]]) %in% splitpoint)
    splitindex[!(levels(data[[xselect]]) %in% lev)] <- NA_integer_
    splitindex <- splitindex - min(splitindex, na.rm = TRUE) + 1L
  }
  if (is.data.frame(data[[xselect]])) {
    # select the column
    all.columns = seq(1, ncol(data[[xselect]]))
    # iterate over each column of the dataframe
    p = sapply(all.columns, function(i) {
      mytestREG(x = data[[
        xselect]][[i]][which(weights == 1)], y = y, R = R)
    })
    
    minp <- min(p[2,], na.rm = TRUE)
    minp <- 1 - (1 - minp) ^ sum(!is.na(p[2,]))
    
    if (minp > alpha)
      return(NULL)
    
    # select the column to split
    cselect <- NA
    if (length(which(p[2,] == min(p[2,], na.rm = T))) > 1) {
      cselect <- which.max(p[1,])
    } else{
      cselect <- which.min(p[2,])
    }
    
    splitindex <-
      s.opt(response, x[which(weights == 1), cselect], rnd.splt)
  }
  if (is.numeric(x)) {
    splitindex <- s.opt(response, x[which(weights == 1)], rnd.splt)
  }
  if (is.list(x)) {
    if (is.fdata(x$fdata.est)) {
      bselect <- 1:dim(x1)[2]
      p1 <- c()
      p1        <-
        sapply(bselect, function(i)
          mytestREG(x1[, i], y, R = R))
      colnames(p1) <- colnames(x1)
      if (length(which(p1[2,] == min(p1[2,], na.rm = T))) > 1) {
        bselect <- which.max(p1[1,])
      } else{
        bselect <- which.min(p1[2,])
      }
      
      splitindex <- s.opt(y = y, X = x1[, bselect], rnd.splt)
    }
  }
  
  # return split as partysplit object
  if (is.numeric(x)) {
    return(partysplit(
      varid = as.integer(xselect),
      breaks = splitindex,
      info = list(p.value = 1 - (1 - p) ^ sum(!is.na(p)))
    ))
  }
  if (is.data.frame(x)) {
    temp = list (
      sp = partysplit(
        varid = as.integer(cselect),
        breaks = splitindex,
        info = list(p.value = 1 - (1 - p) ^ sum(!is.na(p)))
      ),
      varselect = xselect
    )
    return(temp)
  }
  if (is.factor(x)) {
    return(partysplit(
      varid = as.integer(xselect),
      index = splitindex,
      info = list(p.value = 1 - (1 - p) ^ sum(!is.na(p)))
    ))
  }
  if (is.list(x)) {
    if (is.fdata(x$fdata.est)) {
      return(list(
        sp = partysplit(
          varid = as.integer(bselect),
          breaks = splitindex,
          info = list(p.value = 1 - (1 - p[2,]) ^
                        sum(!is.na(p[2,])))
        ),
        varselect = xselect
      ))
    }
  }
}

s.opt <- function(y, X, rnd = T) {
  # find the split minimizing variance
  s  <- sort(X)
  obj <- c()
  for (i in 1:length(s)) {
    data1  <- y[which(X < s[i])]
    data2  <- y[which(X >= s[i])]
    v1 = var(data1)
    v2 = var(data2)
    n1 = length(data1)
    n2 = length(data2)
    n = n1 + n2
    
    obj[i] = (n1 * v1 + n2 * v2) / n
  }
  
  if (all(is.na(obj)))
  {
    splitindex <- length(obj)
  }
  else {
    splitindex <- s[which.min(obj)]
  }
  
  return(splitindex)
}

# convert a list of graphs and/or functions to matrix
list2matrix <- function(lst,
                        struct2copy = NULL,
                        Y = NULL,
                        nb = NULL) {
  m.data <- c()
  if (!is.null(Y)) {
    n.var <- which(names(lst) != names(Y))
  }
  else {
    n.var <- seq(1, length(lst))
  }
  
  for (j in n.var) {
    if (class(lst[[j]]) == "fdata") {
      foo <- min.basis(lst[[j]], numbasis = nb)
      fd3 <-
        fdata2fd(foo$fdata.est,
                 type.basis = "bspline",
                 nbasis = foo$numbasis.opt)
      foo$coef <- t(fd3$coefs)
      foo <- foo$coef
      colnames(foo) <-
        paste(names(lst)[j], colnames(foo), sep = ".")
    }
    else if (class(lst[[j]]) == "data.frame" |
             (class(lst[[j]]) == "list" &
              class(lst[[j]][[1]]) == "igraph")) {
      if (!is.null(struct2copy)) {
        graph.data.train <-
          struct2copy[, grep(paste("V", j, sep = ""), 
                             names(struct2copy))]
        foo <-
          graph.to.shellness.distr.df(lst[[j]], 
                                      ncol(graph.data.train))
        foo.names <- names(graph.data.train)
      }
      else {
        foo <- graph.to.shellness.distr.df(lst[[j]], NULL)
        foo.names <- names(foo)
        # concatenates the variable prefix
        foo.names <- paste(names(lst)[j], foo.names, sep = ".")
      }
      colnames(foo) <- foo.names
    }
    
    # fill in m.data or initialize it
    if (!is.null(m.data)) {
      m.data <- cbind(m.data, foo)
    }
    else {
      m.data <- foo
    }
  }
  if (!is.null(Y)) {
    m.data <- cbind(m.data, Y)
  }
  return(m.data)
}

my.predict <- function(Y = "Y",
                       model,
                       newdata,
                       nb = NULL) {
  m.data <- c()
  if (class(newdata) == "list") {
    m.data <- list2matrix(newdata, 
                          struct2copy = model$data, 
                          Y = NULL, 
                          nb)
  }
  else if (class(newdata) == "data.frame") {
    m.data = newdata
  }
  
  response.pred <-
    partykit::predict.party(model, newdata = m.data)
  return(response.pred)
}