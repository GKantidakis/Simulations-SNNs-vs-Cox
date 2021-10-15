
# A simulation study to compare the predictive performance of survival neural networks
# with Cox models for clinical trial data

# R code for the simulations


######################################################################
######################################################################
# file "the_functions2.R" (necessary to call it when tuning the PLANN
# models or running the simulations)

# contains functions for the estimation of the predictive performance 
# (discrimination, calibration) for all methods. 
######################################################################
######################################################################

## start of the file

# To install packages
# install_packages <- c("survival", "pec", "caret", "e1071")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("survcomp")


# libraries to call
library(survival)
library(pec)
library(caret)
library(e1071)
library(survcomp)

# In this file you will find functions to calculate the Brier and IBS for Cox models and SNN
# the prognostic index of a survival neural network and variable importance


# A function that calculates the quantiles, the mean and the standard deviation
summa <- function(x) {
  
  res1 <- round(quantile(x, probs = c(0, 0.025, 0.5, 0.975, 1), na.rm = TRUE), digits = 3)
  res2 <- round(mean(x, na.rm = TRUE), digits = 3)
  res3 <- round(sd(x, na.rm = TRUE), digits = 3)
  res <- c(res1[1:3], res2, res1[4:5], res3)
  names(res) <- c("Min.", "2.5% Qu.", "Median", "Mean", "97.5% Qu.", "Max.", "Sd.")
  return(res)
}

# Calculates Brier score at specific time point using the risk regression R package

brier_general_new <- function(risk_matrix, data, times) {
  
  f_mod <- Surv(survt_ovr, surv_status) ~ 1
  
  PredError <- riskRegression:::Score(object = list("Cox model" = risk_matrix),
                                      formula = f_mod,
                                      data = data,
                                      metrics = "Brier",
                                      times = times,
                                      null.model = FALSE,
                                      summary = "ibs",
                                      conservative = TRUE, # speed up computation
                                      conf.int = FALSE,
                                      split.method = "none")
  
  # note: we use the IBS at 5 years
  return(list(Time = times, Brier = as.numeric(unlist(PredError$Brier$score[, 3])),
              Int_brier = as.numeric(unlist(PredError$Brier$score[, 4]))[5+1]))
  
}

# EXAMPLE of how to apply the function on training and test data
# cox_all <- coxph(Surv(survt_ovr, surv_status) ~ ., data = training_ovr,
#                  method = "breslow", x = TRUE, y = TRUE)
# times <- c(0, 1, 2, 3, 4, 5, 6, 7, 7.9973)
# risk_cox_all <- riskRegression::predictRisk(cox_all, test_ovr,
#                                              times = times)
# brier_general_new(risk_matrix = risk_cox_all, data = test_ovr,
#                   times = times)



# calculate Brier and Integrated Brier Score (IBS) for the neural networks
# prob matrix is of the form nrow(data) x 6
brier_nn <- function(risk_matrix, data) {
  
  f_nn <- Surv(survt_ovr, surv_status) ~ 1
  # we calculate the Brier score utilizing the pec package
  # we subtract 1 day at 8 years to avoid zero censoring prob weight at this time  
  bsc <- riskRegression:::Score(object = list("Survival Neural Network" = risk_matrix),
                                formula = f_nn,
                                data = data,
                                metrics = "Brier",
                                times = c(0, 1, 2, 3, 4, 4.999),
                                null.model = FALSE,
                                summary = "ibs",
                                conservative = TRUE, # speed up computation
                                conf.int = FALSE,
                                split.method = "none")
  
  # note we use the IBS at 5 years
  return(list(Time = seq(0, 5, length.out = 6), 
              Brier = as.numeric(unlist(bsc$Brier$score[, 3])),
              Int_brier = as.numeric(unlist(bsc$Brier$score[, 4]))[5+1]))
}


# This function calculates the (non-linear) prognostic index

pi_calculator <- function(trained_model, datanew) {
  
  df1 <- data.frame(hazard = predict_proba(trained_model,
                                           as.matrix(datanew[, c(1:5, 12:(11 + max(datanew$interval)))]),
                                           batch_size = 64))
  df1$id <- datanew$id # ids of the patients
  df1$survival <- datanew$survival # survival time in years
  groups <- split(df1, f = df1$id) 
  
  haz_mat <- NULL
  for (id in 1:length(groups)) {
    
    haz_mat <- rbind(haz_mat, groups[[id]]$hazard)
    
  }
  
  pi <- vector(mode = "list", length = ncol(haz_mat))
  
  for (t in 1:length(pi)) {
    
    pi[[t]] <- log(haz_mat[, t] / (1 - haz_mat[, t]), base = exp(1))
    pi[[t]] <- as.vector(scale(pi[[t]], center = TRUE, scale = FALSE)) # centered prognostic index
  }
  
  pi_person <- vector(mode = "numeric", length = ncol(haz_mat))
  pi_total <- vector(mode = "numeric", length = nrow(haz_mat))
  
  
  for (i in 1:nrow(haz_mat)) {
    for (t in 1:length(pi)) {
      
      pi_person[t] <- pi[[t]][i]    
      pi_total[i] <- mean(pi_person, na.rm = TRUE)  
      # if pi total was not calculated set to zero
      if (is.nan(pi_total[i]) | is.infinite(pi_total[i]) ) pi_total[i] <- 0
    }
    
  }
  
  return(pi_total)
  
}



# adjust pi_calculator for library nnet
pi_calculator2 <- function(trained_model, datanew) {
  
  df1 <- data.frame(hazard = predict(trained_model,
                                     as.matrix(datanew[, c(1:5, 12 + max(datanew$interval))])))
  df1$id <- datanew$id # ids of the patients
  df1$survival <- datanew$survival # survival time in years
  groups <- split(df1, f = df1$id) 
  
  haz_mat <- NULL
  for (id in 1:length(groups)) {
    
    haz_mat <- rbind(haz_mat, groups[[id]]$hazard)
    
  }
  
  pi <- vector(mode = "list", length = ncol(haz_mat))
  
  for (t in 1:length(pi)) {
    
    pi[[t]] <- log(haz_mat[, t] / (1 - haz_mat[, t]), base = exp(1))
    pi[[t]] <- as.vector(scale(pi[[t]], center = TRUE, scale = TRUE)) # centered prognostic index
  }
  
  pi_person <- vector(mode = "numeric", length = ncol(haz_mat))
  pi_total <- vector(mode = "numeric", length = nrow(haz_mat))
  
  
  for (i in 1:nrow(haz_mat)) {
    for (t in 1:length(pi)) {
      
      pi_person[t] <- pi[[t]][i]    
      pi_total[i] <- mean(pi_person, na.rm = TRUE)
      # if pi total was not calculated set to zero
      if (is.nan(pi_total[i]) | is.infinite(pi_total[i]) ) pi_total[i] <- 0
    }
    
  }
  
  return(pi_total)
  
}

# this is connection weight method performed in keras

var_imp_keras <- function(model){
  
  list_weights <- get_weights(model)
  # list_weights[[1]] are the weights between input-hidden node
  # list_weights[[2]] are the weights between bias-hidden node
  # list_weights[[3]] are the weights between hidden node-output
  # list_weights[[4]] is the weight between bias-output
  # names input nodes, number input nodes
  
  # n_input <- length(list_weights[[1]]) / length(list_weights[[2]])
  # number hidden nodes , number output nodes
  # n_nodes <- length(list_weights[[3]])
  # n_outputs <- length(list_weights[[4]])
  
  # matrix multiplication input x hidden with hidden x output
  mega_mat <- t(list_weights[[1]] %*% list_weights[[3]])
  
  names <- c("trtReg_DI", "sexmale", "hist_respgood", "excisioncomplete",
             "age", paste0("interval_", 1:(ncol(mega_mat) - 5)))
  
  colnames(mega_mat) <- names
  mega_mat_abs <- abs(mega_mat)
  totals <- sum(mega_mat_abs)
  mega_mat_rel <- as.data.frame(mega_mat_abs/ totals)
  rels_garson <- as.vector(as.numeric(mega_mat_rel)) # relative importance
  rels_olden <- as.vector(mega_mat) # raw connection weights
  
  return(list(Names = names, Rels_Garson = round(rels_garson, 3), 
              Rels_Olden = round(rels_olden, 3)))
}


# variable importance for library nnet
var_imp_nnet <- function(model){
  
  # number of input nodes
  n_inputs <- model$n[1]
  # number of hidden nodes
  n_hidden <- model$n[2]
  # number of output nodes
  n_output <- model$n[3]
  
  # dimensions for weight matrices
  ncol1 <- n_inputs + 1
  ncol2 <- n_hidden + 1
  nrows1 <- n_hidden
  nrows2 <- n_output
  length1 <- ncol1*nrows1
  length2 <- ncol2*nrows2
  
  # selecting weights
  weights <- model$wts
  weights1 <- as.data.frame(matrix(weights[1:length1],
                                   ncol = ncol1, nrow = nrows1,
                                   byrow=TRUE))
  
  names <- c("trtReg_DI", "sexmale", "hist_respgood", "excisioncomplete",
             "age", "interval")
  colnames(weights1) <- c("bias", names)
  
  weights2 <- t(matrix(weights[(length1 + 1):(length1 + length2)],
                       nrow = nrows2, ncol = ncol2 , byrow = TRUE))
  rownames(weights2) <-  c ( "bias", paste("n", seq(1:n_hidden), sep=""))
  
  # calculating variable contribution
  # to do so we use the weights without including the bias nodes
  mega_mat <-  matrix(0, ncol = n_inputs, nrow = n_output)
  for (i in 1:n_inputs) {
    for (j in 1:n_output){
      mega_mat[j, i] <-  sum(weights1[, i+1] * weights2[2:(n_hidden+1), j ])
    }
  }
  
  colnames(mega_mat) <- names
  mega_mat_abs <- abs(mega_mat)
  totals <- sum(mega_mat_abs)
  mega_mat_rel <- as.data.frame(mega_mat_abs/ totals)
  rels_garson <- as.vector(as.numeric(mega_mat_rel)) # relative importance
  rels_olden <- as.vector(mega_mat) # raw connection weights
  
  return(list(Names = names, Rels_Garson = round(rels_garson, 3), 
              Rels_Olden = round(rels_olden, 3)))
}



# calibration statistics
calibration_stats <- function(prob_vector, data, times = 2, cutpoints = 4) {
  
  if (var(prob_vector, na.rm = TRUE) != 0) { # run if not all probabilities are equal
    
    cuts <- unique(quantile(prob_vector, seq(0, 1, length = cutpoints + 1), na.rm = TRUE))
    
    
    data$group <-  as.factor(as.numeric(cut(x = prob_vector, breaks = as.vector(cuts),
                                            include.lowest = TRUE)))
    
    average_probs <- vector(mode = "numeric", length = length(cuts) - 1)
    for (i in 1:(length(cuts) - 1)) {
      average_probs[i] <- mean(prob_vector[which(as.numeric(data$group) == i)], na.rm = TRUE)
    }
    
    km <- summary(survfit(Surv(survt_ovr, surv_status) ~ group,
                          data = data, conf.type = "log-log"),
                  times = times, extend = TRUE)
    km_surv <- km$surv
    
    mse1 <- sum((km_surv - average_probs)^2, na.rm = TRUE) / length(average_probs)
    
  } else { # else if they are equal set to NA
    
    mse1 <- NA  
    
  }
  
  
  # return the calculated mean squared error
  return(list(mse_groups = mse1))
  
}


## end of the file



######################################################################
######################################################################
# file "functions_nn.R" (necessary to call it when tuning the PLANN
# models or running the simulations)

# contains functions for data pre-processing and measure 
# calculation of PLANNs
######################################################################
######################################################################

## start of the file

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("survcomp")

library(survcomp)
###############################################################

data_train_creator <- function(data){
  
  N <- nrow(data)
  # assign survival times to 8 intervals, each is a 1-year period
  data$interval <- cut(data$survt_ovr, breaks = c(0:7, max(data$survt_ovr, 8)), labels = FALSE)
  data$survival <- cut(data$survt_ovr, breaks = c(0:7, max(data$survt_ovr, 8)), labels = FALSE)
  data$id <- 1:N
  n.times <- data$interval
  data_long <- data[rep(seq_len(N), times = n.times), ]
  
  # create the correct intervals
  for(i in unique(data_long$id)) {
    n_length <- length(data_long$interval[data_long$id == i])
    data_long$interval[data_long$id == i] <- 1:n_length
  }
  
  data_long$status_long <- vector(mode = "numeric",
                                  length = nrow(data_long))
  
  # put indication 1 on surv_status at the interval that patient dies
  for (i in 1:nrow(data_long)) {
    if (data_long$surv_status[i] != data_long$status_long[i] &&
        data_long$survival[i] == data_long$interval[i])
      data_long$status_long[i] <- 1
  }
  
  intervals <- dummy_cols(as.factor(data_long$interval))
  colnames(intervals) <- gsub(".data", "interval", colnames(intervals))
  data_long <- data.frame(data_long, intervals[, seq(2, max(n.times) + 1, by = 1)])
  data_long$interval_scaled <- scale(data_long$interval)
  return(data_long)
  
}

# function that creates data in the right long format for
# test set
data_test_creator <- function(data, enforce_intervals = FALSE){
  
  N <- nrow(data)
  # assign survival times to 8 intervals
  data$interval <- max(cut(data$survt_ovr,
                           breaks = c(0:7, max(data$survt_ovr, 8)),
                           labels = FALSE)) # replicate the test data 
  
  
  if (enforce_intervals == TRUE) {
    data$interval <- rep(5, times = length(data$interval))
  } 
  
  # the true interval survival
  data$survival <- cut(data$survt_ovr,
                       breaks = c(0:7, max(data$survt_ovr, 8)),
                       labels = FALSE) 
  
  
  data$id <- 1001:(1000 + N) # define the patient ids abstractly
  n.times <- data$interval
  data_long <- data[rep(seq_len(N), times = n.times), ]
  
  # create the correct intervals
  for(i in unique(data_long$id)) {
    n_length <- length(data_long$interval[data_long$id == i])
    data_long$interval[data_long$id == i] <- 1:n_length
    
  }
  
  data_long$status_long <- vector(mode = "numeric",
                                  length = nrow(data_long))
  
  # put indication 1 on surv_status at the intervals on
  # which a patient has died
  for (i in 1:nrow(data_long)) {
    if (data_long$surv_status[i] == 1 && 
        data_long$survival[i] <= data_long$interval[i]) 
      data_long$status_long[i] <- 1
  }
  
  if (enforce_intervals == TRUE) {
    data_long$survival <- ifelse(data_long$survival > 5, 5, data_long$survival) # administrative survival 5 years
  } 
  
  intervals2 <- dummy_cols(as.factor(data_long$interval))
  colnames(intervals2) <- gsub(".data", "interval", colnames(intervals2))
  data_long <- data.frame(data_long, intervals2[, seq(2, max(n.times) + 1, by = 1)])
  data_long$interval_scaled <- scale(data_long$interval)
  return(data_long)
  
}


# function that calculates the metrics for a neural network from keras  
measures_calculator <- function(trained_model,
                                datashort,
                                datalong) {
  
  real_times <- datashort$survt_ovr
  real_status <- datashort$surv_status
  
  df1 <- data.frame(hazard = predict_proba(trained_model,
                                           as.matrix(datalong[, c(1:5, 12:(11 + max(datalong$interval)))]),
                                           batch_size = 64))
  df1$id <- datalong$id # ids of the patients
  df1$survival <- datalong$survival # survival time in years
  groups <- split(df1, f = df1$id) 
  
  true_surv <- unlist(lapply(groups, function(x) {
    surv_obj <- x$survival
    true_res <- surv_obj[1] 
    return(true_res)}
  ))
  group_probs <- lapply(groups, function(x) {
    x <- cumprod(1 - x$hazard)})
  pred_mat <- do.call("rbind", group_probs)
  
  relative_probs <- vector(mode = "numeric",
                           length = length(group_probs))
  for (i in 1:length(group_probs)){ 
    temp <- group_probs[[i]]
    ind <- which(1:max(datalong$interval) == true_surv[i])
    relative_probs[i] <- temp[ind] 
  }
  N0 <- length(unique(df1$id)) # number of unique persons
  
  # in the data frame create random id numbers 
  # to label the patients
  df2 <- data.frame(relative_probs = relative_probs,
                    survt_ovr = real_times,
                    surv_status = real_status, 
                    id = (1001):(1000 + N0))
  df2$prediction <- 1 - round(df2$relative_probs, digits = 0)
  
  prob_matrix <-  cbind(1, pred_mat[, 1:5]) # prediction matrix up to 5 years
  risk_mat <- 1 - prob_matrix
  
  brier_obj <- brier_nn(risk_matrix = risk_mat, data = df2)
  int_brier <- as.numeric(brier_obj$Int_brier)
  brier_set <- brier_obj$Brier
  
  
  # estimate c-index
  pi_total <- pi_calculator(trained_model, datalong)
  cindex <- concordance.index(x = pi_total,
                              surv.time = real_times,
                              surv.event = real_status)$c.index
  # assign 0.5 to c-indexes when prognostic indexes are the same 
  if (is.na(cindex)) cindex <- 0.5 
  
  all_weights <- get_weights(trained_model)
  nr_weights <- length(unlist(all_weights))
  
  # calculate variable importance
  var_imp <- var_imp_keras(model = trained_model)
  
  # calculate calibration statistics at 2 years
  cal_stats2 <- calibration_stats(prob_vector = prob_matrix[, 2+1], data = datashort, times = 2)
  cal_stats5 <- calibration_stats(prob_vector = prob_matrix[, 5+1], data = datashort, times = 5)
  
  return(list(weights = nr_weights, 
              node_size = ncol(get_weights(trained_model)[[1]]),
              Integrated_brier = int_brier,
              Brier_scores = brier_set,
              Cindex = cindex,
              Rel_imp_Garson = var_imp$Rels_Garson,
              Rel_imp_Olden = var_imp$Rels_Olden,
              Mse_groups2 = cal_stats2$mse_groups,
              Mse_groups5 = cal_stats5$mse_groups))
}

# function that calculates the metrics from a neural network with nnet
measures_calculator2 <- function(trained_model,
                                 datashort,
                                 datalong) {
  
  real_times <- datashort$survt_ovr
  real_status <- datashort$surv_status
  
  df1 <- data.frame(hazard = predict(trained_model,
                                     as.matrix(datalong[, c(1:5, 12 + max(datalong$interval))])))
  df1$id <- datalong$id # ids of the patients
  df1$survival <- datalong$survival # survival time in years
  groups <- split(df1, f = df1$id) 
  
  true_surv <- unlist(lapply(groups, function(x) {
    surv_obj <- x$survival
    true_res <- surv_obj[1] 
    return(true_res)}
  ))
  group_probs <- lapply(groups, function(x) {
    x <- cumprod(1 - x$hazard)})
  pred_mat <- do.call("rbind", group_probs)
  
  relative_probs <- vector(mode = "numeric",
                           length = length(group_probs))
  for (i in 1:length(group_probs)){ 
    temp <- group_probs[[i]]
    ind <- which(1:max(datalong$interval) == true_surv[i])
    relative_probs[i] <- temp[ind] 
  }
  N0 <- length(unique(df1$id)) # number of unique persons
  
  # in the data frame create random id numbers 
  # to label the patients
  df2 <- data.frame(relative_probs = relative_probs,
                    survt_ovr = real_times,
                    surv_status = real_status, 
                    id = (1001):(1000 + N0))
  df2$prediction <- 1 - round(df2$relative_probs, digits = 0)
  
  prob_matrix <-  cbind(1, pred_mat[, 1:5]) # prediction matrix up to 5 years
  risk_mat <- 1 - prob_matrix
  
  brier_obj <- brier_nn(risk_matrix = risk_mat, data = df2)
  int_brier <- as.numeric(brier_obj$Int_brier)
  brier_set <- brier_obj$Brier
  
  # estimate c-index
  pi_total <- pi_calculator2(trained_model, datalong)
  cindex <- concordance.index(x = pi_total,
                              surv.time = real_times,
                              surv.event = real_status)$c.index
  # assign 0.5 to c-indexes when prognostic indexes are the same or NaN
  if (is.na(cindex)) cindex <- 0.5 
  nr_weights <- length(trained_model$wts)
  
  # calculate variable importance
  var_imp <- var_imp_nnet(model = trained_model)
  
  # calculate calibration statistics at 2 years
  cal_stats2 <- calibration_stats(prob_vector = prob_matrix[, 2+1], data = datashort, times = 2)
  cal_stats5 <- calibration_stats(prob_vector = prob_matrix[, 5+1], data = datashort, times = 5)
  
  return(list(weights = nr_weights, 
              node_size = trained_model$n[2],
              Integrated_brier = int_brier,
              Brier_scores = brier_set,
              Cindex = cindex,
              Rel_imp_Garson = var_imp$Rels_Garson,
              Rel_imp_Olden = var_imp$Rels_Olden,
              Mse_groups2 = cal_stats2$mse_groups,
              Mse_groups5 = cal_stats5$mse_groups))
}


## end of the file



######################################################################
######################################################################
# file "create_train.R"
# Splits randomly generated synthetic data ("data.RData", n = 1000 for each scenario)
# into training and test sets. The training set produced 
# ("training_ovr_scaled.RData", n = 500) is used with 5-fold cross-validation
# to tune PLANN original and extended.
######################################################################
######################################################################

## start of the file

# To install packages
# install_packages <- c("survival", "caret")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

library(survival)
library(caret)
load("data.RData") # data with 1000 simulated patients

# create training and test data
# we use 500 synthetic patients for training

dataset <- data1[, 1:7]
# first for bo06
N <- nrow(dataset)
round(table(dataset$surv_status)/ N, 3) # event for 38.7% of the people

set.seed(12345)
index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
training_ovr <- dataset[index, ]
# test_ovr  <- dataset[-index,]
round(table(training_ovr$surv_status) / nrow(training_ovr), 3)
# round(table(test_ovr$surv_status) / nrow(test_ovr), 3)
# str(training_ovr)


# scale the data
# creating training and test ovr scaled
response_pair <- c("survt_ovr", "surv_status")
# for the training set
X_train_scaled <- training_ovr[, !(colnames(training_ovr) %in% response_pair)]

for(i in 1:length(colnames(X_train_scaled))) {
  if(class(X_train_scaled[,i]) == "numeric" || class(X_train_scaled[,i]) == "integer") {
    X_train_scaled[, i] <- as.vector(scale(X_train_scaled[, i])) }
}
X_train_ovr <- model.matrix(~., X_train_scaled)[, -1] 
Y_train_ovr <- training_ovr[, colnames(training_ovr) %in% response_pair] # create matrix Y
training_ovr_scaled <- as.data.frame(cbind(X_train_ovr, Y_train_ovr))

# for the test set
# x_test_scaled <- test_ovr[, !(colnames(test_ovr) %in% response_pair)]
# for(i in 1:length(colnames(x_test_scaled))) {
#   if(class(x_test_scaled[,i]) == "numeric" || class(x_test_scaled[,i]) == "integer") {
#     x_test_scaled[, i] <- as.vector(scale(x_test_scaled[, i])) }
# }
# X_test_ovr <- model.matrix(~., x_test_scaled)[, -1]
# Y_test_ovr <- test_ovr[, colnames(test_ovr) %in% response_pair]
# test_ovr_scaled <- as.data.frame(cbind(X_test_ovr, Y_test_ovr))

save(training_ovr_scaled, file = "training_ovr_scaled.RData")
#save(test_ovr_scaled, file = "test_ovr_scaled.RData")


## end of the file



######################################################################
######################################################################
# file "snn_training_nnet.R" 

# Used to tune the parameters for PLANN original 
# with library nnet for e.g 61% censoring, yearly intervals
# before performing the simulation studies
######################################################################
######################################################################

## start of the file

# To install packages
# install_packages <- c("survival", "nnet", "fastDummies",
#                       "pec", "doParallel", "tictoc")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }


# train on training data (parameters selected in this step)
library(survival)
library(nnet)
library(fastDummies)
library(pec)
library(doParallel)
library(tictoc)


load("training_ovr_scaled.RData")
source("the_functions2.R")
source("functions_nn.R")


# set up the cross-validation

# most popular optimization algorithms used are the Stochastic Gradient Descent (SGD), ADAM and RMSprop
# you need to tune certain parameters such as learning rate or momentum


nfolds <- 5
set.seed(12345)
folds <- createFolds(training_ovr_scaled$survt_ovr, k = 5, list = TRUE)

node_size <- seq(2, 8, by = 1) # grid of node sizes
weight_decay <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3)
combis <- expand.grid(node_size, weight_decay)


# initialize objects
cv_weights <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_intbrier <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_cindex <- matrix(0, nrow = nfolds, ncol = nrow(combis))

# run the cross-validation
tic("Running the cross-validation")

for (i in 1:nfolds) {
  
  cat("Started iteration i = ", i, "\n")
  indices <- folds[[i]]
  cat("Creating the training set ...", "\n")
  train_set <-  data_train_creator(training_ovr_scaled[-indices, ]) # create the train set
  cat("Creating the validation set ...", "\n")
  validation_set <- data_test_creator(data = training_ovr_scaled[indices, ]) # create the validation set
  
  # create the matrices to be used for keras library
  train_x <- as.matrix(train_set[, c(1:5, 12 + max(train_set$interval))]) # predictors: 5 variables + intervals
  dimnames(train_x) <- NULL # the object must have empty dimnames
  train_y <- train_set$status_long
  validation_x <- as.matrix(validation_set[, c(1:5, 12 + max(validation_set$interval))])
  dimnames(validation_x) <- NULL # the object must have empty dimnames
  validation_y <- validation_set$status_long
  
  for (j in 1:nrow(combis)) {
    
    cat("Testing combination number:", j, "of repeat", i, " out of 5", "\n")
    cat("calculating for node size:", combis[j, 1], "and weight decay", combis[j, 2], "...", "\n")
    
    set.seed(12345)  
    
    # start building the model
    fit_nnet <- nnet(x = train_x, y = train_y, size = combis[j, 1],
                     decay = combis[j, 2], maxit = 1000, trace = FALSE, entropy = TRUE)
    
    
    
    # now that the model has run lets calulate the measures
    # the total weights are 13*(node_size + bias_input) + (node_size + bias)*output 
    values <- measures_calculator2(trained_model = fit_nnet,
                                   datashort = training_ovr_scaled[indices, ],
                                   datalong = validation_set)
    
    cv_weights[i, j] <- values$weights
    cv_intbrier[i, j] <- values$Integrated_brier
    cv_cindex[i, j] <- values$Cindex
    cat("The IBS is:", round(values$Integrated_brier, 3),
        "and the C-index is:", round(values$Cindex, 3), "\n")
    
  }
}

time_list <- toc()

cat((time_list$toc - time_list$tic) / 60, "minutes elapsed", "\n") 


df_logistic_ovr <- as.data.frame(cbind(node_size = combis[, 1],
                                       weight_decay = combis[, 2],
                                       weights = colMeans(cv_weights),
                                       integrated_brier = colMeans(cv_intbrier),
                                       cindex = colMeans(cv_cindex))) 

# save results to be easily accessible
save(df_logistic_ovr, file = "results_logistic_ovr.RData")
#save(time_list, file = "time_cv_logistic_ovr.RData")

# best 3 combinations for IBS, C-index
ind <- head(order(df_logistic_ovr$integrated_brier, decreasing = FALSE), 3)
df_logistic_ovr[ind, ]

ind2 <- head(order(df_logistic_ovr$cindex, decreasing = TRUE), 3)
df_logistic_ovr[ind2, ]


## end of the file


######################################################################
######################################################################
# file "snn_training_keras.R" 

# Used to tune the parameters for PLANN extended
# with library keras for e.g 61% censoring, yearly intervals
# and sigmoid  activation function for the input-hidden layer
# # before performing the simulation studies
######################################################################
######################################################################

## start of file

# To install packages
# install_packages <- c("survival", "fastDummies",
#                       "pec", "doParallel", "tictoc")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

# how to set up keras in R
# first install anaconda 3
# in anaconda command prompt use the following:
# pip install tensorflow
# pip install keras

# test by running: import tensorflow as tf

# this is done
# install.packages("reticulate")
# install.packages("tensorflow")
# install.packages("keras")

# library(tensorflow)
# install_tensorflow()
# library(keras)
# install_keras()

# train on training data (parameters selected in this step)

library(survival)
library(keras)
library(fastDummies)
library(pec)
library(doParallel)
library(tictoc)


load("training_ovr_scaled.RData")
source("the_functions2.R")
source("functions_nn.R")


# set up the cross-validation

# most popular optimization algorithms used are the Stochastic Gradient Descent (SGD), ADAM and RMSprop
# you need to tune certain parameters such as learning rate or momentum


nfolds <- 5
set.seed(12345)
folds <- createFolds(training_ovr_scaled$survt_ovr, k = 5, list = TRUE)

node_size <- seq(1, 16, by = 3) # grid of node sizes
dropout_rate <- c(0.1, 0.2, 0.4)
lr <- c(0.1, 0.2, 0.4)
class_weights <- c(0.95, 1, 1.05)
momentum <- c(0.8, 0.9)
combis <- expand.grid(node_size, dropout_rate, lr, class_weights, momentum)


# initialize objects
cv_weights <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_intbrier <- matrix(0, nrow = nfolds, ncol = nrow(combis))
cv_cindex <- matrix(0, nrow = nfolds, ncol = nrow(combis))

# run the cross-validation
tic("Running the cross-validation")

for (i in 1:nfolds) {
  
  cat("Started iteration i = ", i, "\n")
  indices <- folds[[i]]
  cat("Creating the training set ...", "\n")
  train_set <-  data_train_creator(training_ovr_scaled[-indices, ]) # create the train set
  cat("Creating the validation set ...", "\n")
  validation_set <- data_test_creator(data = training_ovr_scaled[indices, ]) # create the validation set
  
  # create the matrices to be used for keras library
  train_x <- as.matrix(train_set[, c(1:5, 12:(11 + max(train_set$interval)))]) # predictors: 5 variables + intervals
  dimnames(train_x) <- NULL # the object must have empty dimnames
  train_y <- train_set$status_long
  validation_x <- as.matrix(validation_set[, c(1:5, 12:(11 + max(validation_set$interval)))])
  dimnames(validation_x) <- NULL # the object must have empty dimnames
  validation_y <- validation_set$status_long
  
  for (j in 1:nrow(combis)) {
    
    cat("Testing combination number:", j, "of repeat", i, " out of 5", "\n")
    cat("calculating for node size:", combis[j, 1], ", dropout rate:", combis[j, 2], "\n",
        "and learning rate", combis[j, 3], "and weak class weight", combis[j, 4],
        "and momentum", combis[j, 5], "...", "\n")
    
    
    options(keras.view_metrics = FALSE)
    
    k_clear_session() # to avoid clutter from old models / layers in cross validation
    
    tensorflow::tf$random$set_seed(12345) # for update to TensorFlow 2.0
    
    # start building the model
    fit_keras <- keras_model_sequential()
    # Add layers to the model
    # here we have logistic activation function for the inputs but also for the outputs
    # we create a densely connected ANN to the output (input shape 15 + 10)
    fit_keras %>%
      layer_dense(units = combis[j, 1], activation = 'sigmoid',
                  input_shape = c(5 + max(train_set$interval))) %>% 
      layer_dropout(rate = combis[j, 2]) %>%
      layer_dense(units = 1, activation = 'sigmoid')
    fit_keras %>% compile(
      loss = 'binary_crossentropy', # for binary class classification problem
      optimizer = optimizer_sgd(lr = combis[j, 3], momentum = combis[j, 5])
    )
    
    early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 5)
    
    result <- fit_keras %>% fit(
      train_x, 
      train_y, 
      epochs = 50, 
      batch_size = 64,
      validation_data = list(validation_x, validation_y),
      class_weight = list("0" = 1, "1" = combis[j, 4]),
      callbacks = c(early_stopping), # to enforce early stopping in case the loss function stops improving
      verbose = 0 # don't display progress bar
    )
    
    # now that the model has run lets calulate the measures
    # the total weights are 13*(node_size + bias_input) + (node_size + bias)*output 
    values <- measures_calculator(trained_model = fit_keras,
                                  datashort = training_ovr_scaled[indices, ],
                                  datalong = validation_set)
    
    cv_weights[i, j] <- values$weights
    cv_intbrier[i, j] <- values$Integrated_brier
    cv_cindex[i, j] <- values$Cindex
    cat("The IBS is:", round(values$Integrated_brier, 3),
        "and the C-index is:", round(values$Cindex, 3), "\n")
    
  }
}

time_list <- toc()

cat((time_list$toc - time_list$tic) / 60, "minutes elapsed", "\n") 


df_sigmoid_ovr <- as.data.frame(cbind(node_size = combis[, 1],
                                      dropout_rate = combis[, 2],
                                      learning_rate = combis[, 3],
                                      momentum = combis[, 5],
                                      weak_weight = combis[, 4],
                                      weights = colMeans(cv_weights),
                                      integrated_brier = colMeans(cv_intbrier),
                                      cindex = colMeans(cv_cindex))) 

# save results to be easily accessible
save(df_sigmoid_ovr, file = "results_sigmoid_ovr.RData")
save(time_list, file = "time_cv_sigmoid_ovr.RData")

# best 3 combinations for IBS, C-index
ind <- head(order(df_sigmoid_ovr$integrated_brier, decreasing = FALSE), 3)
df_sigmoid_ovr[ind, ]

ind2 <- head(order(df_sigmoid_ovr$cindex, decreasing = TRUE), 3)
df_sigmoid_ovr[ind2, ]

## end of the file



######################################################################
######################################################################
# file "simulations1_ibs.R" (core file to run simulations using 
# "R environment.RData" that sets up necessary parameters based on the original BO06 data)

# START OF SIMULATIONS
# censoring times with Weibull distribution as original ~61% censoring
# Simulations for IBS n = 250 (can be adjusted to 1000 by changing 250 with 1000)
# Similarly simulations can be run for C-index regarding the ML techniques
# by choosing optimally tuned parameters for this measure


# you need to create folders results_simulations1 (used for n = 250)
# and results_simulations2 (used for n = 1000)
# and subfolders results_simulations1/snn_ibs, results_simulations2/snn_ibs
# to store the results for all methods in the destination of this R file

######################################################################
######################################################################


## start of file

# the following simulations are based on 4 categorical variables: trt, sex, hist_resp and 
# excision as well as in 1 continuous variable (age). 

# Initialises the R environment (R objects) to generate all data
# during this study based on the original BO06 data (the BO06 dataset is not required)
load("R environment.RData")

# To install packages
# install_packages <- c("dplyr", "data.table", "tictoc", "survival",
#                        "fastDummies", "nnet", "pec", "caret", "ggpubr")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

# how to set up keras in R
# first install anaconda 3
# in anaconda command prompt use the following:
# pip install tensorflow
# pip install keras

# test by running: import tensorflow as tf

# this is done
# install.packages("reticulate")
# install.packages("tensorflow")
# install.packages("keras")

# library(tensorflow)
# install_tensorflow()
# library(keras)
# install_keras()


library(dplyr)
library(data.table)
library(tictoc)
library(survival)
library(fastDummies)
library(nnet)
library(keras)
library(pec)
library(caret)
library(ggpubr)
source("the_functions2.R")
source("functions_nn.R")



# part 1: count unique combinations of categorical variables
# use the mean age for this case to simulate for variable age

# tabel <-  bo06 %>% count(trt, sex, hist_resp, excision)
# tabel$prop <- tabel$n / nrow(bo06)
# 
# options(warn = -1)
# stats <- bo06 %>%
#   group_by(trt, sex, hist_resp, excision) %>%
#   summarise_at(vars(age), funs(mean = "mean", sd = "sd"))
# options(warn = 0)
# 
# tabel$age_mean <- stats$mean
# tabel$age_sd <- stats$sd
# tabel$case <- 1:nrow(tabel)
# 
# # create 16 lists with the different scenarios
# lists <- split(tabel[, c("trt", "sex", "hist_resp", "excision",
#                          "age_mean", "age_sd")], f = tabel$case) 

# to simulate survival and censoring times
# we need to calculate betas, mu and sigma of the original bo06 data
# assumed distribution for the y variable is loglogistic

# logreg <-  survreg(Surv(survt_ovr, surv_status) ~
#                      trt + sex + hist_resp + excision + age,
#                    data = bo06, dist = "lognormal",
#                    robust = TRUE)
# 
# coefs <- matrix(logreg$coefficients) # intercept is parameter mu
# sigma <- as.vector(logreg$scale)

# censoring distribution of bo06 as weibull
# cens_times <- Surv(bo06$survt_ovr, bo06$surv_status == 0)
# plot(cens_times)
# summary(cens_times)
# distr <- survreg(cens_times ~ 1, dist = "weibull")
# cens <- rweibull(n = nrow(bo06), shape = 1 / distr$scale, scale = exp(distr$coefficients))
# summary(cens)

# data simulation procedure

data_simulator <- function(n = 250) {
  
  ind <- sample(x = 1:nrow(tabel), size = n, prob = tabel$prop, replace = TRUE)
  data_list <- vector(mode = "list", length = length(ind))
  
  for (i in 1:length(ind)) {
    
    value <- ind[i]
    line <- lists[[value]]
    line$age <- rnorm(1, mean = line$age_mean, sd = line$age_sd)
    line$age <- ifelse(line$age < 3.600757, 3.600757, line$age) 
    line$age <- ifelse(line$age > 40.85132, 40.85132, line$age) 
    line <- line[, c("trt", "sex", "hist_resp", "excision", "age")] # keep variables of interest
    data_list[[i]] <- line
  }
  
  data_sim <- as.data.frame(rbindlist(data_list))
  
  # create the model matrix for the simulated data
  mat_data <- model.matrix(~., data_sim[, c("trt", "sex", "hist_resp", "excision", "age")])
  
  # calculate the survival time
  time <- exp(mat_data %*% coefs + sigma*rnorm(n = n, 0, 1))
  # shape and scale parameters calculated with a Weibull distribution for censoring times of BO06
  cens <- rweibull(n = n, shape = 2.0287, scale = 5.7177) 
  cens[cens == 0] <- (1 / 365.25)
  surv_status <- as.numeric(time <= cens)
  survt_ovr <- pmin(time, cens)
  data_sim <- cbind(data_sim, survt_ovr, surv_status, time)
  
  return(data_sim)
}

data1 <- data_simulator(n = 250)

# summary(bo06)
# summary(data1)

table(data1$surv_status)[1] / nrow(data1)


################################################################
# repeat the simulation 1000 times
set.seed(12345)
big_list <- replicate(n = 1000, {
  list(data_simulator(n = 250))
})

# censoring for the big list
cens_perc <- vector(mode = "numeric", length = length(big_list))
#surv_time <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)) {
  datas <- big_list[[i]]
  cens_perc[i] <- table(datas$surv_status)[1] / nrow(datas)
  #surv_time[i] <- max(datas$survt_ovr)
}

summary(cens_perc) # mean around 60.9 %


# START OF SIMULATION FUNCTIONS
########################################################
########################################################
########################################################

# function that runs a Cox model and calculates Brier and Integrated Brier score 
# C-index and the calibration measures for all simulated datasets

cox_simulator <- function(data) {
  
  results_cox <- vector(mode = "list", length = length(data))
  
  for (i in 1:length(data)) {
    
    dataset <- data[[i]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    times <- c(0, 1, 2, 3, 4, 4.999)
    # calculating the Brier and the Integrated Brier scores
    cox_all <- coxph(Surv(survt_ovr, surv_status) ~ trt + sex + hist_resp + excision + age,
                     data = training_ovr,
                     method = "breslow", x = TRUE, y = TRUE)
    probs_cox_all <- predictSurvProb(object = cox_all, newdata = test_ovr,
                                     times = times)
    risk_cox_all <- riskRegression::predictRisk(cox_all, test_ovr,
                                                times = times)
    brier_cox <- brier_general_new(risk_matrix = risk_cox_all, data = test_ovr,
                                   times = times)
    
    lp <- predict(cox_all, newdata = test_ovr, type = "lp")
    cindex_cox <- concordance.index(x = lp,
                                    surv.time = test_ovr$survt_ovr,
                                    surv.event = test_ovr$surv_status)$c.index
    
    calib2y <- calibration_stats(prob_vector = probs_cox_all[, 2+1], data = test_ovr, times = 2)
    calib5y <- calibration_stats(prob_vector = probs_cox_all[, 5+1], data = test_ovr, times = 5)
    
    results_cox[[i]] <- list(Brier = brier_cox$Brier,
                             Int_brier = brier_cox$Int_brier,
                             Cindex = cindex_cox,
                             Calib2y = calib2y,
                             Calib5y = calib5y)
    
    if (any(i == seq(50, 1000, by = 50))) print(i) # print working progress
    
  }
  return(results_cox)
}

results_cox <- cox_simulator(data = big_list)

save(results_cox, file = "results_simulations1/simulations_cox.RData")

# enter the list to calculate the measures
brier_cox <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
ibs_cox <- vector(mode = "numeric", length = length(big_list))
cindex_cox <- vector(mode = "numeric", length = length(big_list))
calib_cox2y_groups <- vector(mode = "numeric", length = length(big_list))
calib_cox5y_groups <- vector(mode = "numeric", length = length(big_list))


for (i in 1:length(big_list)){
  brier_cox[i,] <- results_cox[[i]]$Brier
  ibs_cox[i] <- results_cox[[i]]$Int_brier
  cindex_cox[i] <- results_cox[[i]]$Cindex
  calib_cox2y_groups[i] <- results_cox[[i]]$Calib2y$mse_groups
  calib_cox5y_groups[i] <- results_cox[[i]]$Calib5y$mse_groups
}

# brier score at 2 years
summa(brier_cox[, 3])
# brier score at 5 years
summa(brier_cox[, 6])
# ibs at 5 years
summa(ibs_cox)
# C-index
summa(cindex_cox)
# Calibration 2 years
summa(calib_cox2y_groups)
# calibration 5 years
summa(calib_cox5y_groups)


###############################################################
# Cox simulator B = 1000, n = 250, split sample = 1/2
# To apply this wrong method, patients censored before year 2 will be removed
# from the data. Then, the Brier score and the rest of the measures defined are calculated


cox_simulator_wr1 <- function(data) {
  
  results_cox_wr1 <- vector(mode = "list", length = length(data))
  
  for (i in 1:length(data)) {
    
    dataset <- data[[i]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    # remove censored cases before year 2 (for the training data)
    cases_removed1y <- sum(training_ovr$survt_ovr <= 1 & training_ovr$surv_status == 0)
    cases_removed2y <- sum(training_ovr$survt_ovr > 1 &  training_ovr$survt_ovr < 2
                           & training_ovr$surv_status == 0)
    cases_removed <- sum(training_ovr$survt_ovr < 2
                         & training_ovr$surv_status == 0)
    training_ovr <- training_ovr[!(training_ovr$survt_ovr < 2 & training_ovr$surv_status == 0), ]
    
    times <- c(0, 1, 2, 3, 4, 4.999)
    # calculating the Brier and the Integrated Brier scores
    cox_all <- coxph(Surv(survt_ovr, surv_status) ~ trt + sex + hist_resp + excision + age,
                     data = training_ovr,
                     method = "breslow", x = TRUE, y = TRUE)
    probs_cox_all <- predictSurvProb(object = cox_all, newdata = test_ovr,
                                     times = times)
    risk_cox_all <- riskRegression::predictRisk(cox_all, test_ovr,
                                                times = times)
    brier_cox <- brier_general_new(risk_matrix = risk_cox_all, data = test_ovr,
                                   times = times)
    
    lp <- predict(cox_all, newdata = test_ovr, type = "lp")
    cindex_cox <- concordance.index(x = lp,
                                    surv.time = test_ovr$survt_ovr,
                                    surv.event = test_ovr$surv_status)$c.index
    
    calib2y <- calibration_stats(prob_vector = probs_cox_all[, 2+1], data = test_ovr, times = 2)
    calib5y <- calibration_stats(prob_vector = probs_cox_all[, 5+1], data = test_ovr, times = 5)
    
    results_cox_wr1[[i]] <- list(Cases_removed1y = cases_removed1y,
                                 Cases_removed2y = cases_removed2y,
                                 Cases_removed = cases_removed,
                                 Brier = brier_cox$Brier,
                                 Int_brier = brier_cox$Int_brier,
                                 Cindex = cindex_cox,
                                 Calib2y = calib2y,
                                 Calib5y = calib5y)
    
    if (any(i == seq(50, 1000, by = 50))) print(i) # print working progress
    
  }
  return(results_cox_wr1)
}

results_cox_wr1 <- cox_simulator_wr1(data = big_list)

save(results_cox_wr1, file = "results_simulations1/simulations_cox_wr1.RData")


# reach the list to estimate the statistics
brier_cox_wr1 <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
cases_removed1y <- vector(mode = "numeric", length = length(big_list))
cases_removed2y <- vector(mode = "numeric", length = length(big_list))
cases_removed <- vector(mode = "numeric", length = length(big_list))
ibs_cox_wr1 <- vector(mode = "numeric", length = length(big_list))
cindex_cox_wr1 <- vector(mode = "numeric", length = length(big_list))
calib_cox2y_groups_wr1 <- vector(mode = "numeric", length = length(big_list))
calib_cox5y_groups_wr1 <- vector(mode = "numeric", length = length(big_list))


for (i in 1:length(big_list)){
  brier_cox_wr1[i,] <- results_cox_wr1[[i]]$Brier
  cases_removed1y[i] <- results_cox_wr1[[i]]$Cases_removed1y
  cases_removed2y[i] <- results_cox_wr1[[i]]$Cases_removed2y
  cases_removed[i] <- results_cox_wr1[[i]]$Cases_removed
  ibs_cox_wr1[i] <- results_cox_wr1[[i]]$Int_brier
  cindex_cox_wr1[i] <- results_cox_wr1[[i]]$Cindex
  calib_cox2y_groups_wr1[i] <- results_cox_wr1[[i]]$Calib2y$mse_groups
  calib_cox5y_groups_wr1[i] <- results_cox_wr1[[i]]$Calib5y$mse_groups
}

#cases removed year 1
summary(cases_removed1y)
# cases removed year 2
summary(cases_removed2y)
# cases removed
summary(cases_removed)
# brier score at 2 years
summa(brier_cox_wr1[, 3])
# brier score at 5 years
summa(brier_cox_wr1[, 6])
# ibs at 5 years
summa(ibs_cox_wr1)
# C-index
summa(cindex_cox_wr1)
# Calibration 2 years
summa(calib_cox2y_groups_wr1)
# calibration 5 years
summa(calib_cox5y_groups_wr1)


##############################################################################
# Cox simulator B = 1000, n = 250, split sample = 1/2
# patients curtailed at 5 years (administrative censoring) on the training set

cox_simulator_wr2 <- function(data) {
  
  results_cox_wr2 <- vector(mode = "list", length = length(data))
  
  for (i in 1:length(data)) {
    
    dataset <- data[[i]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    # curtail patients at 5 years on the training data
    cases_curtailed <- sum(training_ovr$survt_ovr >= 5)
    training_ovr$surv_status <- training_ovr$surv_status*(training_ovr$survt_ovr <= 5)
    training_ovr$survt_ovr <- pmin(training_ovr$survt_ovr, 5)
    
    times <- c(0, 1, 2, 3, 4, 4.999)
    # calculating the Brier and the Integrated Brier scores
    cox_all <- coxph(Surv(survt_ovr, surv_status) ~ trt + sex + hist_resp + excision + age,
                     data = training_ovr,
                     method = "breslow", x = TRUE, y = TRUE)
    probs_cox_all <- predictSurvProb(object = cox_all, newdata = test_ovr,
                                     times = times)
    risk_cox_all <- riskRegression::predictRisk(cox_all, test_ovr,
                                                times = times)
    brier_cox <- brier_general_new(risk_matrix = risk_cox_all, data = test_ovr,
                                   times = times)
    
    lp <- predict(cox_all, newdata = test_ovr, type = "lp")
    cindex_cox <- concordance.index(x = lp,
                                    surv.time = test_ovr$survt_ovr,
                                    surv.event = test_ovr$surv_status)$c.index
    
    calib2y <- calibration_stats(prob_vector = probs_cox_all[, 2+1], data = test_ovr, times = 2)
    calib5y <- calibration_stats(prob_vector = probs_cox_all[, 5+1], data = test_ovr, times = 5)
    
    results_cox_wr2[[i]] <- list(Cases_curtailed = cases_curtailed,
                                 Brier = brier_cox$Brier,
                                 Int_brier = brier_cox$Int_brier,
                                 Cindex = cindex_cox,
                                 Calib2y = calib2y,
                                 Calib5y = calib5y)
    
    if (any(i == seq(50, 1000, by = 50))) print(i) # print working progress
    
  }
  return(results_cox_wr2)
}

results_cox_wr2 <- cox_simulator_wr2(data = big_list)

save(results_cox_wr2, file = "results_simulations1/simulations_cox_wr2.RData")


# reach the list to estimate the statistics
brier_cox_wr2 <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
cases_curtailed <- vector(mode = "numeric", length = length(big_list))
ibs_cox_wr2 <- vector(mode = "numeric", length = length(big_list))
cindex_cox_wr2 <- vector(mode = "numeric", length = length(big_list))
calib_cox2y_groups_wr2 <- vector(mode = "numeric", length = length(big_list))
calib_cox5y_groups_wr2 <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)){
  brier_cox_wr2[i,] <- results_cox_wr2[[i]]$Brier
  cases_curtailed[i] <- results_cox_wr2[[i]]$Cases_curtailed
  ibs_cox_wr2[i] <- results_cox_wr2[[i]]$Int_brier
  cindex_cox_wr2[i] <- results_cox_wr2[[i]]$Cindex
  calib_cox2y_groups_wr2[i] <- results_cox_wr2[[i]]$Calib2y$mse_groups
  calib_cox5y_groups_wr2[i] <- results_cox_wr2[[i]]$Calib5y$mse_groups
}

# cases curtailed
summary(cases_curtailed)
# brier score at 2 years
summa(brier_cox_wr2[, 3])
# brier score at 5 years
summa(brier_cox_wr2[, 6])
# ibs at 5 years
summa(ibs_cox_wr2)
# C-index
summa(cindex_cox_wr2)
# Calibration 2 years
summa(calib_cox2y_groups_wr2)
# calibration 5 years
summa(calib_cox5y_groups_wr2)



###############################################################
###############################################################
###############################################################
# apply the method of Biganzoli on nnet library
# function that runs a NN and calculates Brier and Integrated Brier score 
# for all simulated datasets

nnet_simulator <- function(data) {
  
  results_nnet <- vector(mode = "list", length = length(data))
  
  for (a in 1:length(data)) {
    
    dataset <- data[[a]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    # scaling the training and test data 
    # creating training and test ovr scaled
    response_pair <- c("survt_ovr", "surv_status")
    # for the training set
    X_train_scaled <- training_ovr[, !(colnames(training_ovr) %in% response_pair)]
    
    for(i in 1:length(colnames(X_train_scaled))) {
      if(class(X_train_scaled[,i]) == "numeric" || class(X_train_scaled[,i]) == "integer") {
        X_train_scaled[, i] <- as.vector(scale(X_train_scaled[, i])) }
    }
    X_train_ovr <- model.matrix(~., X_train_scaled)[, -1] 
    Y_train_ovr <- training_ovr[, colnames(training_ovr) %in% response_pair] # create matrix Y
    training_ovr_scaled <- as.data.frame(cbind(X_train_ovr, Y_train_ovr))
    
    # for the test set
    x_test_scaled <- test_ovr[, !(colnames(test_ovr) %in% response_pair)]
    for(i in 1:length(colnames(x_test_scaled))) {
      if(class(x_test_scaled[,i]) == "numeric" || class(x_test_scaled[,i]) == "integer") {
        x_test_scaled[, i] <- as.vector(scale(x_test_scaled[, i])) }
    }
    X_test_ovr <- model.matrix(~., x_test_scaled)[, -1]
    Y_test_ovr <- test_ovr[, colnames(test_ovr) %in% response_pair]
    test_ovr_scaled <- as.data.frame(cbind(X_test_ovr, Y_test_ovr))
    
    training_ovr_long <- data_train_creator(training_ovr_scaled) # 5 + 8 = 13 prognostic variables
    test_ovr_long <- data_test_creator(test_ovr_scaled)
    
    # create the matrices to be used for keras library
    train_x <- as.matrix(training_ovr_long[, c(1:5, 12 + max(training_ovr_long$interval))])
    dimnames(train_x) <- NULL # the object must have empty dimnames
    train_y <- training_ovr_long$status_long
    
    set.seed(12345)  
    
    # start building the model
    fit_nnet <- nnet(x = train_x, y = train_y, size = 7, decay = 0.3,
                     maxit = 1000, trace = FALSE, entropy = TRUE)
    
    metrics_model <- measures_calculator2(trained_model = fit_nnet,
                                          datashort = test_ovr_scaled,
                                          datalong = test_ovr_long)
    
    if (any(a == seq(20, 1000, by = 20))) {
      cat("The IBS of dataset", a, "is:", round(metrics_model$Integrated_brier, 3), "\n")
    }
    
    results_nnet[[a]] <- list(Time = 0:5, 
                              Brier = metrics_model$Brier_scores,
                              Int_brier = metrics_model$Integrated_brier,
                              Cindex = metrics_model$Cindex,
                              Calib2y_mse_groups = metrics_model$Mse_groups2,
                              Calib5y_mse_groups = metrics_model$Mse_groups5)
  }
  
  return(results_nnet)
  
}

# run the cross-validation
tic("Running simulations for the neural network of Biganzoli in keras")

results_nnet <- nnet_simulator(data = big_list)

time_list <- toc()
cat((time_list$toc - time_list$tic) / 60, "minutes elapsed", "\n") 
save(results_nnet, file = "results_simulations1/snn_ibs/simulations_nnet.RData")
# save(time_list, file = "results_simulations1/snn_ibs/simulations_nnet_time.RData")

# enter the list to calculate the measures
brier_nnet <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
ibs_nnet <- vector(mode = "numeric", length = length(big_list))
cindex_nnet <- vector(mode = "numeric", length = length(big_list))
calib_nnet2y_groups <- vector(mode = "numeric", length = length(big_list))
calib_nnet5y_groups <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)){
  brier_nnet[i,] <- results_nnet[[i]]$Brier
  ibs_nnet[i] <- results_nnet[[i]]$Int_brier
  cindex_nnet[i] <- results_nnet[[i]]$Cindex
  calib_nnet2y_groups[i] <- results_nnet[[i]]$Calib2y_mse_groups
  calib_nnet5y_groups[i] <- results_nnet[[i]]$Calib5y_mse_groups
}

# brier score at 2 years
summa(brier_nnet[, 3])
# brier score at 5 years
summa(brier_nnet[, 6]) 
# ibs at 5 years
summa(ibs_nnet)
# C-index
summa(cindex_nnet)
# Calibration 2 years
summa(calib_nnet2y_groups)
# calibration 5 years
summa(calib_nnet5y_groups)

################################################################################
# nnet simulator wrong 1
# B = 1000, n = 250, split sample = 1/2
# To apply this wrong method, patients censored before year 2 will be removed
# from the data. Then, the Brier score and the rest of the measures defined are calculated

nnet_simulator_wr1 <- function(data) {
  
  results_nnet_wr1 <- vector(mode = "list", length = length(data))
  
  for (a in 1:length(data)) {
    
    dataset <- data[[a]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    # remove  patients censored before 2 years
    cases_removed <- sum(training_ovr$survt_ovr < 2 
                         & training_ovr$surv_status == 0)
    training_ovr <- training_ovr[!(training_ovr$survt_ovr < 2 & training_ovr$surv_status == 0), ]
    
    # scaling the training and test data 
    # creating training and test ovr scaled
    response_pair <- c("survt_ovr", "surv_status")
    # for the training set
    X_train_scaled <- training_ovr[, !(colnames(training_ovr) %in% response_pair)]
    
    for(i in 1:length(colnames(X_train_scaled))) {
      if(class(X_train_scaled[,i]) == "numeric" || class(X_train_scaled[,i]) == "integer") {
        X_train_scaled[, i] <- as.vector(scale(X_train_scaled[, i])) }
    }
    X_train_ovr <- model.matrix(~., X_train_scaled)[, -1] 
    Y_train_ovr <- training_ovr[, colnames(training_ovr) %in% response_pair] # create matrix Y
    training_ovr_scaled <- as.data.frame(cbind(X_train_ovr, Y_train_ovr))
    
    # for the test set
    x_test_scaled <- test_ovr[, !(colnames(test_ovr) %in% response_pair)]
    for(i in 1:length(colnames(x_test_scaled))) {
      if(class(x_test_scaled[,i]) == "numeric" || class(x_test_scaled[,i]) == "integer") {
        x_test_scaled[, i] <- as.vector(scale(x_test_scaled[, i])) }
    }
    X_test_ovr <- model.matrix(~., x_test_scaled)[, -1]
    Y_test_ovr <- test_ovr[, colnames(test_ovr) %in% response_pair]
    test_ovr_scaled <- as.data.frame(cbind(X_test_ovr, Y_test_ovr))
    
    training_ovr_long <- data_train_creator(training_ovr_scaled) # 5 + 8 = 13 prognostic variables
    test_ovr_long <- data_test_creator(test_ovr_scaled)
    
    # create the matrices to be used for keras library
    train_x <- as.matrix(training_ovr_long[, c(1:5, 12 + max(training_ovr_long$interval))])
    dimnames(train_x) <- NULL # the object must have empty dimnames
    train_y <- training_ovr_long$status_long
    
    set.seed(12345)  
    
    # start building the model
    fit_nnet <- nnet(x = train_x, y = train_y, size = 7, decay = 0.3,
                     maxit = 1000, trace = FALSE, entropy = TRUE)
    
    metrics_model <- measures_calculator2(trained_model = fit_nnet,
                                          datashort = test_ovr_scaled,
                                          datalong = test_ovr_long)
    
    if (any(a == seq(20, 1000, by = 20))) {
      cat("The IBS of dataset", a, "is:", round(metrics_model$Integrated_brier, 3), "\n")
    }
    
    results_nnet_wr1[[a]] <- list(Cases_removed = cases_removed, 
                                  Time = 0:5, 
                                  Brier = metrics_model$Brier_scores,
                                  Int_brier = metrics_model$Integrated_brier,
                                  Cindex = metrics_model$Cindex,
                                  Calib2y_mse_groups = metrics_model$Mse_groups2,
                                  Calib5y_mse_groups = metrics_model$Mse_groups5)
  }
  
  return(results_nnet_wr1)
  
}

# run the cross-validation
tic("Running simulations for the neural network of Biganzoli in nnet")

results_nnet_wr1 <- nnet_simulator_wr1(data = big_list)

time_list_wr1 <- toc()
cat((time_list_wr1$toc - time_list_wr1$tic) / 60, "minutes elapsed", "\n") 
save(results_nnet_wr1, file = "results_simulations1/snn_ibs/simulations_nnet_wr1.RData")
#save(time_list_wr1, file = "results_simulations1/snn_ibs/simulations_nnet_time_wr1.RData")

# enter the list to calculate the measures
brier_nnet_wr1 <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
cases_removed <- vector(mode = "numeric", length = length(big_list))
ibs_nnet_wr1 <- vector(mode = "numeric", length = length(big_list))
cindex_nnet_wr1 <- vector(mode = "numeric", length = length(big_list))
calib_nnet2y_groups_wr1 <- vector(mode = "numeric", length = length(big_list))
calib_nnet5y_groups_wr1 <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)){
  brier_nnet_wr1[i,] <- results_nnet_wr1[[i]]$Brier
  cases_removed[i] <- results_nnet_wr1[[i]]$Cases_removed
  ibs_nnet_wr1[i] <- results_nnet_wr1[[i]]$Int_brier
  cindex_nnet_wr1[i] <- results_nnet_wr1[[i]]$Cindex
  calib_nnet2y_groups_wr1[i] <- results_nnet_wr1[[i]]$Calib2y_mse_groups
  calib_nnet5y_groups_wr1[i] <- results_nnet_wr1[[i]]$Calib5y_mse_groups
}

# cases removed before 2 years
summary(cases_removed)
# brier score at 2 years
summa(brier_nnet_wr1[, 3])
# brier score at 5 years
summa(brier_nnet_wr1[, 6]) 
# ibs at 5 years
summa(ibs_nnet_wr1)
# C-index
summa(cindex_nnet_wr1)
# Calibration 2 years
summa(calib_nnet2y_groups_wr1)
# calibration 5 years
summa(calib_nnet5y_groups_wr1)



##############################################################################
# nnet simulator wrong 2 B = 1000, n = 250, split sample = 1/2
# patients curtailed at 5 years (administrative censoring) on the training set

nnet_simulator_wr2 <- function(data) {
  
  results_nnet_wr2 <- vector(mode = "list", length = length(data))
  
  for (a in 1:length(data)) {
    
    dataset <- data[[a]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    # curtail patients at 5 years on the training data
    cases_curtailed <- sum(training_ovr$survt_ovr >= 5)
    training_ovr$surv_status <- training_ovr$surv_status*(training_ovr$survt_ovr <= 5)
    training_ovr$survt_ovr <- pmin(training_ovr$survt_ovr, 5)
    
    # scaling the training and test data 
    # creating training and test ovr scaled
    response_pair <- c("survt_ovr", "surv_status")
    # for the training set
    X_train_scaled <- training_ovr[, !(colnames(training_ovr) %in% response_pair)]
    
    for(i in 1:length(colnames(X_train_scaled))) {
      if(class(X_train_scaled[,i]) == "numeric" || class(X_train_scaled[,i]) == "integer") {
        X_train_scaled[, i] <- as.vector(scale(X_train_scaled[, i])) }
    }
    X_train_ovr <- model.matrix(~., X_train_scaled)[, -1] 
    Y_train_ovr <- training_ovr[, colnames(training_ovr) %in% response_pair] # create matrix Y
    training_ovr_scaled <- as.data.frame(cbind(X_train_ovr, Y_train_ovr))
    
    # for the test set
    x_test_scaled <- test_ovr[, !(colnames(test_ovr) %in% response_pair)]
    for(i in 1:length(colnames(x_test_scaled))) {
      if(class(x_test_scaled[,i]) == "numeric" || class(x_test_scaled[,i]) == "integer") {
        x_test_scaled[, i] <- as.vector(scale(x_test_scaled[, i])) }
    }
    X_test_ovr <- model.matrix(~., x_test_scaled)[, -1]
    Y_test_ovr <- test_ovr[, colnames(test_ovr) %in% response_pair]
    test_ovr_scaled <- as.data.frame(cbind(X_test_ovr, Y_test_ovr))
    
    
    training_ovr_long <- data_train_creator(training_ovr_scaled) # 5 + 8 = 13 prognostic variables
    test_ovr_long <- data_test_creator(test_ovr_scaled, enforce_intervals = TRUE)
    
    # create the matrices to be used for keras library
    train_x <- as.matrix(training_ovr_long[, c(1:5, 12 + max(training_ovr_long$interval))])
    dimnames(train_x) <- NULL # the object must have empty dimnames
    train_y <- training_ovr_long$status_long
    
    set.seed(12345)  
    
    # start building the model
    fit_nnet <- nnet(x = train_x, y = train_y, size = 7, decay = 0.3,
                     maxit = 1000, trace = FALSE, entropy = TRUE)
    
    metrics_model <- measures_calculator2(trained_model = fit_nnet,
                                          datashort = test_ovr_scaled,
                                          datalong = test_ovr_long)
    
    if (any(a == seq(20, 1000, by = 20))) {
      cat("The IBS of dataset", a, "is:", round(metrics_model$Integrated_brier, 3), "\n")
    }
    
    results_nnet_wr2[[a]] <- list(Cases_curtailed = cases_curtailed, 
                                  Time = 0:5, 
                                  Brier = metrics_model$Brier_scores,
                                  Int_brier = metrics_model$Integrated_brier,
                                  Cindex = metrics_model$Cindex,
                                  Calib2y_mse_groups = metrics_model$Mse_groups2,
                                  Calib5y_mse_groups = metrics_model$Mse_groups5)
  }
  
  return(results_nnet_wr2)
  
}

# run the cross-validation
tic("Running simulations for the neural network of Biganzoli in nnet")

results_nnet_wr2 <- nnet_simulator_wr2(data = big_list)

time_list_wr2 <- toc()
cat((time_list_wr2$toc - time_list_wr2$tic) / 60, "minutes elapsed", "\n") 
save(results_nnet_wr2, file = "results_simulations1/snn_ibs/simulations_nnet_wr2.RData")
#save(time_list_wr2, file = "results_simulations1/snn_ibs/simulations_nnet_time_wr2.RData")

# enter the list to calculate the measures
brier_nnet_wr2 <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
cases_curtailed <- vector(mode = "numeric", length = length(big_list))
ibs_nnet_wr2 <- vector(mode = "numeric", length = length(big_list))
cindex_nnet_wr2 <- vector(mode = "numeric", length = length(big_list))
calib_nnet2y_groups_wr2 <- vector(mode = "numeric", length = length(big_list))
calib_nnet5y_groups_wr2 <- vector(mode = "numeric", length = length(big_list))


for (i in 1:length(big_list)){
  brier_nnet_wr2[i,] <- results_nnet_wr2[[i]]$Brier
  cases_curtailed[i] <- results_nnet_wr2[[i]]$Cases_curtailed
  ibs_nnet_wr2[i] <- results_nnet_wr2[[i]]$Int_brier
  cindex_nnet_wr2[i] <- results_nnet_wr2[[i]]$Cindex
  calib_nnet2y_groups_wr2[i] <- results_nnet_wr2[[i]]$Calib2y_mse_groups
  calib_nnet5y_groups_wr2[i] <- results_nnet_wr2[[i]]$Calib5y_mse_groups
}

# cases curtailed at 5 years
summary(cases_curtailed)
# brier score at 2 years
summa(brier_nnet_wr2[, 3])
# brier score at 5 years
summa(brier_nnet_wr2[, 6]) 
# ibs at 5 years
summa(ibs_nnet_wr2)
# C-index
summa(cindex_nnet_wr2)
# Calibration 2 years
summa(calib_nnet2y_groups_wr2)
# calibration 5 years
summa(calib_nnet5y_groups_wr2)


###############################################################
###############################################################
###############################################################

# apply the method of Biganzoli on keras library
# function that runs a NN and calculates Brier and Integrated Brier score 
# for all simulated datasets
# set global default to never show metrics
options(keras.view_metrics = FALSE)

keras_simulator <- function(data) {
  
  results_keras <- vector(mode = "list", length = length(data))
  
  for (a in 1:length(data)) {
    
    dataset <- data[[a]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    # scaling the training and test data 
    # creating training and test ovr scaled
    response_pair <- c("survt_ovr", "surv_status")
    # for the training set
    X_train_scaled <- training_ovr[, !(colnames(training_ovr) %in% response_pair)]
    
    for(i in 1:length(colnames(X_train_scaled))) {
      if(class(X_train_scaled[,i]) == "numeric" || class(X_train_scaled[,i]) == "integer") {
        X_train_scaled[, i] <- as.vector(scale(X_train_scaled[, i])) }
    }
    X_train_ovr <- model.matrix(~., X_train_scaled)[, -1] 
    Y_train_ovr <- training_ovr[, colnames(training_ovr) %in% response_pair] # create matrix Y
    training_ovr_scaled <- as.data.frame(cbind(X_train_ovr, Y_train_ovr))
    
    # for the test set
    x_test_scaled <- test_ovr[, !(colnames(test_ovr) %in% response_pair)]
    for(i in 1:length(colnames(x_test_scaled))) {
      if(class(x_test_scaled[,i]) == "numeric" || class(x_test_scaled[,i]) == "integer") {
        x_test_scaled[, i] <- as.vector(scale(x_test_scaled[, i])) }
    }
    X_test_ovr <- model.matrix(~., x_test_scaled)[, -1]
    Y_test_ovr <- test_ovr[, colnames(test_ovr) %in% response_pair]
    test_ovr_scaled <- as.data.frame(cbind(X_test_ovr, Y_test_ovr))
    
    training_ovr_long <- data_train_creator(training_ovr_scaled) # 5 + 8 = 13 prognostic variables
    validation_ovr_long <- data_test_creator(training_ovr_scaled)
    test_ovr_long <- data_test_creator(test_ovr_scaled)
    
    k_clear_session() # to avoid clutter from old models / layers in cross validation
    
    tensorflow::tf$random$set_seed(12345) # for update to TensorFlow 2.0
    # use_session_with_seed(seed = 12345, disable_gpu = TRUE,
    #                       disable_parallel_cpu = TRUE,
    #                       quiet = TRUE)
    
    # model with one hidden layer
    
    fit_keras1 <- keras_model_sequential()
    # Add layers to the model
    # here we have logistic activation function for the inputs but also for the outputs
    # we create a densely connected ANN to the output
    
    
    fit_keras1 %>%
      layer_dense(units = 16, activation = 'tanh',
                  input_shape = c(5 + max(training_ovr_long$interval))) %>%
      layer_dropout(rate = 0.1) %>%
      layer_dense(units = 1, activation = 'sigmoid')
    fit_keras1 %>% compile(
      loss = 'binary_crossentropy', # for binary class classification problem
      optimizer = optimizer_sgd(lr = 0.2, momentum = 0.8)
    )
    early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 5)
    
    # create the matrices to be used for keras library
    train_x <- as.matrix(training_ovr_long[, c(1:5, 12:(11 + max(training_ovr_long$interval)))])
    dimnames(train_x) <- NULL # the object must have empty dimnames
    train_y <- training_ovr_long$status_long
    validation_x <- as.matrix(validation_ovr_long[, c(1:5, 12:(11 + max(validation_ovr_long$interval)))])
    dimnames(validation_x) <- NULL # the object must have empty dimnames
    validation_y <- validation_ovr_long$status_long
    
    result_keras <- fit_keras1 %>% fit(
      train_x,
      train_y,
      epochs = 50,
      batch_size = 64,
      validation_data = list(validation_x, validation_y),
      class_weight = list("0" = 1, "1" = 1.05),
      callbacks = c(early_stopping),
      verbose = 0 # don't display progress bar
    )
    
    
    metrics_model <- measures_calculator(trained_model = fit_keras1,
                                         datashort = test_ovr_scaled,
                                         datalong = test_ovr_long)
    
    if (any(a == seq(10, 1000, by = 10))) {
      cat("The IBS of dataset", a, "is:", round(metrics_model$Integrated_brier, 3), "\n")
    }
    
    results_keras[[a]] <- list(Time = 0:5, 
                               Brier = metrics_model$Brier_scores,
                               Int_brier = metrics_model$Integrated_brier,
                               Cindex = metrics_model$Cindex,
                               Calib2y_mse_groups = metrics_model$Mse_groups2,
                               Calib5y_mse_groups = metrics_model$Mse_groups5)
  }
  
  return(results_keras)
  
}

# run the cross-validation
tic("Running simulations for the neural network of Biganzoli in keras")

results_keras <- keras_simulator(data = big_list)

time_list <- toc()
cat((time_list$toc - time_list$tic) / 60, "minutes elapsed", "\n") 
save(results_keras, file = "results_simulations1/snn_ibs/simulations_keras.RData")
save(time_list, file = "results_simulations1/snn_ibs/simulations_keras_time.RData")

# enter the list to calculate the measures
brier_keras <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
ibs_keras <- vector(mode = "numeric", length = length(big_list))
cindex_keras <- vector(mode = "numeric", length = length(big_list))
calib_keras2y_groups <- vector(mode = "numeric", length = length(big_list))
calib_keras5y_groups <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)){
  brier_keras[i,] <- results_keras[[i]]$Brier
  ibs_keras[i] <- results_keras[[i]]$Int_brier
  cindex_keras[i] <- results_keras[[i]]$Cindex
  calib_keras2y_groups[i] <- results_keras[[i]]$Calib2y_mse_groups
  calib_keras5y_groups[i] <- results_keras[[i]]$Calib5y_mse_groups
}

# brier score at 2 years
summa(brier_keras[, 3])
# brier score at 5 years
summa(brier_keras[, 6]) 
# ibs at 5 years
summa(ibs_keras)
# C-index
summa(cindex_keras)
# Calibration 2 years
summa(calib_keras2y_groups)
# calibration 5 years
summa(calib_keras5y_groups)

################################################################################
# B = 1000, n = 250, split sample = 1/2
# To apply this wrong method, patients censored before year 2 will be removed
# from the data. Then, the Brier score and the rest of the measures defined are calculated

options(keras.view_metrics = FALSE)

keras_simulator_wr1 <- function(data) {
  
  results_keras_wr1 <- vector(mode = "list", length = length(data))
  
  for (a in 1:length(data)) {
    
    dataset <- data[[a]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    # remove  patients censored before 2 years
    cases_removed <- sum(training_ovr$survt_ovr < 2 
                         & training_ovr$surv_status == 0)
    training_ovr <- training_ovr[!(training_ovr$survt_ovr < 2 & training_ovr$surv_status == 0), ]
    
    # scaling the training and test data 
    # creating training and test ovr scaled
    response_pair <- c("survt_ovr", "surv_status")
    # for the training set
    X_train_scaled <- training_ovr[, !(colnames(training_ovr) %in% response_pair)]
    
    for(i in 1:length(colnames(X_train_scaled))) {
      if(class(X_train_scaled[,i]) == "numeric" || class(X_train_scaled[,i]) == "integer") {
        X_train_scaled[, i] <- as.vector(scale(X_train_scaled[, i])) }
    }
    X_train_ovr <- model.matrix(~., X_train_scaled)[, -1] 
    Y_train_ovr <- training_ovr[, colnames(training_ovr) %in% response_pair] # create matrix Y
    training_ovr_scaled <- as.data.frame(cbind(X_train_ovr, Y_train_ovr))
    
    # for the test set
    x_test_scaled <- test_ovr[, !(colnames(test_ovr) %in% response_pair)]
    for(i in 1:length(colnames(x_test_scaled))) {
      if(class(x_test_scaled[,i]) == "numeric" || class(x_test_scaled[,i]) == "integer") {
        x_test_scaled[, i] <- as.vector(scale(x_test_scaled[, i])) }
    }
    X_test_ovr <- model.matrix(~., x_test_scaled)[, -1]
    Y_test_ovr <- test_ovr[, colnames(test_ovr) %in% response_pair]
    test_ovr_scaled <- as.data.frame(cbind(X_test_ovr, Y_test_ovr))
    
    training_ovr_long <- data_train_creator(training_ovr_scaled) # 5 + 8 = 13 prognostic variables
    validation_ovr_long <- data_test_creator(training_ovr_scaled)
    test_ovr_long <- data_test_creator(test_ovr_scaled)
    
    k_clear_session() # to avoid clutter from old models / layers in cross validation
    
    tensorflow::tf$random$set_seed(12345) # for update to TensorFlow 2.0
    
    # model with one hidden layer
    
    fit_keras1 <- keras_model_sequential()
    # Add layers to the model
    # here we have logistic activation function for the inputs but also for the outputs
    # we create a densely connected ANN to the output
    
    fit_keras1 %>%
      layer_dense(units = 16, activation = 'tanh',
                  input_shape = c(5 + max(training_ovr_long$interval))) %>%
      layer_dropout(rate = 0.1) %>%
      layer_dense(units = 1, activation = 'sigmoid')
    fit_keras1 %>% compile(
      loss = 'binary_crossentropy', # for binary class classification problem
      optimizer = optimizer_sgd(lr = 0.2, momentum = 0.8)
    )
    early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 5)
    
    # create the matrices to be used for keras library
    train_x <- as.matrix(training_ovr_long[, c(1:5, 12:(11 + max(training_ovr_long$interval)))])
    dimnames(train_x) <- NULL # the object must have empty dimnames
    train_y <- training_ovr_long$status_long
    validation_x <- as.matrix(validation_ovr_long[, c(1:5, 12:(11 + max(validation_ovr_long$interval)))])
    dimnames(validation_x) <- NULL # the object must have empty dimnames
    validation_y <- validation_ovr_long$status_long
    
    result_keras <- fit_keras1 %>% fit(
      train_x,
      train_y,
      epochs = 50,
      batch_size = 64,
      validation_data = list(validation_x, validation_y),
      class_weight = list("0" = 1, "1" = 1.05),
      callbacks = c(early_stopping),
      verbose = 0 # don't display progress bar
    )
    
    
    
    metrics_model <- measures_calculator(trained_model = fit_keras1,
                                         datashort = test_ovr_scaled,
                                         datalong = test_ovr_long)
    
    if (any(a == seq(10, 1000, by = 10))) {
      cat("The IBS of dataset", a, "is:", round(metrics_model$Integrated_brier, 3), "\n")
    }
    
    results_keras_wr1[[a]] <- list(Cases_removed = cases_removed, 
                                   Time = 0:5, 
                                   Brier = metrics_model$Brier_scores,
                                   Int_brier = metrics_model$Integrated_brier,
                                   Cindex = metrics_model$Cindex,
                                   Calib2y_mse_groups = metrics_model$Mse_groups2,
                                   Calib5y_mse_groups = metrics_model$Mse_groups5)
  }
  
  return(results_keras_wr1)
  
}

# run the cross-validation
tic("Running simulations for the neural network of Biganzoli in keras")

results_keras_wr1 <- keras_simulator_wr1(data = big_list)

time_list_wr1 <- toc()
cat((time_list_wr1$toc - time_list_wr1$tic) / 60, "minutes elapsed", "\n") 
save(results_keras_wr1, file = "results_simulations1/snn_ibs/simulations_keras_wr1.RData")
#save(time_list_wr1, file = "results_simulations1/snn_ibs/simulations_keras_time_wr1.RData")

# enter the list to calculate the measures
brier_keras_wr1 <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
cases_removed <- vector(mode = "numeric", length = length(big_list))
ibs_keras_wr1 <- vector(mode = "numeric", length = length(big_list))
cindex_keras_wr1 <- vector(mode = "numeric", length = length(big_list))
calib_keras2y_groups_wr1 <- vector(mode = "numeric", length = length(big_list))
calib_keras5y_groups_wr1 <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)){
  brier_keras_wr1[i,] <- results_keras_wr1[[i]]$Brier
  cases_removed[i] <- results_keras_wr1[[i]]$Cases_removed
  ibs_keras_wr1[i] <- results_keras_wr1[[i]]$Int_brier
  cindex_keras_wr1[i] <- results_keras_wr1[[i]]$Cindex
  calib_keras2y_groups_wr1[i] <- results_keras_wr1[[i]]$Calib2y_mse_groups
  calib_keras5y_groups_wr1[i] <- results_keras_wr1[[i]]$Calib5y_mse_groups
}

# cases removed before 2 years
summary(cases_removed)
# brier score at 2 years
summa(brier_keras_wr1[, 3])
# brier score at 5 years
summa(brier_keras_wr1[, 6]) 
# ibs at 5 years
summa(ibs_keras_wr1)
# C-index
summa(cindex_keras_wr1)
# Calibration 2 years
summa(calib_keras2y_groups_wr1)
# calibration 5 years
summa(calib_keras5y_groups_wr1)


##############################################################################
# Keras simulator B = 1000, n = 250, split sample = 1/2
# patients curtailed at 5 years (administrative censoring) on the training set

options(keras.view_metrics = FALSE)

keras_simulator_wr2 <- function(data) {
  
  results_keras_wr2 <- vector(mode = "list", length = length(data))
  
  for (a in 1:length(data)) {
    
    dataset <- data[[a]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    # curtail patients at 5 years on the training data
    cases_curtailed <- sum(training_ovr$survt_ovr >= 5)
    training_ovr$surv_status <- training_ovr$surv_status*(training_ovr$survt_ovr <= 5)
    training_ovr$survt_ovr <- pmin(training_ovr$survt_ovr, 5)
    
    # scaling the training and test data 
    # creating training and test ovr scaled
    response_pair <- c("survt_ovr", "surv_status")
    # for the training set
    X_train_scaled <- training_ovr[, !(colnames(training_ovr) %in% response_pair)]
    
    for(i in 1:length(colnames(X_train_scaled))) {
      if(class(X_train_scaled[,i]) == "numeric" || class(X_train_scaled[,i]) == "integer") {
        X_train_scaled[, i] <- as.vector(scale(X_train_scaled[, i])) }
    }
    X_train_ovr <- model.matrix(~., X_train_scaled)[, -1] 
    Y_train_ovr <- training_ovr[, colnames(training_ovr) %in% response_pair] # create matrix Y
    training_ovr_scaled <- as.data.frame(cbind(X_train_ovr, Y_train_ovr))
    
    # for the test set
    x_test_scaled <- test_ovr[, !(colnames(test_ovr) %in% response_pair)]
    for(i in 1:length(colnames(x_test_scaled))) {
      if(class(x_test_scaled[,i]) == "numeric" || class(x_test_scaled[,i]) == "integer") {
        x_test_scaled[, i] <- as.vector(scale(x_test_scaled[, i])) }
    }
    X_test_ovr <- model.matrix(~., x_test_scaled)[, -1]
    Y_test_ovr <- test_ovr[, colnames(test_ovr) %in% response_pair]
    test_ovr_scaled <- as.data.frame(cbind(X_test_ovr, Y_test_ovr))
    
    
    training_ovr_long <- data_train_creator(training_ovr_scaled) 
    validation_ovr_long <- data_test_creator(training_ovr_scaled)
    test_ovr_long <- data_test_creator(test_ovr_scaled, enforce_intervals = TRUE)
    
    k_clear_session() # to avoid clutter from old models / layers in cross validation
    
    tensorflow::tf$random$set_seed(12345) # for update to TensorFlow 2.0
    
    # model with one hidden layer
    
    fit_keras1 <- keras_model_sequential()
    # Add layers to the model
    # here we have logistic activation function for the inputs but also for the outputs
    # we create a densely connected ANN to the output
    
    fit_keras1 %>%
      layer_dense(units = 16, activation = 'tanh',
                  input_shape = c(5 + max(training_ovr_long$interval))) %>%
      layer_dropout(rate = 0.1) %>%
      layer_dense(units = 1, activation = 'sigmoid')
    fit_keras1 %>% compile(
      loss = 'binary_crossentropy', # for binary class classification problem
      optimizer = optimizer_sgd(lr = 0.2, momentum = 0.8)
    )
    early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 5)
    
    # create the matrices to be used for keras library
    train_x <- as.matrix(training_ovr_long[, c(1:5, 12:(11 + max(training_ovr_long$interval)))])
    dimnames(train_x) <- NULL # the object must have empty dimnames
    train_y <- training_ovr_long$status_long
    validation_x <- as.matrix(validation_ovr_long[, c(1:5, 12:(11 + max(validation_ovr_long$interval)))])
    dimnames(validation_x) <- NULL # the object must have empty dimnames
    validation_y <- validation_ovr_long$status_long
    
    result_keras <- fit_keras1 %>% fit(
      train_x,
      train_y,
      epochs = 50,
      batch_size = 64,
      validation_data = list(validation_x, validation_y),
      class_weight = list("0" = 1, "1" = 1.05),
      callbacks = c(early_stopping),
      verbose = 0 # don't display progress bar
    )
    
    
    
    metrics_model <- measures_calculator(trained_model = fit_keras1,
                                         datashort = test_ovr_scaled,
                                         datalong = test_ovr_long)
    
    if (any(a == seq(10, 1000, by = 10))) {
      cat("The IBS of dataset", a, "is:", round(metrics_model$Integrated_brier, 3), "\n")
    }
    
    results_keras_wr2[[a]] <- list(Cases_curtailed = cases_curtailed, 
                                   Time = 0:5, 
                                   Brier = metrics_model$Brier_scores,
                                   Int_brier = metrics_model$Integrated_brier,
                                   Cindex = metrics_model$Cindex,
                                   Calib2y_mse_groups = metrics_model$Mse_groups2,
                                   Calib5y_mse_groups = metrics_model$Mse_groups5)
  }
  
  return(results_keras_wr2)
  
}

# run the cross-validation
tic("Running simulations for the neural network of Biganzoli in keras")

results_keras_wr2 <- keras_simulator_wr2(data = big_list)

time_list_wr2 <- toc()
cat((time_list_wr2$toc - time_list_wr2$tic) / 60, "minutes elapsed", "\n") 
save(results_keras_wr2, file = "results_simulations1/snn_ibs/simulations_keras_wr2.RData")
#save(time_list_wr2, file = "results_simulations1/snn_ibs/simulations_keras_time_wr2.RData")

# enter the list to calculate the measures
brier_keras_wr2 <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
cases_curtailed <- vector(mode = "numeric", length = length(big_list))
ibs_keras_wr2 <- vector(mode = "numeric", length = length(big_list))
cindex_keras_wr2 <- vector(mode = "numeric", length = length(big_list))
calib_keras2y_groups_wr2 <- vector(mode = "numeric", length = length(big_list))
calib_keras5y_groups_wr2 <- vector(mode = "numeric", length = length(big_list))


for (i in 1:length(big_list)){
  brier_keras_wr2[i,] <- results_keras_wr2[[i]]$Brier
  cases_curtailed[i] <- results_keras_wr2[[i]]$Cases_curtailed
  ibs_keras_wr2[i] <- results_keras_wr2[[i]]$Int_brier
  cindex_keras_wr2[i] <- results_keras_wr2[[i]]$Cindex
  calib_keras2y_groups_wr2[i] <- results_keras_wr2[[i]]$Calib2y_mse_groups
  calib_keras5y_groups_wr2[i] <- results_keras_wr2[[i]]$Calib5y_mse_groups
}

# cases curtailed at 5 years
summary(cases_curtailed)
# brier score at 2 years
summa(brier_keras_wr2[, 3])
# brier score at 5 years
summa(brier_keras_wr2[, 6]) 
# ibs at 5 years
summa(ibs_keras_wr2)
# C-index
summa(cindex_keras_wr2)
# Calibration 2 years
summa(calib_keras2y_groups_wr2)
# calibration 5 years
summa(calib_keras5y_groups_wr2)

## end of the file




######################################################################
######################################################################
# file "simulations1_cindex.R" (core file to run simulations using 
# "R environment.RData" that sets up necessary parameters based on the original BO06 data)

# START OF SIMULATIONS
# censoring times with Weibull distribution as original ~61% censoring
# Simulations for Cindex n = 250 (can be adjusted to 1000 by changing 250 with 1000)
# Similarly simulations can be run for IBS regarding the ML techniques
# by choosing optimally tuned parameters for this measure


# you need to create folders results_simulations1 (used for n = 250)
# and results_simulations2 (used for n = 1000)
# and subfolders results_simulations1/snn_ibs, results_simulations2/snn_ibs
# to store the results for all methods in the destination of this R file

######################################################################
######################################################################

## start of the file

# Initialises the R environment (R objects) to generate all data
# during this study based on the original BO06 data (the BO06 dataset is not required)
load("R environment.RData")

# To install packages
# install_packages <- c("dplyr", "data.table", "tictoc", "survival",
#                        "fastDummies", "nnet", "pec", "caret", "ggpubr")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

# how to set up keras in R
# first install anaconda 3
# in anaconda command prompt use the following:
# pip install tensorflow
# pip install keras

# test by running: import tensorflow as tf

# this is done
# install.packages("reticulate")
# install.packages("tensorflow")
# install.packages("keras")

# library(tensorflow)
# install_tensorflow()
# library(keras)
# install_keras()

library(dplyr)
library(data.table)
library(tictoc)
library(survival)
library(fastDummies)
library(nnet)
library(keras)
library(pec)
library(caret)
library(ggpubr)
source("the_functions2.R")
source("functions_nn.R")



# part 1: count unique combinations of categorical variables
# use the mean age for this case to simulate for variable age

# tabel <-  bo06 %>% count(trt, sex, hist_resp, excision)
# tabel$prop <- tabel$n / nrow(bo06)
# 
# options(warn = -1)
# stats <- bo06 %>%
#   group_by(trt, sex, hist_resp, excision) %>%
#   summarise_at(vars(age), funs(mean = "mean", sd = "sd"))
# options(warn = 0)
# 
# tabel$age_mean <- stats$mean
# tabel$age_sd <- stats$sd
# tabel$case <- 1:nrow(tabel)

# create 16 lists with the different scenarios
# lists <- split(tabel[, c("trt", "sex", "hist_resp", "excision",
#                          "age_mean", "age_sd")], f = tabel$case) 

# to simulate survival and censoring times
# we need to calculate betas, mu and sigma of the original bo06 data
# assumed distribution for the y variable is loglogistic
# logreg <-  survreg(Surv(survt_ovr, surv_status) ~
#                      trt + sex + hist_resp + excision + age,
#                    data = bo06, dist = "lognormal",
#                    robust = TRUE)
# 
# coefs <- matrix(logreg$coefficients) # intercept is parameter mu
# sigma <- as.vector(logreg$scale)


# data simulation procedure

data_simulator <- function(n = 250) {
  
  ind <- sample(x = 1:nrow(tabel), size = n, prob = tabel$prop, replace = TRUE)
  data_list <- vector(mode = "list", length = length(ind))
  
  for (i in 1:length(ind)) {
    
    value <- ind[i]
    line <- lists[[value]]
    line$age <- rnorm(1, mean = line$age_mean, sd = line$age_sd)
    line$age <- ifelse(line$age < 3.600757, 3.600757, line$age) 
    line$age <- ifelse(line$age > 40.85132, 40.85132, line$age) 
    line <- line[, c("trt", "sex", "hist_resp", "excision", "age")] # keep variables of interest
    data_list[[i]] <- line
  }
  
  data_sim <- as.data.frame(rbindlist(data_list))
  
  # create the model matrix for the simulated data
  mat_data <- model.matrix(~., data_sim[, c("trt", "sex", "hist_resp", "excision", "age")])
  
  # calculate the survival time
  time <- exp(mat_data %*% coefs + sigma*rnorm(n = n, 0, 1))
  # shape and scale parameters calculated with a Weibull distribution for censoring times of BO06
  cens <- rweibull(n = n, shape = 2.0287, scale = 5.7177) 
  cens[cens == 0] <- (1 / 365.25)
  surv_status <- as.numeric(time <= cens)
  survt_ovr <- pmin(time, cens)
  data_sim <- cbind(data_sim, survt_ovr, surv_status, time)
  
  return(data_sim)
}


data1 <- data_simulator(n = 250)

# summary(bo06)
summary(data1)

table(data1$surv_status)[1] / nrow(data1)


################################################################
# repeat the simulation 1000 times
set.seed(12345)
big_list <- replicate(n = 1000, {
  list(data_simulator(n = 250))
})

# censoring for the big list
cens_perc <- vector(mode = "numeric", length = length(big_list))
#surv_time <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)) {
  datas <- big_list[[i]]
  cens_perc[i] <- table(datas$surv_status)[1] / nrow(datas)
  #surv_time[i] <- max(datas$survt_ovr)
}

summary(cens_perc) # mean around 60.9 %



# load results 
# load("results_simulations1/snn_cindex/simulations_nnet.RData")
# load("results_simulations1/snn_cindex/simulations_keras.RData")



###############################################################
###############################################################
###############################################################
# apply the method of Biganzoli on nnet library
# function that runs a NN and calculates Brier and Integrated Brier score 
# for all simulated datasets

nnet_simulator <- function(data) {
  
  results_nnet <- vector(mode = "list", length = length(data))
  
  for (a in 1:length(data)) {
    
    dataset <- data[[a]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    # scaling the training and test data 
    # creating training and test ovr scaled
    response_pair <- c("survt_ovr", "surv_status")
    # for the training set
    X_train_scaled <- training_ovr[, !(colnames(training_ovr) %in% response_pair)]
    
    for(i in 1:length(colnames(X_train_scaled))) {
      if(class(X_train_scaled[,i]) == "numeric" || class(X_train_scaled[,i]) == "integer") {
        X_train_scaled[, i] <- as.vector(scale(X_train_scaled[, i])) }
    }
    X_train_ovr <- model.matrix(~., X_train_scaled)[, -1] 
    Y_train_ovr <- training_ovr[, colnames(training_ovr) %in% response_pair] # create matrix Y
    training_ovr_scaled <- as.data.frame(cbind(X_train_ovr, Y_train_ovr))
    
    # for the test set
    x_test_scaled <- test_ovr[, !(colnames(test_ovr) %in% response_pair)]
    for(i in 1:length(colnames(x_test_scaled))) {
      if(class(x_test_scaled[,i]) == "numeric" || class(x_test_scaled[,i]) == "integer") {
        x_test_scaled[, i] <- as.vector(scale(x_test_scaled[, i])) }
    }
    X_test_ovr <- model.matrix(~., x_test_scaled)[, -1]
    Y_test_ovr <- test_ovr[, colnames(test_ovr) %in% response_pair]
    test_ovr_scaled <- as.data.frame(cbind(X_test_ovr, Y_test_ovr))
    
    training_ovr_long <- data_train_creator(training_ovr_scaled) # 5 + 8 = 13 prognostic variables
    test_ovr_long <- data_test_creator(test_ovr_scaled)
    
    # create the matrices to be used for keras library
    train_x <- as.matrix(training_ovr_long[, c(1:5, 12 + max(training_ovr_long$interval))])
    dimnames(train_x) <- NULL # the object must have empty dimnames
    train_y <- training_ovr_long$status_long
    
    set.seed(12345)  
    
    # start building the model
    fit_nnet <- nnet(x = train_x, y = train_y, size = 3, decay = 0.1,
                     maxit = 1000, trace = FALSE, entropy = TRUE)
    
    metrics_model <- measures_calculator2(trained_model = fit_nnet,
                                          datashort = test_ovr_scaled,
                                          datalong = test_ovr_long)
    
    if (any(a == seq(20, 1000, by = 20))) {
      cat("The C-index of dataset", a, "is:", round(metrics_model$Cindex, 3), "\n")
    }
    
    results_nnet[[a]] <- list(Time = 0:5, 
                              Brier = metrics_model$Brier_scores,
                              Int_brier = metrics_model$Integrated_brier,
                              Cindex = metrics_model$Cindex,
                              Calib2y_mse_groups = metrics_model$Mse_groups2,
                              Calib5y_mse_groups = metrics_model$Mse_groups5)
  }
  
  return(results_nnet)
  
}

# run the cross-validation
tic("Running simulations for the neural network of Biganzoli in keras")

results_nnet <- nnet_simulator(data = big_list)

time_list <- toc()
cat((time_list$toc - time_list$tic) / 60, "minutes elapsed", "\n") 
save(results_nnet, file = "results_simulations1/snn_cindex/simulations_nnet.RData")
save(time_list, file = "results_simulations1/snn_cindex/simulations_nnet_time.RData")

# enter the list to calculate the measures
brier_nnet <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
ibs_nnet <- vector(mode = "numeric", length = length(big_list))
cindex_nnet <- vector(mode = "numeric", length = length(big_list))
calib_nnet2y_groups <- vector(mode = "numeric", length = length(big_list))
calib_nnet5y_groups <- vector(mode = "numeric", length = length(big_list))


for (i in 1:length(big_list)){
  brier_nnet[i,] <- results_nnet[[i]]$Brier
  ibs_nnet[i] <- results_nnet[[i]]$Int_brier
  cindex_nnet[i] <- results_nnet[[i]]$Cindex
  calib_nnet2y_groups[i] <- results_nnet[[i]]$Calib2y_mse_groups
  calib_nnet5y_groups[i] <- results_nnet[[i]]$Calib5y_mse_groups
}

# brier score at 2 years
summa(brier_nnet[, 3])
# brier score at 5 years
summa(brier_nnet[, 6]) 
# ibs at 5 years
summa(ibs_nnet)
# C-index
summa(cindex_nnet)
# Calibration 2 years
summa(calib_nnet2y_groups)
# calibration 5 years
summa(calib_nnet5y_groups)




###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# apply the method of Biganzoli on keras library
# function that runs a NN and calculates Brier and Integrated Brier score 
# for all simulated datasets
# set global default to never show metrics
options(keras.view_metrics = FALSE)

keras_simulator <- function(data) {
  
  results_keras <- vector(mode = "list", length = length(data))
  
  for (a in 1:length(data)) {
    
    dataset <- data[[a]][, 1:7]
    N <- nrow(dataset)
    set.seed(12345)
    index <- createDataPartition(dataset$surv_status, p = 1/2, list = FALSE)
    training_ovr <- dataset[index, ]
    test_ovr  <- dataset[-index,]
    
    # scaling the training and test data 
    # creating training and test ovr scaled
    response_pair <- c("survt_ovr", "surv_status")
    # for the training set
    X_train_scaled <- training_ovr[, !(colnames(training_ovr) %in% response_pair)]
    
    for(i in 1:length(colnames(X_train_scaled))) {
      if(class(X_train_scaled[,i]) == "numeric" || class(X_train_scaled[,i]) == "integer") {
        X_train_scaled[, i] <- as.vector(scale(X_train_scaled[, i])) }
    }
    X_train_ovr <- model.matrix(~., X_train_scaled)[, -1] 
    Y_train_ovr <- training_ovr[, colnames(training_ovr) %in% response_pair] # create matrix Y
    training_ovr_scaled <- as.data.frame(cbind(X_train_ovr, Y_train_ovr))
    
    # for the test set
    x_test_scaled <- test_ovr[, !(colnames(test_ovr) %in% response_pair)]
    for(i in 1:length(colnames(x_test_scaled))) {
      if(class(x_test_scaled[,i]) == "numeric" || class(x_test_scaled[,i]) == "integer") {
        x_test_scaled[, i] <- as.vector(scale(x_test_scaled[, i])) }
    }
    X_test_ovr <- model.matrix(~., x_test_scaled)[, -1]
    Y_test_ovr <- test_ovr[, colnames(test_ovr) %in% response_pair]
    test_ovr_scaled <- as.data.frame(cbind(X_test_ovr, Y_test_ovr))
    
    training_ovr_long <- data_train_creator(training_ovr_scaled) # 5 + 8 = 13 prognostic variables
    validation_ovr_long <- data_test_creator(training_ovr_scaled)
    test_ovr_long <- data_test_creator(test_ovr_scaled)
    
    k_clear_session() # to avoid clutter from old models / layers in cross validation
    
    tensorflow::tf$random$set_seed(12345) # for update to TensorFlow 2.0
    
    fit_keras1 <- keras_model_sequential()
    # Add layers to the model
    # here we have logistic activation function for the inputs but also for the outputs
    # we create a densely connected ANN to the output
    
    
    fit_keras1 %>%
      layer_dense(units = 13, activation = 'relu',
                  input_shape = c(5 + max(training_ovr_long$interval))) %>%
      layer_dropout(rate = 0.1) %>%
      layer_dense(units = 1, activation = 'sigmoid')
    fit_keras1 %>% compile(
      loss = 'binary_crossentropy', # for binary class classification problem
      optimizer = optimizer_sgd(lr = 0.2, momentum = 0.9)
    )
    early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 5)
    
    # create the matrices to be used for keras library
    train_x <- as.matrix(training_ovr_long[, c(1:5, 12:(11 + max(training_ovr_long$interval)))])
    dimnames(train_x) <- NULL # the object must have empty dimnames
    train_y <- training_ovr_long$status_long
    validation_x <- as.matrix(validation_ovr_long[, c(1:5, 12:(11 + max(validation_ovr_long$interval)))])
    dimnames(validation_x) <- NULL # the object must have empty dimnames
    validation_y <- validation_ovr_long$status_long
    
    result_keras <- fit_keras1 %>% fit(
      train_x,
      train_y,
      epochs = 50,
      batch_size = 64,
      validation_data = list(validation_x, validation_y),
      class_weight = list("0" = 1, "1" = 1),
      callbacks = c(early_stopping),
      verbose = 0 # don't display progress bar
    )
    
    
    metrics_model <- measures_calculator(trained_model = fit_keras1,
                                         datashort = test_ovr_scaled,
                                         datalong = test_ovr_long)
    
    if (any(a == seq(10, 1000, by = 10))) {
      cat("The C-index of dataset", a, "is:", round(metrics_model$Cindex, 3), "\n")
    }
    
    results_keras[[a]] <- list(Time = 0:5, 
                               Brier = metrics_model$Brier_scores,
                               Int_brier = metrics_model$Integrated_brier,
                               Cindex = metrics_model$Cindex,
                               Calib2y_mse_groups = metrics_model$Mse_groups2,
                               Calib5y_mse_groups = metrics_model$Mse_groups5)
  }
  
  return(results_keras)
  
}

# run the cross-validation
tic("Running simulations for the neural network of Biganzoli in keras")

results_keras <- keras_simulator(data = big_list)

time_list <- toc()
cat((time_list$toc - time_list$tic) / 60, "minutes elapsed", "\n") 
save(results_keras, file = "results_simulations1/snn_cindex/simulations_keras.RData")
save(time_list, file = "results_simulations1/snn_cindex/simulations_keras_time.RData")

# enter the list to calculate the measures
brier_keras <- matrix(data = 0, nrow = length(big_list), ncol = 6) # Brier score until 5 years
ibs_keras <- vector(mode = "numeric", length = length(big_list))
cindex_keras <- vector(mode = "numeric", length = length(big_list))
calib_keras2y_groups <- vector(mode = "numeric", length = length(big_list))
calib_keras5y_groups <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)){
  brier_keras[i,] <- results_keras[[i]]$Brier
  ibs_keras[i] <- results_keras[[i]]$Int_brier
  cindex_keras[i] <- results_keras[[i]]$Cindex
  calib_keras2y_groups[i] <- results_keras[[i]]$Calib2y_mse_groups
  calib_keras5y_groups[i] <- results_keras[[i]]$Calib5y_mse_groups
}

# brier score at 2 years
summa(brier_keras[, 3])
# brier score at 5 years
summa(brier_keras[, 6]) 
# ibs at 5 years
summa(ibs_keras)
# C-index
summa(cindex_keras)
# Calibration 2 years
summa(calib_keras2y_groups)
# calibration 5 years
summa(calib_keras5y_groups)

## end of the file






######################################################################
######################################################################
# file "example_20perc.R"
# censoring times with Weibull distribution  ~20% censoring
# Simulations for n = 1000 

######################################################################
######################################################################

## start of the file

# Initialises the R environment (R objects) to generate all data
# during this study based on the original BO06 data (the BO06 dataset is not required)
load("R environment.RData")

# To install packages
# install_packages <- c("dplyr", "data.table", "tictoc", "survival",
#                        "fastDummies", "nnet", "pec", "caret", "ggpubr")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

# how to set up keras in R
# first install anaconda 3
# in anaconda command prompt use the following:
# pip install tensorflow
# pip install keras

# test by running: import tensorflow as tf

# this is done
# install.packages("reticulate")
# install.packages("tensorflow")
# install.packages("keras")

# library(tensorflow)
# install_tensorflow()
# library(keras)
# install_keras()

library(dplyr)
library(data.table)
library(tictoc)
library(survival)
library(fastDummies)
library(nnet)
library(keras)
library(pec)
library(caret)
library(ggpubr)
source("the_functions2.R")
source("functions_nn.R")


# data simulation procedure

data_simulator <- function(n = 1000) {
  
  ind <- sample(x = 1:nrow(tabel), size = n, prob = tabel$prop, replace = TRUE)
  data_list <- vector(mode = "list", length = length(ind))
  
  for (i in 1:length(ind)) {
    
    value <- ind[i]
    line <- lists[[value]]
    line$age <- rnorm(1, mean = line$age_mean, sd = line$age_sd)
    line$age <- ifelse(line$age < 3.600757, 3.600757, line$age) 
    line$age <- ifelse(line$age > 40.85132, 40.85132, line$age) 
    line <- line[, c("trt", "sex", "hist_resp", "excision", "age")] # keep variables of interest
    data_list[[i]] <- line
  }
  
  data_sim <- as.data.frame(rbindlist(data_list))
  
  # create the model matrix for the simulated data
  mat_data <- model.matrix(~., data_sim[, c("trt", "sex", "hist_resp", "excision", "age")])
  
  # calculate the survival time
  time <- exp(mat_data %*% coefs + sigma*rnorm(n = n, 0, 1))
  # shape and scale parameters calculated with a Weibull distribution for censoring times of BO06
  cens <- rweibull(n = n, shape = 0.75, scale = 76) 
  cens[cens == 0] <- (1 / 365.25)
  surv_status <- as.numeric(time <= cens)
  survt_ovr <- pmin(time, cens)
  data_sim <- cbind(data_sim, survt_ovr, surv_status, time)
  
  return(data_sim)
}


data1 <- data_simulator(n = 1000)

# summary(bo06)
summary(data1)

table(data1$surv_status)[1] / nrow(data1)


################################################################
# repeat the simulation 1000 times
set.seed(12345)
big_list <- replicate(n = 1000, {
  list(data_simulator(n = 1000))
})

# censoring for the big list
cens_perc <- vector(mode = "numeric", length = length(big_list))
#surv_time <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)) {
  datas <- big_list[[i]]
  cens_perc[i] <- table(datas$surv_status)[1] / nrow(datas)
  #surv_time[i] <- max(datas$survt_ovr)
}

summary(cens_perc) # mean around 20 %

# etc for the functions of the methods as for 61% as original...


## end of the file





######################################################################
######################################################################
# file "example_40perc.R"
# censoring times with Weibull distribution  ~40% censoring
# Simulations for n = 1000 

######################################################################
######################################################################


## start of the file

# Initialises the R environment (R objects) to generate all data
# during this study based on the original BO06 data (the BO06 dataset is not required)
load("R environment.RData")

# To install packages
# install_packages <- c("dplyr", "data.table", "tictoc", "survival",
#                        "fastDummies", "nnet", "pec", "caret", "ggpubr")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

# how to set up keras in R
# first install anaconda 3
# in anaconda command prompt use the following:
# pip install tensorflow
# pip install keras

# test by running: import tensorflow as tf

# this is done
# install.packages("reticulate")
# install.packages("tensorflow")
# install.packages("keras")

# library(tensorflow)
# install_tensorflow()
# library(keras)
# install_keras()

library(dplyr)
library(data.table)
library(tictoc)
library(survival)
library(fastDummies)
library(nnet)
library(keras)
library(pec)
library(caret)
library(ggpubr)
source("the_functions2.R")
source("functions_nn.R")

# data simulation procedure

data_simulator <- function(n = 1000) {
  
  ind <- sample(x = 1:nrow(tabel), size = n, prob = tabel$prop, replace = TRUE)
  data_list <- vector(mode = "list", length = length(ind))
  
  for (i in 1:length(ind)) {
    
    value <- ind[i]
    line <- lists[[value]]
    line$age <- rnorm(1, mean = line$age_mean, sd = line$age_sd)
    line$age <- ifelse(line$age < 3.600757, 3.600757, line$age) 
    line$age <- ifelse(line$age > 40.85132, 40.85132, line$age) 
    line <- line[, c("trt", "sex", "hist_resp", "excision", "age")] # keep variables of interest
    data_list[[i]] <- line
  }
  
  data_sim <- as.data.frame(rbindlist(data_list))
  
  # create the model matrix for the simulated data
  mat_data <- model.matrix(~., data_sim[, c("trt", "sex", "hist_resp", "excision", "age")])
  
  # calculate the survival time
  time <- exp(mat_data %*% coefs + sigma*rnorm(n = n, 0, 1))
  # shape and scale parameters calculated with a Weibull distribution for censoring times of BO06
  cens <- rweibull(n = n, shape = 0.75, scale = 20.5) 
  cens[cens == 0] <- (1 / 365.25)
  surv_status <- as.numeric(time <= cens)
  survt_ovr <- pmin(time, cens)
  data_sim <- cbind(data_sim, survt_ovr, surv_status, time)
  
  return(data_sim)
}


data1 <- data_simulator(n = 1000)

# summary(bo06)
summary(data1)

table(data1$surv_status)[1] / nrow(data1)


################################################################
# repeat the simulation 1000 times
set.seed(12345)
big_list <- replicate(n = 1000, {
  list(data_simulator(n = 1000))
})

# censoring for the big list
cens_perc <- vector(mode = "numeric", length = length(big_list))
#surv_time <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)) {
  datas <- big_list[[i]]
  cens_perc[i] <- table(datas$surv_status)[1] / nrow(datas)
  #surv_time[i] <- max(datas$survt_ovr)
}

summary(cens_perc) # mean around 40 %

# etc for the functions of the methods as for 61% as original...


## end of the file




######################################################################
######################################################################
# file "example_61perc_user.R"
# censoring times with Weibull distribution  ~61% censoring (user defined)
# Simulations for n = 1000 

######################################################################
######################################################################

## start of the file

# Initialises the R environment (R objects) to generate all data
# during this study based on the original BO06 data (the BO06 dataset is not required)
load("R environment.RData")

# To install packages
# install_packages <- c("dplyr", "data.table", "tictoc", "survival",
#                        "fastDummies", "nnet", "pec", "caret", "ggpubr")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

# how to set up keras in R
# first install anaconda 3
# in anaconda command prompt use the following:
# pip install tensorflow
# pip install keras

# test by running: import tensorflow as tf

# this is done
# install.packages("reticulate")
# install.packages("tensorflow")
# install.packages("keras")

# library(tensorflow)
# install_tensorflow()
# library(keras)
# install_keras()

library(dplyr)
library(data.table)
library(tictoc)
library(survival)
library(fastDummies)
library(nnet)
library(keras)
library(pec)
library(caret)
library(ggpubr)
source("the_functions2.R")
source("functions_nn.R")

# data simulation procedure

data_simulator <- function(n = 1000) {
  
  ind <- sample(x = 1:nrow(tabel), size = n, prob = tabel$prop, replace = TRUE)
  data_list <- vector(mode = "list", length = length(ind))
  
  for (i in 1:length(ind)) {
    
    value <- ind[i]
    line <- lists[[value]]
    line$age <- rnorm(1, mean = line$age_mean, sd = line$age_sd)
    line$age <- ifelse(line$age < 3.600757, 3.600757, line$age) 
    line$age <- ifelse(line$age > 40.85132, 40.85132, line$age) 
    line <- line[, c("trt", "sex", "hist_resp", "excision", "age")] # keep variables of interest
    data_list[[i]] <- line
  }
  
  data_sim <- as.data.frame(rbindlist(data_list))
  
  # create the model matrix for the simulated data
  mat_data <- model.matrix(~., data_sim[, c("trt", "sex", "hist_resp", "excision", "age")])
  
  # calculate the survival time
  time <- exp(mat_data %*% coefs + sigma*rnorm(n = n, 0, 1))
  # shape and scale parameters calculated with a Weibull distribution for censoring times of BO06
  cens <- rweibull(n = n, shape = 0.75, scale = 6.8) 
  cens[cens == 0] <- (1 / 365.25)
  surv_status <- as.numeric(time <= cens)
  survt_ovr <- pmin(time, cens)
  data_sim <- cbind(data_sim, survt_ovr, surv_status, time)
  
  return(data_sim)
}


data1 <- data_simulator(n = 1000)

# summary(bo06)
summary(data1)

table(data1$surv_status)[1] / nrow(data1)


################################################################
# repeat the simulation 1000 times
set.seed(12345)
big_list <- replicate(n = 1000, {
  list(data_simulator(n = 1000))
})

# censoring for the big list
cens_perc <- vector(mode = "numeric", length = length(big_list))
#surv_time <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)) {
  datas <- big_list[[i]]
  cens_perc[i] <- table(datas$surv_status)[1] / nrow(datas)
  #surv_time[i] <- max(datas$survt_ovr)
}

summary(cens_perc) # mean around 61 %

# etc for the functions of the methods as for 61% as original...


## end of the file




######################################################################
######################################################################
# file "example_80perc.R"
# censoring times with Weibull distribution  ~80% censoring
# Simulations for n = 1000 

######################################################################
######################################################################

## start of the file

# Initialises the R environment (R objects) to generate all data
# during this study based on the original BO06 data (the BO06 dataset is not required)
load("R environment.RData")

# To install packages
# install_packages <- c("dplyr", "data.table", "tictoc", "survival",
#                        "fastDummies", "nnet", "pec", "caret", "ggpubr")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

# how to set up keras in R
# first install anaconda 3
# in anaconda command prompt use the following:
# pip install tensorflow
# pip install keras

# test by running: import tensorflow as tf

# this is done
# install.packages("reticulate")
# install.packages("tensorflow")
# install.packages("keras")

# library(tensorflow)
# install_tensorflow()
# library(keras)
# install_keras()

library(dplyr)
library(data.table)
library(tictoc)
library(survival)
library(fastDummies)
library(nnet)
library(keras)
library(pec)
library(caret)
library(ggpubr)
source("the_functions2.R")
source("functions_nn.R")

# data simulation procedure

data_simulator <- function(n = 1000) {
  
  ind <- sample(x = 1:nrow(tabel), size = n, prob = tabel$prop, replace = TRUE)
  data_list <- vector(mode = "list", length = length(ind))
  
  for (i in 1:length(ind)) {
    
    value <- ind[i]
    line <- lists[[value]]
    line$age <- rnorm(1, mean = line$age_mean, sd = line$age_sd)
    line$age <- ifelse(line$age < 3.600757, 3.600757, line$age) 
    line$age <- ifelse(line$age > 40.85132, 40.85132, line$age) 
    line <- line[, c("trt", "sex", "hist_resp", "excision", "age")] # keep variables of interest
    data_list[[i]] <- line
  }
  
  data_sim <- as.data.frame(rbindlist(data_list))
  
  # create the model matrix for the simulated data
  mat_data <- model.matrix(~., data_sim[, c("trt", "sex", "hist_resp", "excision", "age")])
  
  # calculate the survival time
  time <- exp(mat_data %*% coefs + sigma*rnorm(n = n, 0, 1))
  # shape and scale parameters calculated with a Weibull distribution for censoring times of BO06
  cens <- rweibull(n = n, shape = 0.75, scale = 2.4) 
  cens[cens == 0] <- (1 / 365.25)
  surv_status <- as.numeric(time <= cens)
  survt_ovr <- pmin(time, cens)
  data_sim <- cbind(data_sim, survt_ovr, surv_status, time)
  
  return(data_sim)
}


data1 <- data_simulator(n = 1000)

# summary(bo06)
summary(data1)

table(data1$surv_status)[1] / nrow(data1)


################################################################
# repeat the simulation 1000 times
set.seed(12345)
big_list <- replicate(n = 1000, {
  list(data_simulator(n = 1000))
})

# censoring for the big list
cens_perc <- vector(mode = "numeric", length = length(big_list))
#surv_time <- vector(mode = "numeric", length = length(big_list))

for (i in 1:length(big_list)) {
  datas <- big_list[[i]]
  cens_perc[i] <- table(datas$surv_status)[1] / nrow(datas)
  #surv_time[i] <- max(datas$survt_ovr)
}

summary(cens_perc) # mean around 80 %

# etc for the functions of the methods as for 61% as original...


## end of the file




######################################################################
######################################################################
# file "lineplots.R" (after analysis for visualisation of the results)

# creates lineplots based on stored results in destination folders
# called results_simulations1 (used for n = 250) and
# results_simulations2 (used for n = 1000)
######################################################################
######################################################################

## start of the file 

# To install packages
# install_packages <- c("ggplot2", "gridExtra")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

library(ggplot2)
library(gridExtra)


# For n = 1000

# load results 
load("results_simulations2/simulations_cox.RData")
load("results_simulations2/snn_ibs/simulations_nnet.RData")
load("results_simulations2/snn_ibs/simulations_keras.RData")

summa <- function(x) {
  
  res1 <- round(quantile(x, probs = c(0, 0.025, 0.5, 0.975, 1), na.rm = TRUE), digits = 3)
  res2 <- round(mean(x, na.rm = TRUE), digits = 3)
  res3 <- round(sd(x, na.rm = TRUE), digits = 3)
  res <- c(res1[1:3], res2, res1[4:5], res3)
  names(res) <- c("Min.", "2.5% Qu.", "Median", "Mean", "97.5% Qu.", "Max.", "Sd.")
  return(res)
}


quant_low <- function(x) {
  res <- mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE)
  return(res)
}

quant_high <- function(x) {
  res <- mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)
  return(res)
}


# for the Cox models

# enter the list to calculate the measures
brier_cox <- matrix(data = 0, nrow = 1000, ncol = 6) # Brier score until 5 years
ibs_cox <- vector(mode = "numeric", length = 1000)
cindex_cox <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  brier_cox[i,] <- results_cox[[i]]$Brier
  ibs_cox[i] <- results_cox[[i]]$Int_brier
  cindex_cox[i] <- results_cox[[i]]$Cindex
}


line1 <- data.frame(times = 0:5,
                    brier = colMeans(brier_cox, na.rm = TRUE),
                    ci_left = apply(brier_cox, 2, quant_low),
                    ci_right = apply(brier_cox, 2, quant_high),
                    Model = "Cox")

# for nnet

brier_nnet <- matrix(data = 0, nrow = 1000, ncol = 6) # Brier score until 5 years
ibs_nnet <- vector(mode = "numeric", length = 1000)
cindex_nnet <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  brier_nnet[i,] <- results_nnet[[i]]$Brier
  ibs_nnet[i] <- results_nnet[[i]]$Int_brier
  cindex_nnet[i] <- results_nnet[[i]]$Cindex

}




line2 <- data.frame(times = 0:5,
                    brier = colMeans(brier_nnet, na.rm = TRUE),
                    ci_left = apply(brier_nnet, 2, quant_low),
                    ci_right = apply(brier_nnet, 2, quant_high),
                    Model = "PLANN original")


# for keras
# enter the list to calculate the measures
brier_keras <- matrix(data = 0, nrow = 1000, ncol = 6) # Brier score until 5 years
ibs_keras <- vector(mode = "numeric", length = 1000)
cindex_keras <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  brier_keras[i,] <- results_keras[[i]]$Brier
  ibs_keras[i] <- results_keras[[i]]$Int_brier
  cindex_keras[i] <- results_keras[[i]]$Cindex
}



line3 <- data.frame(times = 0:5,
                    brier = colMeans(brier_keras, na.rm = TRUE),
                    ci_left = apply(brier_keras, 2, quant_low),
                    ci_right = apply(brier_keras, 2, quant_high),
                    Model = "PLANN extended")



df <- rbind(line1, line2, line3)

df$Model <- as.factor(df$Model)
df$Model <- factor(df$Model, levels = c("Cox", "PLANN original", "PLANN extended"))


brier_models1000 <- ggplot(df, aes(x = times, y = brier, color = Model)) + 
  geom_line(size = 0.9) + 
  geom_point() + 
  xlab("Time in years since surgery") + 
  ylab("Brier score") + 
  scale_x_continuous(breaks = 0:5, limits = c(0, 5)) + 
  ylim(c(0, 0.30)) + 
  theme_classic() + 
  ggtitle("N2 = 1000") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_ribbon(aes(ymin = ci_left, ymax = ci_right), linetype = 2, alpha = 0.05) + 
  # scale_color_manual(values=c("red", "blue", "green"), 
  #                    breaks=c("Cox","PLANN original","PLANN extended")) + 
  theme(legend.position = "top")


###############################################################################################


# for n = 250

# load results 
load("results_simulations1/simulations_cox.RData")
load("results_simulations1/snn_ibs/simulations_nnet.RData")
load("results_simulations1/snn_ibs/simulations_keras.RData")


# for the Cox models

# enter the list to calculate the measures
brier_cox <- matrix(data = 0, nrow = 1000, ncol = 6) # Brier score until 5 years
ibs_cox <- vector(mode = "numeric", length = 1000)
cindex_cox <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  brier_cox[i,] <- results_cox[[i]]$Brier
  ibs_cox[i] <- results_cox[[i]]$Int_brier
  cindex_cox[i] <- results_cox[[i]]$Cindex
}


line4 <- data.frame(times = 0:5,
                    brier = colMeans(brier_cox, na.rm = TRUE),
                    ci_left = apply(brier_cox, 2, quant_low),
                    ci_right = apply(brier_cox, 2, quant_high),
                    Model = "Cox")

# for nnet

brier_nnet <- matrix(data = 0, nrow = 1000, ncol = 6) # Brier score until 5 years
ibs_nnet <- vector(mode = "numeric", length = 1000)
cindex_nnet <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  brier_nnet[i,] <- results_nnet[[i]]$Brier
  ibs_nnet[i] <- results_nnet[[i]]$Int_brier
  cindex_nnet[i] <- results_nnet[[i]]$Cindex
}


line5 <- data.frame(times = 0:5,
                    brier = colMeans(brier_nnet, na.rm = TRUE),
                    ci_left = apply(brier_nnet, 2, quant_low),
                    ci_right = apply(brier_nnet, 2, quant_high),
                    Model = "PLANN original")


# for keras
# enter the list to calculate the measures
brier_keras <- matrix(data = 0, nrow = 1000, ncol = 6) # Brier score until 5 years
ibs_keras <- vector(mode = "numeric", length = 1000)
cindex_keras <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  brier_keras[i,] <- results_keras[[i]]$Brier
  ibs_keras[i] <- results_keras[[i]]$Int_brier
  cindex_keras[i] <- results_keras[[i]]$Cindex
}



line6 <- data.frame(times = 0:5,
                    brier = colMeans(brier_keras, na.rm = TRUE),
                    ci_left = apply(brier_keras, 2, quant_low),
                    ci_right = apply(brier_keras, 2, quant_high),
                    Model = "PLANN extended")



df2 <- rbind(line4, line5, line6)

df2$Model <- as.factor(df2$Model)
df2$Model <- factor(df2$Model, levels = c("Cox", "PLANN original", "PLANN extended"))


brier_models250 <- ggplot(df2, aes(x = times, y = brier, color = Model)) + 
  geom_line(size = 0.9) + 
  geom_point() + 
  xlab("Time in years since surgery") + 
  ylab("Brier score") + 
  scale_x_continuous(breaks = 0:5, limits = c(0, 5)) + 
  ylim(c(0, 0.30)) + 
  theme_classic() + 
  ggtitle("N1 = 250") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_ribbon(aes(ymin = ci_left, ymax = ci_right), linetype = 2, alpha = 0.05) + 
  # scale_color_manual(values=c("red", "blue", "green"), 
  #                    breaks=c("Cox","PLANN original","PLANN extended")) + 
  theme(legend.position = "top")


# figure 1 of the manuscript
gridExtra::grid.arrange(brier_models250, brier_models1000, ncol = 2)

## end of the file




######################################################################
######################################################################
# file "boxplots.R" (after analysis for visualisation of the results)

# creates boxplots based on stored results in destination folders
# called results_simulations1 (used for n = 250) and
# results_simulations2 (used for n = 1000)
######################################################################
######################################################################

## start of the file

# To install packages
# install_packages <- c("ggplot2", "gridExtra")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

library(ggplot2)
library(gridExtra)


# For n = 1000

# load results 
load("results_simulations2/simulations_cox.RData")
load("results_simulations2/snn_ibs/simulations_nnet.RData")
load("results_simulations2/snn_ibs/simulations_keras.RData")


quant_low <- function(x) {
  res <- mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE)
  return(res)
}

quant_high <- function(x) {
  res <- mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)
  return(res)
}


# for the Cox models

# enter the list to calculate the measures
calib_cox2y_groups <- vector(mode = "numeric", length = 1000)
calib_cox5y_groups <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  calib_cox2y_groups[i] <- results_cox[[i]]$Calib2y$mse_groups
  calib_cox5y_groups[i] <- results_cox[[i]]$Calib5y$mse_groups
  
}



df1 <- data.frame(Calib2y_groups = calib_cox2y_groups,
                  Calib5y_groups = calib_cox5y_groups,
                  Model = "Cox",
                  Sample_size = "1000")

# for nnet

calib_nnet2y_groups <- vector(mode = "numeric", length = 1000)
calib_nnet5y_groups <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  calib_nnet2y_groups[i] <- results_nnet[[i]]$Calib2y_mse_groups
  calib_nnet5y_groups[i] <- results_nnet[[i]]$Calib5y_mse_groups
}

df2 <- data.frame(Calib2y_groups = calib_nnet2y_groups,
                  Calib5y_groups = calib_nnet5y_groups,
                  Model = "PLANN original",
                  Sample_size = "1000")

# for keras
# enter the list to calculate the measures
calib_keras2y_groups <- vector(mode = "numeric", length = 1000)
calib_keras5y_groups <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  calib_keras2y_groups[i] <- results_keras[[i]]$Calib2y_mse_groups
  calib_keras5y_groups[i] <- results_keras[[i]]$Calib5y_mse_groups
}

df3 <- data.frame(Calib2y_groups = calib_keras2y_groups,
                  Calib5y_groups = calib_keras5y_groups,
                  Model = "PLANN extended",
                  Sample_size = "1000")



# for n = 250

# load results 
load("results_simulations1/simulations_cox.RData")
load("results_simulations1/snn_ibs/simulations_nnet.RData")
load("results_simulations1/snn_ibs/simulations_keras.RData")



# for the Cox models

# enter the list to calculate the measures
calib_cox2y_groups <- vector(mode = "numeric", length = 1000)
calib_cox5y_groups <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  calib_cox2y_groups[i] <- results_cox[[i]]$Calib2y$mse_groups
  calib_cox5y_groups[i] <- results_cox[[i]]$Calib5y$mse_groups
}

df4 <- data.frame(Calib2y_groups = calib_cox2y_groups,
                  Calib5y_groups = calib_cox5y_groups,
                  Model = "Cox",
                  Sample_size = "250")


# for nnet

calib_nnet2y_groups <- vector(mode = "numeric", length = 1000)
calib_nnet5y_groups <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  calib_nnet2y_groups[i] <- results_nnet[[i]]$Calib2y_mse_groups
  calib_nnet5y_groups[i] <- results_nnet[[i]]$Calib5y_mse_groups
}

df5 <- data.frame(Calib2y_groups = calib_nnet2y_groups,
                  Calib5y_groups = calib_nnet5y_groups,
                  Model = "PLANN original",
                  Sample_size = "250")


# for keras
# enter the list to calculate the measures
calib_keras2y_groups <- vector(mode = "numeric", length = 1000)
calib_keras5y_groups <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  calib_keras2y_groups[i] <- results_keras[[i]]$Calib2y_mse_groups
  calib_keras5y_groups[i] <- results_keras[[i]]$Calib5y_mse_groups
}


df6 <- data.frame(Calib2y_groups = calib_keras2y_groups,
                  Calib5y_groups = calib_keras5y_groups,
                  Model = "PLANN extended",
                  Sample_size = "250")


df <- rbind(df4, df5, df6, df1, df2, df3)
df$Model <- as.factor(df$Model)
df$Model <- factor(df$Model, levels = c("Cox", "PLANN original", "PLANN extended"))
df$Sample_size <- as.factor(df$Sample_size)
df$Sample_size <- factor(df$Sample_size, levels = c("250", "1000"))





####################################################################################

# boxplots for miscalibration


calib1 <- ggplot(df, aes(x=Sample_size, y=Calib2y_groups, fill = Model)) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_brewer(palette="Reds") +  #, breaks = c("250", "1000")) +
  ylim(c(0, 0.05)) + 
  labs(title = "t = 2 years",x = "Sample size", y = "Miscalibration (MSE for 4 groups)", fill = " ") +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5)) # Remove legend by setting to "none" 


calib2 <- ggplot(df, aes(x=Sample_size, y=Calib5y_groups, fill = Model)) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_brewer(palette="Reds") +  #, breaks = c("250", "1000")) +
  ylim(c(0, 0.05)) + 
  labs(title = "t = 5 years",x = "Sample size", y = "Miscalibration (MSE for 4 groups)", fill = " ") +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5)) # Remove legend by setting to "none" 


# figure 3 of the manuscript
gridExtra::grid.arrange(calib1, calib2, ncol = 2)



####################################################################################

# boxplots for miscalibration


calib1 <- ggplot(df, aes(x=Sample_size, y=Calib2y_groups, fill = Model)) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_brewer(palette="Reds") +  #, breaks = c("250", "1000")) +
  ylim(c(0, 0.05)) + 
  labs(title = "t = 2 years",x = "Sample size", y = "Miscalibration (MSE for 4 groups)", fill = " ") +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5)) # Remove legend by setting to "none" 


calib2 <- ggplot(df, aes(x=Sample_size, y=Calib5y_groups, fill = Model)) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_brewer(palette="Reds") +  #, breaks = c("250", "1000")) +
  ylim(c(0, 0.05)) + 
  labs(title = "t = 5 years",x = "Sample size", y = "Miscalibration (MSE for 4 groups)", fill = " ") +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5)) # Remove legend by setting to "none" 



gridExtra::grid.arrange(calib1, calib2, ncol = 2)


## end of the file




######################################################################
######################################################################
# file "bargraphs.R" (after analysis for visualisation of the results)

# creates bargraphs based on stored results in destination folders
# called results_simulations1 (used for n = 250) and
# results_simulations2 (used for n = 1000)
######################################################################
######################################################################

## start of the file

# To install packages
# install_packages <- c("ggplot2", "gridExtra")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }


library(ggplot2)
library(gridExtra)


quant_low <- function(x) {
  res <- mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE)
  return(res)
}

quant_high <- function(x) {
  res <- mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)
  return(res)
}




# for n = 250

# load results 
load("results_simulations1/simulations_cox.RData")
load("results_simulations1/snn_ibs/simulations_nnet.RData")
load("results_simulations1/snn_ibs/simulations_keras.RData")



# for the Cox models

# enter the list to calculate the measures
ibs_cox <- vector(mode = "numeric", length = 1000)
cindex_cox <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  ibs_cox[i] <- results_cox[[i]]$Int_brier
  cindex_cox[i] <- results_cox[[i]]$Cindex
}



line1 <- data.frame(ibs = mean(ibs_cox, na.rm = TRUE),
                    ibs_left = quant_low(ibs_cox),
                    ibs_right = quant_high(ibs_cox),
                    cindex = mean(cindex_cox, na.rm = TRUE),
                    cindex_left = quant_low(cindex_cox),
                    cindex_right = quant_high(cindex_cox),
                    Model = "Cox",
                    Sample_size = "250")



# for nnet

ibs_nnet <- vector(mode = "numeric", length = 1000)
cindex_nnet <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_nnet[i] <- results_nnet[[i]]$Int_brier
  cindex_nnet[i] <- results_nnet[[i]]$Cindex
}



line2 <- data.frame(ibs = mean(ibs_nnet, na.rm = TRUE),
                    ibs_left = quant_low(ibs_nnet),
                    ibs_right = quant_high(ibs_nnet),
                    cindex = mean(cindex_nnet, na.rm = TRUE),
                    cindex_left = quant_low(cindex_nnet),
                    cindex_right = quant_high(cindex_nnet),
                    Model = "PLANN original",
                    Sample_size = "250")


# for keras
# enter the list to calculate the measures
ibs_keras <- vector(mode = "numeric", length = 1000)
cindex_keras <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  ibs_keras[i] <- results_keras[[i]]$Int_brier
  cindex_keras[i] <- results_keras[[i]]$Cindex
}


line3 <- data.frame(ibs = mean(ibs_keras, na.rm = TRUE),
                    ibs_left = quant_low(ibs_keras),
                    ibs_right = quant_high(ibs_keras),
                    cindex = mean(cindex_keras, na.rm = TRUE),
                    cindex_left = quant_low(cindex_keras),
                    cindex_right = quant_high(cindex_keras),
                    Model = "PLANN extended",
                    Sample_size = "250")


# For n = 1000

# load results 
load("results_simulations2/simulations_cox.RData")
load("results_simulations2/snn_ibs/simulations_nnet.RData")
load("results_simulations2/snn_ibs/simulations_keras.RData")


# enter the list to calculate the measures
ibs_cox <- vector(mode = "numeric", length = 1000)
cindex_cox <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  ibs_cox[i] <- results_cox[[i]]$Int_brier
  cindex_cox[i] <- results_cox[[i]]$Cindex
}


line4 <- data.frame(ibs = mean(ibs_cox, na.rm = TRUE),
                    ibs_left = quant_low(ibs_cox),
                    ibs_right = quant_high(ibs_cox),
                    cindex = mean(cindex_cox, na.rm = TRUE),
                    cindex_left = quant_low(cindex_cox),
                    cindex_right = quant_high(cindex_cox),
                    Model = "Cox",
                    Sample_size = "1000")

# for nnet

ibs_nnet <- vector(mode = "numeric", length = 1000)
cindex_nnet <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_nnet[i] <- results_nnet[[i]]$Int_brier
  cindex_nnet[i] <- results_nnet[[i]]$Cindex
}


line5 <- data.frame(ibs = mean(ibs_nnet, na.rm = TRUE),
                    ibs_left = quant_low(ibs_nnet),
                    ibs_right = quant_high(ibs_nnet),
                    cindex = mean(cindex_nnet, na.rm = TRUE),
                    cindex_left = quant_low(cindex_nnet),
                    cindex_right = quant_high(cindex_nnet),
                    Model = "PLANN original",
                    Sample_size = "1000")


# for keras
# enter the list to calculate the measures
ibs_keras <- vector(mode = "numeric", length = 1000)
cindex_keras <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  ibs_keras[i] <- results_keras[[i]]$Int_brier
  cindex_keras[i] <- results_keras[[i]]$Cindex
}


line6 <- data.frame(ibs = mean(ibs_keras, na.rm = TRUE),
                    ibs_left = quant_low(ibs_keras),
                    ibs_right = quant_high(ibs_keras),
                    cindex = mean(cindex_keras, na.rm = TRUE),
                    cindex_left = quant_low(cindex_keras),
                    cindex_right = quant_high(cindex_keras),
                    Model = "PLANN extended",
                    Sample_size = "1000")


df <- rbind(line1, line2, line3, line4, line5, line6)

df$Model <- as.factor(df$Model)
df$Model <- factor(df$Model, levels = c("Cox", "PLANN original", "PLANN extended"))
df$Sample_size <- as.factor(df$Sample_size)
df$Sample_size <- factor(df$Sample_size, levels = c("250", "1000"))


# For the C-index
# Use 95% confidence intervals 
cindex61 <- ggplot(df, aes(x=Sample_size, y=cindex, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  labs(title = " ",x = "Sample size", y = "Concordance index", fill = "Model") +
  geom_errorbar(aes(ymin=cindex_left, ymax=cindex_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 

# For the IBS
# Use 95% confidence intervals 
ibs61 <- ggplot(df, aes(x=Sample_size, y=ibs, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  labs(title = " ",x = "Sample size", y = "Integrated Brier Score at 5 years", fill = "Model") +
  geom_errorbar(aes(ymin=ibs_left, ymax=ibs_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 

# figure 2 of the manuscript
gridExtra::grid.arrange(cindex61, ibs61, ncol = 2)


######################################################################################
######################################################################################
# Bargraphs for Cox models

# To install packages
# install_packages <- c("ggplot2", "gridExtra")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

library(ggplot2)
library(gridExtra)



quant_low <- function(x) {
  res <- mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE)
  return(res)
}

quant_high <- function(x) {
  res <- mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)
  return(res)
}




# for n = 250

# load results 
load("results_simulations1/simulations_cox.RData")
load("results_simulations1/simulations_cox_wr1.RData")
load("results_simulations1/simulations_cox_wr2.RData")


# for the Cox models

# enter the list to calculate the measures
#brier_cox <- matrix(data = 0, nrow = 1000, ncol = 6) # Brier score until 5 years
ibs_cox <- vector(mode = "numeric", length = 1000)
cindex_cox <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_cox[i] <- results_cox[[i]]$Int_brier
  cindex_cox[i] <- results_cox[[i]]$Cindex
}



line1 <- data.frame(ibs = mean(ibs_cox, na.rm = TRUE),
                    ibs_left = quant_low(ibs_cox),
                    ibs_right = quant_high(ibs_cox),
                    cindex = mean(cindex_cox, na.rm = TRUE),
                    cindex_left = quant_low(cindex_cox),
                    cindex_right = quant_high(cindex_cox),
                    Model = "Cox",
                    Sample_size = "250")


# for the Cox models with patients censored removed at 2 years (training data)

# enter the list to calculate the measures
ibs_cox_wr1 <- vector(mode = "numeric", length = 1000)
cindex_cox_wr1 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_cox_wr1[i] <- results_cox_wr1[[i]]$Int_brier
  cindex_cox_wr1[i] <- results_cox_wr1[[i]]$Cindex
}



line2 <- data.frame(ibs = mean(ibs_cox_wr1, na.rm = TRUE),
                    ibs_left = quant_low(ibs_cox_wr1),
                    ibs_right = quant_high(ibs_cox_wr1),
                    cindex = mean(cindex_cox_wr1, na.rm = TRUE),
                    cindex_left = quant_low(cindex_cox_wr1),
                    cindex_right = quant_high(cindex_cox_wr1),
                    Model = "Cox: removed before 2 years",
                    Sample_size = "250")


# for the Cox models with patients curtailed at 5 years (training data)

# enter the list to calculate the measures
ibs_cox_wr2 <- vector(mode = "numeric", length = 1000)
cindex_cox_wr2 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_cox_wr2[i] <- results_cox_wr2[[i]]$Int_brier
  cindex_cox_wr2[i] <- results_cox_wr2[[i]]$Cindex
}



line3 <- data.frame(ibs = mean(ibs_cox_wr2, na.rm = TRUE),
                    ibs_left = quant_low(ibs_cox_wr2),
                    ibs_right = quant_high(ibs_cox_wr2),
                    cindex = mean(cindex_cox_wr2, na.rm = TRUE),
                    cindex_left = quant_low(cindex_cox_wr2),
                    cindex_right = quant_high(cindex_cox_wr2),
                    Model = "Cox: curtailed at 5 years",
                    Sample_size = "250")


# For n = 1000

# load results 
load("results_simulations2/simulations_cox.RData")
load("results_simulations2/simulations_cox_wr1.RData")
load("results_simulations2/simulations_cox_wr2.RData")



# enter the list to calculate the measures
ibs_cox <- vector(mode = "numeric", length = 1000)
cindex_cox <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  ibs_cox[i] <- results_cox[[i]]$Int_brier
  cindex_cox[i] <- results_cox[[i]]$Cindex
}


line4 <- data.frame(ibs = mean(ibs_cox, na.rm = TRUE),
                    ibs_left = quant_low(ibs_cox),
                    ibs_right = quant_high(ibs_cox),
                    cindex = mean(cindex_cox, na.rm = TRUE),
                    cindex_left = quant_low(cindex_cox),
                    cindex_right = quant_high(cindex_cox),
                    Model = "Cox",
                    Sample_size = "1000")


# for the Cox models with patients censored removed at 2 years (training data)

# enter the list to calculate the measures
ibs_cox_wr1 <- vector(mode = "numeric", length = 1000)
cindex_cox_wr1 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_cox_wr1[i] <- results_cox_wr1[[i]]$Int_brier
  cindex_cox_wr1[i] <- results_cox_wr1[[i]]$Cindex
}



line5 <- data.frame(ibs = mean(ibs_cox_wr1, na.rm = TRUE),
                    ibs_left = quant_low(ibs_cox_wr1),
                    ibs_right = quant_high(ibs_cox_wr1),
                    cindex = mean(cindex_cox_wr1, na.rm = TRUE),
                    cindex_left = quant_low(cindex_cox_wr1),
                    cindex_right = quant_high(cindex_cox_wr1),
                    Model = "Cox: removed before 2 years",
                    Sample_size = "1000")


# for the Cox models with patients curtailed at 5 years (training data)

# enter the list to calculate the measures
ibs_cox_wr2 <- vector(mode = "numeric", length = 1000)
cindex_cox_wr2 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_cox_wr2[i] <- results_cox_wr2[[i]]$Int_brier
  cindex_cox_wr2[i] <- results_cox_wr2[[i]]$Cindex
}



line6 <- data.frame(ibs = mean(ibs_cox_wr2, na.rm = TRUE),
                    ibs_left = quant_low(ibs_cox_wr2),
                    ibs_right = quant_high(ibs_cox_wr2),
                    cindex = mean(cindex_cox_wr2, na.rm = TRUE),
                    cindex_left = quant_low(cindex_cox_wr2),
                    cindex_right = quant_high(cindex_cox_wr2),
                    Model = "Cox: curtailed at 5 years",
                    Sample_size = "1000")


df <- rbind(line1, line2, line3, line4, line5, line6)

df$Model <- as.factor(df$Model)
df$Model <- factor(df$Model, levels = c("Cox", "Cox: removed before 2 years", "Cox: curtailed at 5 years"))

df$Sample_size <- as.factor(df$Sample_size)
df$Sample_size <- factor(df$Sample_size, levels = c("250", "1000"))


# For the C-index
# Use 95% confidence intervals 
cox_scenario1 <- ggplot(df, aes(x=Sample_size, y=cindex, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_brewer(palette="Reds") +  #, breaks = c("250", "1000")) +
  labs(title = " ",x = "Sample size", y = "Concordance index", fill = "Model") +
  geom_errorbar(aes(ymin=cindex_left, ymax=cindex_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 


# For the IBS
# Use 95% confidence intervals 
cox_scenario2 <- ggplot(df, aes(x=Sample_size, y=ibs, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  labs(title = " ",x = "Sample size", y = "Integrated Brier Score at 5 years", fill = "Model") +
  scale_fill_brewer(palette="Reds") +  #, breaks = c("250", "1000")) +
  geom_errorbar(aes(ymin=ibs_left, ymax=ibs_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 

gridExtra::grid.arrange(cox_scenario1, cox_scenario2, ncol = 2)


######################################################################################
######################################################################################
# Bargraphs for PLANN original

# To install packages
# install_packages <- c("ggplot2", "gridExtra")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

library(ggplot2)
library(gridExtra)



quant_low <- function(x) {
  res <- mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE)
  return(res)
}

quant_high <- function(x) {
  res <- mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)
  return(res)
}


load("results_simulations1/snn_ibs/simulations_nnet.RData")
load("results_simulations1/snn_ibs/simulations_nnet_wr1.RData")
load("results_simulations1/snn_ibs/simulations_nnet_wr2.RData")

# for nnet

ibs_nnet <- vector(mode = "numeric", length = 1000)
cindex_nnet <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_nnet[i] <- results_nnet[[i]]$Int_brier
  cindex_nnet[i] <- results_nnet[[i]]$Cindex
}



line1b <- data.frame(ibs = mean(ibs_nnet, na.rm = TRUE),
                     ibs_left = quant_low(ibs_nnet),
                     ibs_right = quant_high(ibs_nnet),
                     cindex = mean(cindex_nnet, na.rm = TRUE),
                     cindex_left = quant_low(cindex_nnet),
                     cindex_right = quant_high(cindex_nnet),
                     Model = "PLANN original",
                     Sample_size = "250")

# nnet: patients censored before 2 years removed

ibs_nnet_wr1 <- vector(mode = "numeric", length = 1000)
cindex_nnet_wr1 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_nnet_wr1[i] <- results_nnet_wr1[[i]]$Int_brier
  cindex_nnet_wr1[i] <- results_nnet_wr1[[i]]$Cindex
}



line2b <- data.frame(ibs = mean(ibs_nnet_wr1, na.rm = TRUE),
                     ibs_left = quant_low(ibs_nnet_wr1),
                     ibs_right = quant_high(ibs_nnet_wr1),
                     cindex = mean(cindex_nnet_wr1, na.rm = TRUE),
                     cindex_left = quant_low(cindex_nnet_wr1),
                     cindex_right = quant_high(cindex_nnet_wr1),
                     Model = "PLANN: removed before 2 years",
                     Sample_size = "250")


# nnet: patients curtailed at 5 years

ibs_nnet_wr2 <- vector(mode = "numeric", length = 1000)
cindex_nnet_wr2 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_nnet_wr2[i] <- results_nnet_wr2[[i]]$Int_brier
  cindex_nnet_wr2[i] <- results_nnet_wr2[[i]]$Cindex
}



line3b <- data.frame(ibs = mean(ibs_nnet_wr2, na.rm = TRUE),
                     ibs_left = quant_low(ibs_nnet_wr2),
                     ibs_right = quant_high(ibs_nnet_wr2),
                     cindex = mean(cindex_nnet_wr2, na.rm = TRUE),
                     cindex_left = quant_low(cindex_nnet_wr2),
                     cindex_right = quant_high(cindex_nnet_wr2),
                     Model = "PLANN: curtailed at 5 years",
                     Sample_size = "250")




load("results_simulations2/snn_ibs/simulations_nnet.RData")
load("results_simulations2/snn_ibs/simulations_nnet_wr1.RData")
load("results_simulations2/snn_ibs/simulations_nnet_wr2.RData")


# for nnet

ibs_nnet <- vector(mode = "numeric", length = 1000)
cindex_nnet <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_nnet[i] <- results_nnet[[i]]$Int_brier
  cindex_nnet[i] <- results_nnet[[i]]$Cindex
}



line4b <- data.frame(ibs = mean(ibs_nnet, na.rm = TRUE),
                     ibs_left = quant_low(ibs_nnet),
                     ibs_right = quant_high(ibs_nnet),
                     cindex = mean(cindex_nnet, na.rm = TRUE),
                     cindex_left = quant_low(cindex_nnet),
                     cindex_right = quant_high(cindex_nnet),
                     Model = "PLANN original",
                     Sample_size = "1000")

# nnet: patients censored before 2 years removed

ibs_nnet_wr1 <- vector(mode = "numeric", length = 1000)
cindex_nnet_wr1 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_nnet_wr1[i] <- results_nnet_wr1[[i]]$Int_brier
  cindex_nnet_wr1[i] <- results_nnet_wr1[[i]]$Cindex
}



line5b <- data.frame(ibs = mean(ibs_nnet_wr1, na.rm = TRUE),
                     ibs_left = quant_low(ibs_nnet_wr1),
                     ibs_right = quant_high(ibs_nnet_wr1),
                     cindex = mean(cindex_nnet_wr1, na.rm = TRUE),
                     cindex_left = quant_low(cindex_nnet_wr1),
                     cindex_right = quant_high(cindex_nnet_wr1),
                     Model = "PLANN: removed before 2 years",
                     Sample_size = "1000")


# nnet: patients curtailed at 5 years

ibs_nnet_wr2 <- vector(mode = "numeric", length = 1000)
cindex_nnet_wr2 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_nnet_wr2[i] <- results_nnet_wr2[[i]]$Int_brier
  cindex_nnet_wr2[i] <- results_nnet_wr2[[i]]$Cindex
}



line6b <- data.frame(ibs = mean(ibs_nnet_wr2, na.rm = TRUE),
                     ibs_left = quant_low(ibs_nnet_wr2),
                     ibs_right = quant_high(ibs_nnet_wr2),
                     cindex = mean(cindex_nnet_wr2, na.rm = TRUE),
                     cindex_left = quant_low(cindex_nnet_wr2),
                     cindex_right = quant_high(cindex_nnet_wr2),
                     Model = "PLANN: curtailed at 5 years",
                     Sample_size = "1000")


df2 <- rbind(line1b, line2b, line3b, line4b, line5b, line6b)

df2$Model <- as.factor(df2$Model)
df2$Model <- factor(df2$Model, levels = c("PLANN original",
                                          "PLANN: removed before 2 years",
                                          "PLANN: curtailed at 5 years"))

df2$Sample_size <- as.factor(df2$Sample_size)
df2$Sample_size <- factor(df2$Sample_size, levels = c("250", "1000"))


# For the C-index
# Use 95% confidence intervals 
nnet_scenario1 <- ggplot(df2, aes(x=Sample_size, y=cindex, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_brewer(palette="Greens") +  #, breaks = c("250", "1000")) +
  labs(title = " ",x = "Sample size", y = "Concordance index", fill = "Model") +
  geom_errorbar(aes(ymin=cindex_left, ymax=cindex_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 


# For the IBS
# Use 95% confidence intervals 
nnet_scenario2 <- ggplot(df2, aes(x=Sample_size, y=ibs, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  labs(title = " ",x = "Sample size", y = "Integrated Brier Score at 5 years", fill = "Model") +
  scale_fill_brewer(palette="Greens") +  #, breaks = c("250", "1000")) +
  geom_errorbar(aes(ymin=ibs_left, ymax=ibs_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 

# figure 4 of the manuscript
gridExtra::grid.arrange(nnet_scenario1, nnet_scenario2, ncol = 2)



######################################################################################
######################################################################################
# Bargraphs for PLANN extended

# To install packages
# install_packages <- c("ggplot2", "gridExtra")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }

library(ggplot2)
library(gridExtra)


quant_low <- function(x) {
  res <- mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE)
  return(res)
}

quant_high <- function(x) {
  res <- mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE)
  return(res)
}


load("results_simulations1/snn_ibs/simulations_keras.RData")
load("results_simulations1/snn_ibs/simulations_keras_wr1.RData")
load("results_simulations1/snn_ibs/simulations_keras_wr2.RData")

# for nnet

ibs_keras <- vector(mode = "numeric", length = 1000)
cindex_keras <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_keras[i] <- results_keras[[i]]$Int_brier
  cindex_keras[i] <- results_keras[[i]]$Cindex
}



line1c <- data.frame(ibs = mean(ibs_keras, na.rm = TRUE),
                     ibs_left = quant_low(ibs_keras),
                     ibs_right = quant_high(ibs_keras),
                     cindex = mean(cindex_keras, na.rm = TRUE),
                     cindex_left = quant_low(cindex_keras),
                     cindex_right = quant_high(cindex_keras),
                     Model = "PLANN extended",
                     Sample_size = "250")

# nnet: patients censored before 2 years removed

ibs_keras_wr1 <- vector(mode = "numeric", length = 1000)
cindex_keras_wr1 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_keras_wr1[i] <- results_keras_wr1[[i]]$Int_brier
  cindex_keras_wr1[i] <- results_keras_wr1[[i]]$Cindex
}



line2c <- data.frame(ibs = mean(ibs_keras_wr1, na.rm = TRUE),
                     ibs_left = quant_low(ibs_keras_wr1),
                     ibs_right = quant_high(ibs_keras_wr1),
                     cindex = mean(cindex_keras_wr1, na.rm = TRUE),
                     cindex_left = quant_low(cindex_keras_wr1),
                     cindex_right = quant_high(cindex_keras_wr1),
                     Model = "PLANN: removed before 2 years",
                     Sample_size = "250")


# nnet: patients curtailed at 5 years

ibs_keras_wr2 <- vector(mode = "numeric", length = 1000)
cindex_keras_wr2 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_keras_wr2[i] <- results_keras_wr2[[i]]$Int_brier
  cindex_keras_wr2[i] <- results_keras_wr2[[i]]$Cindex
}



line3c <- data.frame(ibs = mean(ibs_keras_wr2, na.rm = TRUE),
                     ibs_left = quant_low(ibs_keras_wr2),
                     ibs_right = quant_high(ibs_keras_wr2),
                     cindex = mean(cindex_keras_wr2, na.rm = TRUE),
                     cindex_left = quant_low(cindex_keras_wr2),
                     cindex_right = quant_high(cindex_keras_wr2),
                     Model = "PLANN: curtailed at 5 years",
                     Sample_size = "250")




load("results_simulations2/snn_ibs/simulations_keras.RData")
load("results_simulations2/snn_ibs/simulations_keras_wr1.RData")
load("results_simulations2/snn_ibs/simulations_keras_wr2.RData")


# for nnet

ibs_keras <- vector(mode = "numeric", length = 1000)
cindex_keras <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_keras[i] <- results_keras[[i]]$Int_brier
  cindex_keras[i] <- results_keras[[i]]$Cindex
}



line4c <- data.frame(ibs = mean(ibs_keras, na.rm = TRUE),
                     ibs_left = quant_low(ibs_keras),
                     ibs_right = quant_high(ibs_keras),
                     cindex = mean(cindex_keras, na.rm = TRUE),
                     cindex_left = quant_low(cindex_keras),
                     cindex_right = quant_high(cindex_keras),
                     Model = "PLANN extended",
                     Sample_size = "1000")

# nnet: patients censored before 2 years removed

ibs_keras_wr1 <- vector(mode = "numeric", length = 1000)
cindex_keras_wr1 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_keras_wr1[i] <- results_keras_wr1[[i]]$Int_brier
  cindex_keras_wr1[i] <- results_keras_wr1[[i]]$Cindex
}



line5c <- data.frame(ibs = mean(ibs_keras_wr1, na.rm = TRUE),
                     ibs_left = quant_low(ibs_keras_wr1),
                     ibs_right = quant_high(ibs_keras_wr1),
                     cindex = mean(cindex_keras_wr1, na.rm = TRUE),
                     cindex_left = quant_low(cindex_keras_wr1),
                     cindex_right = quant_high(cindex_keras_wr1),
                     Model = "PLANN: removed before 2 years",
                     Sample_size = "1000")


# nnet: patients curtailed at 5 years

ibs_keras_wr2 <- vector(mode = "numeric", length = 1000)
cindex_keras_wr2 <- vector(mode = "numeric", length = 1000)


for (i in 1:1000){
  ibs_keras_wr2[i] <- results_keras_wr2[[i]]$Int_brier
  cindex_keras_wr2[i] <- results_keras_wr2[[i]]$Cindex
}



line6c <- data.frame(ibs = mean(ibs_keras_wr2, na.rm = TRUE),
                     ibs_left = quant_low(ibs_keras_wr2),
                     ibs_right = quant_high(ibs_keras_wr2),
                     cindex = mean(cindex_keras_wr2, na.rm = TRUE),
                     cindex_left = quant_low(cindex_keras_wr2),
                     cindex_right = quant_high(cindex_keras_wr2),
                     Model = "PLANN: curtailed at 5 years",
                     Sample_size = "1000")


df3 <- rbind(line1c, line2c, line3c, line4c, line5c, line6c)

df3$Model <- as.factor(df3$Model)
df3$Model <- factor(df3$Model, levels = c("PLANN extended",
                                          "PLANN: removed before 2 years",
                                          "PLANN: curtailed at 5 years"))

df3$Sample_size <- as.factor(df3$Sample_size)
df3$Sample_size <- factor(df3$Sample_size, levels = c("250", "1000"))


# For the C-index
# Use 95% confidence intervals 
keras_scenario1 <- ggplot(df3, aes(x=Sample_size, y=cindex, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_brewer(palette="Blues") +  #, breaks = c("250", "1000")) +
  labs(title = " ",x = "Sample size", y = "Concordance index", fill = "Model") +
  geom_errorbar(aes(ymin=cindex_left, ymax=cindex_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 


# For the IBS
# Use 95% confidence intervals 
keras_scenario2 <- ggplot(df3, aes(x=Sample_size, y=ibs, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  labs(title = " ",x = "Sample size", y = "Integrated Brier Score at 5 years", fill = "Model") +
  scale_fill_brewer(palette="Blues") +  #, breaks = c("250", "1000")) +
  geom_errorbar(aes(ymin=ibs_left, ymax=ibs_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 


# figure 5 of the manuscript
gridExtra::grid.arrange(keras_scenario1, keras_scenario2, ncol = 2)


## end of the file


