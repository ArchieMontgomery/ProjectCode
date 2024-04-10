# Function to generate any data set with m subject, n individuals and r proportion
# of missingness in the last measurement
sim_params = c()
dropout_time_n = function(n_measure, n_indiv, miss.proportion)   {
  #R = matrix(rnorm(n_measure^2), nrow = n_measure, ncol = n_measure, byrow = TRUE)
  R = genPositiveDefMat(n_measure)
  #sim_params = as.vector(R$Sigma)
  init = sample(1:100, 1)
  final = sample(1:100, 1)
  mu = seq(from = init, to = final, length.out = n_measure)
  range = final - init
  for (i in 1:(n_measure)) {
    mu[i] = sample(rnorm(100, mu[i], 1.75), 1) 
  }
  sim_params = c(mu, autocorr.mat(n_measure, 0.9))
  y1 = as.data.frame(MASS::mvrnorm(n_indiv, mu = mu, Sigma = autocorr.mat(n_measure, 0.9)))
  
  for(i in 1:n_measure){
    colnames(y1)[i] = paste("Y",i,sep="")
  } 
  y.miss = prodNA(as.data.frame(y1[,n_measure]), noNA = miss.proportion)
  colnames(y.miss) = paste("Y", n_measure, sep = "")
  y1 = y1[,-n_measure]
  y1 = cbind(y1, y.miss)
  colnames(y1) = paste("Y", 1:n_measure, sep = "")
  y1$ID = c(1:n_indiv)
  y_long = reshape(y1, 
                   timevar = "Y", 
                   varying = c(paste("Y",1:n_measure,sep="")), 
                   idvar = "ID",
                   direction ="long",
                   sep = "")
  y_long = y_long[order(y_long$ID),]
  y_long$time = c(1:n_measure)
  y_long = y_long[,c(1,3,2)]
  y_wide = reshape(y_long,
                   timevar = "time",
                   idvar = "ID",
                   direction = "wide",
                   sep = "")
  params_y = list(sim_params, y_long, y_wide)
  return(params_y)
}
vals_n_data = dropout_time_n(5, 100, 0.2)
true_vals = vals_n_data[[1]]
sample_d_long = vals_n_data[[2]]
ggplot(sample_d_long, aes(time, Y)) +
  geom_line(aes(col = factor(ID))) +
  geom_point() +
  theme(legend.position = "none")

model.marginal.lmer = lmer(Y ~ time + (1|ID), data = sample_d_long)
  model.marginal.lme = lme(fixed = Y ~ time, random = ~ 1|ID, 
                           data = sample_d_long, subset = complete.cases(Y))
  
  # Design matrix
  X = model.matrix(model.marginal.lmer)
  colnames(X) = c("Intercept", "Time")
  
  # Regression coefficients
  beta = summary(model.marginal.lmer)$coefficients[,1]
  
  # X Beta
  X_beta = X%*%beta
  
  # lmer estimates
  lmer_estimates = X_beta[1:3]
  
  # EE (lmer)
  diff.lmer = abs(lmer_estimates - true_vals[1:3])
  EE.lmer[i] = mean(diff.lmer[1:3])
EE.lmer

