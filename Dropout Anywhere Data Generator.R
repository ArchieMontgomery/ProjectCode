# Function to generate any data set with m subject, n individuals and r proportion
# of missingness in all measurements but the first one

autocorr.mat <- function(p, rho) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}

autocorr.mat(5, 0.9)
sim_params = c()
missing.dataset.generator = function(n_measure, n_indiv, miss.proportion)   {
  #R = matrix(rnorm(n_measure^2), nrow = n_measure, ncol = n_measure, byrow = TRUE)
  R = genPositiveDefMat(n_measure, rangeVar = c(1,2), covMethod = "unifcorrmat")
  #sim_params = as.vector(R$Sigma)
  init = sample(1:100, 1)
  final = sample(1:100, 1)
  mu = seq(from = init, to = final, length.out = n_measure)
  range = final - init
  for (i in 1:(n_measure)) {
    mu[i] = sample(rnorm(100, mu[i], 1), 1) 
  }
  mu = mu[order(mu)]
  sim_params = c(mu, R$Sigma)
  y1 = as.data.frame(MASS::mvrnorm(n_indiv, mu = mu, Sigma = as.matrix(autocorr.mat(n_measure, 0.9))))
  for(i in 1:n_measure){
    colnames(y1)[i] = paste("Y",i,sep="")
  }
  y.miss = prodNA(as.data.frame(y1[,2]), noNA = miss.proportion)
  colnames(y.miss) = paste("Y", 2, sep = "")
  y1 = y1[,-2]
  y1 = cbind(y1, y.miss)
  y1 = y1[,c(1,n_measure,2:(n_measure-1))]
  y1 = y1[order(rowSums(is.na(y1))), ]
  for (i in 1:n_indiv) {
    if (sum(is.na(y1[i,]) == TRUE) > 0) {
      y1[i,2:n_measure] = NA 
    } 
  }
  
  y.misses = list()
  for (i in 1:(n_measure - 2)) {
    y.misses[[i]] = prodNA(as.data.frame(y1[,(i+2)]), noNA = miss.proportion)
    colnames(y.misses[[i]]) = paste("Y", (i+2), sep = "")
    y1 = y1[,-(i+2)]
    y1 = add_column(y1, y.misses[[i]], .after = (i+1))
    #y1 = cbind(y1, y.misses[[i]])
    for (j in 1:n_indiv) {
      if (sum(is.na(y1[j,]) == TRUE) > 0) {
        y1[j,(i+2):n_measure] = NA 
      } 
    }
  }
  y1 = y1[order(rowSums(is.na(y1))), ]
  y1$ID = c(1:n_indiv)
  y1_long = reshape(y1, 
                    timevar = "Y", 
                    varying = c(paste("Y",1:n_measure,sep="")), 
                    idvar = "ID",
                    direction ="long",
                    sep = "")
  y1_long = y1_long[order(y1_long$ID),]
  y1_long$time = c(1:n_measure)
  params_dropout = list(sim_params, y1_long, y1)
  return(params_dropout)
}
A = missing.dataset.generator(10, 500, 0.05)
trues = A[[1]]
y_long = A[[2]]
model.marginal.lmer = lmer(Y ~ time + (1|ID), data = y_long)
model.marginal.lme = lme(fixed = Y ~ time, random = ~ 1|ID, 
                         data = y_long, subset = complete.cases(Y))

ggplot(y_long, aes(time, Y, col = factor(ID))) +
  geom_point() +
  geom_line() +
  theme(legend.position = "none")

y_long.complete = na.omit(y_long)
y_long.complete$ID = as.factor(y_long.complete$ID)
y_long.complete$time = as.numeric(y_long.complete$time)

oxboys.lmm.pred = predict(model.marginal.lmer)
ggplot(data = y_long.complete, aes(x = time, y = Y)) +
  geom_point(aes(col=ID)) +
  geom_line(aes(y=oxboys.lmm.pred, col = ID))+
  labs(title = "Height vs. Age", subtitle="with subject-specific intercepts")+
  theme(legend.position = "none")

# Design matrix
X = model.matrix(model.marginal.lmer)
colnames(X) = c("Intercept", "Time")

# Regression coefficients
beta = summary(model.marginal.lmer)$coefficients[,1]

# X Beta
X_beta = X%*%beta

# lmer mu estimates
lmer_mu_estimates = X_beta[1:10]

# lmer sigma estimates
Sigma = getVarCov(model.marginal.lme, type = "marginal")
lmer_sigma_estimates = as.vector(unlist(Sigma[1]))

lmer_estimates = c(lmer_mu_estimates, lmer_sigma_estimates)
mean(abs(lmer_estimates - trues))

