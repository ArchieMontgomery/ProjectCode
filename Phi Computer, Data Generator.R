rm(list = ls())
# Loading packages
install.packages("aplore3")
install.packages("optimx")
install.packages("numDeriv")
require(aplore3)
require(optimx)
require(numDeriv)
require(dplyr)

# Loading dataset
missing.dataset.generator = function(n_measure, n_indiv, miss.proportion)   {
  #R = matrix(rnorm(n_measure^2), nrow = n_measure, ncol = n_measure, byrow = TRUE)
  R = genPositiveDefMat(n_measure)
  #sim_params = as.vector(R$Sigma)
  for (i in 1:n_measure){
    mu = sample(1:100, n_measure)
  }
  sim_params = c(mu, R$Sigma)
  y1 = as.data.frame(MASS::mvrnorm(n_indiv, mu = mu, Sigma = R$Sigma))
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
  
  y1$R = +(rowMeans(is.na(y1)) == 0)
  #return(y1)
  # Define log-likelihood function for logistic regression model
  negll = function(par){
    #Extract guesses for alpha and beta1
    alpha = par[1]
    betas = matrix(c(par[2:(n_measure+1)]), nrow = 1, ncol = n_measure)
    #Define dependent and independent variables
    y = y1$R
    Y = c()
    for (i in 1:n_measure) {
      Y[i] = na.omit(y1[,i])
    }
    Y = as.matrix(Y)
    #Calculate pi and xb
    xb = alpha + betas%*%Y
    pi = exp(xb) / (1 + exp(xb))
    #Set high values for 0 < pi < 1
    if(any(na.omit(pi)) > 1 || any(na.omit(pi)) < 0) {
      val = 1e+200
    } else {
      val = -sum(y * log(pi) + (1 - y) * log(1 - pi))
    }
    val
  }
  # 4. Define fradient function for logistic regression model
  # (see applied logistic regression)
  negll.grad = function(par){
    #Extract guesses for alpha and beta1
    alpha = par[1]
    betas = matrix(c(par[2:(n_measure+1)]), nrow = 1, ncol = n_measure)
    #Define dependent and independent variables
    y = y1$R
    Y = c()
    for (i in 1:n_measure) {
      Y[i] = na.omit(y1[,i])
    }
    Y = as.matrix(Y)
    #Create output vector
    n = length(par[1])
    gg = as.vector(rep(0, n))
    #Calculate pi and xb
    xb = alpha + betas%*%Y
    pi = exp(xb) / (1 + exp(xb))
    #Calculate gradients for alpha and beta1
    gg = c()
    gg[1] = - sum(y - pi)
    for (i in 1:n_measure) {
      gg[i+1] = - sum(Y[i] * (y - pi))
    }
    return(gg)
  }
  
  mygrad = negll.grad(c(rep(0, (n_measure + 1))))
  numgrad = grad(x = c(rep(0, (n_measure + 1))), func = negll)
  all.equal(mygrad, numgrad)
  
  opt = optimx(par = c(rep(0, (n_measure + 1))), negll,
               gr = negll.grad,
               control = list(trace = 0, all.methods = TRUE))
  ests_and_data = list(opt, y1)
  logit_first = exp(sum(ests_and_data[[1]][5,1:(n_measure+1)]*cbind(1, ests_and_data[[2]][1,1:n_measure])))
  p = logit_first / (1 + logit_first)
  final_results = list(opt, y1, p)
  return(final_results)
}

y_wide = missing.dataset.generator(4, 50, 0.25)
phis = as.data.frame(y_wide[1])
phis = as.vector(unlist(phis[1,1:5]))

# print results of optimisation
# remove not needed information for purpose of presentation
X = summary(y_wide[[1]], order = "convcode") %>%
  select(-value, -niter, -gevals, -fevals)
print(xtable(X, type = "latex"), file = "filename2.tex")
