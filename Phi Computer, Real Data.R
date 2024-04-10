rm(list = ls())
# Loading packages
install.packages("aplore3")
install.packages("optimx")
install.packages("HLMdiag")
install.packages("numDeriv")
require(aplore3)
require(optimx)
require(HLMdiag)
require(numDeriv)
require(dplyr)


phi_estimator_real = function(data) {
  n_measure = ncol(data) - 1
  data$R = +(rowMeans(is.na(data)) == 0)
  # 3. Define log-likelihood function for logistic regression model -------------
  # (see applied logistic regression)
  negll = function(par){
    #Extract guesses for alpha and beta1
    alpha = par[1]
    betas = matrix(c(par[2:(n_measure+1)]), nrow = 1, ncol = n_measure)
    #Define dependent and independent variables
    y = data$R
    Y = c()
    for (i in 1:n_measure) {
      Y[i] = na.omit(data[,(i+1)])
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
  negll.grad <- function(par){
    #Extract guesses for alpha and beta1
    alpha = par[1]
    betas = matrix(c(par[2:(n_measure+1)]), nrow = 1, ncol = n_measure)
    #Define dependent and independent variables
    y = data$R
    Y = c()
    for (i in 1:n_measure) {
      Y[i] = na.omit(data[,(i+1)])
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
  ests_and_data = list(opt, data)
  logit_first = exp(sum(ests_and_data[[1]][1,1:(n_measure+1)]*cbind(1, ests_and_data[[2]][1,2:(n_measure+1)])))
  p = logit_first / (1 + logit_first)
  final_results = list(opt, data, p)
  return(final_results)
}

B = phi_estimator_real(vals_n_data[[3]])
init_phis = as.vector(unlist(B[[1]][1,1:5]))

#y_wide = missing.dataset.generator(4, 50, 0.25)
exp(y_wide[1,1:5])
# print results of optimisation
# remove not needed information for purpose of presentation
X = summary(y_wide[[1]], order = "convcode") %>%
  select(-value, -niter, -gevals, -fevals)
print(xtable(X, type = "latex"), file = "filename2.tex")
