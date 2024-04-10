
oo <- options(repos = "https://cran.r-project.org/")
install.packages("Matrix")
install.packages("lme4")
install.packages("maxLik")
install.packages("matlib") #optimx, nlmeanb
install.packages("MASS")
install.packages("missMethods")
install.packages("broom.mixed")
install.packages("missForest")
install.packages("dplyr")
install.packages("gridExtra")
install.packages("nlme")
install.packages("clusterGeneration")
install.packages("optimx")
install.packages("mice")
library(clusterGeneration)
library(tibble)
library(broom.mixed)
library(dplyr)
library(mice)
library(optimx)
library(nlme)
library(gridExtra)
library(missForest)
library(missMethods)
library(MASS)
library(matlib)
library(maxLik)
library(lme4)
library(Matrix)
options(oo)

# Creating Orthodont data frame
Orthodont.nofactor = Orthodont[,-4]
Orthodont.nofactor_wide = reshape(Orthodont.nofactor,
                                  idvar = "Subject", 
                                  timevar = "age", 
                                  direction='wide',
                                  sep = "")
colnames(Orthodont.nofactor_wide) = c("Subject", "Dist8", "Dist10", "Dist12", "Dist14")
UMI_data_wide = impute_mean(Orthodont.nofactor_wide, type = "columnwise", convert_tibble = TRUE)
UMI_data_long = reshape(UMI_data_wide, 
                        timevar = "time", 
                        varying = c("Dist8", "Dist10", "Dist12", "Dist14"), 
                        idvar = "ID",
                        direction ="long",
                        sep = "")
UMI_data_long = UMI_data_long[order(UMI_data_long$ID),]
UMI_data_long = UMI_data_long[,c(4,2,3)]
# Option 1: Calculate mean of each column
UMI_means = c()
for (j in 1:4) {
  UMI_means[j] = mean(UMI_data_wide[,j])
}
############################################

###                 Function for n repeated measures of a normal model, 
###                             returning MLEs of parameters
NMAR_selection_n_with_real(Orthodont.nofactor, age, Subject, distance)

NMAR_selection_n_with_real = function(data, timevar, idvar, resvar){
  if(missing(data)|missing(timevar)|missing(idvar)){
    
    return('ERROR: ONE OR MORE FUNCTION ARGUMENTS MISSING')
    
  } else{
    
    ## Missingness Matrix
    
    missingmatrix = function(data){
      
      data = as.data.frame(data)
      n = dim(data)
      M = matrix(NA, nrow=n[1], ncol=n[2])
      for (i in 1:n[1]){
        for (j in 1:n[2]){
          if (is.na(data[i,j])){
            M[i,j]= 0
          }else{
            M[i,j] =1
          }
        }
      }
      return(M)
    }
    
    timevarstr = deparse(substitute(timevar))
    idvarstr = deparse(substitute(idvar))
    resvarstr = deparse(substitute(resvar))
    #factorstr = deparse(substitute(factor))
    
    data = data[,c(resvarstr, timevarstr, idvarstr)]
    names(data) = c("Y", "time", "ID")
    data = data[order(data$ID),]
    
    # number of subjects
    nsub = length(unique(data$ID))
    
    # number of repeated measurements
    nrepeat = length(unique(data$time))
    
    # produce missingness matrix
    fullmissing = missingmatrix(data)
    
    ## we only want the missingness of the final response value in each subject
    finalmissing = matrix(NA, nrow=1, ncol=nsub)
    
    for (i in 1:nsub){
      
      finalmissing[i] = fullmissing[nrepeat*i, 1]
      
    }
    
    finalmissing = as.data.frame(t(finalmissing))
    
    ## converting to a wide frame so we can get the GLM working
    datawide = reshape(data, 
                       idvar = "ID", 
                       timevar = "time", 
                       direction='wide',
                       sep = "")
    
    newdata = cbind(finalmissing, datawide)
    names(newdata)[1] = "R"
    
    # number of complete cases
    comp.cases = sum(newdata$R)
    
    # number of complete cases
    incomp.cases = nsub - comp.cases
    
    ## removing the subject column as this messes with glm
    newdata = newdata[,-2]
    
    #Y = data[which(data$time==0),1]
    #for (z in 2:nrepeat){
    #  Y = cbind(Y, data[which(data$time==z),1])
    #}
    #m = min(which(is.na(Y[,nrepeat])==TRUE))
    
    ## f(R_i | Y_i)
    model.dropout = glm(R~., family=binomial(link=logit), data=newdata)
    
    ## f(Y_i)
    model.marginal = lme(fixed = Y ~ time, random = ~1|ID, 
                         subset = complete.cases(Y),
                         data = data)
    models = list(model.marginal, model.dropout)
    
    # ZDZ^T + sigma^2 (Variance of marginal model)
    V_i = getVarCov(model.marginal,
                    type = "marginal")
    
    #margVar = V_i[1]
    #margVar = matrix(unlist(margVar), ncol = nrepeat, byrow = TRUE)
    
    #detVar = det(margVar)
    
    #sigma_11 = margVar[1,1]
    
    # complete cases data set
    complete.data = newdata[newdata$R == '1',]
    complete.data = complete.data[,-1]
    # Incomplete cases data set
    incomplete.data = newdata[newdata$R == '0',]
    incomplete.data = incomplete.data[,-1]
    #return(incomplete.data)
    
    
    
    ### Estimates ###
    
    ## Auxillary estimates ##
    ###  grabbing some initial values for the y's
    ### fist grabbing a mu hat, just by doing a complete case analysis and taking the mean of each column
    meanvec = matrix(NA, nrow=1, ncol=nrepeat)
    for (i in 1:nrepeat){
      meanvec[i] = mean(complete.data[,i])
    }
    meanvec = as.vector(meanvec)
    
    ##now need to estimate a covariance matrix
    mumat = matrix(meanvec, nrow=comp.cases, ncol=nrepeat, byrow=TRUE)
    difmu = as.matrix(complete.data - mumat)
    sum1 = t(difmu)%*%difmu
    sampcov = as.matrix((1/(comp.cases-1))*sum1)
    
    phis = cbind(matrix(0.01, nrow=1, ncol=nrepeat),0.05)
    phis = as.vector(phis)
    
    #newdata_wide = newdata[,-1]
    #newdata_wide$ID = c(1:nsub)
    #init_phis = as.vector(unlist(phi_estimator_real(newdata_wide)[[1]][1,1:(nrepeat+1)]))
    
    init_vals = list(meanvec, sampcov, phis)
    #return(init_vals)
    
    
    
    ### Likelihoods ###
    
    # optim (nonignorable) (Method 1)
    
    logLik_nonignorable_optim_min_n = function(theta) {
      mu = matrix(c(theta[1:nrepeat]),nrow=1, ncol=nrepeat)
      Sigma = matrix(c(theta[(nrepeat+1):(nrepeat+(nrepeat*nrepeat))]), 
                     nrow=nrepeat, ncol=nrepeat, byrow=TRUE)
      phis = matrix(c(theta[(nrepeat+(nrepeat*nrepeat)+1):length(theta)]), 
                    nrow=1, ncol=(nrepeat+1))
      detsig = det(Sigma)
      invsig = as.matrix(solve(Sigma))
      
      ll1 = 0
      ## Part 1 ##
      for (i in 1:comp.cases) {
        Ys = as.matrix(cbind(1, complete.data[i,]))
        logits = phis%*%t(Ys)        
        p = exp(logits)/(1 + exp(logits))
        diff = as.matrix(complete.data[i,] - mu)
        mat1 = diff%*%invsig
        mat1 = mat1%*%t(diff)
        mat1 = as.numeric(mat1)
        ll1 = ll1 + (-1/2)*log(detsig) - (1/2)*mat1 + log(p)
      }
      
      ll2 = 0
      ## Part 2 ##
      for (i in 1:incomp.cases) {
        gfunc = function(ql){
          ######### this is the issue here need to include an x
          Ys = cbind(1,incomplete.data[i,1:(nrepeat-1)],ql)
          logits = phis%*%t(Ys)
          p = exp(logits)/(1+exp(logits))
          diff = as.matrix(cbind(incomplete.data[i,1:(nrepeat-1)], ql) - mu)
          mat1 = diff%*%invsig
          mat1 = mat1%*%t(diff)
          mat1 = as.numeric(mat1)
          return(exp(-(1/2)*mat1)*(1-p))
        }
        hfunc = function(x){
          return((1/sqrt(2*pi))*exp(((x-mu[nrepeat])^2)*(-1/2))) # have changed the distribution to a mean mu[n_obs] sd 1 dist to see if it changes anything
        }
        ## here 100 is number of importance samples we do
        impsample = matrix(NA, nrow=1, ncol=100)
        for (q in 1:100){
          r = rnorm(1, mean=mu[nrepeat], sd=1)
          impsample[q] = (gfunc(r)/hfunc(r))
        }
        impsample = as.vector(impsample)
        intapprox = mean(impsample)
        covminn = Sigma[-nrepeat,-nrepeat]
        detcovminn = abs(det(as.matrix(covminn)))
        d = (1/2)*log(detcovminn)
        ll2 = ll2 - d + log(intapprox)
      }
      
      ll = ll1 + ll2
      llf = -ll
    }
    
    param_list = list()
    wrapper_function = function(x){
      #cat(x)
      ## paramatrix = rbind(paramatrix, x)
      param_list[[length(param_list)+1]] <<-x
      return(logLik_nonignorable_optim_min_n(x))
    }
    
    ests = optim(unlist(init_vals), fn = wrapper_function, 
                 method = "SANN", control=list(maxit=50))
    #estim_and_true = list(head(ests$par, - (nrepeat + 1)), params_y[[1]])
    mus_and_data = list(head(ests$par, - (nrepeat + 1)), head(unlist(init_vals), -(nrepeat + 1)), newdata)
    return(mus_and_data)
  }
}

sim_ours = NMAR_selection_n_with_real(vals_n_data[[2]], time, ID, Y)
mean(abs(sim_ours[[1]] - vals_n_data[[1]]))
#mean(sim_ours[[1]] - sim_ours[[2]])
trues = vals_n_data[[1]]
y_long = vals_n_data[[2]]
model.marginal.lmer = lmer(Y ~ time + (1|ID), data = y_long)
model.marginal.lme = lme(fixed = Y ~ time, random = ~ 1|ID, 
                         data = y_long, subset = complete.cases(Y))

ggplot(y_long, aes(time, Y, col = factor(ID))) +
  geom_point() +
  geom_line() +
  theme(legend.position = "none")

str(y_long.complete)
y_long.complete = na.omit(y_long)
y_long.complete$ID = as.factor(y_long.complete$ID)
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
lmer_mu_estimates = X_beta[1:4]

# lmer sigma estimates
Sigma = getVarCov(model.marginal.lme, type = "marginal")
lmer_sigma_estimates = as.vector(unlist(Sigma[1]))

lmer_estimates = c(lmer_mu_estimates, lmer_sigma_estimates)
mean(abs(lmer_estimates - trues))

##############################################################################

# Function simulating data n times and obtaining estimation Errors for each 
# simulation, and comparing with UMI imputation for last measurement

EE.estimator.real = function(n_sim, miss.prop) {
  start.time = Sys.time()
  diff.ours = NULL
  diff.UMI = NULL
  diff.lmer = NULL
  diff.MI = NULL
  EE.ours = numeric(n_sim)
  EE.UMI = numeric(n_sim)
  EE.lmer = numeric(n_sim)
  EE.MI = numeric(n_sim)
  Dist8_mu = mean(Orthodont.nofactor_wide$Dist8)
  Dist10_mu = mean(Orthodont.nofactor_wide$Dist10)
  Dist12_mu = mean(Orthodont.nofactor_wide$Dist12)
  Dist14_mu = mean(Orthodont.nofactor_wide$Dist14)
  true_means = c(Dist8_mu, Dist10_mu, Dist12_mu, Dist14_mu)
  for (i in 1:n_sim){
    # Data set up
    sample_d = Orthodont.nofactor_wide[sample(1:nrow(Orthodont.nofactor_wide), 
                                              nrow(Orthodont.nofactor_wide), 
                                              replace = TRUE), ]
    sample_d = distinct(sample_d, Subject, .keep_all = TRUE)
    Orth.mis = prodNA(as.data.frame(sample_d[,5]), noNA = miss.prop)
    sample_d = sample_d[,1:4]
    sample_d = cbind(sample_d, Orth.mis)
    colnames(sample_d)[5] = c("Dist14")
    sample_d.miss = reshape(data = sample_d, 
                            varying = c("Dist8", "Dist10", "Dist12", "Dist14"),
                            timevar = "Age",
                            idvar = "Subject",
                            direction = "long",
                            sep = "")
    sample_d.miss = sample_d.miss[order(sample_d.miss$Subject),]
    colnames(sample_d.miss) = c("ID", "time", "Y")
    
    ############################ LMER ###############################
    
    # f(Y_i)
    model.marginal.lmer = lmer(Y ~ time + (1|ID), data = sample_d.miss)
    model.marginal.lme = lme(fixed = Y ~ time, random = ~ 1|ID, 
                             data = sample_d.miss, subset = complete.cases(Y))
    
    # Design matrix
    X = model.matrix(model.marginal.lmer)
    colnames(X) = c("Intercept", "Time")
    
    # Regression coefficients
    beta = summary(model.marginal.lmer)$coefficients[,1]
    
    # X Beta
    X_beta = X%*%beta
    
    # lmer estimates
    lmer_estimates = X_beta[1:4]
    
    # EE (lmer)
    diff.lmer = c(diff.lmer, abs(lmer_estimates - true_means))
    EE.lmer[i] = mean(diff.lmer[(1+(i-1)*4):(4*i)])
    
    ############################ UMI ###############################
    
    # Mean of each Response
    UMI_1 = mean(sample_d[,2], na.rm = TRUE)
    UMI_2 = mean(sample_d[,3], na.rm = TRUE)
    UMI_3 = mean(sample_d[,4], na.rm = TRUE)
    UMI_4 = mean(sample_d[,5], na.rm = TRUE)
    UMIs = c(UMI_1, UMI_2, UMI_3, UMI_4)
    
    # Mean of each Response
    UMI_data_wide = impute_mean(sample_d, type = "columnwise", convert_tibble = TRUE)
    
    # Option 1: Calculate mean of each column
    #UMI_means = c()
    #for (j in 1:4) {
    #  UMI_means[j] = mean(UMI_data_wide[,(j+1)])
    #}
    
    #Option 2: Use lmer
    UMI_data_long = reshape(UMI_data_wide, 
                            timevar = "time", 
                            varying = c("Dist8", "Dist10", "Dist12", "Dist14"), 
                            idvar = "ID",
                            direction ="long",
                            sep = "")
    UMI_data_long = UMI_data_long[order(UMI_data_long$ID),]
    UMI_data_long = UMI_data_long[,c(4,2,3)]
    colnames(UMI_data_long) = c("ID", "time", "Y")
    
    model.marginal.UMI = lmer(Y ~ time + (1|ID), data = UMI_data_long)
    model.marginal.lme.UMI = lme(fixed = Y ~ time, random = ~ 1|ID, 
                                 data = UMI_data_long)
    
    # Design matrix
    X_UMI = model.matrix(model.marginal.UMI)
    colnames(X_UMI) = c("Intercept", "Time")
    
    # Regression coefficients
    beta_UMI = summary(model.marginal.UMI)$coefficients[,1]
    
    # X Beta
    X_beta_UMI = X_UMI%*%beta_UMI
    
    # UMI mu estimates
    UMI_mu_estimates = X_beta_UMI[1:4]
    
    # Difference between estimates and true
    diff.UMI = c(diff.UMI, abs(UMI_mu_estimates - true_means))
    
    # EE (UMI)
    EE.UMI[i] = mean(diff.UMI[(1+(i-1)*4):(4*i)])
    
    ############################ MI ###############################
    
    ini = mice(sample_d.miss, maxit=0) # Long Imputation
    pred1 = ini$predictorMatrix
    pred1[,"ID"] = -2 # set ID as class variable for 2l.norm
    meth1 = ini$method
    meth1[which(meth1 == "pmm")] = "pmm"
    sample_d.miss = sample_d.miss %>% mutate(ID = as.integer(factor(ID)))
    imp.long = mice(sample_d.miss, m=2, method = meth1, predictorMatrix = pred1)
    #stripplot(imp.long, Distance, xlab = "Imputation number", col=mdc(1:2), pch=20, cex = 2.8)
    fit = with(imp.long, lmer(Y ~ time + (1|ID)))
    Orthodont.imp.summary = summary(pool(fit))
    
    # Design matrix (MI)
    X_MI = model.matrix(fit$analyses[[1]])
    colnames(X_MI) = c("Intercept", "Time")
    
    # Regression coefficients (MI)
    beta_MI = Orthodont.imp.summary[,2]
    
    # X Beta (MI)
    X_beta_MI = X_MI%*%beta_MI
    
    # MI estimates
    MI_estimates = X_beta_MI[1:4]
    
    # EE (MI)
    diff.MI = c(diff.MI, abs(MI_estimates - true_means))
    EE.MI[i] = mean(diff.MI[(1+(i-1)*4):(4*i)])
    
    ############################ OURS ###############################
    
    # Our function
    sim_ours = NMAR_selection_n_with_real(sample_d.miss, time, ID, Y)
    
    # EE (ours)
    diff.ours = c(diff.ours, abs(sim_ours[[1]][1:4] - true_means))
    EE.ours[i] = mean(diff.ours[(1+(i-1)*4):(4*i)])
  }
  end.time = Sys.time()
  time.taken = round(end.time - start.time, 2)
  outputs = list(time.taken, EE.ours, EE.lmer, EE.UMI, EE.MI)
  return(outputs)
}
EE.estimates = EE.estimator.real(50, 0.4)
EE.estimates[[2]] = (EE.estimates[[2]][!(EE.estimates[[2]] > 2)])
mean(EE.estimates[[2]])
mean(EE.estimates[[3]])
mean(EE.estimates[[4]])
mean(EE.estimates[[5]])

EE.estimates.10 = EE.estimates
EE.estimates.20 = EE.estimates
EE.estimates.30 = EE.estimates
EE.estimates.40 = EE.estimates
EE.estimates.full = c(EE.estimates.10, EE.estimates.20, EE.estimates.30, EE.estimates.40)

##############################################################################

# Create data frame of estimates for each method
our_estimate = as.data.frame(cbind(EE.estimates.10[[2]], "Ours"))
LMER_estimate = as.data.frame(cbind(EE.estimates.10[[3]], "CCA"))
UMI_estimate = as.data.frame(cbind(EE.estimates.10[[4]], "UMI"))
MI_estimate = as.data.frame(cbind(EE.estimates.10[[5]], "MI"))

colnames(our_estimate) = c("EE", "Method")
colnames(LMER_estimate) = c("EE", "Method")
colnames(UMI_estimate) = c("EE", "Method")
colnames(MI_estimate) = c("EE", "Method")

EE_dataframe = rbind(our_estimate, LMER_estimate, UMI_estimate, MI_estimate)
EE_dataframe = EE_dataframe %>% mutate(Method=factor(Method, 
                                                     levels=c("Ours", "CCA", "UMI", "MI")))

# Violin Plot
p_10 = ggplot(EE_dataframe, aes(x = Method, y = as.numeric(EE), fill = Method))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x = "Method", y = "EE (10%)") +
  theme_bw() +
  theme(axis.text.x=element_text(size=15, colour = "black"),
        axis.text.y=element_text(size=15, colour = "black"),
        axis.title=element_text(size=18)) +
  theme(legend.key.size = unit(1.5,"line"))

##############################################################################

# Create data frame of estimates for each method
our_estimate = as.data.frame(cbind(EE.estimates.20[[2]], "Ours"))
LMER_estimate = as.data.frame(cbind(EE.estimates.20[[3]], "CCA"))
UMI_estimate = as.data.frame(cbind(EE.estimates.20[[4]], "UMI"))
MI_estimate = as.data.frame(cbind(EE.estimates.20[[5]], "MI"))

colnames(our_estimate) = c("EE", "Method")
colnames(LMER_estimate) = c("EE", "Method")
colnames(UMI_estimate) = c("EE", "Method")
colnames(MI_estimate) = c("EE", "Method")

EE_dataframe = rbind(our_estimate, LMER_estimate, UMI_estimate, MI_estimate)
EE_dataframe = EE_dataframe %>% mutate(Method=factor(Method, 
                                                     levels=c("Ours", "CCA", "UMI", "MI")))

# Violin Plot
p_20 = ggplot(EE_dataframe, aes(x = Method, y = as.numeric(EE), fill = Method))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x = "Method", y = "EE (20%)") +
  theme_bw() +
  theme(axis.text.x=element_text(size=15, colour = "black"),
        axis.text.y=element_text(size=15, colour = "black"),
        axis.title=element_text(size=18)) +
  theme(legend.key.size = unit(1.5,"line"))

##############################################################################

# Create data frame of estimates for each method
our_estimate = as.data.frame(cbind(EE.estimates.30[[2]], "Ours"))
LMER_estimate = as.data.frame(cbind(EE.estimates.30[[3]], "CCA"))
UMI_estimate = as.data.frame(cbind(EE.estimates.30[[4]], "UMI"))
MI_estimate = as.data.frame(cbind(EE.estimates.30[[5]], "MI"))

colnames(our_estimate) = c("EE", "Method")
colnames(LMER_estimate) = c("EE", "Method")
colnames(UMI_estimate) = c("EE", "Method")
colnames(MI_estimate) = c("EE", "Method")

EE_dataframe = rbind(our_estimate, LMER_estimate, UMI_estimate, MI_estimate)
EE_dataframe = EE_dataframe %>% mutate(Method=factor(Method, 
                                                     levels=c("Ours", "CCA", "UMI", "MI")))

# Violin Plot
p_30 = ggplot(EE_dataframe, aes(x = Method, y = as.numeric(EE), fill = Method))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x = "Method", y = "EE (30%)") +
  theme_bw() +
  theme(axis.text.x=element_text(size=15, colour = "black"),
        axis.text.y=element_text(size=15, colour = "black"),
        axis.title=element_text(size=18)) +
  theme(legend.key.size = unit(1.5,"line"))

#############################################################################

# Create data frame of estimates for each method
our_estimate = as.data.frame(cbind(EE.estimates.40[[2]], "Ours"))
LMER_estimate = as.data.frame(cbind(EE.estimates.40[[3]], "CCA"))
UMI_estimate = as.data.frame(cbind(EE.estimates.40[[4]], "UMI"))
MI_estimate = as.data.frame(cbind(EE.estimates.40[[5]], "MI"))

colnames(our_estimate) = c("EE", "Method")
colnames(LMER_estimate) = c("EE", "Method")
colnames(UMI_estimate) = c("EE", "Method")
colnames(MI_estimate) = c("EE", "Method")

EE_dataframe = rbind(our_estimate, LMER_estimate, UMI_estimate, MI_estimate)
EE_dataframe = EE_dataframe %>% mutate(Method=factor(Method, 
                                                     levels=c("Ours", "CCA", "UMI", "MI")))

# Violin Plot
p_40 = ggplot(EE_dataframe, aes(x = Method, y = as.numeric(EE), fill = Method))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x = "Method", y = "EE (40%)") +
  theme_bw() +
  theme(axis.text.x=element_text(size=15, colour = "black"),
        axis.text.y=element_text(size=15, colour = "black"),
        axis.title=element_text(size=18)) +
  theme(legend.key.size = unit(1.5,"line"))

gridExtra::grid.arrange(p_10, p_20, p_30, p_40)

#############################################################################

## Convergence plots

paramatrix = do.call(rbind, param_list)
print(paramatrix)
for (i in 1:nrepeat){
  plot(paramatrix[,i], type='l', main=paste('Convergence of Mean', i, sep=" "))
}
for (i in ((nrepeat+1)+(nrepeat^2)):((2*nrepeat)+1+(nrepeat^2))){
  plot(paramatrix[,i], type='l', main=paste('Convergence of Phi', (i-((nrepeat+1)+(nrepeat^2))), sep=" "))
}

#estimates = NMAR_selection_n_with_real(y_long[[2]], time, ID, Y)
#PE = 1/20*(sum(head(estimates, -5) - y_long[[1]]))
#results1 = NMAR_selection_n_with_real(y_long, time, ID, Y)

# List of Parameters

parlist =  list()
for (i in 1:5){
  selecttt = NMAR_selection_n_with_real(y_long, time, ID, Y)
  parlist = c(parlist,list(selecttt$par[1:6]))
}
parlist[[1]]
meannss = matrix(NA, nrow=1, ncol=6)
for (i in 1:6){
  meannss[i] = mean(parlist[[1]][i], parlist[[2]][i], 
                    parlist[[3]][i], parlist[[4]][i], 
                    parlist[[5]][i])
}
longdata = function(Y){
  M = matrix(NA, nrow=((dim(Y)[1])*(dim(Y)[2])), ncol=dim(Y)[2])
  for (i in 1:dim(Y)[1]){
    M[(((dim(Y)[2])*(i-1)+1):((dim(Y)[2])*i)),3]=i
    for (j in 1:(dim(Y)[2])){
      M[(((dim(Y)[2])*(i-1)+j)),2]=j
      M[(((dim(Y)[2])*(i-1)+j)),1]=Y[i,j]
    }
  }
  M = as.data.frame(M)
  names(M) = c('Y', 'Time', 'ID')
  return(M)
}
data5 = longdata(Y)[,1:3]
cbind(t(Y[1,(1:2)]),1)



sampcov = matrix(NA, nrow=4, ncol=4)
for (w in 1:4){
  for (a in 1:4){
    veccc = matrix(NA, nrow=1, ncol=nsub)
    for (o in 1:nsub){
      veccc[o] = (Y[o,w]-meanvec[w])*(Y[o,a]-meanvec[a])
    }
    sampcov[w,a] = (1/(nsub-1))*sum(veccc)
    print(sampcov[w,a])
  }
}

meanvec = c(2,3,4,5,6)
cov = matrix(c(0.5,0.25,0.25,0.25,0.25,0.25,0.50,0.25,0.25,0.25,0.25,0.25,0.50,0.25,0.25,0.25,0.25,0.25,0.50,0.25,0.25,0.25,0.25,0.25,0.50),nrow=5, byrow=TRUE)
Y = rmvnorm(n=100, mean=meanvec, sigma=cov)
Y[60:100,5]=NA

NMAR_selection_n_with_real(longdata(Y)[1,3], 'Time', 'ID', 'Y')

plot_optim_convergence = function(opt_result) {
  if (!inherits(opt_result, "list")) {
    stop("Input must be an 'optim' object.")
  }
  
  par_names = paste("Parameter", 1:length(opt_result$par))
  par_colors = rainbow(length(opt_result$par))
  
  plot(1:length(opt_result$value), opt_result$value, type = "l",
       xlab = "Iteration", ylab = "Objective Value",
       main = "Convergence of Objective Function")
  legend("topright", legend = par_names, col = par_colors, lty = 1)
  
  for (i in 1:length(opt_result$par)) {
    lines(1:length(opt_result$value), opt_result$par_history[,i], col = par_colors[i])
  }
}
plot_optim_convergence(results1)
str(results1)
