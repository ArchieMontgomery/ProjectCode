## function returning table of summaries for 3 repeated measurements MAR model
dropout = function(data, timev, idv, factor) {
  ntimes = count(data, deparse(substitute(idv)))$freq[1]
  
  #reshape to wide format
  data_wide = reshape(
    data = data,
    direction = "wide",
    idvar = c(deparse(substitute(idv)),deparse(substitute(factor))), 
    timevar = deparse(substitute(timev)),
    sep = ""
  )
  
  #return(ntimes)
  
  
  count(ahd, "subject")$freq[1]
  
  data_wide = data_wide[,c(3,4,5,6,1,2)]
  
  #for (i in 1:data_wide)
  data_wide$R1 = as.numeric(!is.na(data_wide[,1]))
  data_wide$R2 = as.numeric(!is.na(data_wide[,2]))
  data_wide$R3 = as.numeric(!is.na(data_wide[,3]))
  data_wide$R4 = as.numeric(!is.na(data_wide[,4]))
  
  pi1 = glm(R2~data_wide[,1], data = data_wide, family = binomial(link = logit))
  pi2 = glm(R3~data_wide[,1]+data_wide[,2], data = data_wide, family = binomial(link = logit))
  pi3 = glm(R4~data_wide[,1]+data_wide[,2]+data_wide[,3], data = data_wide, family = binomial(link = logit))
  
  results = list(summary(pi1), summary(pi2), summary(pi3))
  return(results)
  
}

# Autocorrelation function
autocorr.mat <- function(p, rho) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}


# Function to generate any data set with m subject, n individuals and r proportion
# of missingness in all measurements but the first one
sim_params = c()
missing.dataset.generator = function(n_measure, n_indiv, miss.proportion)   {
  #R = matrix(rnorm(n_measure^2), nrow = n_measure, ncol = n_measure, byrow = TRUE)
  R = genPositiveDefMat(n_measure)
  #sim_params = as.vector(R$Sigma)
  init = sample(1:100, 1)
  final = sample(1:100, 1)
  mu = seq(from = init, to = final, length.out = n_measure)
  range = final - init
  for (i in 1:n_measure) {
    mu[i] = sample(rnorm(100, mu[i], 1.5), 1) 
  }
  #a = sample(init:final, n_measure-2, replace=TRUE)
  #mu = c(init, a, final)
  #mu = mu[order(mu)]
  sim_params = c(mu, autocorr.mat(n_measure, 0.9))
  y1 = as.data.frame(MASS::mvrnorm(n_indiv, mu = mu, Sigma = autocorr.mat(n_measure, 0.9)))
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
  #return(params_dropout)
  #missing.dataset.generator(4, 50, 0.25)
  
  data_separator = function(data, timevar, idvar, resvar){
    
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
      
      # number of incomplete cases at each point
      incomp.cases = numeric(nrepeat-1)
      for (i in 1:(nrepeat-1)) {
        incomp.cases[i] = sum(is.na(newdata[,(i+3)]))
      }
      incomp.cases = c(0, incomp.cases)
      
      #Difference between number of NAs at each time point
      diff.incomp.cases = incomp.cases[-1] - incomp.cases[-length(incomp.cases)]
      diff.incomp.cases = rev(c(diff.incomp.cases, 0))
      diff.incomp.cases.rev = rev(diff.incomp.cases)
      
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
      
      # Design matrix for fixed effects, X
      X = model.matrix(model.marginal, data = data)
      
      # Fixed effects coefficients, Beta
      Beta = c(summary(model.marginal)$coefficients$fixed)
      
      # X Beta
      XBeta = X%*%Beta
      XBeta_i = unique(XBeta)
      # ZDZ^T + sigma^2 (Variance of marginal model)
      V_i = getVarCov(model.marginal,
                      type = "marginal",
                      individual = 2)
      
      # Convert to vector
      V_i_vector = c(as.matrix(unlist(V_i)))
      
      # Create vector of arbitrary phis, (exchange for function to compute phis)
      phis = cbind(matrix(0.05, nrow=1, ncol=nrepeat),0.1)
      phis = as.vector(phis)
      
      # List of starting parameters
      init_vals = list(t(XBeta_i), V_i_vector, phis)
      #return(init_vals)
      
      
      # complete cases data set
      complete.data = newdata[newdata$R == '1',]
      complete.data = complete.data[,-1]
      # Incomplete cases data set
      incomplete.data = newdata[newdata$R == '0',]
      incomplete.data = incomplete.data[,-1]
      newdata = newdata[,-1]
      
      
      
      #incomplete.data[is.na(incomplete.data)] = 0
      
      # Order data by number of NAs in each row
      incomplete.data = incomplete.data[order(rowSums(is.na(incomplete.data))), ]
      newdata = newdata[order(rowSums(is.na(newdata))), ]
      
      # Separate data into (nrepeat - 1) data sets of structure in Fitzmaurice
      separated_data = list()
      separated_data[[1]] = newdata[,1:2]
      for (i in 1:(nrepeat - 2)) {
        separated_data[[i+1]] = newdata[(1:(nrow(newdata) - sum(diff.incomp.cases.rev[1:(i)]))),1:(i+2)]
      }
      sep_data_long = list()
      for (i in 1:(nrepeat - 1)) {
        sep_data_long[[i]] = separated_data[[i]]
        sep_data_long[[i]]$ID = c(1:nrow(sep_data_long[[i]]))
        sep_data_long[[i]] = reshape(sep_data_long[[i]], 
                                     timevar = "Y", 
                                     varying = c(paste("Y",(1:(i+1)),sep="")), 
                                     idvar = "ID",
                                     direction ="long",
                                     sep = "")
        sep_data_long[[i]] = sep_data_long[[i]][order(sep_data_long[[i]]$ID),]
        sep_data_long[[i]]$time = c(1:(i+1))
      }
      return(sep_data_long)
    }
  }
  
  sep_data = data_separator(params_dropout[[2]], time, ID, Y)
  
  selection.two = function(data, timevar, idvar, resvar){
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
                   control=list(maxit=100), 
                   method = c("SANN"))
      #estim_and_true = list(head(ests$par, - (nrepeat + 1)), params_y[[1]])
      mus_and_data = list(ests, init_vals)
      return(mus_and_data)
    }
  }
  
  EE.miss.everywhere = function(n_measures) {
    start.time = Sys.time()
    # Parameter estimates for each of the separated data sets
    parameter_ests = list()
    for (i in 1:(n_measures-1)) {
      parameter_ests[[i]] = head(selection.two(sep_data[[i]], time, ID, Y)[[1]]$par, - (i+2))
    }
    
    # Length of the last separated data set
    max_length = length(tail(parameter_ests, n=1)[[1]])
    
    # First two mean parameter estimates, as we assume that the first measurement
    # is observed for all subjects
    mu_1 = mean(sapply(parameter_ests, `[[`, 1))
    mu_2 = mean(sapply(parameter_ests, `[[`, 2))
    mean_mu_params = c(mu_1, mu_2)
    mean_sigma_params = c()
    
    # Adds each subsequent mean parameter to a vector
    for (i in 1:(length(parameter_ests)-1)) {
      mean_mu_params[i+2] = mean(sapply(parameter_ests[-(1:i)], `[[`, (i+2)))
    }
    
    # Makes a new vector that contains only the sigma estimates from each data set
    sigma_parameter_ests = list()
    for (i in 1:length(parameter_ests)) {
      sigma_parameter_ests[[i]] = parameter_ests[[i]][-(1:(i+1))] 
    }
    
    # Adds NA values to each set of estimates to create lists of same length
    sigma_parameter_ests = lapply(sigma_parameter_ests, `length<-`, max(lengths(sigma_parameter_ests)))
    
    # Extracts the sigma estimates and adds them to a vector
    for (i in 1:length(tail(sigma_parameter_ests, n=1)[[1]])) {
      mean_sigma_params[i] = mean(sapply(sigma_parameter_ests, `[[`, i), na.rm = TRUE)
    }
    
    # Vector of final parameter estimates
    fin_paramter_ests = c(mean_mu_params, mean_sigma_params)
    EE = mean(abs(fin_paramter_ests - params_dropout[[1]]))
    end.time = Sys.time()
    time.taken = round(end.time - start.time, 2)
    fin_results = list(fin_paramter_ests, params_dropout[[1]], 
                       params_dropout[[2]], params_dropout[[3]], EE, 
                       time.taken)
    return(fin_results)
  }
  EE.miss.everywhere(n_measure)
}
missing.dataset.generator(3, 20, 0.1)
#################################################################

EE.estimator.data = function(n_sim, n_rm, n_subs, miss_prop) {
  start.time = Sys.time()
  diff.ours = NULL
  diff.UMI = NULL
  diff.lmer = NULL
  diff.MI = NULL
  EE.ours = numeric(n_sim)
  EE.UMI = numeric(n_sim)
  EE.lmer = numeric(n_sim)
  EE.MI = numeric(n_sim)
  for (i in 1:n_sim){
    sim = missing.dataset.generator(n_rm, n_subs, miss_prop)
    our_ests = sim[[1]]
    true_params = sim[[2]]
    sample_d_long = sim[[3]]
    sample_d_wide = sim[[4]]
    
    ############################ LMER ###############################
    
    # f(Y_i)
    model.marginal.lmer = lmer(Y ~ time + (1|ID), data = sample_d_long)
    model.marginal.lme = lme(fixed = Y ~ time, random = ~ 1|ID, 
                             data = sample_d_long, subset = complete.cases(Y),
                             na.action=na.exclude)
    
    # Design matrix
    X = model.matrix(model.marginal.lmer)
    colnames(X) = c("Intercept", "Time")
    
    # Regression coefficients
    beta = summary(model.marginal.lmer)$coefficients[,1]
    
    # X Beta
    X_beta = X%*%beta
    
    # lmer mu estimates
    lmer_mu_estimates = X_beta[1:n_rm]
    
    # lmer sigma estimates
    Sigma = getVarCov(model.marginal.lme, type = "marginal")
    lmer_sigma_estimates = as.vector(unlist(Sigma[1]))
    
    lmer_estimates = c(lmer_mu_estimates, lmer_sigma_estimates)
    
    # EE (lmer)
    diff.lmer = c(diff.lmer, abs(lmer_estimates - true_params))
    EE.lmer[i] = mean(diff.lmer[(1+(i-1)*(n_rm + (n_rm^2))):((n_rm + (n_rm^2))*i)])
    
    ############################ UMI ###############################
    
    # Mean of each Response
    UMI_data_wide = impute_mean(sample_d_wide, type = "columnwise", convert_tibble = TRUE)
    
    UMI_data_long = reshape(UMI_data_wide, 
                            timevar = "Y", 
                            varying = c(paste("Y",1:n_rm,sep="")), 
                            idvar = "ID",
                            direction ="long",
                            sep = "")
    UMI_data_long = UMI_data_long[order(UMI_data_long$ID),]
    UMI_data_long$time = c(1:n_rm)
    #UMI_data_long = UMI_data_long[,c(1,3,2)]
    
    model.marginal.UMI = lmer(Y ~ time + (1|ID), data = UMI_data_long)
    model.marginal.lme.UMI = lme(fixed = Y ~ time, random = ~ 1|ID, 
                                 data = UMI_data_long, na.action=na.exclude)
    
    # Design matrix
    X_UMI = model.matrix(model.marginal.UMI)
    colnames(X_UMI) = c("Intercept", "Time")
    
    # Regression coefficients
    beta_UMI = summary(model.marginal.UMI)$coefficients[,1]
    
    # X Beta
    X_beta_UMI = X_UMI%*%beta_UMI
    
    # UMI mu estimates
    UMI_mu_estimates = X_beta_UMI[1:n_rm]
    
    # UMI sigma estimates
    Sigma_UMI = getVarCov(model.marginal.lme.UMI, type = "marginal")
    UMI_sigma_estimates = as.vector(unlist(Sigma_UMI[1]))
    
    # UMI estimates
    UMI_estimates = c(UMI_mu_estimates, UMI_sigma_estimates)
    
    # Difference between estimates and true
    diff.UMI = c(diff.UMI, abs(UMI_estimates - true_params))
    
    # EE (UMI)
    EE.UMI[i] = mean(diff.UMI[(1+(i-1)*(n_rm + (n_rm^2))):((n_rm + (n_rm^2))*i)])
    
    ############################ MI ###############################
    
    ini = mice(sample_d_long, maxit=0) # Long Imputation
    pred1 = ini$predictorMatrix
    pred1[,"ID"] = -2 # set ID as class variable for 2l.norm
    meth1 = ini$method
    meth1[which(meth1 == "pmm")] = "pmm"
    sample_d_long = sample_d_long %>% mutate(ID = as.integer(factor(ID)))
    imp.long = mice(sample_d_long, m=2, method = meth1, predictorMatrix = pred1)
    #stripplot(imp.long, Distance, xlab = "Imputation number", col=mdc(1:2), pch=20, cex = 2.8)
    fit = with(imp.long, lme(fixed = Y ~ time, random = ~1|ID, na.action=na.exclude))
    Orthodont.imp.summary = summary(pool(fit))
    
    # Design matrix (MI)
    X_MI = model.matrix(fit$analyses[[1]])
    colnames(X_MI) = c("Intercept", "Time")
    
    # Regression coefficients (MI)
    beta_MI = Orthodont.imp.summary[,2]
    
    # X Beta (MI)
    X_beta_MI = X_MI%*%beta_MI
    
    # MI estimates
    MI_mu_estimates = X_beta_MI[1:n_rm]
    
    # MI sigma estimates
    Sigma_MI = getVarCov(fit$analyses[[1]], type = "conditional")
    MI_sigma_estimates = as.vector(unlist(Sigma_MI[1]))
    
    # MI estimates
    MI_estimates = c(MI_mu_estimates, MI_sigma_estimates)
    
    # EE (MI)
    diff.MI = c(diff.MI, abs(MI_estimates - true_params))
    EE.MI[i] = mean(diff.MI[(1+(i-1)*(n_rm + (n_rm^2))):((n_rm + (n_rm^2))*i)])
    
    ############################ OURS ###############################
    
    # EE (ours)
    diff.ours = c(diff.ours, abs(our_ests - true_params))
    EE.ours[i] = mean(diff.ours[(1+(i-1)*(n_rm + (n_rm^2))):((n_rm + (n_rm^2))*i)])
  }
  end.time = Sys.time()
  time.taken = round(end.time - start.time, 2)
  outputs = list(time.taken, EE.lmer, EE.UMI, EE.MI, EE.ours, lmer_estimates, 
                 UMI_estimates, MI_estimates, our_ests, sim)
  return(outputs)
}

X = EE.estimator.data(10, 5, 50, 0.3)

# Remove outliers

X[[2]] = (X[[2]][!(X[[2]] > 2)])
X[[3]] = (X[[3]][!(X[[3]] > 2)])
X[[4]] = (X[[4]][!(X[[4]] > 2)])
X[[5]] = (X[[5]][!(X[[5]] > 2)])


# Calculate means

round(mean(X[[2]]), 3)
round(mean(X[[3]]), 3)
round(mean(X[[4]]), 3)
round(mean(X[[5]]), 3)
c(round(mean(X[[2]]), 3), round(mean(X[[3]]), 3), 
  round(mean(X[[4]]), 3), round(mean(X[[5]]), 3))

# Set results equal to corresponding variable name
# EE.estimates.anywhere.nrm.nsub.missprop
EE.estimates.anywhere.3.50.1 = X
EE.estimates.anywhere.3.50.2 = X
EE.estimates.anywhere.3.50.3 = X
EE.estimates.anywhere.5.50.1 = X
EE.estimates.anywhere.5.50.2 = X
#################################################################

# Simulate certain number of missingness anywhere data sets

EEs = NULL
simulating_EE = function(n_sim, n_rm, n_subs, miss_perc) {
  for (k in 1:n_sim) {
    sim = missing.dataset.generator(n_rm, n_subs, miss_perc)
    EEs = c(EEs, sim)
  }
  mean_EE = mean(EEs)
  EEs_and_mean = list(EEs, mean_EE)
  return(EEs_and_mean)
}
simulating_EE(15, 3, 20, 0.2)
