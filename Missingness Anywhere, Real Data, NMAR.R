# Function to generate any data set with m subject, n individuals and r proportion
# of missingness in all measurements but the first one
sim_params = c()
Orthodont.nofactor_wide
missing.dataset.generator.real = function(data, miss.proportion) {
  n_measure = ncol(data) - 1
  n_indiv = nrow(data)
  sim_params = c()
  for (i in 1:n_measure) {
    sim_params = c(sim_params, mean(data[,(i+1)]))
  }
  names(data)[1] = "ID"
  for (i in 1:n_measure) {
    names(data)[i+1] = paste("Y",i,sep = "")
  }
  y1 = data
  y1 = y1[-1]
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
  #missing.dataset.generator.real(Orthodont.nofactor_wide, 0.15)
  
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
        #llf = -ll
      }
      
      param_list = list()
      wrapper_function = function(x){
        #cat(x)
        ## paramatrix = rbind(paramatrix, x)
        param_list[[length(param_list)+1]] <<-x
        return(logLik_nonignorable_optim_min_n(x))
      }
      
      ests = optim(unlist(init_vals), fn = wrapper_function, control=list(maxit=20))
      #estim_and_true = list(head(ests$par, - (nrepeat + 1)), params_y[[1]])
      mus_and_data = list(ests, init_vals)
      return(mus_and_data)
    }
  }
  
  PE.Miss.everywhere = function(n_measures) {
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
    EE = mean(fin_paramter_ests[1:n_measures] - params_dropout[[1]])
    fin_results = list(fin_paramter_ests, params_dropout[[1]], 
                       params_dropout[[2]], params_dropout[[3]], EE)
    return(fin_results)
  }
  PE.Miss.everywhere(n_measure)
}

X = missing.dataset.generator.real(Orthodont.nofactor_wide, 0.1)


#################################################################

PE.estimator.miss.anywhere.real = function(n_sim, data, miss_prop) {
  start.time = Sys.time()
  diff.ours = NULL
  diff.UMI = NULL
  diff.lmer = NULL
  diff.MI = NULL
  PE.ours = numeric(n_sim)
  PE.UMI = numeric(n_sim)
  PE.lmer = numeric(n_sim)
  PE.MI = numeric(n_sim)
  
  # Number of repeated measures
  n_measure = ncol(data) - 1
  
  # Number of individuals
  n_indiv = nrow(data)
  
  # Means of Data
  true_means = c()
  for (i in 1:n_measure) {
    true_means = c(true_means, mean(data[,(i+1)]))
  }
  
  for (i in 1:n_sim){
    sim = missing.dataset.generator.real(data, miss_prop)
    our_ests = sim[[1]]
    true_params = sim[[2]]
    sample_d_long = sim[[3]]
    sample_d_wide = sim[[4]]
    
    ############################ LMER ###############################
    
    # f(Y_i)
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
    
    # lmer mu estimates
    lmer_mu_estimates = X_beta[1:n_measure]
    
    # EE (lmer)
    diff.lmer = c(diff.lmer, lmer_mu_estimates - true_means)
    PE.lmer[i] = mean(diff.lmer[(1+(i-1)*n_measure):(n_measure*i)])
    
    ############################ UMI ###############################
    
    # Mean of each Response
    UMI_data_wide = impute_mean(sample_d_wide, type = "columnwise", convert_tibble = TRUE)
    
    # Option 1: Calculate mean of each column
    UMI_means = c()
    for (j in 1:n_measure) {
      UMI_means[j] = mean(UMI_data_wide[,j])
    }
    
    # Option 2: Use lmer
    #UMI_data_long = reshape(UMI_data_wide, 
    #                        timevar = "Y", 
    #                        varying = c(paste("Y",1:n_measure,sep="")), 
    #                        idvar = "ID",
    #                        direction ="long",
    #                        sep = "")
    #UMI_data_long = UMI_data_long[order(UMI_data_long$ID),]
    #UMI_data_long$time = c(1:n_measure)
    #UMI_data_long = UMI_data_long[,c(1,3,2)]
    
    #model.marginal.UMI = lmer(Y ~ time + (1|ID), data = UMI_data_long)
    #model.marginal.lme.UMI = lme(fixed = Y ~ time, random = ~ 1|ID, 
    #                             data = UMI_data_long)
    
    # Design matrix
    #X_UMI = model.matrix(model.marginal.UMI)
    #colnames(X_UMI) = c("Intercept", "Time")
    
    # Regression coefficients
    #beta_UMI = summary(model.marginal.UMI)$coefficients[,1]
    
    # X Beta
    #X_beta_UMI = X_UMI%*%beta_UMI
    
    # UMI mu estimates
    #UMI_mu_estimates = X_beta_UMI[1:n_measure]
    
    # Difference between estimates and true
    diff.UMI = c(diff.UMI, (UMI_means - true_means))
    #diff.UMI = c(diff.UMI, (UMI_mu_estimates - true_means))
    
    # EE (UMI)
    PE.UMI[i] = mean(diff.UMI[(1+(i-1)*n_measure):(n_measure*i)], na.rm = TRUE)
    
    ############################ MI ###############################
    
    ini = mice(sample_d_long, maxit=0) # Long Imputation
    pred1 = ini$predictorMatrix
    pred1[,"ID"] = -2 # set ID as class variable for 2l.norm
    meth1 = ini$method
    meth1[which(meth1 == "pmm")] = "pmm"
    sample_d_long = sample_d_long %>% mutate(ID = as.integer(factor(ID)))
    imp.long = mice(sample_d_long, m=2, method = meth1, predictorMatrix = pred1)
    #stripplot(imp.long, Distance, xlab = "Imputation number", col=mdc(1:2), pch=20, cex = 2.8)
    fit = with(imp.long, lme(fixed = Y ~ time, random = ~1|ID))
    Orthodont.imp.summary = summary(pool(fit))
    
    # Design matrix (MI)
    X_MI = model.matrix(fit$analyses[[1]])
    colnames(X_MI) = c("Intercept", "Time")
    
    # Regression coefficients (MI)
    beta_MI = Orthodont.imp.summary[,2]
    
    # X Beta (MI)
    X_beta_MI = X_MI%*%beta_MI
    
    # MI estimates
    MI_mu_estimates = X_beta_MI[1:n_measure]
    
    # EE (MI)
    diff.MI = c(diff.MI, (MI_mu_estimates - true_means))
    PE.MI[i] = mean(diff.MI[(1+(i-1)*n_measure):(n_measure*i)])
    
    ############################ OURS ###############################
    
    # EE (ours)
    diff.ours = c(diff.ours, our_ests[1:n_measure] - true_means)
    PE.ours[i] = mean(diff.ours[(1+(i-1)*n_measure):(n_measure*i)])
  }
  end.time = Sys.time()
  time.taken = round(end.time - start.time, 2)
  outputs = list(time.taken, PE.lmer, PE.UMI, PE.MI, PE.ours)
  return(outputs)
}

PEs = PE.estimator.miss.anywhere.real(2, Orthodont.nofactor_wide, 0.3)
mean(PEs[5:9])
lmer_PEs = PEs[[2]]
UMI_PEs = PEs[[3]]
MI_PEs = PEs[[4]]
our_PEs = PEs[[5]]
##################################################################

# Create data frame of estimates for each method
LMER_estimate = as.data.frame(cbind(lmer_PEs, "LMER"))
UMI_estimate = as.data.frame(cbind(UMI_PEs, "UMI"))
MI_estimate = as.data.frame(cbind(MI_PEs, "MI"))
our_estimate = as.data.frame(cbind(our_PEs, "Ours"))

colnames(our_estimate) = c("PE", "Method")
colnames(LMER_estimate) = c("PE", "Method")
colnames(UMI_estimate) = c("PE", "Method")
colnames(MI_estimate) = c("PE", "Method")

EE_dataframe = rbind(LMER_estimate, UMI_estimate, MI_estimate, our_estimate)
EE_dataframe = EE_dataframe %>% mutate(Method=factor(Method, 
                                                     levels=c("LMER", "UMI", "MI", "Ours")))

# Violin Plot
p_30 = ggplot(EE_dataframe, aes(x = Method, y = as.numeric(PE), fill = Method))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x = "Method", y = "EE (30%)") +
  theme_bw() +
  theme(axis.text.x=element_text(size=15, colour = "black"),
        axis.text.y=element_text(size=15, colour = "black"),
        axis.title=element_text(size=18)) +
  theme(legend.key.size = unit(1.5,"line"))

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