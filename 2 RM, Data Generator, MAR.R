
oo <- options(repos = "https://cran.r-project.org/")
install.packages("Matrix")
install.packages("lme4")
install.packages("maxLik")
install.packages("matlib") #optimx, nlmeanb
install.packages("MASS")
install.packages("mice")
install.packages("clusterGeneration")
install.packages("missForest")
install.packages("tibble")
install.packages("nlme")
install.packages("missMethods")
install.packages("optimx")
library(optimx)
library(clusterGeneration)
library(mice)
library(missMethods)
library(tibble)
library(nlme)
library(missForest)
library(MASS)
library(matlib)
library(maxLik)
library(lme4)
library(Matrix)
options(oo)


### 2 repeated measures model ###



# Function to generate any data set with m subject, n individuals and r proportion
# of missingness in the last measurement
sim_params = c()
one_dataset = function(n_measure, n_indiv, miss.proportion)   {
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
y_long = one_dataset(9, 45, 0.05)[[2]]
#mean1 = mean(which(y_long$time == 1))
#params = dataset(3, 100, 0.15)[[1]]


###                 Function for n repeated measures of a normal model, 
###                             returning MLEs of parameters
#NMAR.selection.n.with.data(y_long, time, ID, Y)

##################################################################
##################################################################
##################################################################

sim_params = c()
NMAR_selection_n_with_data = function(n_measure, n_indiv, miss.proportion)   {
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
  
  NMAR.selection.n.with.data = function(data, timevar, idvar, resvar){
    
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
        mu = matrix(c(theta[1:2]),nrow=1, ncol=2)
        Sigma = matrix(c(theta[3:6]), nrow=nrepeat, ncol=nrepeat, byrow=TRUE)
        detsig = det(Sigma)
        invsig = as.matrix(solve(Sigma))
        
        ll1 = 0
        ## Part 1 ##
        for (i in 1:comp.cases) {
          #Ys = as.matrix(cbind(1, complete.data[i,]))
          #logits = phis%*%t(Ys)        
          #p = exp(logits)/(1 + exp(logits))
          diff = as.matrix(complete.data[i,] - mu)
          mat1 = diff%*%invsig
          mat1 = mat1%*%t(diff)
          mat1 = as.numeric(mat1)
          ll1 = ll1 + (-1/2)*log(detsig) - (1/2)*mat1
        }
        
        ll2 = 0
        ## Part 2 ##
        for (i in 1:incomp.cases) {
          #Ys = cbind(1,incomplete.data[i,1:(nrepeat-1)],ql)
          #logits = phis%*%t(Ys)
          #p = exp(logits)/(1+exp(logits))
          diff = as.matrix(incomplete.data[i,1] - mu[1])
          mat1 = diff^2/(Sigma[1,1])
          ll2 = ll2 - 0.5*log(Sigma[1,1]) - 0.5*mat1
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
                   method = "SANN", control=list(maxit=500))
      estim_and_true = list(head(ests$par, - (nrepeat + 1)), params_y[[1]])
      mus_and_data = list(estim_and_true, init_vals, params_y[[2]])
      return(mus_and_data)
    }
  }
  NMAR.selection.n.with.data(y_long, time, ID, Y)
}
X = NMAR_selection_n_with_data(2, 100, 0.2)
X1 = abs(X[[1]][[1]][1] - X[[1]][[1]][2])
X2 = abs(X[[1]][[2]][1] - X[[1]][[2]][2])
abs(X2 - X1)
##############################################################################

# Function simulating data n times and obtaining Prediction Errors for each simulation,
# and comparing with UMI imputation for last measurement

PE.estimator = function(n_sim, n_measures, n_subjects, miss.prop) {
  start.time = Sys.time()
  mu_diff_ours = NULL
  mu_diff_true = NULL
  mu_diff_est = numeric(n_sim)
  #PE.ours = numeric(n_sim)
  for (i in 1:length(mu_diff_est)){
    # Our function
    sim_ours = NMAR_selection_n_with_data(n_measures, n_subjects, miss.prop)
    mu_diff_ours = c(mu_diff_ours, abs(sim_ours[[1]][[1]][1] - sim_ours[[1]][[1]][2]))
    mu_diff_true = c(mu_diff_true, abs(sim_ours[[1]][[2]][1] - sim_ours[[1]][[2]][2]))
    mu_diff_est = abs(mu_diff_true - mu_diff_ours)
    #PE.ours[i] = mean(diff.ours[i])
  }
  end.time = Sys.time()
  time.taken = round(end.time - start.time,2)
  outputs = list(time.taken, mu_diff_ours, mu_diff_true, mu_diff_est)
  return(outputs)
}

# n_sim = number of simulations to perform (always use 50)
# n_measures = number of measurements in data (3, 5, 7, 9)
# n_subjects = number of subjects (20, 50, 100)
# miss_prop = proportion of missingness in the measurements (0.1, 0.2, 0.3)

# run this next line excatly, only changing n_measures, using 
# n_measures = 3
# n_measures = 5
# n_measures = 7
# n_measures = 9
X = PE.estimator(n_sim = 10, n_measures = 2, n_subjects = 500, miss.prop = 0.2)
X = unlist(X[[4]])
X = as.vector(X)
X = X[!(X > 3)]

mean(X)
