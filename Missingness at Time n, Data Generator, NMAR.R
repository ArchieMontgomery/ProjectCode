
oo <- options(repos = "https://cran.r-project.org/")
install.packages("Matrix")
install.packages("lme4")
install.packages("maxLik")
install.packages("matlib") #optimx, nlmeanb
install.packages("MASS")
install.packages("mice")
install.packages("dplyr")
install.packages("clusterGeneration")
install.packages("ggchicklet", repos = "https://cinc.rud.is")
install.packages("missForest")
install.packages("broom")
install.packages("tibble")
install.packages("broom.mixed")
install.packages("nlme")
install.packages("missMethods")
install.packages("optimx")
library(optimx)
library(clusterGeneration)
library(broom.mixed)
library(mice)
library(broom)
library(ggchicklet)
library(missMethods)
library(tibble)
library(nlme)
library(dplyr)
library(missForest)
library(MASS)
library(matlib)
library(maxLik)
library(lme4)
library(Matrix)
options(oo)


### 2 repeated measures model ###


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
  #for (i in 1:n_measure){
  #  mu = sample(1:100, n_measure)
  #}
  init = sample(1:100, 1)
  final = sample(1:100, 1)
  mu = seq(from = init, to = final, length.out = n_measure)
  range = final - init
  for (i in 1:(n_measure)) {
    mu[i] = sample(rnorm(100, mu[i], 1), 1) 
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
                   method = "SANN", control=list(maxit=200))
      estim_and_true = list(head(ests$par, - (nrepeat + 1)), params_y[[1]])
      mus_and_data = list(estim_and_true, init_vals, params_y[[2]])
      return(mus_and_data)
    }
  }
  NMAR.selection.n.with.data(y_long, time, ID, Y)
}
NMAR_selection_n_with_data(5, 100, 0.1)

##############################################################################

# Function simulating data n times and obtaining estimation Errors for each 
# simulation, and comparing with UMI imputation for last measurement

EE.estimator.real = function(n_sim, n_measure, n_subjects, miss.prop) {
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
    # Data set up
    vals_n_data = NMAR_selection_n_with_data(n_measure, n_subjects, miss.prop)
    true_vals = vals_n_data[[1]][[2]]
    sample_d_long = vals_n_data[[3]]
    
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
    
    # lmer estimates
    lmer_estimates = X_beta[1:n_measure]
    
    # EE (lmer)
    diff.lmer = c(diff.lmer, abs(lmer_estimates - true_vals[1:n_measure]))
    EE.lmer[i] = mean(diff.lmer[(1+(i-1)*n_measure):(n_measure*i)])
    
    ############################ UMI ###############################
    
    # Mean of each Response
    #UMI_1 = mean(sample_d[,2], na.rm = TRUE)
    #UMI_2 = mean(sample_d[,3], na.rm = TRUE)
    #UMI_3 = mean(sample_d[,4], na.rm = TRUE)
    #UMI_4 = mean(sample_d[,5], na.rm = TRUE)
    #UMIs = c(UMI_1, UMI_2, UMI_3, UMI_4)
    
    # Mean of each Response
    UMI_data_wide = reshape(sample_d_long,
                     timevar = "time",
                     idvar = "ID",
                     direction = "wide",
                     sep = "")
    UMI_data_wide = impute_mean(UMI_data_wide, type = "columnwise", 
                                convert_tibble = TRUE)
    
    # Option 1: Calculate mean of each column
    #UMI_means = c()
    #for (j in 1:4) {
    #  UMI_means[j] = mean(UMI_data_wide[,(j+1)])
    #}
    
    #Option 2: Use lmer
    UMI_data_long = reshape(UMI_data_wide, 
                            timevar = "Y", 
                            varying = c(paste("Y",1:n_measure,sep="")), 
                            idvar = "ID",
                            direction ="long",
                            sep = "")
    UMI_data_long = UMI_data_long[order(UMI_data_long$ID),]
    UMI_data_long$time = c(1:n_measure)
    UMI_data_long = UMI_data_long[,c(1,3,2)]
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
    UMI_mu_estimates = X_beta_UMI[1:n_measure]
    
    # Difference between estimates and true
    diff.UMI = c(diff.UMI, abs(UMI_mu_estimates - true_vals[1:n_measure]))
    
    # EE (UMI)
    EE.UMI[i] = mean(diff.UMI[(1+(i-1)*n_measure):(n_measure*i)])
    
    ############################ MI ###############################
    
    ini = mice(sample_d_long, maxit=0) # Long Imputation
    pred1 = ini$predictorMatrix
    pred1[,"ID"] = -2 # set ID as class variable for 2l.norm
    meth1 = ini$method
    meth1[which(meth1 == "pmm")] = "pmm"
    sample_d_long = sample_d_long %>% mutate(ID = as.integer(factor(ID)))
    imp.long = mice(sample_d_long, m=2, method = meth1, predictorMatrix = pred1)
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
    MI_estimates = X_beta_MI[1:n_measure]
    
    # EE (MI)
    diff.MI = c(diff.MI, abs(MI_estimates - true_vals[1:n_measure]))
    EE.MI[i] = mean(diff.MI[(1+(i-1)*n_measure):(n_measure*i)])
    
    ############################ OURS ###############################
    
    # Our function
    sim_ours = vals_n_data[[1]][[1]]
    
    # EE (ours)
    diff.ours = c(diff.ours, abs(sim_ours[1:n_measure] - true_vals[1:n_measure]))
    EE.ours[i] = mean(diff.ours[(1+(i-1)*n_measure):(n_measure*i)])
  }
  end.time = Sys.time()
  time.taken = round(end.time - start.time, 2)
  outputs = list(time.taken, EE.ours, EE.lmer, EE.UMI, EE.MI)
  return(outputs)
}
EE.estimates = EE.estimator.real(2, 3, 50, 0.1)

# Remove outliers

EE.estimates[[2]] = (EE.estimates[[2]][!(EE.estimates[[2]] > 2)])
EE.estimates[[3]] = (EE.estimates[[3]][!(EE.estimates[[3]] > 2)])
EE.estimates[[4]] = (EE.estimates[[4]][!(EE.estimates[[4]] > 2)])
EE.estimates[[5]] = (EE.estimates[[5]][!(EE.estimates[[5]] > 2)])

EE.estimates.7.50.0.3[[2]] = (EE.estimates.7.50.0.1[[2]][!(EE.estimates.7.50.0.1[[2]] > 2)])
EE.estimates.7.50.0.1[[3]] = (EE.estimates.7.50.0.1[[3]][!(EE.estimates.7.50.0.1[[3]] > 2)])
EE.estimates.7.50.0.1[[4]] = (EE.estimates.7.50.0.1[[4]][!(EE.estimates.7.50.0.1[[4]] > 2)])
EE.estimates.7.50.0.1[[5]] = (EE.estimates.7.50.0.1[[5]][!(EE.estimates.7.50.0.1[[5]] > 2)])

# Calculate means

round(mean(EE.estimates[[2]]), 3)
round(mean(EE.estimates[[3]]), 3)
round(mean(EE.estimates[[4]]), 3)
round(mean(EE.estimates[[5]]), 3)


# 20 subjects variations
EE.estimates.3.20.0.1 = EE.estimates
EE.estimates.3.20.0.2 = EE.estimates
EE.estimates.3.20.0.3 = EE.estimates
EE.estimates.5.20.0.1 = EE.estimates
EE.estimates.5.20.0.2 = EE.estimates
EE.estimates.5.20.0.3 = EE.estimates
EE.estimates.7.20.0.1 = EE.estimates
EE.estimates.7.20.0.2 = EE.estimates
EE.estimates.7.20.0.3 = EE.estimates

# 50 subjects variations
EE.estimates.3.50.0.1 = EE.estimates
EE.estimates.3.50.0.2 = EE.estimates
EE.estimates.3.50.0.3 
#EE.estimates.5.50.0.1 = EE.estimates
EE.estimates.5.50.0.2 = EE.estimates
EE.estimates.5.50.0.3 = EE.estimates
#EE.estimates.7.50.0.1 = EE.estimates

EE.estimates.full = c(EE.estimates.10, EE.estimates.20, EE.estimates.30, EE.estimates.40)

##############################################################################

# Create data frame of estimates for each method

## Bar Chart
our_estimate.3.20.10 = as.data.frame(cbind(EE.estimates.3.20.0.1[[2]], "Ours", "10%"))
LMER_estimate.3.20.10 = as.data.frame(cbind(EE.estimates.3.20.0.1[[3]], "CCA", "10%"))
UMI_estimate.3.20.10 = as.data.frame(cbind(EE.estimates.3.20.0.1[[4]], "UMI", "10%"))
MI_estimate.3.20.10 = as.data.frame(cbind(EE.estimates.3.20.0.1[[5]], "MI", "10%"))

our_estimate.3.20.20 = as.data.frame(cbind(EE.estimates.3.20.0.2[[2]], "Ours", "20%"))
LMER_estimate.3.20.20 = as.data.frame(cbind(EE.estimates.3.20.0.2[[3]], "CCA", "20%"))
UMI_estimate.3.20.20 = as.data.frame(cbind(EE.estimates.3.20.0.2[[4]], "UMI", "20%"))
MI_estimate.3.20.20 = as.data.frame(cbind(EE.estimates.3.20.0.2[[5]], "MI", "20%"))

our_estimate.3.20.30 = as.data.frame(cbind(EE.estimates.3.20.0.3[[2]], "Ours", "30%"))
LMER_estimate.3.20.30 = as.data.frame(cbind(EE.estimates.3.20.0.3[[3]], "CCA", "30%"))
UMI_estimate.3.20.30 = as.data.frame(cbind(EE.estimates.3.20.0.3[[4]], "UMI", "30%"))
MI_estimate.3.20.30 = as.data.frame(cbind(EE.estimates.3.20.0.3[[5]], "MI", "30%"))

EE_dataframe = rbind(our_estimate.3.20.10, LMER_estimate.3.20.10, 
                     UMI_estimate.3.20.10, MI_estimate.3.20.10,
                     our_estimate.3.20.20, LMER_estimate.3.20.20, 
                     UMI_estimate.3.20.20, MI_estimate.3.20.20,
                     our_estimate.3.20.30, LMER_estimate.3.20.30, 
                     UMI_estimate.3.20.30, MI_estimate.3.20.30)

colnames(EE_dataframe) = c("EE", "Method", "Missingness")
EE_dataframe$EE = as.numeric(EE_dataframe$EE)
new_df = subset(EE_dataframe, EE<2) 

EE_dataframe = EE_dataframe %>% mutate(Method=factor(Method, 
                                                     levels=c("Ours", "CCA", "UMI", "MI")))

data_summary <- function(data, varname, groupnames){
  require(dplyr)
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum = rename(data_sum, c("mean" = varname))
  return(data_sum)
}
EE_dataframe = data_summary(EE_dataframe, varname="EE", 
                    groupnames=c("Missingness", "Method"))

library(RColorBrewer)
install.packages("wesanderson")
library(wesanderson)
display.brewer.all(colorblindFriendly = TRUE)

ggplot(EE_dataframe, aes(x=Method, y=as.numeric(EE), fill=Missingness)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_hline(yintercept = mean(EE_dataframe$EE), 
             linetype = "dashed", size = 1) +
  #coord_polar() +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3)) +
  #geom_chicklet(radius = grid::unit(10, "mm")) +
  #geom_rect(aes(xmin = 0.5, xmax = 3.5, ymin = 0, ymax = 1.95), fill = "gray95") +
  labs(x = "Method", y = "EE") +
  #geom_text(aes(label = EE), position = position_dodge(0.9), 
  #          vjust = 2, size = 4, color = "#ffffff") + 
  theme_bw() +
  theme(axis.text.x=element_text(size=15, colour = "black"),
        axis.text.y=element_text(size=15, colour = "black"),
        axis.title=element_text(size=18),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        #axis.text.x = element_text(face = "bold", size = 18),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.key.size = unit(2.5, "line"),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12))
  #geom_errorbar(aes(ymin=EE-sd, ymax=EE+sd), width=.2,
  #              position=position_dodge(.5)) 

## 3 RM, 50 Subjects, 10% missingness
our_estimate = as.data.frame(cbind(EE.estimates.3.50.0.1[[2]], "Ours"))
LMER_estimate = as.data.frame(cbind(EE.estimates.3.50.0.1[[3]], "CCA"))
UMI_estimate = as.data.frame(cbind(EE.estimates.3.50.0.1[[4]], "UMI"))
MI_estimate = as.data.frame(cbind(EE.estimates.3.50.0.1[[5]], "MI"))

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
        axis.title=element_text(size=18),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.key.size = unit(1.5,"line"))

p <- ggplot(df2, aes(x=dose, y=len, fill=supp)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(.9)) 

##############################################################################

# Create data frame of estimates for each method
our_estimate = as.data.frame(cbind(EE.estimates.5.50.0.1[[2]], "Ours"))
LMER_estimate = as.data.frame(cbind(EE.estimates.5.50.0.1[[3]], "CCA"))
UMI_estimate = as.data.frame(cbind(EE.estimates.5.50.0.1[[4]], "UMI"))
MI_estimate = as.data.frame(cbind(EE.estimates.5.50.0.1[[5]], "MI"))

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
our_estimate = as.data.frame(cbind(EE.estimates.7.50.0.1[[2]], "Ours"))
LMER_estimate = as.data.frame(cbind(EE.estimates.7.50.0.1[[3]], "CCA"))
UMI_estimate = as.data.frame(cbind(EE.estimates.7.50.0.1[[4]], "UMI"))
MI_estimate = as.data.frame(cbind(EE.estimates.7.50.0.1[[5]], "MI"))

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

# Extra Work
EE.estimator = function(n_sim, n_measures, n_subjects, miss.prop) {
  start.time = Sys.time()
  diff.ours = NULL
  PE.ours = numeric(n_sim)
  for (i in 1:length(PE.ours)){
    # Our function
    sim_ours = NMAR_selection_n_with_data(n_measures, n_subjects, miss.prop)
    diff.ours = c(diff.ours, abs(sim_ours[[1]][[1]] - sim_ours[[1]][[2]]))
    PE.ours[i] = mean(diff.ours[i])
  }
  end.time = Sys.time()
  time.taken = round(end.time - start.time,2)
  outputs = list(time.taken, PE.ours)
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
X = EE.estimator(n_sim = 50, n_measures = 2, n_subjects = 200, miss.prop = 0.2)
X = unlist(X[[2]])
X = as.vector(X)
X2 = X[!(X > 4)]
length(X[(X > 4)])
X1 = X[!(X > 3)]
length(X[(X > 3)])


mean(X)
mean(X2)
mean(X1)
