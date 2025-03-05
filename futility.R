#This file runs the simulations of the paper "Improving the (approximate) sequential probability ratio 
#test by avoiding overshoot" that include a stop for futility

library(nloptr)

#Load the boosting functions
source("boosting_functions.R")

#Set the simulation parameters
betas=c(0.01,0.05,0.1,0.15, 0.2, 0.25, 0.3)
alpha=0.05
mu_N=0
mu_A=0.3
rho=0.583
m=10000
n=10000

#Initialize vectors for stopping times and power
mean_stop_sprt=rep(0,length(betas))
mean_stop_sprt_cons=rep(0,length(betas))
mean_stop_sprt_ds=rep(0,length(betas))
mean_stop_boosted=rep(0,length(betas))
power_sprt=rep(0,length(betas))
power_sprt_cons=rep(0,length(betas))
power_sprt_ds=rep(0,length(betas))
power_boosted=rep(0,length(betas))
mean_type_I_sprt=rep(0,length(betas))
mean_type_I_sprt_cons=rep(0,length(betas))
mean_type_I_sprt_ds=rep(0,length(betas))
mean_type_I_boosted=rep(0,length(betas))

#Set seed for reproducibility
set.seed(123)

count=1
for(beta in betas){

  stop_sprt=rep(n,m)                #Stopping times for a specific beta (SPRT approx. thresholds)
  stop_sprt_cons=rep(n,m)           #Stopping times for a specific beta (SPRT cons. thresholds)
  stop_sprt_ds=rep(n,m)             #Stopping times for a specific beta (SPRT Siegmund's thresholds)
  stop_boosted=rep(n,m)             #Stopping times for a specific beta (boosted SPRT)
  decision_sprt=rep(0,m)            #decisions for a specific beta (SPRT approx. thresholds)
  decision_sprt_cons=rep(0,m)       #decisions for a specific beta (SPRT cons. thresholds)
  decision_sprt_ds=rep(0,m)         #decisions for a specific beta (SPRT Siegmund's thresholds)
  decision_boosted=rep(0,m)         #decisions for a specific beta (boosted SPRT)
  type_I_sprt=rep(0,m)              #Type I errors for a specific beta (SPRT approx. thresholds)
  type_I_sprt_cons=rep(0,m)         #Type I errors for a specific beta (SPRT cons. thresholds)
  type_I_sprt_ds=rep(0,m)           #Type I errors for a specific beta (SPRT Siegmund's thresholds)
  type_I_boosted=rep(0,m)           #Type I errors for a specific beta (boosted SPRT)
  
  for(j in 1:m){

    data=rnorm(n,mean=mu_A,sd=1)                #Create normally distributed data with mean mu_A and variance 1
    
    #SPRT
    lr=dnorm(data,mu_A,1)/dnorm(data,mu_N,1)    #Likelihood ratio
    LR=cumprod(lr)                              #Cumulative likelihood ratio
    stop_sprt[j]=min(c(which(LR>=(1-beta)/alpha | LR<=(beta)/(1-alpha)), n))    #Calculate stop for SPRT
    decision_sprt[j]=(LR[stop_sprt[j]]>=(1-beta)/alpha)                         #Calculate decision for SPRT
    
    stop_sprt_cons[j]=min(c(which(LR>=1/alpha | LR<=beta), n))    #Calculate stop for cons. SPRT
    decision_sprt_cons[j]=(LR[stop_sprt_cons[j]]>=1/alpha)        #Calculate decision for cons. SPRT
    
    stop_sprt_ds[j]=min(c(which(LR>=(1-beta)/(alpha*exp(mu_A*rho)) | LR<=exp(mu_A*rho)*beta/(1-alpha)), n)) #Calculate stop for Siegmund's SPRT
    decision_sprt_ds[j]=(LR[stop_sprt_ds[j]]>=(1-beta)/(alpha*exp(mu_A*rho)))                               #Calculate decision for Siegmund's SPRT
    
    #Parameters for boosted SPRT
    b=rep(1,n)                      
    lr_boosted=lr 
    lr_inv=1/lr
    lr_inv_boosted=lr_inv
    b_inv=rep(1,n) 
    nu_t=beta
    
    for(i in 1:n){
      
      #Calculate the boosting factors
      b_full=e_boosted_type2(prod(lr_boosted[1:(i-1)]), prod(lr_inv_boosted[1:(i-1)]), prod(b[1:(i-1)]), prod(b_inv[1:(i-1)]), alpha, beta, mu_A, -mu_A)
      
      b_inv[i]=b_full[2]
      b[i]=b_full[1]
      nu_t=b_full[3]
      
      #Calculate boosted LR and boosted inverse LR
      lr_inv_boosted[i]=b_inv[i]*lr_inv_boosted[i]
      lr_boosted[i]=b[i]*lr_boosted[i]
      
      
      #Set stopping time and decision for boosted SPRT
      if(prod(lr_boosted[1:i])>=(1/alpha)|prod(lr_boosted[1:i])<=nu_t){
        stop_boosted[j]=i
        decision_boosted[j]=(prod(lr_boosted[1:i])>=(1/alpha))
        break
      }
  
    }
    #Calculate the estimate for the type I error in each trial
    type_I_sprt[j]=(1/prod(lr[1:stop_sprt[j]]))*decision_sprt[j]
    type_I_sprt_cons[j]=(1/prod(lr[1:stop_sprt_cons[j]]))*decision_sprt_cons[j]
    type_I_sprt_ds[j]=(1/prod(lr[1:stop_sprt_ds[j]]))*decision_sprt_ds[j]
    type_I_boosted[j]=(1/prod(lr[1:stop_boosted[j]]))*decision_boosted[j]
  }
  #Save the mean power and stops for each mu_A
  power_boosted[count]=mean(decision_boosted)
  power_sprt[count]=mean(decision_sprt)
  power_sprt_cons[count]=mean(decision_sprt_cons)
  power_sprt_ds[count]=mean(decision_sprt_ds)
  mean_stop_boosted[count]=mean(stop_boosted)
  mean_stop_sprt[count]=mean(stop_sprt)
  mean_stop_sprt_cons[count]=mean(stop_sprt_cons)
  mean_stop_sprt_ds[count]=mean(stop_sprt_ds)
  mean_type_I_boosted[count]=mean(type_I_boosted)
  mean_type_I_sprt[count]=mean(type_I_sprt)
  mean_type_I_sprt_ds[count]=mean(type_I_sprt_ds)
  count=count+1
}


#Save the results
results_df=data.frame(idx=betas, mean_stop_boosted, mean_stop_sprt, power_boosted, power_sprt, 
                      mean_type_I_boosted, mean_type_I_sprt, mean_stop_sprt_cons, power_sprt_cons, mean_type_I_sprt_cons,
                      mean_stop_sprt_ds, power_sprt_ds, mean_type_I_sprt_ds)
save(results_df, file="results/futility.rda")


###############################same simulations for alpha=0.01

alpha=0.01

#Initialize vectors for stopping times and power
mean_stop_sprt=rep(0,length(betas))
mean_stop_sprt_cons=rep(0,length(betas))
mean_stop_sprt_ds=rep(0,length(betas))
mean_stop_boosted=rep(0,length(betas))
power_sprt=rep(0,length(betas))
power_sprt_cons=rep(0,length(betas))
power_sprt_ds=rep(0,length(betas))
power_boosted=rep(0,length(betas))
mean_type_I_sprt=rep(0,length(betas))
mean_type_I_sprt_cons=rep(0,length(betas))
mean_type_I_sprt_ds=rep(0,length(betas))
mean_type_I_boosted=rep(0,length(betas))

#Set seed for reproducibility
set.seed(123)

count=1
for(beta in betas){
  
  stop_sprt=rep(n,m)                #Stopping times for a specific beta (SPRT approx. thresholds)
  stop_sprt_cons=rep(n,m)           #Stopping times for a specific beta (SPRT cons. thresholds)
  stop_sprt_ds=rep(n,m)             #Stopping times for a specific beta (SPRT Siegmund's thresholds)
  stop_boosted=rep(n,m)             #Stopping times for a specific beta (boosted SPRT)
  decision_sprt=rep(0,m)            #decisions for a specific beta (SPRT approx. thresholds)
  decision_sprt_cons=rep(0,m)       #decisions for a specific beta (SPRT cons. thresholds)
  decision_sprt_ds=rep(0,m)         #decisions for a specific beta (SPRT Siegmund's thresholds)
  decision_boosted=rep(0,m)         #decisions for a specific beta (boosted SPRT)
  type_I_sprt=rep(0,m)              #Type I errors for a specific beta (SPRT approx. thresholds)
  type_I_sprt_cons=rep(0,m)         #Type I errors for a specific beta (SPRT cons. thresholds)
  type_I_sprt_ds=rep(0,m)           #Type I errors for a specific beta (SPRT Siegmund's thresholds)
  type_I_boosted=rep(0,m)           #Type I errors for a specific beta (boosted SPRT)
  
  for(j in 1:m){
    
    data=rnorm(n,mean=mu_A,sd=1)                #Create normally distributed data with mean mu_A and variance 1
    
    #SPRT
    lr=dnorm(data,mu_A,1)/dnorm(data,mu_N,1)    #Likelihood ratio
    LR=cumprod(lr)                              #Cumulative likelihood ratio
    stop_sprt[j]=min(c(which(LR>=(1-beta)/alpha | LR<=(beta)/(1-alpha)), n))    #Calculate stop for SPRT
    decision_sprt[j]=(LR[stop_sprt[j]]>=(1-beta)/alpha)                         #Calculate decision for SPRT
    
    stop_sprt_cons[j]=min(c(which(LR>=1/alpha | LR<=beta), n))    #Calculate stop for cons. SPRT
    decision_sprt_cons[j]=(LR[stop_sprt_cons[j]]>=1/alpha)        #Calculate decision for cons. SPRT
    
    stop_sprt_ds[j]=min(c(which(LR>=(1-beta)/(alpha*exp(mu_A*rho)) | LR<=exp(mu_A*rho)*beta/(1-alpha)), n)) #Calculate stop for Siegmund's SPRT
    decision_sprt_ds[j]=(LR[stop_sprt_ds[j]]>=(1-beta)/(alpha*exp(mu_A*rho)))                               #Calculate decision for Siegmund's SPRT
    
    #Parameters for boosted SPRT
    b=rep(1,n)                      
    lr_boosted=lr 
    lr_inv=1/lr
    lr_inv_boosted=lr_inv
    b_inv=rep(1,n) 
    nu_t=beta
    
    for(i in 1:n){
      
      #Calculate the boosting factors
      b_full=e_boosted_type2(prod(lr_boosted[1:(i-1)]), prod(lr_inv_boosted[1:(i-1)]), prod(b[1:(i-1)]), prod(b_inv[1:(i-1)]), alpha, beta, mu_A, -mu_A)
      
      b_inv[i]=b_full[2]
      b[i]=b_full[1]
      nu_t=b_full[3]
      
      #Calculate boosted LR and boosted inverse LR
      lr_inv_boosted[i]=b_inv[i]*lr_inv_boosted[i]
      lr_boosted[i]=b[i]*lr_boosted[i]
      
      
      #Set stopping time and decision for boosted SPRT
      if(prod(lr_boosted[1:i])>=(1/alpha)|prod(lr_boosted[1:i])<=nu_t){
        stop_boosted[j]=i
        decision_boosted[j]=(prod(lr_boosted[1:i])>=(1/alpha))
        break
      }
      
    }
    #Calculate the estimate for the type I error in each trial
    type_I_sprt[j]=(1/prod(lr[1:stop_sprt[j]]))*decision_sprt[j]
    type_I_sprt_cons[j]=(1/prod(lr[1:stop_sprt_cons[j]]))*decision_sprt_cons[j]
    type_I_sprt_ds[j]=(1/prod(lr[1:stop_sprt_ds[j]]))*decision_sprt_ds[j]
    type_I_boosted[j]=(1/prod(lr[1:stop_boosted[j]]))*decision_boosted[j]
  }
  #Save the mean power and stops for each mu_A
  power_boosted[count]=mean(decision_boosted)
  power_sprt[count]=mean(decision_sprt)
  power_sprt_cons[count]=mean(decision_sprt_cons)
  power_sprt_ds[count]=mean(decision_sprt_ds)
  mean_stop_boosted[count]=mean(stop_boosted)
  mean_stop_sprt[count]=mean(stop_sprt)
  mean_stop_sprt_cons[count]=mean(stop_sprt_cons)
  mean_stop_sprt_ds[count]=mean(stop_sprt_ds)
  mean_type_I_boosted[count]=mean(type_I_boosted)
  mean_type_I_sprt[count]=mean(type_I_sprt)
  mean_type_I_sprt_ds[count]=mean(type_I_sprt_ds)
  count=count+1
}


#Save the results
results_df=data.frame(idx=betas, mean_stop_boosted, mean_stop_sprt, power_boosted, power_sprt, 
                      mean_type_I_boosted, mean_type_I_sprt, mean_stop_sprt_cons, power_sprt_cons, mean_type_I_sprt_cons,
                      mean_stop_sprt_ds, power_sprt_ds, mean_type_I_sprt_ds)
save(results_df, file="results/futility_alpha001.rda")