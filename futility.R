#This file runs the simulations shown in Section 5 of the paper 
#"Improving the (approximate) sequential probability ratio test by avoiding overshoot"

library(nleqslv)
library(BB)

#Load the boosting functions
source("boosting_functions.R")

#Set the simulation parameters
betas=c(0.01,0.05,0.1,0.15, 0.2, 0.25, 0.3)
alpha=0.05
mu_N=0
mu_A=0.3
m=10000
n=10000

#Initialize vectors for stopping times and power
mean_stop_sprt=rep(0,length(betas))
mean_stop_boosted=rep(0,length(betas))
power_sprt=rep(0,length(betas))
power_boosted=rep(0,length(betas))
mean_type_I_sprt=rep(0,length(betas))
mean_type_I_boosted=rep(0,length(betas))
options(warn=2)
#Set seed for reproducibility
set.seed(123)

count=1
for(beta in betas){

  stop_sprt=rep(n,m)          #Stopping times for a specific beta
  stop_boosted=rep(n,m)       #Stopping times for a specific beta
  decision_sprt=rep(0,m)      #decisions for a specific beta
  decision_boosted=rep(0,m)   #decisions for a specific beta
  type_I_sprt=rep(0,m)        #Type I errors for a specific beta
  type_I_boosted=rep(0,m)     #Type I errors for a specific beta
  
  for(j in 1:m){

    data=rnorm(n,mean=mu_A,sd=1)                #Create normally distributed data with mean mu_A and variance 1
    
    #SPRT
    lr=dnorm(data,mu_A,1)/dnorm(data,mu_N,1)    #Likelihood ratio
    LR=cumprod(lr)                              #Cumulative likelihood ratio
    stop_sprt[j]=min(c(which(LR>=(1-beta)/alpha | LR<=(beta)/(1-alpha)), n))    #Calculate stop for SPRT
    decision_sprt[j]=(LR[stop_sprt[j]]>=(1-beta)/alpha)                         #Calculate decision for SPRT
    
    #Boosted SPRT
    b=rep(1,n)                      #boosting factors
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
    type_I_boosted[j]=(1/prod(lr[1:stop_boosted[j]]))*decision_boosted[j]
  }
  #Save the mean power and stops for each mu_A
  power_boosted[count]=mean(decision_boosted)
  power_sprt[count]=mean(decision_sprt)
  mean_stop_boosted[count]=mean(stop_boosted)
  mean_stop_sprt[count]=mean(stop_sprt)
  mean_type_I_boosted[count]=mean(type_I_boosted)
  mean_type_I_sprt[count]=mean(type_I_sprt)
  count=count+1
}

#Save the results
results_df=data.frame(idx=betas, mean_stop_boosted, mean_stop_sprt, power_boosted, power_sprt, 
                      mean_type_I_boosted, mean_type_I_sprt)
save(results_df, file="results/futility.rda")


