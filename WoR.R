#This file runs the simulations shown in Section 4.2 of the paper 
#"Improving the (approximate) sequential probability ratio test by avoiding overshoot"

#Set simulation parameters
mu_N=0.5
mu_A=0.55
alpha=0.01

m=1000
ns=c(1000, 2000, 5000, 10000, 20000, 50000, 100000)

#Initialize vectors for stopping times and power
mean_stop_sprt=rep(0,length(ns))
mean_stop_boosted=rep(0,length(ns))
power_sprt=rep(0,length(ns))
power_boosted=rep(0,length(ns))

#Set seed for reproducibility
set.seed(123)

count=1
for(n in ns){
  
  stop_sprt=rep(n,m)          #Stopping times for a specific n
  stop_boosted=rep(n,m)       #Stopping times for a specific n
  decision_sprt=rep(0,m)      #decisions for a specific n
  decision_boosted=rep(0,m)   #decisions for a specific n
  for(j in 1:m){
    data=sample(c(rep(1,n*mu_A),rep(0,n*(1-mu_A))))       #Create binary data of length n with mu_A*n ones
    m_t=1/alpha                                           #Used for boosting
    mart_fac=rep(1,n)                                     #Martingale factors
    mart_fac_boosted=rep(1,n)                             #Boosted martingale factors
    mu_A_pred=rep(0,n)                                    #Conditional mean 
    lambda=rep(0,n)                                       #Hyperparameter used for RiLACS 
    boost_factor=rep(1,n)                                 #Obtained boosting factors
    
    mu_A_pred[1]=mu_N
    
    
    for(i in 1:(n-1)){
      #Calculate conditional mean
      if(i>1){
        mu_A_pred[i]=(n*mu_N-sum(data[1:(i-1)]))/(n-i+1)
      }
      #Calculate hyperparameter lambda
      if(mu_A_pred[i]==0){
        lambda[i]=Inf                                
      }else{
        lambda[i]=min(1/mu_A_pred[i], 2*(2*mu_A-1))
      }
      #Calculate boosting factor
      if(i>1 & ((1+lambda[i]*(1-mu_A_pred[i]))*prod(mart_fac_boosted[1:(i-1)])>=1/alpha) & stop_boosted[j]==n & m_t>1 & lambda[i]<Inf){
        boost_factor[i]=(1-m_t*mu_A_pred[i])/((1-mu_A_pred[i])*(1-lambda[i]*mu_A_pred[i]))
      }
      #Calculate martingale factor
      if(stop_sprt[j]==n){
        if(data[i]==0 & lambda[i]==Inf){
          mart_fac[i]=1
        }else{
          mart_fac[i]=(1+lambda[i]*(data[i]-mu_A_pred[i]))
        }
      }
      #Calculate boosted martingale factor
      if(stop_boosted[j]==n){
        if(data[i]==0 & lambda[i]==Inf){
          mart_fac_boosted[i]=1
        }else{
          mart_fac_boosted[i]=boost_factor[i]*mart_fac[i]
        }
      }
      
      #Calculate m_t for the next step (used for boosting)
      m_t=1/(alpha*prod(mart_fac_boosted[1:i]))
      
      #Set stopping times and decisions
      if(prod(mart_fac_boosted[1:i])>=1/alpha & stop_boosted[j]==n){
        stop_boosted[j]=i
        decision_boosted[j]=(prod(mart_fac_boosted[1:i])>=1/alpha)
      }
      if(prod(mart_fac[1:i])>=1/alpha & stop_sprt[j]==n){
        stop_sprt[j]=i
        decision_sprt[j]=(prod(mart_fac[1:i])>=1/alpha)
        break
      }
      
    }
  }
  #Save the mean power and stops for each n
  power_boosted[count]=mean(decision_boosted)
  power_sprt[count]=mean(decision_sprt)
  mean_stop_boosted[count]=mean(stop_boosted)
  mean_stop_sprt[count]=mean(stop_sprt)
  count=count+1
}

#Save the results
results_df=data.frame(idx=ns, power_boosted, power_sprt, mean_stop_boosted, mean_stop_sprt)
save(results_df, file="results/WoR.rda")

