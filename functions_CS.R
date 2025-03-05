#This file contains auxiliary functions for calculating a boosted CS as described in Section S.1 of the paper 
#"Improving the (approximate) sequential probability ratio test by avoiding overshoot"

###Auxiliary function for the construction of boosted CSs.

##Input:
#mu_N:          Mean value for which to be checked whether H_0:mu=mu_N is rejected  
#alpha:         Individual significance level
#n_ind:         Time for which it is checked whether H_0 is rejected
#delta:         Parameter for the CS
#data:          Data for which CS is to be computed

###Output:      Wealth - 1/\alpha at step n_ind  

reject_boosted=function(mu_N, alpha, n_ind, delta, data){
  Z=data-mu_N
  m_t=1/alpha
  wealth=1
  for(i in 1:n_ind){
      b=e_boosted(m_t,delta)
      wealth=wealth*max(1,b)*exp(Z[i]*delta-delta^2/2)
      m_t=max(1/(alpha*wealth),1)
  }
  return(min(wealth-1/alpha,10))
}

###Auxiliary function for the construction of boosted CSs.

##Input:
#alpha:         Individual significance level
#n:             Time up to which CS is to be constructed
#delta:         Parameter for the CS
#data:          Data for which CS is to be computed

###Output:      #Values of mu_N such that reject_boosted(mu_N, alpha, t, data) equals exactly 0 for all steps
#               from 1 to n. These are the smallest \mu_N such that H_0:mu=mu_N is not rejected.

boosted_CS=function(alpha, n, delta, data){
  CS=rep(0,n)
  for(j in 1:n){
    CS[j]=uniroot(function(x) reject_boosted(x, alpha, j, delta, data), lower=-10, upper=10)$root
  }
  return(cummax(CS))
}



