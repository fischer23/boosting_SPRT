#This file contains the functions for boosting the martingale factors in the paper 
#"Improving the (approximate) sequential probability ratio test by avoiding overshoot"

###Function to compute boosting factors for Gaussian testing problems
##Input:
#m_t:           cutoff value >=1 for the truncation function (m_t=1/(M*\alpha)).  
#delta:         Real parameter used for boosting the martingale factors. In case of the simple null H_i:X_i~N(mu_N,1)
#               vs. alternative H_i^A:X_i~N(mu_A,1) delta should be set to delta=mu_A-mu_N.

###Output:      boosting factor.

e_boosted=function(m_t, delta){
  b_factor=function(b){
    return(b*(1-pnorm(delta/2-log(m_t/b)/delta))+m_t*(1-pnorm((log(m_t/b)+delta^2/2)/delta))-1)
  }
  return(uniroot(b_factor, lower=0, upper=100000000000)$root)
}

###Function to compute boosting for Vovk's conformal martingales
##Input:
#m_t:           cutoff value for the truncation function (m_t=1/(M*alpha)).  
#kappa:         Parameter in (0,1) used for Vovk's calibrator (f_{\kappa}(u)=\kappa u^{\kappa-1})

###Output:      boosting factor.

e_boosted_calib=function(m_t, kappa){
  b_factor=function(b){
    return(b*max((1-(m_t/(b*kappa))^(kappa/(kappa-1))),0)+m_t*min((m_t/(b*kappa))^(1/(kappa-1)),1)-1)
  }
  return(uniroot(b_factor, lower=0, upper=100000000000)$root)
}


###Function to compute boosting factors for Gaussian testing problems with stop for futility.
##Input:
#m_t:           cutoff value >=1 for the truncation function (m_t=1/(M*\alpha)).  
#delta:         Real parameter used for boosting the martingale factors. In case of the simple null H_i:X_i~N(mu_N,1)
#               vs. alternative H_i^A:X_i~N(mu_A,1) delta should be set to delta=mu_A-mu_N.
#l_t:           lower cutoff value <1 for the truncation function (l_t=\nu_t/M).

###Output:      boosting factor.


e_boosted_futility=function(m_t, l_t, delta){
  b_factor=function(b){
    return(b*(pnorm(delta/2-log(l_t/b)/delta)-pnorm(delta/2-log(m_t/b)/delta))+m_t*(1-pnorm((log(m_t/b)+delta^2/2)/delta))-1)
  }
  return(uniroot(b_factor, lower=0, upper=100000000000)$root)
}