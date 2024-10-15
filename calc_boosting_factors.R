#This file calculates the boosting factors reported in the Sections 3.1, 4.3 and 5 of the paper 
#"Improving the (approximate) sequential probability ratio test by avoiding overshoot"

source("boosting_functions.R")


###Calculations for Section 3.1

deltas=c(0.1,0.5,1,2,3)
LRs=c(0.5,1,2,4,10)
alpha=0.05

boosting_factors=matrix(, nrow=length(deltas), ncol=length(LRs))

for(delta in 1:length(deltas)){
  for(LR in 1:length(LRs)){
   m_t=1/(alpha*LRs[LR]) 
    boosting_factors[delta, LR]=e_boosted(m_t, deltas[delta])
  }
}

boosting_factors

#Calculations for Section 4.3
e_boosted_calib(20, 0.5)
e_boosted_calib(20, 0.1)


#Calculations for Section 5
deltas=c(0.1,0.5,1,2,3)
LRs=c(0.5,1,2,4,10)
alpha=0.05
nu=0.4

boosting_factors=matrix(, nrow=length(deltas), ncol=length(LRs))

for(delta in 1:length(deltas)){
  for(LR in 1:length(LRs)){
    m_t=1/(alpha*LRs[LR]) 
    boosting_factors[delta, LR]=e_boosted_futility(m_t, nu/LRs[LR], deltas[delta])
  }
}

boosting_factors

