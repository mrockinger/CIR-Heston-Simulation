# Scheme of Tse/Wan (2013)
# Note: Note vectorized, yet


IPZScheme = function(v0, kappa, theta, sigma, T, n, M){
  
  v = matrix(, M, n + 1)
  v[,1] = v0

  Delta = T / n
  
  emkt = exp(- kappa * Delta)
  h = 4 * kappa * emkt / (sigma^2 * (1 - emkt)) 
  d = 4 * kappa * theta / sigma^2
  
  N_u = 2^(15 + ceiling(log2(n))) + 1 # this seems to big / is very slow.
  u = seq(0, 1, len = N_u)
  
  
  # Precompute quantiles of v_t|v_s via nearest-neighbor interpolation (Algorithm 2) if necessary. 
  # These quantiles will only be used, if the Poisson variate m_p is 0
  # This will take a while because of the for loop! p is generated quite fast even for the original N_u value. Can this be vectorized?
  
  v_min = 0.0001
  v_max = 8 * sigma
  
  v.v = seq(v_min, v_max, len = N_u)
  
  # Here, I use shape = d / 2 since we will only use this approximation if the Possion variate in eq. (21) is 0
  
  p = pgamma(h / (2 * emkt) * v.v, shape = d / 2) 
  q = rep(0, N_u)
  
  for(i in 1:N_u){
    
    #print(i)
    
    if(u[i] < p[1]){
      
      q[i] = 0
      
    }else{
      
      q[i] = v.v[which.min(abs(p-u[i]))] # which.min is not a vectorized function
      
    }
    
  }
  
  
  # Now, simulate v_t | v_s with Algorithm 1
  
  
  for(j in 1:M){
    
    for(i in 1:n){
      
      if(v[j, i] == 0){
        
        # print("IPZ")
        
        U = runif(1)
        v[j, i + 1] = q[which.min(abs(U-u))]
        
      }else if(v[j, i] > 0){
        
        m_p = rpois(1, v[j, i] * h / 2)
        
        if(m_p == 0){
          
          # print("IPZ")
          
          U = runif(1)
          v[j, i + 1] = q[which.min(abs(U-u))]
          
        }else if(m_p > 0){
          
          v[j, i + 1] = 2 * emkt / h * rgamma(1, shape = m_p + d  / 2) 
          
        } 
        
      }
      
    }
    
  }
  
  return(v)
  
}

# system.time({test = IPZScheme(0.1, 3, 0.19, 0.4, 1, 1/12)})
# MomentsCIR(p = 1:4, kappa = 3, theta = 0.19, sigma = 0.4, v0 = 0.1, t = 1 / 12, conditional = F)
# system.time({test1 = replicate(1000, TseWanCIR(0.1, 3, 0.19, 0.4, 1, 1 / 12))[2,]})
# MomentsCIR(p = 1:4, kappa = 3, theta = 0.19, sigma = 0.4, v0 = 0.1, t = 1 / 12)

