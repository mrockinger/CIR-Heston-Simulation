# compute Realized Skewness from non-centered moments   

getSkFromMoments = function(m1,m2,m3) { 
  
  Sk = (m3 - 3*m2*m1 + 2*m1^3) / (m2-m1^2)^1.5 
  
}
 
