# compute Realized Skewness from non-centered moments   

getKuFromMoments = function(m1, m2, m3, m4){ 
  
  Ku = (m4 - 4*m3*m1 + 6*m2*m1^2 - 3*m1^4 ) / (m2-m1^2)^2

}