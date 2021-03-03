function v = fTWV(M,n,T,kappa,theta,sigma,v0)
% Calculate the variance via the IPZ scheme of TW

D = T/n;

v = zeros(n + 1, M);
v(1,:) = v0.*ones(1,M);

emkt = exp(- kappa * D);
h = 4 * kappa * emkt / (sigma^2 * (1 - emkt));
d = 4 * kappa * theta / sigma^2;

N_u = 2^(15 + ceil(log2(n))) + 1; % this seems to big / is very slow.

v_min = 0.0001;
v_max = 8 * sigma;
  
grid = linspace(v_min, v_max, N_u);

% Calculate the quantiles (Algorithm 2) for algorithm 1 (IPZ)
p = gamcdf(h / (2 * emkt) * grid, d / 2);
u = linspace(0, 1, N_u);
q = zeros(1, N_u);
    
for j=1:N_u
      
    if u(j) < p(1)
        
        q(j) = 0;
        
    else
        
        [~,index] = min(abs(p - u(j)));
        q(j) = grid(index);
        
        %q(j) = grid(find(abs(p - u(j)) == min(abs(p - u(j))))); gives
        %the same result but is hard to read

    end
end

 for j=1:M
    
    for k=1:n
      
      if(v(k, j) == 0)
        
        U = unifrnd(0,1);
        [~,index] = min(abs(U-u));
        v(k + 1, j) = q(index);
        
      elseif(v(k, j) > 0)
        
        m_p = poissrnd(v(k, j) * h / 2);
        
        if(m_p == 0)
          
          U = unifrnd(0,1);
          [~,index] = min(abs(U-u));
          v(k + 1, j) = q(index);
          
        elseif(m_p > 0)
          
          v(k + 1, j) = 2 * emkt / h * gamrnd(m_p + d  / 2, 1);
          
        end 
        
      end
      
    end
    
 end
 
 