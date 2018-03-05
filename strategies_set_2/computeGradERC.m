function gval = computeGradERC (x)

global Q
  
  n = size(Q,1) ; 
  
  if(size(x,1)==1)
     x = x';
  end

  Qw = Q * x;
  y = x .* (Q*x) ; 
  
  gval = [] ; 
  for h = 1:n
      cumulative_grad = 0;
      
      for i = 1:n        
          
        for j = i + 1:n
          
          if i == h
              gradI = Qw(i) + Q(i, i) * x(i);
          else 
              gradI = Q(i, h) * x(i);
          end
              
          if j == h
              gradJ = Qw(j) + Q(j, j) * x(j);
          else 
              gradJ = Q(j, h) * x(j);
          end
          
          cumulative_grad  = cumulative_grad + ((y(i) - y(j)) * (gradI - gradJ));           
        end 
      end
      gval(h) = cumulative_grad;
  end

  gval = 4*gval;

  %  finite_diff
    grad = zeros(n, 1);
    interval = 1e-6;
    
    for i = 1:n
        x_minus = x;
        x_minus(i) = x_minus(i) - interval;
        x_plus = x;
        x_plus(i) = x_plus(i) + interval;
        result = (computeObjERC(x_plus) - computeObjERC(x_minus)) / (2 * interval);
        
        grad(i) = result;
    end
    
    if (sum(abs(gval - grad')) / sum(abs(grad))) > 1e-2
        disp('gradients mismatched');
    end
    
        
    %gval = grad;
end
