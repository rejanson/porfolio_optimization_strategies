function gval = computeGradERC (x)

global Q
  
  n = size(Q,1) ; 
  
%   std = sqrt(x * (Q*(x'))); 

  if(size(x,1)==1)
     x = x';
  end

  y = (Q*x) ; 
  
  gval = [] ; 
  for z = 1:n
      xij = 0;
      for i = 1:n        
        for j = i+1:n
          if i == z
            gradI = y(z);
            %gradI = gradtemp(z);
            gradJ = Q(j,z)*x(z);
          elseif j == z
            gradI = Q(i,z)*x(z);
            gradJ = y(z);
            %gradJ = gradtemp(z);
          else
            gradI = Q(i,z)*x(z);
            gradJ = Q(j,z)*x(z);
          end
          xij  = xij + ((y(i) - y(j)) * (gradI - gradJ));           
        end 
      end
      gval(z) = xij;
  end

  gval = 4*gval;


  % Insert your gradiant computations here
  % You can use finite differences to check the gradient

  % gval = zeros(n,1); 
  
  %  finite_diff
    grad = zeros(n, 1);
    diff = 1e-6;
    for i = 1:n
        x_minus = x;
        x_minus(i) = x_minus(i) - diff;
        x_plus = x;
        x_plus(i) = x_plus(i) + diff;
        result = (computeObjERC(x_plus) - computeObjERC(x_minus)) / 2 * diff;
        
        grad(i) = result;
    end
    gval = grad;
end
