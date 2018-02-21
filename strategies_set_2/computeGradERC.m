function gval = computeGradERC (x)

global Q
  
  n = size(Q,1) ;  

  if(size(x,1)==1)
     x = x';
  end
  
  % Insert your gradiant computations here
  % You can use finite differences to check the gradient

  gval = zeros(n,1); 
 
end
