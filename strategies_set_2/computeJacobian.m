function J = computeJacobian (x)

global A_eq A_ineq;
  
  J = sparse([ A_eq; A_ineq ]);
