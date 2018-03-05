function [x_optimal, cash_optimal, proportion_optimal] = strat_robust_optim(x_init, cash_init, mu, Q, cur_prices )
%STRAT_ROBUST_OPTIM Summary of this function goes here
%   Detailed explanation goes here
% Random data for 10 stocks

portfolio_value =  cur_prices * x_init + cash_init;

N = size(mu);
N = N(1);

w_minVar = cplexqp(Q, zeros(N,1) , [], [],  ones(1, N), [1;], zeros(N, 1), ones(N, 1));
ret_minVar = w_minVar' * mu; %r_rf / 6;


% Initial portfolio ("equally weighted" or "1/n")
w0 = ones(N,1) ./ N;
% Bounds on variables
lb_rMV = zeros(N,1); 
ub_rMV = inf*ones(N,1);
% Target portfolio return estimation error
var_matr = diag(diag(Q));
rob_init = w0' * var_matr * w0; % r.est.err. of 1/n portf
rob_bnd = rob_init; % target return estimation error
% Compute minimum variance portfolio (MVP)
% Target portfolio return = return of MVP
Portf_Retn = ret_minVar;



% Formulate and solve robust mean-variance problem
f_rMV = zeros(N,1); % objective function
% Constraints
A_rMV = sparse([ mu'; ones(1,N)]);
lhs_rMV = [Portf_Retn; 1]; rhs_rMV = [inf; 1];
% Create CPLEX model
cplex_rMV = Cplex('Robust_MV');
cplex_rMV.addCols(f_rMV, [], lb_rMV, ub_rMV);
cplex_rMV.addRows(lhs_rMV, A_rMV, rhs_rMV);
% Add quadratic objective
cplex_rMV.Model.Q = 2*Q;
% Add quadratic constraint on return estimation error (robustness constraint)
cplex_rMV.addQCs(zeros(size(f_rMV)), var_matr, 'L', rob_bnd, {'qc_robust'});
% Solve
cplex_rMV.DisplayFunc = [];
cplex_rMV.solve();

if cplex_rMV.Solution.status == 3
    proportion_optimal = x_init .* cur_prices' / portfolio_value;
    x_optimal = x_init;
else
    proportion_optimal = cplex_rMV.Solution.x;
    x_optimal = round(proportion_optimal * portfolio_value ./ (cur_prices' + 1e-9));
end


transaction_fees = 5e-3 * sum(abs(x_optimal - x_init)' * cur_prices');

cash_optimal = cash_init + ...
    (sum(cur_prices * x_init) - transaction_fees - cur_prices * x_optimal);

%if cash is neg, adjust shares by -1 until cash is pos
j = 0;
cash_decrement = 1e2;
while cash_optimal < 0
    j = j + 1;
    x_optimal = round(proportion_optimal * (portfolio_value -  cash_decrement * j) ...
        ./ (cur_prices' + 1e-2));
    transaction_fees = 5e-3 * sum(abs(x_optimal - x_init)' * cur_prices');
    cash_optimal = cash_init + ...
        (sum(cur_prices * x_init) - transaction_fees - cur_prices * x_optimal);
end

end

