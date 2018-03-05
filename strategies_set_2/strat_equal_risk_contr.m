function [ x_optimal, cash_optimal, proportion_optimal] = strat_equal_risk_contr(x_init, cash_init, mu, Q, cur_prices )
%STRAT_EQUAL_RISK_CONTR Summary of this function goes here
%   Detailed explanation goes here

global A_ineq A_eq

portfolio_value =  cur_prices * x_init + cash_init;

N = size(mu);
N = N(1);

A_eq = ones(1,N);
b_eq = 1;

% Inequality constraints
A_ineq = [];
b_ineql = [];
b_inequ = [];

%w0 = repmat(1.0/N, N, 1);
w0 = cur_prices' .* x_init / (cur_prices * x_init + cash_init);

options.lb = zeros(1,N);       % lower bounds on variables
options.lu = ones (1,N);       % upper bounds on variables
options.cl = [b_eq' b_ineql']; % lower bounds on constraints
options.cu = [b_eq' b_inequ']; % upper bounds on constraints

% Set the IPOPT options
options.ipopt.jac_c_constant        = 'yes';
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.tol                   = 1e-10;
options.ipopt.print_level = 1;

% The callback functions
funcs.objective         = @computeObjERC;
funcs.constraints       = @computeConstraints;
funcs.gradient          = @computeGradERC;
funcs.jacobian          = @computeJacobian;
funcs.jacobianstructure = @computeJacobian;

%ipopt 
[proportion_optimal info] = ipopt(w0',funcs,options);
proportion_optimal = proportion_optimal';
x_optimal = round(proportion_optimal * portfolio_value ./ (cur_prices' + 1e-2));

%transaction cost 

transaction_fees = 5e-3 * sum(abs(x_optimal - x_init)' * cur_prices');

cash_optimal = cash_init + ...
    (sum(cur_prices * x_init) - transaction_fees - cur_prices * x_optimal);

%adjust transaction cost
j = 0;
cash_decrement = 1e2;
while cash_optimal < 0
%     for i = 1:N
%         if x_optimal(i) > 0
%             x_optimal(i) = x_optimal(i) - 1;
%         endﬂ
%     end
    j = j + 1;
    x_optimal = round(proportion_optimal * (portfolio_value -  cash_decrement * j) ...
        ./ (cur_prices' + 1e-2));
    transaction_fees = 5e-3 * sum(abs(x_optimal - x_init)' * cur_prices');
    cash_optimal = cash_init + ...
        (sum(cur_prices * x_init) - transaction_fees - cur_prices * x_optimal);
end

end

