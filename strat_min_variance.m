function [ x_optimal, cash_optimal ] = strat_min_variance(  x_init, cash_init, mu, Q, cur_prices )
%STRAT_MIN_VARIANCE Summary of this function goes here
%   Detailed explanation goes here

portfolio_value =  cur_prices * x_init + cash_init;

N = size(mu);
N = N(1);

Aeq = ones(1, N);
beq = 1;
lb = zeros(N, 1);
ub = [];
ret = zeros(N, 1);

proportion_optimal = quadprog(Q, ret, [], [], Aeq, beq, lb, ub);

x_optimal = round(proportion_optimal * portfolio_value ./ (cur_prices' + 1e-2));

transaction_fees = 5e-3 * sum(abs(x_optimal - x_init)' * cur_prices');

cash_optimal = cash_init + ...
    (sum(cur_prices * x_init) - transaction_fees - cur_prices * x_optimal);

%if cash is neg, adjust shares by -1 until cash is pos
while cash_optimal < 0
    for i = 1:N
        if x_optimal(i) > 0
            x_optimal(i) = x_optimal(i) - 1;
        end
    end
    transaction_fees = 5e-3 * sum(abs(x_optimal - x_init)' * cur_prices');
    cash_optimal = cash_init + ...
        (sum(cur_prices * x_init) - transaction_fees - cur_prices * x_optimal);
end


end

