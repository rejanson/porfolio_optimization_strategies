function [ x_optimal, cash_optimal, proportion_optimal ] = strat_equally_weighted( x_init, cash_init, mu, Q, cur_prices )
%STRAT_EQUALLY_WEIGHTED Summary of this function goes here
% “Equally weighted” (also known as “1/n”) portfolio strategy: asset weights are 
% selected as wti = 1/n, where n is the number of assets. 
% You may need to re-balance your posrtfolio in each period as the number of 
% shares xti changes even when wti = 1/n stays the same in each period. 
% The strategy should be implemented in the function strat_equally_weighted.

%x_init - curr_position 20x1, number of shares
%cash_init - curr_cash float
%mu - 20x1
%Q - 20x2
%curr_prices - 1x20

N = size(mu);
N = N(1);

x_optimal = zeros(N, 1);

market_cap = x_init .* cur_prices';
cap_for_each = (sum(market_cap) + cash_init) / 20;

%compute shares based on market cap for each stock 
for i = 1:N
    x_optimal(i) = round(cap_for_each / cur_prices(i));
end

transaction_fees = 5e-3 * sum(abs(x_optimal - x_init)' * cur_prices');

cash_optimal = cash_init + ...
    (sum(market_cap) - transaction_fees - cur_prices * x_optimal);

%if cash is neg, adjust capitalization for each stock by -$1000 until cash is pos
j = 0;
cash_decrement = 1e3;
while cash_optimal < 0
    
    j = j + 1;
    market_cap = x_init .* cur_prices';
    cap_for_each = (sum(market_cap) + cash_init - j * cash_decrement) / 20;

    for i = 1:N
        x_optimal(i) = round(cap_for_each / cur_prices(i));
    end
    
    transaction_fees = 5e-3 * sum(abs(x_optimal - x_init)' * cur_prices');
    cash_optimal = cash_init + ...
        (sum(market_cap) - transaction_fees - cur_prices * x_optimal);
end

proportion_optimal = ones(N, 1) * 1/20 ;

end

