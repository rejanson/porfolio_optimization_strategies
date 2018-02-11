function  [x_optimal, cash_optimal, proportion_optimal] = strat_buy_and_hold(x_init, cash_init, mu, Q, cur_prices)
   x_optimal = x_init;
   cash_optimal = cash_init;
   proportion_optimal = x_optimal / sum(x_optimal);
end