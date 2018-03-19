function [callMC_European_Price_multi_step, putMC_European_Price_multi_step] = ...
    MC_european_price(S0, K, T, r, mu, sigma, numSteps, numPaths)
%MC_EUROPEAN_PRICE Summary of this function goes here
%   Detailed explanation goes here

S = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);

callMC_European_Price_multi_step = max(0, S(end, :) - K);
putMC_European_Price_multi_step = max(0, K - S(end, :));

callMC_European_Price_multi_step = mean(callMC_European_Price_multi_step)*exp(-r*T);
putMC_European_Price_multi_step = mean(putMC_European_Price_multi_step)*exp(-r*T);

end

