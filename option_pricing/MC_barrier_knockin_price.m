function [callMC_Barrier_Knockin_Price_multi_step, ...
    putMC_Barrier_Knockin_Price_multi_step] ...
    = MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

S = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);

callMC_European_Price = max(0, S(end, :) - K);
putMC_European_Price = max(0, K - S(end, :));

l = zeros(1, numPaths);
for i = 1:numPaths
    if max(S(:, i)) <= Sb
       l(i) = 0;
    else
        l(i) = 1;
    end
end

callMC_Barrier_Knockin_Price_multi_step = mean( l .* callMC_European_Price)*exp(-r*T);
putMC_Barrier_Knockin_Price_multi_step = mean( l .* putMC_European_Price)*exp(-r*T);

end

