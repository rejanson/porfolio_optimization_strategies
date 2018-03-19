function [call_BS_European_Price, putBS_European_Price] = ...
    BS_european_price(S0, K, T, r, sigma)
%BS_EUROPEAN_PRICE Summary of this function goes here
%   Detailed explanation goes here

% C(S, t) = N (d1)S ? N (d2)Ke^(?r(T ?t))
% P(S, t) = N (?d2)Ke^(?r(T ?t)) ? N (?d1)S,
d1 = log(S0 / K) + (r + sigma^2 / 2) * T;
d1 = d1 / (sigma * sqrt(T));
d2 = d1 - sigma * sqrt(T);

call_BS_European_Price = normcdf(d1) * S0 - normcdf(d2) * K * exp(-r * T);

putBS_European_Price = normcdf(-d2) * K * exp(-r * T) - normcdf(-d1) * S0;

end

