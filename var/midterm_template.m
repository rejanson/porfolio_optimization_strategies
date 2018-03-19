clc;
clear all;
format long

addpath('/Applications/CPLEX_Studio128/cplex/matlab/x86-64_osx');


% CSV file with price data
input_file_prices  = 'Daily_closing_prices.csv';

% Read daily prices
if(exist(input_file_prices,'file'))
  fprintf('\nReading daily prices datafile - %s\n', input_file_prices)
  fid = fopen(input_file_prices);
     % Read instrument tickers
     hheader  = textscan(fid, '%s', 1, 'delimiter', '\n');
     headers = textscan(char(hheader{:}), '%q', 'delimiter', ',');
     tickers = headers{1}(2:end);
     % Read time periods
     vheader = textscan(fid, '%[^,]%*[^\n]');
     dates = vheader{1}(1:end);
  fclose(fid);
  data_prices = dlmread(input_file_prices, ',', 1, 1);
else
  error('Daily prices datafile does not exist')
end

% Convert dates into array [year month day]
format_date = 'mm/dd/yyyy';
dates_array = datevec(dates, format_date);
dates_array = dates_array(:,1:3);

% Remove datapoints for year 2014
day_ind_start0 = 1;
day_ind_end0 = length(find(dates_array(:,1)==2014));
data_prices = data_prices(day_ind_end0+1:end,:);
dates_array = dates_array(day_ind_end0+1:end,:);
dates = dates(day_ind_end0+1:end,:);

% Compute means and covariances for Question 2
day_ind_start = 1;
day_ind_end = 39;
cur_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
mu = mean(cur_returns)';  % Expected returns for Question 2
Q = cov(cur_returns);     % Covariances for Question 2

% Number of assets in universe
Na = size(data_prices,2);

% Number of historical scenarios
Ns = size(data_prices,1);

%% Question 1
fprintf('Question 1 Part 1 \n');

% Specify quantile level for VaR/CVaR
alf = 0.95;

% Positions in the portfolio
positions = [100 0 0 0 0 0 0 0 200 500 0 0 0 0 0 0 0 0 0 0]';


%%%%% Insert your code here 

%Historical
loss1 = (data_prices(2:end,:) - data_prices(1:end - 1,:)) * positions;
loss10 = (data_prices(10:end,:) - data_prices(1:end - 9,:)) * positions;
loss1sorted = sort(-loss1)';
loss10sorted = sort(-loss10)';

VaR1 = loss1sorted(ceil((Ns-1) * alf));
VaR10 = loss10sorted(ceil((Ns-1) * alf));
CVaR1 = (1/((Ns-1)*(1-alf))) * ...
    ((ceil((Ns-1)*alf)-(Ns-1)*alf) * VaR1 + sum(loss1sorted(ceil((Ns-1)*alf)+1:Ns-1)));
CVaR10 = (1/((Ns-9)*(1-alf))) *  ...
    ((ceil((Ns-9)*alf)-(Ns-9)*alf) * VaR10 + sum(loss10sorted(ceil((Ns-9)*alf)+1:Ns-9)));   


% Normal

% Compute Normal 1-day VaR from the data
VaR1n = mean(loss1sorted) + norminv(alf,0,1)*std(loss1sorted);
% Compute Normal 1-day CVaR from the data
CVaR1n = mean(loss1sorted) + (normpdf(norminv(alf,0,1))/(1-alf))*std(loss1sorted);

% Compute Normal 10-day VaR from the data
VaR10n = mean(loss10sorted) + norminv(alf,0,1)*std(loss10sorted);
% Compute Normal 10-day CVaR from the data
CVaR10n = mean(loss10sorted) + (normpdf(norminv(1 - alf,0,1))/(1-alf))*std(loss10sorted);

fprintf('Historical 1-day VaR %4.1f%% = $%6.2f,   Historical 1-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR1, 100*alf, CVaR1)
fprintf('    Normal 1-day VaR %4.1f%% = $%6.2f,       Normal 1-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR1n, 100*alf, CVaR1n)
fprintf('Historical 10-day VaR %4.1f%% = $%6.2f,   Historical 10-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR10, 100*alf, CVaR10)
fprintf('    Normal 10-day VaR %4.1f%% = $%6.2f,       Normal 10-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR10n, 100*alf, CVaR10n)


% Plot a histogram of the distribution of losses in portfolio value for 1 day 
figure(1)
[frequencyCounts, binLocations] = hist(loss1 * -1, 100);
hold on
bar(binLocations, frequencyCounts);
line([VaR1 VaR1], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
normfreq = ( 1/(std(loss1sorted)*sqrt(2*pi)) )  ...
    * exp( -0.5*((binLocations-mean(loss1sorted))/std(loss1sorted)).^2 );
normfreq = normfreq * sum(frequencyCounts)/sum(normfreq);
plot(binLocations, normfreq, 'r', 'LineWidth', 3);
line([VaR1n VaR1n], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.');
hold off;

text(1.05*VaR1, max(frequencyCounts)/1.95, 'VaR1')
text(0.7*VaR1n, max(frequencyCounts)/1.95, 'VaR1n')
title('Frequency of Losses for a 1-day window')
xlabel('1-day loss in $ value on 1 unit of stock')
ylabel('Frequency')


% Plot a histogram of the distribution of losses in portfolio value for 10 days
figure(2)
[frequencyCounts, binLocations] = hist(loss10 * -1, 100);
hold on 
bar(binLocations, frequencyCounts);
line([VaR10 VaR10], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
normfreq = ( 1/(std(loss10sorted)*sqrt(2*pi)) )  ...
    * exp( -0.5*((binLocations-mean(loss10sorted))/std(loss10sorted)).^2 );
normfreq = normfreq * sum(frequencyCounts)/sum(normfreq);
plot(binLocations, normfreq, 'r', 'LineWidth', 3);
line([VaR10n VaR10n], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.');
hold off;

text(1.05 * VaR10, max(frequencyCounts)/1.9, 'VaR10')
text(0.7 * VaR10n, max(frequencyCounts)/1.9, 'VaR10n')
title('Frequency of Losses for a 10-day window')
xlabel('1-day loss in $ value on 1 unit of stock')
ylabel('Frequency')



%% Question 1 Part 2
% 1-day 95% VaR for 
alf = 0.95;
position2 = [100 zeros(1,19)];
position3 = [zeros(1,8) 200 zeros(1,11)];
position4 = [zeros(1,9) 500 zeros(1,10)];


position_sequence = [position2; position3; position4];
fprintf('Question 1 Part 2 \n');
for h = 1:size(position_sequence, 1)
    pos = position_sequence(h, :);
    loss1_p2 = (data_prices(2:end,:) - data_prices(1:end - 1,:)) * pos' ;

    %historical 
    loss1sorted_p2 = sort(-loss1_p2)';

    VaR1_p2 = loss1sorted_p2(ceil((Ns-1) * alf));
    fprintf('Historical 1-day VaR %4.1f%% = $%6.2f ', 100*alf, VaR1_p2);
    
    %normal 
    VaR1n_p2 = mean(loss1sorted_p2) + norminv(alf,0,1)*std(loss1sorted_p2);
    fprintf('Normal 1-day VaR %4.1f%% = $%6.2f \n', 100*alf, VaR1n_p2)


end




%% Question 2
% Annual risk-free rate for years 2015-2016 is 2.5%
r_rf = 0.025;

% Initial portfolio weights
init_positions = [5000 950 2000 0 0 0 0 2000 3000 1500 0 0 0 0 0 0 1001 0 0 0]';
init_value = data_prices(day_ind_end+1,:) * init_positions;

w_init = (data_prices(day_ind_end+1,:) .* init_positions')' / init_value;

% Max Sharpe Ratio portfolio weights
w_Sharpe = [ 0 0 0 0 0 0 0 0.385948690661642 0.172970428625544 0 0 0 0 0 0.003409676869715 0.260942060896445 0 0.185966939781285 0 0]';

% Equal Risk Contribution portfolio weights
w_ERC = [0.049946771209069 0.049951626261681 0.049955739901370 0.049998404150207 0.050000297368719 0.050004255546315 0.050006307026730 0.050007308995726 0.050010525832832 0.050013840015521 0.050014404492514 0.050015932843104 0.050016630302524 0.050017212457105 0.050017600497611 0.050017998351827 0.050018997074443 0.050019598350121 0.050019778113513 0.049946771209069]';

w_leveraged_ERC = 2 * w_ERC;

w_minVar = cplexqp(Q, zeros(Na,1) , [], [],  ones(1, Na), [1;], zeros(Na, 1), ones(Na, 1));

w_maxRet = cplexqp(zeros(Na,Na), -1 * mu , [], [],  ones(1, Na), [1;], zeros(Na, 1), ones(Na, 1));

w_equallyWeighted = 1 / Na * ones(Na, 1);

%%%%% Insert your code here 
% % ? Efficient frontier of risky assets under no-short-sales constraint;
% % ? Minimum variance portfolio of risky assets; 
% %? Maximum return portfolio of risky assets;
% %? Equally-weighted (1/N) portfolio of risky assets;
% %? Risk-free asset;
% % ? Leveraged equal risk contribution portfolio;
% ? Efficient frontier of all assets including risk-free asset, if shorting of risk-free asset is
efficient_frontier = [];
rets = [0:0.00005:0.0089];
for targetReturn = rets
    [x,fval] = cplexqp(Q, zeros(Na,1) , [-1 * mu'], [-1 * targetReturn],  [ones(1, Na)], [1;], ...
        zeros(Na, 1), ones(Na, 1));
    efficient_frontier = [efficient_frontier x'*Q*x];
end

% Plot for Question 2, Part 1
figure(3);
plot(sqrt(efficient_frontier), rets);
hold on;
text(sqrt(w_init' * Q * w_init) , mu' * w_init, '\bullet Init');
text(sqrt(w_Sharpe' * Q * w_Sharpe) , mu' * w_Sharpe, '\leftarrow Max Sharpe');
text(sqrt(w_ERC' * Q * w_ERC) , mu' * w_ERC, '\bullet ERC');
text(sqrt(w_leveraged_ERC' * Q * w_leveraged_ERC) , mu' * w_leveraged_ERC, '\bullet Leveraged ERC');
text(sqrt(w_minVar' * Q * w_minVar) , mu' * w_minVar, '\bullet minVar');
text(sqrt(w_maxRet' * Q * w_maxRet) , mu' * w_maxRet, '\bullet maxRet');
text(sqrt(w_equallyWeighted' * Q * w_equallyWeighted) , mu' * w_equallyWeighted, '\bullet ERC');
text(0, r_rf/252, '\bullet Risk Free');
axis([-1e-3 4e-2 -1e-3 1e-2])
%TODO: tangent 
slope = 0.445;
x = -1e-3:1e-4:2.5e-2;
plot(x, slope * x + r_rf/252)

legend('Efficient Frontier', 'Capital Market Line');

title('Efficient Frontier and Capital Market Line')
xlabel('STDEV')
ylabel('Expected Return')

hold off

% Plot for Question 2, Part 2
figure(4);
plot(sqrt(efficient_frontier), rets);
hold on
scatter(sqrt(diag(Q)), mu);


num_port = 1000;
var_random = [];
ret_random = [];
for i = 1:1000
    w_rand = rand(20,1);
    w_rand = w_rand / sum(w_rand);
    %scatter(sqrt(w_rand' * Q * w_rand), mu' * w_rand);
    var_random = [var_random sqrt(w_rand' * Q * w_rand)];
    ret_random = [ret_random mu' * w_rand];
end

scatter(var_random, ret_random);

legend('Efficient Frontier', 'Individual Stocks', 'Random Portfolios');
title('Efficient Frontier vs Individual Stocks vs Random Portfolio')
xlabel('STDEV')
ylabel('Expected Return')


