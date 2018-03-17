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
gains1 = (data_prices(2:end,:) - data_prices(1:end - 1,:)) * positions;
gains10 = (data_prices(11:end,:) - data_prices(1:end - 10,:)) * positions;

VaR1 = prctile(gains1, 100 - 100 * alf);
VaR10 = prctile(gains10, 100 - 100 * alf);

CVaR1 = zeros(size(VaR1, 1), size(VaR1, 2));
CVaR10 = zeros(size(VaR10, 1), size(VaR10, 2));
for i = 1:size(gains1, 2)
    sum_ret1 = [];
    sum_ret10 = [];
    
    for j = 1:size(gains1, 1)
        if gains1(j, i) <= VaR1(1, i)
            sum_ret1 = [gains1(j, i) sum_ret1];
        end
    end
    CVaR1(1, i) = mean(sum_ret1);
    
    for j = 1:size(gains10, 1)
        if gains10(j, i) <= VaR10(1, i)
            sum_ret10 = [gains10(j, i) sum_ret10];
        end
    end
    CVaR10(1, i) = mean(sum_ret10);
end

% Normal

% Compute Normal 1-day VaR from the data
VaR1n = mean(gains1) + norminv(1 - alf,0,1)*std(gains1);
% Compute Normal 1-day CVaR from the data
CVaR1n = mean(gains1) + (normpdf(norminv(1 - alf,0,1))/(1 - alf))*std(gains1);
CVaR1n = CVaR1n * -1;

% Compute Normal 10-day VaR from the data
VaR10n = mean(gains10) + norminv(1 - alf,0,1)*std(gains10);
% Compute Normal 10-day CVaR from the data
CVaR10n = mean(gains10) + (normpdf(norminv(1 - alf,0,1))/(1-alf))*std(gains10);
CVaR10n = CVaR10n * -1;

fprintf('Historical 1-day VaR %4.1f%% = $%6.2f,   Historical 1-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR1, 100*alf, CVaR1)
fprintf('    Normal 1-day VaR %4.1f%% = $%6.2f,       Normal 1-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR1n, 100*alf, CVaR1n)
fprintf('Historical 10-day VaR %4.1f%% = $%6.2f,   Historical 10-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR10, 100*alf, CVaR10)
fprintf('    Normal 10-day VaR %4.1f%% = $%6.2f,       Normal 10-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR10n, 100*alf, CVaR10n)


% Plot a histogram of the distribution of losses in portfolio value for 1 day 
figure(1)
histogram(gains1, 'NumBins', 100);

% Plot a histogram of the distribution of losses in portfolio value for 10 days
figure(2)
histogram(gains10, 'NumBins', 100);


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
    gains1_p2 = (data_prices(2:end,:) - data_prices(1:end - 1,:)) * pos' ;

    %historical 
    VaR1 = prctile(gains1_p2, 100 - 100 * alf);
    fprintf('Historical 1-day VaR %4.1f%% = $%6.2f ', 100*alf, VaR1);
    
    %normal 
    VaR1n = mean(gains1_p2) + norminv(1 - alf,0,1)*std(gains1_p2);
    fprintf('Normal 1-day VaR %4.1f%% = $%6.2f \n', 100*alf, VaR1n)


end




%% Question 2
clc;
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

w_maxRet = cplexqp(zeros(Na,Na), mu , [], [],  ones(1, Na), [1;], zeros(Na, 1), ones(Na, 1));

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
rets = [0.0001:0.00005:0.02];
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
text(0, r_rf/6, '\bullet Risk Free');
axis([-1e-3 2.5e-2 -1e-3 1e-2])
%TODO: tangent 



% Plot for Question 2, Part 2
% figure(4);
