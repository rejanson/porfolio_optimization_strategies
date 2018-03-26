clear all;
clc
format long;

Nout  = 100000; % number of out-of-sample scenarios
Nin   = 5000;   % number of in-sample scenarios
Ns    = 5;      % number of idiosyncratic scenarios for each systemic

C = 8;          % number of credit states

% Filename to save out-of-sample scenarios
filename_save_out  = 'scen_out';

% Read and parse instrument data
instr_data = dlmread('instrum_data.csv', ',');
instr_id   = instr_data(:,1);           % ID
driver     = instr_data(:,2);           % credit driver
beta       = instr_data(:,3);           % beta (sensitivity to credit driver)
recov_rate = instr_data(:,4);           % expected recovery rate
value      = instr_data(:,5);           % value
prob       = instr_data(:,6:6+C-1);     % credit-state migration probabilities (default to A)
exposure   = instr_data(:,6+C:6+2*C-1); % credit-state migration exposures (default to A)
retn       = instr_data(:,6+2*C);       % market returns

K = size(instr_data, 1); % number of  counterparties

% Read matrix of correlations for credit drivers
rho = dlmread('credit_driver_corr.csv', '\t');
sqrt_rho = (chol(rho))'; % Cholesky decomp of rho (for generating correlated Normal random numbers)

disp('======= Credit Risk Model with Credit-State Migrations =======')
disp('============== Monte Carlo Scenario Generation ===============')
disp(' ')
disp(' ')
disp([' Number of out-of-sample Monte Carlo scenarios = ' int2str(Nout)])
disp([' Number of in-sample Monte Carlo scenarios = ' int2str(Nin)])
disp([' Number of counterparties = ' int2str(K)])
disp(' ')

% Find credit-state for each counterparty
% 8 = AAA, 7 = AA, 6 = A, 5 = BBB, 4 = BB, 3 = B, 2 = CCC, 1 = default
[Ltemp, CS] = max(prob, [], 2);
clear Ltemp

% Account for default recoveries
exposure(:, 1) = (1-recov_rate) .* exposure(:, 1);

% Compute credit-state boundaries
CS_Bdry = norminv( cumsum(prob(:,1:C-1), 2) );


if(~exist('scenarios_out.mat','file'))
    
    % TODO: -------- Insert your code here -------- %
    Losses_out = zeros(Nout, K);
    
    
    parfor s = 1:Nout
        if mod(s, 1000) == 0
           s 
        end
        loss = [];
        zscores = mvnrnd(zeros(1,50) , rho);
        
        for i = 1:size(beta, 1)
            sigma_i = sqrt(1 - beta(i)^2);
            beta_i = beta(i);
            z_stat = beta_i * zscores(driver(i)) + sigma_i * mvnrnd(0, 1);
            
            for credit_state_i = 1:size(CS_Bdry, 2)
                current_boundary = CS_Bdry(i, credit_state_i);
                if z_stat < current_boundary
                    loss = [loss exposure(i, credit_state_i)];
                    break;
                end
                
                if credit_state_i >= 7
                    loss = [loss exposure(i, credit_state_i)];
                end
            end
            
        end

        Losses_out(s, :) = loss;

        % TODO: -------- Insert your code here -------- %
    end

    % Calculated out-of-sample losses (100000 x 100)
    % Losses_out

    save('scenarios_out', 'Losses_out')
else
    load('scenarios_out', 'Losses_out')
end

%% Normal approximation computed from out-of-sample scenarios
mu_l = mean(Losses_out)';
var_l = cov(Losses_out);

% Compute portfolio weights
portf_v = sum(value);     % portfolio value
w0{1} = value / portf_v;  % asset weights (portfolio 1)
w0{2} = ones(K, 1) / K;   % asset weights (portfolio 2)
x0{1} = (portf_v ./ value) .* w0{1};  % asset units (portfolio 1)
x0{2} = (portf_v ./ value) .* w0{2};  % asset units (portfolio 2)

% Quantile levels (99%, 99.9%)
alphas = [0.99 0.999];

VaRout = zeros(2, length(alphas));
VaRinN = zeros(2, length(alphas));
CVaRout = zeros(2, length(alphas));
CVaRinN = zeros(2, length(alphas));

% Compute VaR and CVaR (non-Normal and Normal) for 100000 scenarios
for(portN = 1:2)
    losses = Losses_out * cell2mat(x0(portN));
    losses = sort(losses);
    for(q=1:length(alphas))
        alf = alphas(q);
        
        % -------- Insert your code here -------- %
        VaRout(portN,q)  = losses(ceil((Nout-1) * alf));
        VaRinN(portN,q)  = mu_l' * cell2mat(x0(portN)) + norminv(alf,0,1)*std(losses);
        CVaRout(portN,q) = (1/((Nout-9)*(1-alf))) *  ...
                ((ceil((Nout-9)*alf)-(Nout-9)*alf) * VaRout(portN,q) + ...
                sum(losses(ceil((Nout-9)*alf)+1:Nout-9)));   
        CVaRinN(portN,q) = mu_l' * cell2mat(x0(portN)) + ...
            (normpdf(norminv(alf,0,1))/(1-alf))*std(losses);
        % -------- Insert your code here -------- %        
 end
end


%% Perform 100 trials
N_trials = 10;

VaRinMC1 = {};
VaRinMC2 = {};
VaRinN1 = {};
VaRinN2 = {};
CVaRinMC1 = {};
CVaRinMC2 = {};
CVaRinN1 = {};
CVaRinN2 = {};

Cumulative_MC1 = [];
Cumulative_MC2 = [];

for(tr=1:N_trials)
    % Monte Carlo approximation 1

    % -------- Insert your code here -------- %
    Losses_inMC1 = zeros(Nin, K);
    
    for s = 1:ceil(Nin/Ns) % systemic scenarios
        % -------- Insert your code here -------- %
        
        zscores = mvnrnd(zeros(1,50) , rho);

        for si = 1:Ns % idiosyncratic scenarios for each systemic
            % -------- Insert your code here -------- %
            loss = [];
            
            for i = 1:size(beta, 1)
                sigma_i = sqrt(1 - beta(i)^2);
                beta_i = beta(i);
                z_stat = beta_i * zscores(driver(i)) + sigma_i * mvnrnd(0, 1);

                for credit_state_i = 1:size(CS_Bdry, 2)
                    current_boundary = CS_Bdry(i, credit_state_i);
                    if z_stat < current_boundary
                        loss = [loss exposure(i, credit_state_i)];
                        break;
                    end

                    if credit_state_i >= 7
                        loss = [loss exposure(i, credit_state_i)];
                    end
                end

            end

            
            Losses_inMC1((s - 1) * si + si, :) = loss;
        end
        

    end
    
    % Calculated losses for MC1 approximation (5000 x 100)
    % Losses_inMC1
    
    % Monte Carlo approximation 2
    
    % -------- Insert your code here -------- %
    Losses_inMC2 = zeros(Nin, K);
    
    parfor s = 1:Nin % systemic scenarios (1 idiosyncratic scenario for each systemic)
        % -------- Insert your code here -------- %
        zscores = mvnrnd(zeros(1,50) , rho);
        
        loss = [];
            
        for i = 1:size(beta, 1)
            sigma_i = sqrt(1 - beta(i)^2);
            beta_i = beta(i);
            z_stat = beta_i * zscores(driver(i)) + sigma_i * mvnrnd(0, 1);

            for credit_state_i = 1:size(CS_Bdry, 2)
                current_boundary = CS_Bdry(i, credit_state_i);
                if z_stat < current_boundary
                    loss = [loss exposure(i, credit_state_i)];
                    break;
                end

                if credit_state_i >= 7
                    loss = [loss exposure(i, credit_state_i)];
                end
            end

        end


        Losses_inMC2(s, :) = loss;

    end
        
    % Calculated losses for MC2 approximation (5000 x 100)
    % Losses_inMC2
    
    % Compute VaR and CVaR
    for(portN = 1:2)
        for(q=1:length(alphas))
            alf = alphas(q);
            % -------- Insert your code here -------- %            
            % Compute portfolio loss 
            portf_loss_inMC1 = sort(Losses_inMC1 * cell2mat(x0(portN)));
            portf_loss_inMC2 = sort(Losses_inMC2 * cell2mat(x0(portN)));
            mu_MCl = mean(Losses_inMC1)';
            var_MCl = cov(Losses_inMC1);
            mu_MC2 = mean(Losses_inMC2)';
            var_MC2 = cov(Losses_inMC2);
            % Compute portfolio mean loss mu_p_MC1 and portfolio standard deviation of losses sigma_p_MC1
            % Compute portfolio mean loss mu_p_MC2 and portfolio standard deviation of losses sigma_p_MC2
            % Compute VaR and CVaR for the current trial
            VaRinMC1{portN,q}(tr) = portf_loss_inMC1(ceil((Nin-1) * alf));
            VaRinMC2{portN,q}(tr) = portf_loss_inMC2(ceil((Nin-1) * alf));
            VaRinN1{portN,q}(tr) = mu_MCl' * cell2mat(x0(portN)) ...
                    + norminv(alf,0,1)*std(portf_loss_inMC1);
            VaRinN2{portN,q}(tr) =  mu_MC2' * cell2mat(x0(portN)) ...
                    + norminv(alf,0,1)*std(portf_loss_inMC2);
            CVaRinMC1{portN,q}(tr) = (1/((Nin-9)*(1-alf))) *  ...
                ((ceil((Nin-9)*alf)-(Nin-9)*alf) * VaRinMC1{portN,q}(tr) + ...
                sum(portf_loss_inMC1(ceil((Nin-9)*alf)+1:Nin-9)));  
            CVaRinMC2{portN,q}(tr) = (1/((Nin-9)*(1-alf))) *  ...
                ((ceil((Nin-9)*alf)-(Nin-9)*alf) * VaRinMC2{portN,q}(tr) + ...
                sum(portf_loss_inMC2(ceil((Nin-9)*alf)+1:Nin-9)));  
            CVaRinN1{portN,q}(tr) = mu_MCl' * cell2mat(x0(portN)) + ...
                (normpdf(norminv(alf,0,1))/(1-alf))*std(portf_loss_inMC1);
            CVaRinN2{portN,q}(tr) = mu_MC2' * cell2mat(x0(portN)) + ...
                (normpdf(norminv(alf,0,1))/(1-alf))*std(portf_loss_inMC2);
            % -------- Insert your code here -------- %
        end
    end
    
    Cumulative_MC1 = [Cumulative_MC1; portf_loss_inMC1];
    Cumulative_MC2 = [Cumulative_MC2; portf_loss_inMC2];
end



%% Display portfolio VaR and CVaR
for(portN = 1:2)
fprintf('\nPortfolio %d:\n\n', portN)    
 for(q=1:length(alphas))
    alf = alphas(q);
    fprintf('Out-of-sample: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRout(portN,q), 100*alf, CVaRout(portN,q))
    fprintf('In-sample MC1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC1{portN,q}), 100*alf, mean(CVaRinMC1{portN,q}))
    fprintf('In-sample MC2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC2{portN,q}), 100*alf, mean(CVaRinMC2{portN,q}))
    fprintf(' In-sample No: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRinN(portN,q), 100*alf, CVaRinN(portN,q))
    fprintf(' In-sample N1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinN1{portN,q}), 100*alf, mean(CVaRinN1{portN,q}))
    fprintf(' In-sample N2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n\n', 100*alf, mean(VaRinN2{portN,q}), 100*alf, mean(CVaRinN2{portN,q}))
 end
end

% Plot results
% figure(1);
% -------- Insert your code here -------- %
% figure(2);
% -------- Insert your code here -------- %