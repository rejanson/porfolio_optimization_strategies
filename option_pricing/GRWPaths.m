function paths = GRWPaths(initPrice, mu, sigma, T, numSteps, numPaths)
    % Computes numPaths random paths for a geometric random walk
    % mu is the annual drift, sigma the annual volatility
    % T is the total length of time for the path (in years)
    % dT is the time increment (in years)
       
    paths = zeros(numSteps+1, numPaths);
    dT = T/numSteps;
    
    % Vector of paths will store realizations of the asset price
    % First asset price is the initial price
    paths(1,:) = initPrice;
 
    % Generate paths
    for iPath = 1:numPaths
        for iStep = 1:numSteps
            paths(iStep+1, iPath) = paths(iStep, iPath) * exp( (mu - 0.5*sigma^2)*dT + sigma*sqrt(dT)*normrnd(0,1) );
        end
    end 
            
    % Plot paths
%     figure;
%     set(gcf, 'color', 'white');
%     plot(0:numSteps, paths', 'Linewidth', 2);
%     title('Geometric Random Walk Paths', 'FontWeight', 'bold');

end