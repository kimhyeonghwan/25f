function [realizedUtilityA,realizedUtilityB,realizedUtility2A,realizedUtility2B,phat] = welfareOnLambda1(params,lambdaGrid,NSim,nSteps,nt)
% welfareOnLambda1 - calculates realized utilities
% the function is embedded into dependenceOnLambda3 but very slow with
% required number of simulations of 1e6-1e7

% Version 3: 12-May-2017
% Calculates realized utility multiplied by Radon-Nikodym derivative (for
% welfare analysis) using function
% Uses the same sequence of shocks for all labmdas (which reduces lumpiness across lambda)
% Can be computed in parallel

% nt - number of steps used to solve for phat


nGrid = length(lambdaGrid);
phat = nan(nGrid,1);
N = NSim; % number of paths to simulate
n = nSteps;  % number of points in discretization scheme 

realizedUtilityA1 = nan(N,nGrid); % realized utility of agent A in complete markets
realizedUtilityA2 = nan(N,nGrid); % realized utility of agent A with circuit breakers
realizedUtilityB1 = nan(N,nGrid); % realized utility of agent B scaled by RN derivative in complete markets
realizedUtilityB2 = nan(N,nGrid); % realized utility of agent B scaled by RN derivative with circuit breakers

realizedUtility2A1 = nan(N,nGrid); % realized utility of agent A in complete markets when both agents are rational
realizedUtility2A2 = nan(N,nGrid); % realized utility of agent A in complete markets when B is irrational
realizedUtility2A3 = nan(N,nGrid); % realized utility of agent A in incomplete markets when B is irrational
realizedUtility2B1 = nan(N,nGrid); % realized utility of agent B in complete markets when both agents are rational
realizedUtility2B2 = nan(N,nGrid); % realized utility of agent B notscaled by RN derivative in complete markets when B is irrational
realizedUtility2B3 = nan(N,nGrid); % realized utility of agent B notscaled by RN derivative with circuit breakers when B is irrational
parfor i = 1:nGrid
    tic
    iParams = params;
    iParams.lambda = lambdaGrid(i);
    [phat(i),relativeError] = defineThresholdPrice(iParams,nt);
    fprintf('Threshold price is %4.3f with relative error %5.5f \n', phat(i),relativeError)
    toc
    % simulate paths
    rng(1); % set seed for random number generator
    tic
    [iRealizedUtilityA,iRealizedUtilityB,iRealizedUtility2A,iRealizedUtility2B] = computeUtility1(iParams, N, n, phat(i));
    realizedUtilityA1(:,i) = iRealizedUtilityA{1}'; % realized utility of agent A in complete markets
    realizedUtilityA2(:,i) = iRealizedUtilityA{2}'; % realized utility of agent A with circuit breakers
    realizedUtilityB1(:,i) = iRealizedUtilityB{1}'; % realized utility of agent B scaled by RN derivative in complete markets
    realizedUtilityB2(:,i) = iRealizedUtilityB{2}'; % realized utility of agent B scaled by RN derivative with circuit breakers
    
    realizedUtility2A1(:,i) = iRealizedUtility2A{1}'; % realized utility of agent A in complete markets when both agents are rational
    realizedUtility2A2(:,i) = iRealizedUtility2A{2}'; % realized utility of agent A in complete markets when B is irrational
    realizedUtility2A3(:,i) = iRealizedUtility2A{3}'; % realized utility of agent A in incomplete markets when B is irrational
    realizedUtility2B1(:,i) = iRealizedUtility2B{1}'; % realized utility of agent B in complete markets when both agents are rational
    realizedUtility2B2(:,i) = iRealizedUtility2B{2}'; % realized utility of agent B notscaled by RN derivative in complete markets when B is irrational
    realizedUtility2B3(:,i) = iRealizedUtility2B{3}'; % realized utility of agent B notscaled by RN derivative with circuit breakers when B is irrational
    toc
end
realizedUtilityA{1} = realizedUtilityA1;
realizedUtilityA{2} = realizedUtilityA2;
realizedUtilityB{1} = realizedUtilityB1;
realizedUtilityB{2} = realizedUtilityB2;

realizedUtility2A{1} = realizedUtility2A1;
realizedUtility2A{2} = realizedUtility2A2;
realizedUtility2A{3} = realizedUtility2A3;
realizedUtility2B{1} = realizedUtility2B1;
realizedUtility2B{2} = realizedUtility2B2;
realizedUtility2B{3} = realizedUtility2B3;
end
