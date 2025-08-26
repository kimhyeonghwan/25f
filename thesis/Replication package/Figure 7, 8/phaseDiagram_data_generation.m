% This script calculates variables necessary to plot a phase diagram
% differs from previous version by number of lambda points used
%cd('/mnt/nfs6/petukhov/Matlab/Trading hault/Stochastic beliefs grid/Markov chain bump/Phase diagram')
%addpath(genpath('/mnt/nfs6/petukhov/Matlab/Trading hault/Stochastic beliefs grid/Markov chain bump'))

cd('C:/Users/haoxing/Desktop/Work/circuit_breaker/Replication package/Figure 8')

parpool('local',16)

% load numerical results
%% Do analysis using loop

nlambda = 200;

%aGrid = [0.01;0.05;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1;1.2;1.4;1.6;2;2.5;3;4];
aGrid = [0.0001;linspace(0.01,1,200)'];
na = length(aGrid);
lambdaSet = nan(na,nlambda);

phatSet = nan(na,nlambda);
phatError = nan(na,nlambda);
sumNonNan = nan(na,nlambda);

sumCbHigherEpsilon1 = nan(na,nlambda);
sumCbHigherEpsilon2 = nan(na,nlambda);
sumCmHigherEpsilon1 = nan(na,nlambda);
sumCmHigherEpsilon2 = nan(na,nlambda);

sumPDNonNan = nan(na,nlambda);
sumCbPDHigherEpsilon1 = nan(na,nlambda);
sumCbPDHigherEpsilon2 = nan(na,nlambda);
sumCmPDHigherEpsilon1 = nan(na,nlambda);
sumCmPDHigherEpsilon2 = nan(na,nlambda);

cmRV = nan(na,nlambda);
cmRVError = nan(na,nlambda);
cmRVScaled = nan(na,nlambda);
cmRVScaledError = nan(na,nlambda);
cmPD = nan(na,nlambda);
cmPDError = nan(na,nlambda);

cbRV = nan(na,nlambda);
cbRVError = nan(na,nlambda);
cbRVScaled = nan(na,nlambda);
cbRVScaledError = nan(na,nlambda);
cbPD = nan(na,nlambda);
cbPDError = nan(na,nlambda);


cmS = nan(na,nlambda);
cbS = nan(na,nlambda);
cmW1 = nan(na,nlambda);
cbW1 = nan(na,nlambda);
cmLambda = nan(na,nlambda);

%% Parameters

mu = 0.1/250; % daily mean
sigma = 0.03; % daily vol
X_0 = 1; % starting fundamenta value
lambda = 0.9; % share of stock held by rational agent
T = 1; % one day
threshold = 0.05; % price drop
nu = 2*sigma;
D_0 = 1;

nOmegaGrid = 100; % number of points on dynamic grid
nt = 101; % number of time steps

params = struct('mu', mu, 'nu', nu, 'sigma', sigma, 'X_0', X_0, 'lambda', lambda, 'T', T, 'threshold', threshold, 'D_0', D_0); 


%%

directory = 'C:/Users/haoxing/Desktop/Work/circuit_breaker/Replication package/Figure 8';
for i = 1:na
    tic
    params.a = aGrid(i);
    [gridSolution.omegaGrid,gridSolution.lambdaGrid,gridSolution.stopPriceGrid,thetaA3,gridSolution.margUA,gridSolution.stVertices,gridSolution.pW,gridSolution.wGrid] =...
        stopPrices6(params,nt,nOmegaGrid);
                     
    lambdamin = max(gridSolution.lambdaGrid(:,1))
    lambdamax = min(gridSolution.lambdaGrid(:,end))
    iLambdaGrid = linspace(lambdamin,lambdamax,nlambda);
    lambdaSet(i,:) = iLambdaGrid;
	toc
	tic
    parfor j = 1:nlambda
        jlambda = iLambdaGrid(j);
        [phatStruct,uniformComparison,cmVolSummary,cbVolSummary,solutionSummary] = compareVolSd2(gridSolution,params,threshold,jlambda,nt);

        phatSet(i,j) = phatStruct.phat;
        phatError(i,j) = phatStruct.relativeError;

        sumNonNan(i,j) = uniformComparison.sumNonNan;
        sumCbHigherEpsilon1(i,j) = uniformComparison.sumCbHigherEpsilon1;
        sumCbHigherEpsilon2(i,j) = uniformComparison.sumCbHigherEpsilon2;
        sumCmHigherEpsilon1(i,j) = uniformComparison.sumCmHigherEpsilon1;
        sumCmHigherEpsilon2(i,j) = uniformComparison.sumCmHigherEpsilon2;
        
        sumPDNonNan(i,j) = uniformComparison.sumNonNanPD;
        sumCbPDHigherEpsilon1(i,j) = uniformComparison.sumCbPDHigherEpsilon1;
        sumCbPDHigherEpsilon2(i,j) = uniformComparison.sumCbPDHigherEpsilon2;
        sumCmPDHigherEpsilon1(i,j) = uniformComparison.sumCmPDHigherEpsilon1;
        sumCmPDHigherEpsilon2(i,j) = uniformComparison.sumCmPDHigherEpsilon2;

        cmRV(i,j) = cmVolSummary.cmRV;
        cmRVError(i,j) = cmVolSummary.cmRVError;
        cmRVScaled(i,j) = cmVolSummary.cmRVScaled;
        cmRVScaledError(i,j) = cmVolSummary.cmRVScaledError;
        cmPD(i,j) = cmVolSummary.cmPD;
        cmPDError(i,j) = cmVolSummary.cmPDError;

        cbRV(i,j) = cbVolSummary.cbRV;
        cbRVError(i,j) = cbVolSummary.cbRVError;
        cbRVScaled(i,j) = cbVolSummary.cbRVScaled;
        cbRVScaledError(i,j) = cbVolSummary.cbRVScaledError; 
        cbPD(i,j) = cbVolSummary.cbPD;
        cbPDError(i,j) = cbVolSummary.cbPDError;
        
        cmS(i,j)  = solutionSummary.cmS;
        cbS(i,j)  = solutionSummary.cbS;
        cmW1(i,j) = solutionSummary.cmW1;
        cbW1(i,j) = solutionSummary.cbW1;
        cmLambda(i,j) = solutionSummary.cmLambda;
    end
    toc
    %clear gridSolution
    %save('phaseDiagram.mat')
end