function [realizedUtilityA,realizedUtilityB,realizedUtility2A,realizedUtility2B] = computeUtility1(params, N, n, phat)
% computeUtility1 - simulates realized utility of 2 agents in CB and no CB
% cases

% Version 5: 11-May-2017
% computes realized utilities of both agents that are later used in
% welfare analysis

% N - number of paths to simulate
% n - number of steps in discretization scheme
% phat - threshold value for price

mu = params.mu;
sigma = params.sigma;
lambda = params.lambda;
T = params.T;
X_0 = params.D_0;
nu = params.nu;
tGrid = linspace(0,T,n)';
dt = tGrid(2)-tGrid(1);

drift = (mu-sigma^2/2)*T/n*ones(n-1,1);

zBoundary = (log(phat/X_0)-(mu-sigma^2/2)*tGrid+(tGrid-T)*(mu-sigma^2/2)+(exp(-2*nu/sigma*(tGrid-T))-1)*sigma^3/(4*nu))./(sigma*exp(nu/sigma*(T-tGrid)));
Xhat = X_0*exp((mu-sigma^2/2)*tGrid+sigma*zBoundary);

X_T = nan(N,1);
X_tau = nan(N,1);
Z_T = nan(N,1);
Z_tau = nan(N,1);
W_T = nan(N,1);
W_tau = nan(N,1);
S_tau = nan(N,1);
tau = T*ones(N,1);
for i = 1:N
    iShocks = sigma*randn(n-1,1)*sqrt(T/n); % diffusion component of log
    iIncrements = [0; drift+iShocks];
    iLogX = cumsum(iIncrements);
    iZ = (iLogX-(mu-sigma^2/2)*tGrid)/sigma;
    iX = exp(iLogX);
    iW = cumsum(iZ.^2*dt);
    j = find(iX<Xhat,1);
    X_T(i) = iX(end);
    X_tau(i) = iX(end);
    Z_tau(i) = iZ(end);
    Z_T(i) = iZ(end);
    W_tau(i) = iW(end);
    W_T(i) = iW(end);
    S_tau(i) = iX(end);
    if ~isempty(j)
        tau(i) = tGrid(j);
        X_tau(i) = iX(j);
        Z_tau(i) = iZ(j);
        S_tau(i) = phat;
        W_tau(i) = iW(j);
    end
end

xi_T = exp(nu/(2*sigma)*(Z_T.^2-T)-nu^2/(2*sigma^2)*W_T);
xi_tau = exp(nu/(2*sigma)*(Z_tau.^2-tau)-nu^2/(2*sigma^2)*W_tau);

realizedUtilityB{1} = xi_T.*log((1-lambda)*xi_T.*X_T./(lambda+(1-lambda)*xi_T)); % realized utility of agent B scaled by RN derivative in complete markets
realizedUtilityB{2} = xi_tau.*log((1-lambda)*xi_tau.*S_tau./(lambda+(1-lambda)*xi_tau))+xi_tau.*log(X_tau./S_tau)+(mu+nu*Z_tau-sigma^2/2).*xi_tau.*(T-tau); % realized utility of agent B scaled by RN derivative in incomplete markets
realizedUtilityA{1} = log(lambda*X_T./(lambda+(1-lambda)*xi_T)); % realized utility of agent A in complete markets
realizedUtilityA{2} = log(lambda*S_tau./(lambda+(1-lambda)*xi_tau))+(mu-sigma^2/2)*(T-tau)+log(X_tau./S_tau); % realized utility of agent B scaled by RN derivative in incomplete markets

% Calculate realized utility under iobjective beliefs
realizedUtility2A{1} = log(lambda*X_T); % realized utility of agent A in complete markets if both agents are rational
realizedUtility2B{1} = log((1-lambda)*X_T); % realized utility of agent B in complete markets if both agents are rational
realizedUtility2A{2} = log(lambda*X_T./(lambda+(1-lambda)*xi_T)); % realized utility of agent A in complete markets when agent B is irrational
realizedUtility2B{2} = log((1-lambda)*xi_T.*X_T./(lambda+(1-lambda)*xi_T)); % realized utility of agent B in complete markets when agent B is irrational
realizedUtility2A{3} = log(lambda*S_tau./(lambda+(1-lambda)*xi_tau))+log(X_tau./S_tau)+(mu-sigma^2/2)*(T-tau); % realized utility of agent A in incomplete markets when agent B is irrational
realizedUtility2B{3} = log((1-lambda)*xi_tau.*S_tau./(lambda+(1-lambda)*xi_tau))+log(X_tau./S_tau)+(mu-sigma^2/2)*(T-tau); % realized utility of agent B scaled by RN derivative in incomplete markets

end




