function [SGrid,margU1,margU2,theta1Grid,theta2Grid,b1Grid,b2Grid] = equilibriumGrid2Taylor4(params,n1,n2,n3)
% equilibriumGrid1 construct grids for pareto weights, stock prices,
% marginal utilities and portfolios without numerical solution of
% equilibrium

% Version 2: 13-July-2017
% Fixes the error that assumes normal distibution under agent 2 beliefs


epsilon3 = 1e-10; % this value should be determined by maximal necessary lambda
epsilon2 = 1e-3; % this value should be determined by minimal necessary lambda


a = params.a;
sigma = params.sigma;
mu = params.mu;
nu = params.nu;
X_t = params.X_t;
t = params.t;
T = params.T;
mu1 = (mu-sigma^2/2)*(T-t);
sigma1 = sigma*sqrt(T-t);
Z_t = (log(X_t/params.X_0)-(mu-sigma^2/2)*t)/sigma;
deltat = nu*Z_t; % instantaneous disagreement

%% Define grid for scenario 3, agent 1 is marginal, agent 2 is constrained
b1Grid3 = linspace(-a,a-epsilon3,n3)';
b2Grid3 = -b1Grid3;
theta1Grid3 = ones(n3,1);
theta2Grid3 = 1-theta1Grid3;
margUW2Grid3 = 1./(a+b2Grid3);

btilda = X_t*exp(mu1)./(a+b1Grid3);
c = log(btilda)/(sqrt(2)*sigma1);
dU1dB = expecFunction3(sqrt(2)*sigma1,c)./(sqrt(pi)*(a+b1Grid3));
btilda = (a+b1Grid3)/(X_t*exp(mu1));
c = log(btilda)/(sqrt(2)*sigma1);
dU1dTheta = expecFunction3(sqrt(2)*sigma1,c)./(sqrt(pi));
dU1dB(1) = exp(-mu1+sigma1^2/2)/X_t;

SGrid3 = dU1dTheta./dU1dB;
margUW1Grid3 = dU1dB;
%% Define grid for scenario 2, agent 2 is wealthy and marginal, agent 1 is constrained
b1Grid2 = -a*ones(n2,1);
b2Grid2 = -b1Grid2;
%theta1Grid2 = linspace(0+epsilon2,1,n2)';
theta1Grid2 = exp(linspace(log(0+epsilon2),log(1),n2))';
theta2Grid2 = 1-theta1Grid2;

% use Taylor expansion up to the third order
phi = X_t*theta2Grid2/(2*a);
a0 = 1./(1+phi);
a1 = -a0.*phi./(1+phi);
a2 = -a1.*phi./(1+phi);
a3 = -a2.*phi./(1+phi);
a4 = -a3.*phi./(1+phi);
a5 = -a4.*phi./(1+phi);
A = @(k)((mu*k-sigma^2*k/2)*(T-t)+k^2*sigma^3/(4*nu)*(exp(nu/sigma*(T-t)*2)-1));
B = @(k)(k*sigma/nu*(exp(nu/sigma*(T-t))-1));
H = @(k)exp(A(k)+deltat*B(k));
dU2dB = 1/(2*a)*((a0-a1+a2-a3+a4-a5)+(a1-2*a2+3*a3-4*a4+5*a5)*H(1)+(a2-3*a3+6*a4-10*a5)*H(2)+(a3-4*a4+10*a5)*H(3)+(a4-5*a5)*H(4)+a5*H(5));

phi = (2*a)./(X_t*theta2Grid2);
a0 = 1./(1+phi);
a1 = -a0.*phi./(1+phi);
a2 = -a1.*phi./(1+phi);
a3 = -a2.*phi./(1+phi);
a4 = -a3.*phi./(1+phi);
a5 = -a4.*phi./(1+phi);
dU2dTheta = 1./theta2Grid2.*((a0-a1+a2-a3+a4-a5)+(a1-2*a2+3*a3-4*a4+5*a5)*H(-1)+(a2-3*a3+6*a4-10*a5)*H(-2)+(a3-4*a4+10*a5)*H(-3)+(a4-5*a5)*H(-4)+a5*H(-5));
dU2dTheta(end) = X_t/(2*a)*H(1);

SGrid2 = dU2dTheta./dU2dB;
margUW2Grid2 = dU2dB; % marginal utility of agent 2 in scenario 2
margUW1Grid2 = 1./(SGrid2.*theta1Grid2); % marginal utility of agent 2 in scenario 2

%% Define grid for scenario 1: both agents are constrained
SGrid1 = linspace(SGrid2(end),SGrid3(1),n1+2)';
b1Grid1 = -a*ones(n1+2,1);
b2Grid1 = -b1Grid1;
theta1Grid1 = ones(n1+2,1);
theta2Grid1 = 1-theta1Grid1;

margUW2Grid1 = 1/(2*a)*ones(n1+2,1); % marginal utility of agent 2 in scenario 2
margUW1Grid1 = 1./SGrid1; % marginal utility of agent 2 in scenario 2

%% Merge all three scenarios
SGrid = [SGrid2;SGrid1(2:end-1);SGrid3];
theta1Grid = [theta1Grid2;theta1Grid1(2:end-1);theta1Grid3];
theta2Grid = [theta2Grid2;theta2Grid1(2:end-1);theta2Grid3];
b1Grid = [b1Grid2;b1Grid1(2:end-1);b1Grid3];
b2Grid = -b1Grid;

margU1 = [margUW1Grid2;margUW1Grid1(2:end-1);margUW1Grid3];
margU2 = [margUW2Grid2;margUW2Grid1(2:end-1);margUW2Grid3];

end

