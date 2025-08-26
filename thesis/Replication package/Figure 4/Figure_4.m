
% Version 3 9-May-2017 utilizes uploadFunctionsCircuitBreaker4 that uses
% refinement factor to make grid denser
% Version 2 updates: comparative statics of vol wrt lambda is calculated
% using spline differentiation, volatility function uses spline
% differentiation as well

% parameters
mu = 0.1/250; % daily mean
sigma = 0.03; % daily vol
X_0 = 1; % starting fundamenta value
lambda = 0.9; % share of stock held by rational agent
T = 1; % one day
threshold = 0.05; % price drop
nu = sigma;

params = struct('mu', mu, 'nu', nu, 'sigma', sigma, 'X_0', X_0, 'lambda', lambda, 'T', T, 'threshold', threshold);
[StBase,diffusLogStBase,muRtBase,muR2tBase,~,~,theta1tBase,b1tBase,theta2tBase,b2tBase] = uploadFunctionsCompleteMarketsEta4(params);

%% Define stock price, diffusion and other variables for the circuit breakers case
nt = 101; % use nt times refinement points on t axis in discretization scheme
tic
[St,diffusLogSt,muRt,muR2t,theta1t,b1t,phat,tGrid,xBoundary] = uploadFunctionsCircuitBreaker4(params,nt,20);
toc

params.phat = phat;
[StRepAgent,diffLogStRepAgent,~,~,~] = uploadFunctionsCompleteMarketsEta2(params.X_0,params.mu,...
                                    params.sigma,0,params.nu,params.T);  

% plot welfare analysis on lambda
% takes nLamba*20 seconds to compute with 1e6 paths
lambdaGrid = linspace(0.001,0.999,50)';
nSteps = 250;
NSim = 1000000;
nt = 101;
[realizedUtilityA,realizedUtilityB,realizedUtility2A,realizedUtility2B,~] = welfareOnLambda1(params,lambdaGrid,NSim,nSteps,nt);
[welfareLossA,welfareLossB,totalWelfareLoss] = plotWelfare(lambdaGrid,realizedUtilityA,realizedUtilityB,true,false,'$\omega$');

[welfareLossAObj,welfareLossBObj,totalWelfareLossObj] = plotWelfareObj(lambdaGrid,realizedUtility2A,realizedUtility2B,false,'$\omega$');

figure
subplot(1,2,1)
hold on
plot(lambdaGrid,-welfareLossA*100,'-.','color',my_color('pale_blue'),'linewidth',2)
plot(lambdaGrid,-welfareLossB*100,':','color',my_color('green'),'linewidth',2)
plot(lambdaGrid, -totalWelfareLoss*100,'-','color',my_color('maroon'),'linewidth',2)
xlabel('$\omega$','interpreter','latex')
ylabel('Fractional C.E. gain (\%)','interpreter','latex')
box on
legend('Agent A', 'Agent B','Total','location','southwest')

subplot(1,2,2)
hold on
h1=plot(lambdaGrid,-welfareLossAObj*100,'-.','color',my_color('pale_blue'),'linewidth',2);
h2=plot(lambdaGrid,-welfareLossBObj*100,':','color',my_color('green'),'linewidth',2);
h3=plot(lambdaGrid,-totalWelfareLossObj*100,'-','color',my_color('maroon'),'linewidth',2);
h4=line(xlim, [0,0], 'Color', uint8([17 17 17]), 'LineWidth', 1);
xlabel('$\omega$','interpreter','latex')
ylabel('Fractional C.E. gain (\%)','interpreter','latex')
box on
legend([h1 h2 h3], {'Agent A', 'Agnet B', 'Total'}, 'location','southeast')

