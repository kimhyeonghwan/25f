
% parameter values

mu = 0.1/250; % daily mean
sigma = 0.03; % daily vol
D_0 = 1; % starting fundamenta value
T = 1; % one day
threshold = 0.05; % price drop
nu = 2*sigma;

nSteps = 500; % number of time step in simulation
alpha_max = 0.31;
alpha_min = 0.001;
alpha_img = linspace(alpha_min^(1/3),alpha_max^(1/3),40);
alphaGrid = (alpha_img.^3)';
lambdaGrid = 0.5;
lambdaBGrid = linspace(0.3,0.5,50)';
NSim = 1000000;



totalWelfareTwoTypes_alpha = nan(length(alphaGrid), length(lambdaBGrid));

for alpha_iter = 1:length(alphaGrid)

threshold = alphaGrid(alpha_iter);
disp(['alpha= ', num2str(threshold)])

nGrid = length(lambdaGrid);
phat = nan(nGrid,1);

realizedUtilityA1 = nan(NSim,nGrid); % realized utility of agent A in complete markets
realizedUtilityA2 = nan(NSim,nGrid); % realized utility of agent A with circuit breakers
realizedUtilityB1 = nan(NSim,nGrid); % realized utility of agent B scaled by RN derivative in complete markets
realizedUtilityB2 = nan(NSim,nGrid); % realized utility of agent B scaled by RN derivative with circuit breakers

realizedUtility2A1 = nan(NSim,nGrid); % realized utility of agent A in complete markets when both agents are rational
realizedUtility2A2 = nan(NSim,nGrid); % realized utility of agent A in complete markets when B is irrational
realizedUtility2A3 = nan(NSim,nGrid); % realized utility of agent A in incomplete markets when B is irrational
realizedUtility2B1 = nan(NSim,nGrid); % realized utility of agent B in complete markets when both agents are rational
realizedUtility2B2 = nan(NSim,nGrid); % realized utility of agent B notscaled by RN derivative in complete markets when B is irrational
realizedUtility2B3 = nan(NSim,nGrid); % realized utility of agent B notscaled by RN derivative with circuit breakers when B is irrational

parfor i = 1:nGrid
    tic
    iParams = struct('mu', mu, 'nu', nu, 'sigma', sigma, 'D_0', D_0, 'lambda', lambdaGrid(i), 'T', T, 'threshold', threshold);
    
    [phat(i),rel_err] = FindThresholdPrice(iParams);
    disp(['Circuit Breaker Threshold = ', num2str(phat(i)), ' Relative Error (%): ', num2str(rel_err)])
    
    toc
    % simulate paths
    rng(1); % set seed for random number generator
    tic
    [iRealizedUtilityA,iRealizedUtilityB,iRealizedUtility2A,iRealizedUtility2B] = computeUtility1(iParams, NSim, nSteps, phat(i));
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

for j=1:length(lambdaBGrid)
    totalWelfareTwoTypes_alpha(alpha_iter,j) = 1-exp(mean(realizedUtilityA{2}(:,1))*lambdaGrid(1)-mean(realizedUtilityA{1}(:,1))*lambdaGrid(1)...
                                      +mean(realizedUtilityB{2}(:,1))*(1-lambdaGrid(1))*lambdaBGrid(j)-mean(realizedUtilityB{1}(:,1))*(1-lambdaGrid(1))*lambdaBGrid(j)...
                                      +mean(realizedUtility2B{3}(:,1))*(1-lambdaGrid(1))*(1-lambdaBGrid(j))-mean(realizedUtility2B{2}(:,1))*(1-lambdaGrid(1))*(1-lambdaBGrid(j)));
end

end

% plot figure

figure
plot(alphaGrid, -totalWelfareTwoTypes_alpha(:,28)*100, 'LineWidth',3);
xlabel('$\alpha$','interpreter','latex')
ylabel('Fractional C.E. gain (\%)','interpreter','latex')

rawWelfare = -totalWelfareTwoTypes_alpha(:,28)*100;
alpha_smooth = linspace(0,0.3,100);
p = polyfit(alphaGrid,rawWelfare,6);
Welfare_smooth =  polyval(p,alpha_smooth);

figure
hold on
plot(alpha_smooth, Welfare_smooth, 'LineWidth',3);
plot(alpha_smooth, zeros(1,length(alpha_smooth)), '-k');
xlabel('$\alpha$','interpreter','latex')
ylabel('Fractional C.E. gain (\%)','interpreter','latex')

%{
destination = 'welfare_nu2_omega5_smooth_new';
print(destination,'-dpdf')  
axoptions ={'scaled ticks=false, x tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, y tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, xticklabel shift={.1cm}, ylabel style={yshift=0.1cm}'}; % first option turns off scientific labeling of axes, 2nd fixes label styles across subplots, third insures that tick labels do not overlap
matlab2tikz([destination,'.tikz'], 'height', '\figureheight', 'width', '\figurewidth', 'extraAxisOptions',axoptions);
%}

