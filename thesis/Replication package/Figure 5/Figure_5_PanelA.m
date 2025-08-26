
%% Parameters

mu = 0.1/250; % daily mean
sigma = 0.03; % daily vol
X_0 = 1; % starting fundamenta value
lambda = 0.9; % share of stock held by rational agent
T = 1; % one day
threshold = 0.05; % price drop
nu = 2*sigma;
D_0 = 1;


% Utilities A, B

lambdaGrid = linspace(0.001,0.999,30)';
lambdaBGrid = linspace(0.001,0.999,30)';
nSteps = 250;
NSim = 1000000;
nt = 101;

totalWelfareTwoTypes_omega = nan(length(lambdaGrid), length(lambdaBGrid)); 

params = struct('mu', mu, 'nu', nu, 'sigma', sigma, 'X_0', X_0, 'lambda', lambda, 'T', T, 'threshold', threshold, 'D_0', D_0);

[realizedUtilityA,realizedUtilityB,realizedUtility2A,realizedUtility2B,~] = welfareOnLambda1(params,lambdaGrid,NSim,nSteps,nt);


totalWelfareLoss = 1-exp(mean(realizedUtilityB{2}).*(1-lambdaGrid')-mean(realizedUtilityB{1}).*(1-lambdaGrid')+...
                          mean(realizedUtilityA{2}).*lambdaGrid'-mean(realizedUtilityA{1}).*lambdaGrid');


totalWelfareLossObj = 1-exp(mean(realizedUtility2B{3}).*(1-lambdaGrid')-mean(realizedUtility2B{2}).*(1-lambdaGrid')+...
                          mean(realizedUtility2A{3}).*lambdaGrid'-mean(realizedUtility2A{2}).*lambdaGrid');  

for i=1:length(lambdaGrid)
for j=1:length(lambdaBGrid)
    totalWelfareTwoTypes_omega(i,j) = 1-exp(mean(realizedUtilityA{2}(:,i))*lambdaGrid(i)-mean(realizedUtilityA{1}(:,i))*lambdaGrid(i)...
                                      +mean(realizedUtilityB{2}(:,i))*(1-lambdaGrid(i))*lambdaBGrid(j)-mean(realizedUtilityB{1}(:,i))*(1-lambdaGrid(i))*lambdaBGrid(j)...
                                      +mean(realizedUtility2B{3}(:,i))*(1-lambdaGrid(i))*(1-lambdaBGrid(j))-mean(realizedUtility2B{2}(:,i))*(1-lambdaGrid(i))*(1-lambdaBGrid(j)));
end
end
         


%% Plot figure


lambdaGrid = linspace(0.001,0.999,30)';

figure

cs=surfc(lambdaBGrid,lambdaGrid, -totalWelfareTwoTypes_omega*100);
cs(2).ShowText = 'on';
xlabel('$\omega^B$','interpreter','latex')
ylabel('$\omega$','interpreter','latex')
zlabel('Fractional C.E. gain (\%)','interpreter','latex')
xlim([0 1])
ylim([0 1])
zlim([-15 10])
%colormap('gray');

%{
mainPath = '/Users/haoxing/Dropbox/CB_Revision/Replication package/Figure 5';
addpath('/MATLAB/TIKZ Figures/src')
destination = 'Fig5a_bw';
print(destination,'-dpdf')  
axoptions ={'scaled ticks=false, x tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, y tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, xticklabel shift={.1cm}, ylabel style={yshift=0.1cm}'}; % first option turns off scientific labeling of axes, 2nd fixes label styles across subplots, third insures that tick labels do not overlap
matlab2tikz([destination,'.tikz'], 'height', '\figureheight', 'width', '\figurewidth', 'extraAxisOptions',axoptions);
%}
