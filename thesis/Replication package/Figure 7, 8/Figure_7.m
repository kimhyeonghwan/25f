% This script plots graph of price depending on share of wealth of the 2
% agents

addpath('C:/Users/haoxing/Desktop/Work/circuit_breaker/Markov chain bump')

% Parameters

mu = 0.1/250; % daily mean
sigma = 0.03; % daily vol
X_0 = 1; % starting fundamenta value
lambda = 0.9; % share of stock held by rational agent
T = 1; % one day
threshold = 0.05; % price drop
nu = 2*sigma;
D_0 = 1;

params = struct('mu', mu, 'nu', nu, 'sigma', sigma, 'X_0', X_0, 'lambda', lambda, 'T', T, 'threshold', threshold, 'D_0', D_0); 

nOmegaGrid = 100; % number of points on dynamic grid
nt = 101; % number of time steps

a = 0.085;
params.a = a;
tic
[gridSolution.omegaGrid,gridSolution.lambdaGrid,gridSolution.stopPriceGrid,thetaA3,gridSolution.margUA,gridSolution.stVertices] =...
        stopPrices6(params,nt,nOmegaGrid);     
toc
[cmTree.S,~,~,cmTree.xGridArray,cmTree.zGridArray,~,cmTree.wGridArray,vol2,driftA2,~,cmTree.W1] = solveViaMarkovChainCM2(params,nt);

xGrid = cmTree.xGridArray(:,1,1);
zGrid = cmTree.zGridArray(:,1,1);
wGrid = cmTree.wGridArray(1,:,1)';
[~,ix] = min(abs(xGrid-0.97)); % find fundamental grid point closest to 0.97
t = 0.25;
it = (nt-1)*t+1; % time index corresponding to t=0.25

[~,iW] = min(abs(wGrid-zGrid(ix)^2/3*t)); % pick index for w closest to average conditional on x

linearInd = sub2ind(size(cmTree.xGridArray),ix,iW,it);
iCB = find(linearInd==gridSolution.stVertices);
sCBGrid = gridSolution.stopPriceGrid(iCB,:)';
omegaCBGrid = ((gridSolution.omegaGrid(iCB,:).*gridSolution.stopPriceGrid(iCB,:)+params.a)./(gridSolution.stopPriceGrid(iCB,:)+2*params.a))';

   
lambdaGrid = linspace(0.000001,0.999999,nOmegaGrid)';
sCMGrid = nan(nOmegaGrid,1);
omegaCMGrid = nan(nOmegaGrid,1); 
iparams = params;
for i = 1:nOmegaGrid
    iparams.lambda = lambdaGrid(i);
    [iS,~,~,~,~,~,~,~,~,~,iW1] = solveViaMarkovChainCM2(iparams,nt); 
    sCMGrid(i) = iS(ix,iW,it);
    omegaCMGrid(i) =  (iW1(ix,iW,it)+iparams.a)/(iS(ix,iW,it)+2*iparams.a);
end
% Plot static price
figure
plot(omegaCBGrid,sCBGrid,'color',my_color('pale_blue'),'linewidth',2)
hold on
plot(omegaCMGrid,sCMGrid,':','color',my_color('maroon'),'linewidth',2)
xlabel('Share of wealth of agent A: $\omega_{\tau}$','interpreter','latex')
ylabel('$S_{\tau}$','interpreter','latex')
%title(['$\tau=',num2str(t),',\quad \sigma=',num2str(params.vol*100),'\%,\quad a=',num2str(params.a),'$'],'interpreter','latex')
%destination = [directory,'staticPrice'];
%print(destination,'-dpdf')  
%axoptions ={'scaled ticks=false, y tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, xticklabel shift={.1cm}, ylabel style={yshift=0.1cm}'}; % first option turns off scientific labeling of axes, 2nd fixes label styles across subplots, third insures that tick labels do not overlap
%matlab2tikz([destination,'.tikz'],'height', '\figureheight', 'width', '\figurewidth', 'extraAxisOptions',axoptions)

