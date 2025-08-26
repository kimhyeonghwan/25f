%% This script plots heatmap for realized volatility
% reuires uploading mat file phaseDiagram with results of computations
%addpath(genpath('/Users/petukhov/Dropbox (MIT)/Research/Trading hault/NOtes/Computations/symbolic/code/Stochastic beliefs 1/Simulations/Stochastic beliefs grid/Markov chain bump'))
%addpath('C:/Users/haoxing/Desktop/Work/circuit_breaker/Markov chain bump')
load("phaseDiagram.mat")

%% Plot heat map for realized volatility
% create a regular spaced grid
cut = 1;
data = cbRVScaled./cmRVScaled;
data = data(:,1:end-cut)-1;

wShareSet = (cmW1+repmat(aGrid,[1 size(cmW1,2)]))./(cmS+2*repmat(aGrid,[1 size(cmW1,2)]));
wShareSet2 = wShareSet(:,1:end-cut);
%wShareSet2 = [wShareSet2, ones(size(aGrid))];
%data = [data, [data(1,end);zeros(length(aGrid)-1,1)]];
aSet = repmat(aGrid,[1 size(wShareSet2,2)]);

% create uniform grids for wealth share and a
naunif = 300;
nbunif = 300;
aLim = 0.75;
aUniformGrid = linspace(min(aGrid),aLim,naunif)';
wUniformGrid = linspace(0.001,1,nbunif);

aUnifromSet = repmat(aUniformGrid,[1 nbunif]);
wUnifromSet = repmat(wUniformGrid,[naunif 1]);

vq = griddata(aSet(:),wShareSet2(:),data(:),aUnifromSet,wUnifromSet);
colormap = redblue(1000);
%colormap = gray(1000);


hmo = HeatMap(vq','Colormap',colormap);

close
ax = hmo.plot; % 'ax' will be a handle to a standard MATLAB axes.

%colorbar('Peer', ax); % Turn the colorbar on
colorbar(ax); % Turn the colorbar on

%ax.Box = 'off'
nytick = 5;
ax.YTick = linspace(ax.YLim(1),ax.YLim(2),nytick);
ax.YTickLabel = strread(num2str(linspace(0,1,nytick)),'%s');
ylabel(ax,'Agent A wealth share','interpreter','latex')

nxtick = 6;
ax.XTick= linspace(ax.YLim(1),ax.YLim(2),nxtick);
ax.XTickLabel = strread(num2str(linspace(0,2*aLim,nxtick)),'%s')
ax.XTickLabelRotation = 0;
xlabel(ax,'Net bond supply','interpreter','latex')
