function [omegaGrid,lambdaGrid,stopPriceGrid,thetaA,margUA,idxCandid,pW,wGrid] = stopPrices6(params,nt,nOmegaGrid)
%  stopPrices solves for stopping prices and Pareto weights on grid for t, Z
%  and W

n1 = 10;
n2 = floor((nOmegaGrid-n1)/2);
n3 = nOmegaGrid-n1-n2;
if n3 <= 1
    error('Few points on grid')
end

T = params.T;
X_0 = params.X_0;
mu = params.mu;
sigma = params.sigma;
nu = params.nu;

a = params.a;
%% construct grids
% grid for t
tGrid = linspace(0,T,nt);
dt = tGrid(2)-tGrid(1);
h = sqrt(dt);
% grid for z
zGrid = linspace(-h*(nt-1),h*(nt-1),2*nt-1)';
u = max(abs(zGrid))^2*h^2/3; % this parameter is subject to experiments
wGrid = linspace(0,(nt-1)*u,nt)';
nw = length(wGrid);
nz = length(zGrid);
pW = min(zGrid.^2*h^2/u,1);

zGridArray = repmat(zGrid,[1,nw,nt]);
wGridArray = permute(repmat(wGrid,[1,nz,nt]),[2 1 3]);
tGridArray = permute(repmat(tGrid',[1,nz,nw]),[2 3 1]);
xGridArray = X_0*exp((mu-sigma^2/2)*tGridArray+sigma*zGridArray);
etaGridArray = exp(nu/(2*sigma)*(zGridArray.^2-tGridArray)-nu^2/(2*sigma^2)*wGridArray);
%% Calculate SPD and stopping prices conditional on trading halt

% Restrict candidate stopping vertices: this saves memory, accelerates code
% and makes sure the analytical results used in portfolio solution are
% valid
xLow = 0.93; % lowest possible x realization considered (will work for thresholds up to 7%)
% condition on page 4 for solution to be valid: mu2+sigma2^2/2 <
% mu1-sigma1^2/2
mu1 = (mu-sigma^2/2)*(T-tGridArray);
sigma12 = sigma^2*(T-tGridArray);
mu2 = (mu-sigma^2+nu*zGridArray).*(T-tGridArray);
sigma22 = nu^2*(T-tGridArray).^3/3+sigma^2*(T-tGridArray)+(T-tGridArray).^2*nu*sigma;

% remove from the subset vertices that will never be visited
vVisit = false(size(xGridArray));
for i = 1:nt
    vVisit(nt-i+1:nt+i-1,1:i,i) = true; 
end

idxCandid = find((xGridArray > xLow) & (mu2+sigma22/2 < mu1-sigma12/2) & (zGridArray < 0) & vVisit); % linear indices of candidate vertices

nc = length(idxCandid);

% define map objects to store solutions
omegaGrid = nan(nc,nOmegaGrid);
lambdaGrid = nan(nc,nOmegaGrid); % endogenous grid of lambda (corresponding to omegaGrid) for which no trading equilibirum is solved
stopPriceGrid = nan(nc,nOmegaGrid); % no trade prices on lambdaGrid
thetaA = nan(nc,nOmegaGrid); % stock holding of agent A on lambda grid in no trade eq-m
margUA= nan(nc,nOmegaGrid);  % marginal utility of agent A on lambda grid in no trade eq-m 

iparams = params;
for i = 1:nc % solve for stopping prices and 
    iStopNode = idxCandid(i);
    iEta = etaGridArray(iStopNode);
    iX = xGridArray(iStopNode);  
    it = tGridArray(iStopNode);  
    iparams.X_t = iX;
    iparams.t = it;
    [iSGrid,iMargUA,iMargUB,iThetaA,~,ibA,~] = equilibriumGrid2Taylor4(iparams,n1,n2,n3);
    iLambdaGrid = iMargUB*iEta./(iMargUA+iMargUB*iEta);
    iOmegaGrid = (iSGrid.*iThetaA+ibA)./iSGrid;
    omegaGrid(i,:) = iOmegaGrid';
    if iLambdaGrid(end)<0.9999 % this is the case when eta is extemely small and agent B owns share of wealth very close to 0, in this case one can effectively say that for lambda>lambdaMax marginal utility of agent A doesn't change as well as stop price
        iLambdaGrid(end) = 0.9999;
    end
    lambdaGrid(i,:) = iLambdaGrid'; % endogenous grid of lambda (corresponding to omegaGrid) for which no trading equilibirum is solved
    stopPriceGrid(i,:) = iSGrid'; % no trade prices on lambdaGrid
    thetaA(i,:) = iThetaA'; % stock holding of agent A on lambda grid in no trade eq-m
    margUA(i,:) = iMargUA';  % marginal utility of agent A on lambda grid in no trade eq-m 
end

