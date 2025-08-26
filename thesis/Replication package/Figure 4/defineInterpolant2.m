function [StCB,diffLogStCB,muRtCB,muRBtCB,thetaAtCB,phiAtCB,xBoundary] = defineInterpolant2(params,spd,zGridArray,etaGridArray,wGridArray,xGridArray,tGridArray,volatility,driftA,driftB)
% defineInterpolant2 defines interpolant to calculate model characteristics
% for given t, X, eta

% differs from previous version by accepting arrays with volatility and
% drifts computed using transition probabilities

% upload complete market functions for out-of-boundary values
[St,diffLogSt,muRt,muRBt,~,~,thetaAt,phiAt,~,~] = uploadFunctionsCompleteMarketsEta4(params);



T = params.T;
phat = params.phat;
X_0 = params.X_0;
mu = params.mu;
sigma = params.sigma;
nu = params.nu;
lambda = params.lambda;


% replace first out-of-boundary point with boundary point
tGrid = squeeze(tGridArray(1,1,:));
nt = length(tGrid);
zBoundary = (log(phat/X_0)-(mu-sigma^2/2)*tGrid+(tGrid-T)*(mu-sigma^2/2)+(exp(-2*nu/sigma*(tGrid-T))-1)*sigma^3/(4*nu))./(sigma*exp(nu/sigma*(T-tGrid)));
xBoundary = X_0*exp((mu-sigma^2/2)*tGrid+sigma*zBoundary);
zGrid = zGridArray(:,1,1);
for it = 1:nt
    iJ = find(zGrid < zBoundary(it),1,'last');
    if (zBoundary(it)-zGrid(iJ)) > (zGrid(iJ+1)-zBoundary(it)) % this redifinition allows better stability in vol calculation in the left end
        iJ = iJ+1;
    end
    zGridArray(iJ,:,it) = zBoundary(it);
    etaGridArray(iJ,:,it) = exp(nu/(2*sigma)*(zBoundary(it).^2-tGridArray(iJ,:,it))-nu^2/(2*sigma^2)*wGridArray(iJ,:,it));
    spd(iJ,:,it) = (lambda+(1-lambda)*etaGridArray(iJ,:,it))./(phat);
    xGridArray(iJ,:,it) = xBoundary(it);
end

S = (lambda+(1-lambda)*etaGridArray)./spd;
StCB = @(x,eta,tIndex)defineStockPriceCircuitBreaker(x,eta,tIndex,S,zGridArray,wGridArray,tGrid,St,xBoundary,params);
% all variables that depend on diffusion are calculated using right
% derivative
spdDiffusion = nan(size(spd));
spdDiffusion(1:end-1,:,:) = (spd(2:end,:,:)-spd(1:end-1,:,:))./(zGridArray(2:end,:,:)-zGridArray(1:end-1,:,:));

%volatility = nan(size(spd));
%volatility(1:end-1,:,:) = ((S(2:end,:,:)-S(1:end-1,:,:))./(zGridArray(2:end,:,:)-zGridArray(1:end-1,:,:)))./S(1:end-1,:,:);
diffLogStCB = @(x,eta,tIndex)defineStockPriceCircuitBreaker(x,eta,tIndex,volatility,zGridArray,wGridArray,tGrid,diffLogSt,xBoundary,params);
           
etaDiffusionNumerical = nan(size(etaGridArray));
etaDiffusionNumerical(1:end-1,:,:) = (etaGridArray(2:end,:,:)-etaGridArray(1:end-1,:,:))./(zGridArray(2:end,:,:)-zGridArray(1:end-1,:,:));

%driftA = (spdDiffusion./spd).^2-((1-lambda)*etaDiffusionNumerical.*spdDiffusion)./(spd.*(lambda+(1-lambda)*etaGridArray));
muRtCB = @(x,eta,tIndex)defineStockPriceCircuitBreaker(x,eta,tIndex,driftA,zGridArray,wGridArray,tGrid,muRt,xBoundary,params);

%driftB = driftA+volatility.*zGridArray*nu/sigma;
muRBtCB = @(x,eta,tIndex)defineStockPriceCircuitBreaker(x,eta,tIndex,driftB,zGridArray,wGridArray,tGrid,muRBt,xBoundary,params);

thetaA = 1./((lambda+(1-lambda)*etaGridArray)/lambda-(1-lambda)*spd.*etaDiffusionNumerical./(lambda*spdDiffusion));
thetaAtCB = @(x,eta,tIndex)defineStockPriceCircuitBreaker(x,eta,tIndex,thetaA,zGridArray,wGridArray,tGrid,thetaAt,xBoundary,params);

phiA = S.*(lambda./(lambda+(1-lambda)*etaGridArray)-thetaA);
phiAtCB = @(x,eta,tIndex)defineStockPriceCircuitBreaker(x,eta,tIndex,phiA,zGridArray,wGridArray,tGrid,phiAt,xBoundary,params);

end

