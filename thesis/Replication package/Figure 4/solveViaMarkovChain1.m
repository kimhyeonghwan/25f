function [S,spd,tGridArray,xGridArray,zGridArray,etaGridArray,wGridArray] = solveViaMarkovChain1(params,nt)
%  solveViaMarkovChain1 solves for price and volatility on grid for t, Z
%  and W via locally consistent Markov Chains

% acceleration ideas: set boundary for W as max(zCeil^2*T,zBoundary^2*T) (this might help to avoid huge matrices)

zCeil = 3;
T = params.T;
phat = params.phat;
X_0 = params.X_0;
mu = params.mu;
sigma = params.sigma;
nu = params.nu;
lambda = params.lambda;
%% construct grids
% grid for t
tGrid = linspace(0,T,nt);
dt = tGrid(2)-tGrid(1);
h = sqrt(dt);
% grid for z
zBoundary = (log(phat/X_0)-(mu-sigma^2/2)*tGrid+(tGrid-T)*(mu-sigma^2/2)+(exp(-2*nu/sigma*(tGrid-T))-1)*sigma^3/(4*nu))./(sigma*exp(nu/sigma*(T-tGrid)));
zMaxGrid = linspace(-h*(nt-1),h*(nt-1),2*nt-1)';
minIndex  = find(zMaxGrid <= min(zBoundary),1,'last');
maxIndex  = find(zMaxGrid >= zCeil,1);
zGrid = zMaxGrid(minIndex:maxIndex);
u = max(abs(zGrid))^2*h^2; % this parameter is subject to experiments
wGrid = linspace(0,(nt-1)*u,nt)';
nw = length(wGrid);
nz = length(zGrid);
% calculate transition probabilities
pW = zGrid.^2*h^2/u;
pZUp = 1/2;
pZDown = 1-pZUp;
if max(pW)>1
    error('Probability greater than 1: tune grids required \n');
end
%xBoundary = exp((mu-sigma^2/2)*tGrid+sigma*zBoundary); % boundary values of X

zGridArray = repmat(zGrid,[1,nw,nt]);
wGridArray = permute(repmat(wGrid,[1,nz,nt]),[2 1 3]);
tGridArray = permute(repmat(tGrid',[1,nz,nw]),[2 3 1]);
xGridArray = X_0*exp((mu-sigma^2/2)*tGridArray+sigma*zGridArray);
etaGridArray = exp(nu/(2*sigma)*(zGridArray.^2-tGridArray)-nu^2/(2*sigma^2)*wGridArray);
%% Calculate SPD on the grid
spd = nan(nz,nw,nt);
spd(:,:,end) = (lambda+(1-lambda)*etaGridArray(:,:,end))./xGridArray(:,:,end);
spd(nz,:,:) = spdCompleteMarkets(mu,sigma,nu,lambda,T,xGridArray(nz,:,:),...
                           etaGridArray(nz,:,:),nu*zGridArray(nz,:,:),tGridArray(nz,:,:));
                       
startIndex = nan(nt,1); % points to the first z coordinate where I need to compute SPD using Markov chain
for i = 1:nt-1
    iJ = find(zGrid <= zBoundary(i),1,'last');
    spd(1:iJ,:,i) = (lambda+(1-lambda)*etaGridArray(1:iJ,:,i))./(xGridArray(1:iJ,:,i).*exp(-(mu-sigma^2/2)*(tGrid(i)-T)-...
                 sigma^3/(4*nu)*(exp(2*nu/sigma*(T-tGrid(i)))-1)-sigma*(1-exp(nu/sigma*(T-tGrid(i))))*zGridArray(1:iJ,:,i))); % 
    startIndex(i) = iJ+1;
end

pWArray = repmat(pW,[1 nw]);
for iReverse = 1:nt-1
    i = nt-iReverse;
    % NOTE: below efficiency can be improved by omitting w and z where MC never comes 
    izStart = startIndex(i);
    izEnd = nz-1;
    iWStart = 1;
    iWEnd = nt-1;
    spd(izStart:izEnd,iWStart:iWEnd,i) = pZDown*(1-pWArray(izStart:izEnd,iWStart:iWEnd)).*spd(izStart-1:izEnd-1,iWStart:iWEnd,i+1) +...
                                      pZUp*(1-pWArray(izStart:izEnd,iWStart:iWEnd)).*spd(izStart+1:izEnd+1,iWStart:iWEnd,i+1)+...
                                      pZDown*pWArray(izStart:izEnd,iWStart:iWEnd).*spd(izStart-1:izEnd-1,iWStart+1:iWEnd+1,i+1)+...
                                      pZUp*pWArray(izStart:izEnd,iWStart:iWEnd).*spd(izStart+1:izEnd+1,iWStart+1:iWEnd+1,i+1);
    spd(izStart:izEnd,nt,i) = pZDown*spd(izStart-1:izEnd-1,nt,i+1) +...
                                      pZUp*spd(izStart+1:izEnd+1,nt,i+1);                            
    ispd = spd(:,:,i);
    if sum(isnan(ispd(:)))>0
       5; 
    end
end
S = (lambda+(1-lambda)*etaGridArray)./spd;

%{
% calculate diffusion of spd via first difference as left derivative
spdDiffusion = nan(size(spd));
spdDiffusion(2:end,:,:) = (spd(2:end,:,:)-spd(1:end-1,:,:))./(zGridArray(2:end,:,:)-zGridArray(1:end-1,:,:));




spdCM = spdCompleteMarkets(mu,sigma,nu,lambda,T,xGridArray,...
                           etaGridArray,nu*zGridArray,tGridArray);
SCM = (lambda+(1-lambda)*etaGridArray)./spdCM;
it = 50;
iw = floor(it/2+1);
iw = 40;
plot(xGridArray(startIndex(it)-1:end,iw,it),S(startIndex(it)-1:end,iw,it)./xGridArray(startIndex(it)-1:end,iw,it),'o')
hold on
plot(xGridArray(startIndex(it)-1:end,iw,it),SCM(startIndex(it)-1:end,iw,it)./xGridArray(startIndex(it)-1:end,iw,it),'x')
%}


