function [S,spd,tGridArray,xGridArray,zGridArray,etaGridArray,wGridArray,volatility,driftA,driftB,wA,pW,wGrid] = solveViaMarkovChainCM2(params,nt)
%  solveViaMarkovChainCM2 solves for price, volatility, drift on grid for t, Z
%  and W via locally consistent Markov Chains in complete markets case

% acceleration ideas: set boundary for W as max(zCeil^2*T,zBoundary^2*T) (this might help to avoid huge matrices)
% differs from version 1 and 2 by computing volatility and drfit on a grid
% using Markov Chain transition probabilities

% Version 2: 19-07-2017
% Differs from version 1 by refining grid for W


T = params.T;
X_0 = params.X_0;
mu = params.mu;
sigma = params.sigma;
nu = params.nu;
lambda = params.lambda;
a = params.a;
%% construct grids
% grid for t
tGrid = linspace(0,T,nt);
dt = tGrid(2)-tGrid(1);
h = sqrt(dt);
% grid for z
zGrid = linspace(-h*(nt-1),h*(nt-1),2*nt-1)';
% u = max(abs(zGrid))^2*h^2; --Version 1
u = max(abs(zGrid))^2*h^2/3; % this parameter is subject to experiments
wGrid = linspace(0,(nt-1)*u,nt)';
nw = length(wGrid);
nz = length(zGrid);
% calculate transition probabilities
pW = min(zGrid.^2*h^2/u,1);
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
spd = nan(nz,nw,nt); % spd,martingale
sXspd = nan(nz,nw,nt); % stock price times spd,martingale
spd(:,:,end)   = (lambda+(1-lambda)*etaGridArray(:,:,end))./(2*a+xGridArray(:,:,end));
sXspd(:,:,end) = (lambda+(1-lambda)*etaGridArray(:,:,end)).*xGridArray(:,:,end)./(2*a+xGridArray(:,:,end));

waXspd = nan(nz,nw,nt); % wealth of agent A times spd,martingale
waXspd(:,:,end) = (lambda*(xGridArray(:,:,end)+a)-(1-lambda)*etaGridArray(:,:,end)*a)./(2*a+xGridArray(:,:,end));
 
pWArray = repmat(pW,[1 nw]);
for iReverse = 1:nt-1
    i = nt-iReverse;
    izStart = 2;
    izEnd = nz-1;
    iWStart = 1;
    iWEnd = nt-1;
    spd(izStart:izEnd,iWStart:iWEnd,i) = pZDown*(1-pWArray(izStart:izEnd,iWStart:iWEnd)).*spd(izStart-1:izEnd-1,iWStart:iWEnd,i+1) +...
                                      pZUp*(1-pWArray(izStart:izEnd,iWStart:iWEnd)).*spd(izStart+1:izEnd+1,iWStart:iWEnd,i+1)+...
                                      pZDown*pWArray(izStart:izEnd,iWStart:iWEnd).*spd(izStart-1:izEnd-1,iWStart+1:iWEnd+1,i+1)+...
                                      pZUp*pWArray(izStart:izEnd,iWStart:iWEnd).*spd(izStart+1:izEnd+1,iWStart+1:iWEnd+1,i+1); 
                                  
    sXspd(izStart:izEnd,iWStart:iWEnd,i) = pZDown*(1-pWArray(izStart:izEnd,iWStart:iWEnd)).*sXspd(izStart-1:izEnd-1,iWStart:iWEnd,i+1) +...
                                      pZUp*(1-pWArray(izStart:izEnd,iWStart:iWEnd)).*sXspd(izStart+1:izEnd+1,iWStart:iWEnd,i+1)+...
                                      pZDown*pWArray(izStart:izEnd,iWStart:iWEnd).*sXspd(izStart-1:izEnd-1,iWStart+1:iWEnd+1,i+1)+...
                                      pZUp*pWArray(izStart:izEnd,iWStart:iWEnd).*sXspd(izStart+1:izEnd+1,iWStart+1:iWEnd+1,i+1); 
                                  
    waXspd(izStart:izEnd,iWStart:iWEnd,i) = pZDown*(1-pWArray(izStart:izEnd,iWStart:iWEnd)).*waXspd(izStart-1:izEnd-1,iWStart:iWEnd,i+1) +...
                                      pZUp*(1-pWArray(izStart:izEnd,iWStart:iWEnd)).*waXspd(izStart+1:izEnd+1,iWStart:iWEnd,i+1)+...
                                      pZDown*pWArray(izStart:izEnd,iWStart:iWEnd).*waXspd(izStart-1:izEnd-1,iWStart+1:iWEnd+1,i+1)+...
                                      pZUp*pWArray(izStart:izEnd,iWStart:iWEnd).*waXspd(izStart+1:izEnd+1,iWStart+1:iWEnd+1,i+1); 
end
S = sXspd./spd;
wA = waXspd./spd; % wealth of agent A

% Calculate stock price drift and volatility using transition probabilities
volatility = nan(size(spd));
driftLogA = nan(size(spd));
logS = log(S);
for iReverse = 1:nt-1
    i = nt-iReverse;
    
    izStart = 2;
    izEnd = nz-1;
    iWStart = 1;
    iWEnd = nt-1;
    p1 = pZDown*(1-pWArray(izStart:izEnd,iWStart:iWEnd));
    p2 = pZUp*(1-pWArray(izStart:izEnd,iWStart:iWEnd));
    p3 = pZDown*pWArray(izStart:izEnd,iWStart:iWEnd);
    p4 = pZUp*pWArray(izStart:izEnd,iWStart:iWEnd);
    driftLogA(izStart:izEnd,iWStart:iWEnd,i) = p1.*logS(izStart-1:izEnd-1,iWStart:iWEnd,i+1) +...
                                      p2.*logS(izStart+1:izEnd+1,iWStart:iWEnd,i+1)+...
                                      p3.*logS(izStart-1:izEnd-1,iWStart+1:iWEnd+1,i+1)+...
                                      p4.*logS(izStart+1:izEnd+1,iWStart+1:iWEnd+1,i+1)-logS(izStart:izEnd,iWStart:iWEnd,i);
    driftLogA(izStart:izEnd,nt,i) = pZDown*logS(izStart-1:izEnd-1,nt,i+1) +...
                                      pZUp*logS(izStart+1:izEnd+1,nt,i+1) - logS(izStart:izEnd,nt,i); 
    
    iLogS = logS(izStart:izEnd,iWStart:iWEnd,i);
    volatility(izStart:izEnd,iWStart:iWEnd,i) = sqrt(p1.*(logS(izStart-1:izEnd-1,iWStart:iWEnd,i+1)-iLogS).^2 +...
                                      p2.*(logS(izStart+1:izEnd+1,iWStart:iWEnd,i+1)-iLogS).^2+...
                                      p3.*(logS(izStart-1:izEnd-1,iWStart+1:iWEnd+1,i+1)-iLogS).^2+...
                                      p4.*(logS(izStart+1:izEnd+1,iWStart+1:iWEnd+1,i+1)-iLogS).^2-...
                                      driftLogA(izStart:izEnd,iWStart:iWEnd,i).^2);
    iLogS2 = logS(izStart:izEnd,nt,i);
    volatility(izStart:izEnd,nt,i) = sqrt(pZDown*(logS(izStart-1:izEnd-1,nt,i+1)-iLogS2).^2 +...
                                      pZUp*(logS(izStart+1:izEnd+1,nt,i+1)-iLogS2).^2 -...
                                      driftLogA(izStart:izEnd,nt,i).^2);     
end
driftA = (driftLogA+volatility.^2/2)/dt;
volatility = volatility/sqrt(dt);
driftB = driftA+volatility.*zGridArray*nu/sigma;

