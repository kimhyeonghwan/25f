function [S,spd,tGridArray,xGridArray,zGridArray,etaGridArray,wGridArray,volatility,driftA,driftB] = solveViaMarkovChain4(params,nt,refinement)
%  solveViaMarkovChain4 solves for price and volatility on grid for t, Z
%  and W via locally consistent Markov Chains

% acceleration ideas: set boundary for W as max(zCeil^2*T,zBoundary^2*T) (this might help to avoid huge matrices)
% differs from version 3 using refinement parameter that allows computing
% variables on much denser grids without storing data
% nt - number of time steps that will be left in solution
% refiniment - multiply nt to obtain steps in solution process

zCeil = 3; % sets maximal value for Z grid
T = params.T;
phat = params.phat;
X_0 = params.X_0;
mu = params.mu;
sigma = params.sigma;
nu = params.nu;
lambda = params.lambda;
%% construct grids
% grid for t
ntDense = (nt-1)*refinement+1; % number of points on a dense grid
tGridDense = linspace(0,T,ntDense);
dt = tGridDense(2)-tGridDense(1);
h = sqrt(dt);

% grid for z
zBoundary = (log(phat/X_0)-(mu-sigma^2/2)*tGridDense+(tGridDense-T)*(mu-sigma^2/2)+(exp(-2*nu/sigma*(tGridDense-T))-1)*sigma^3/(4*nu))./(sigma*exp(nu/sigma*(T-tGridDense)));
zMaxGrid = linspace(-h*(ntDense-1),h*(ntDense-1),2*ntDense-1)';
minIndex  = find(zMaxGrid <= min(zBoundary),1,'last');
maxIndex  = find(zMaxGrid >= zCeil,1);
zGrid = zMaxGrid(minIndex:maxIndex);
u = max(abs(zGrid))^2*h^2; % this parameter is subject to experiments
wGridDense = linspace(0,(ntDense-1)*u,ntDense)';
nwDense = length(wGridDense);
nz = length(zGrid);
% calculate transition probabilities
pW = zGrid.^2*h^2/u;
pZUp = 1/2;
pZDown = 1-pZUp;
if max(pW)>1
    error('Probability greater than 1: tune grids required \n');
end
%xBoundary = exp((mu-sigma^2/2)*tGrid+sigma*zBoundary); % boundary values of X

%zGridArray = repmat(zGrid,[1,nw,nt]);
%wGridArray = permute(repmat(wGridDense,[1,nz,nt]),[2 1 3]);
%tGridArray = permute(repmat(tGrid',[1,nz,nw]),[2 3 1]);
%xGridArray = X_0*exp((mu-sigma^2/2)*tGridArray+sigma*zGridArray);
%etaGridArray = exp(nu/(2*sigma)*(zGridArray.^2-tGridArray)-nu^2/(2*sigma^2)*wGridArray);
%% Calculate SPD on the grid
% initialize SPD for the maximal t
zMatrix = repmat(zGrid,[1,nwDense]);
wMatrix = permute(repmat(wGridDense,[1,nz]),[2 1]);
etaMatrix1 = exp(nu/(2*sigma)*(zMatrix.^2-T)-nu^2/(2*sigma^2)*wMatrix);
xMatrix1 = X_0*exp((mu-sigma^2/2)*T+sigma*zMatrix);
spdMatrix1 = (lambda+(1-lambda)*etaMatrix1)./xMatrix1;

S1 = xMatrix1;

                       
pWMatrix = repmat(pW,[1 nwDense]);
spdMatrix0 = nan(size(spdMatrix1));

subIndexTime = (1:refinement:ntDense);
subIndexW = (1:refinement:ntDense);

if nt~=length(subIndexTime)
    error('Indexes are not aligned')
end
nw = length(subIndexW);

S          = nan(nz,nw,nt); 
S(:,:,nt) = S1(:,subIndexW);
spd        = nan(nz,nw,nt);
spd(:,:,nt) = spdMatrix1(:,subIndexW);
xGridArray = nan(nz,nw,nt);
xGridArray(:,:,nt) = xMatrix1(:,subIndexW);
etaGridArray = nan(nz,nw,nt);
etaGridArray(:,:,nt) = etaMatrix1(:,subIndexW);
zGridArray = nan(nz,nw,nt);
zGridArray(:,:,nt) = zMatrix(:,subIndexW);
wGridArray = nan(nz,nw,nt);
wGridArray(:,:,nt) = wMatrix(:,subIndexW);
volatility = nan(nz,nw,nt);
driftA     = nan(nz,nw,nt);
tGridArray = permute(repmat(tGridDense(subIndexTime)',[1 nz nw]),[2 3 1]);


k = subIndexTime(end-1);
kIndex = length(subIndexTime)-1;
for iReverse = 1:ntDense-1
    volatility0 = nan(nz,nwDense);
    driftLogA0 = nan(nz,nwDense);
    i  = ntDense-iReverse;
    it = tGridDense(i);
    % find first value on grid for Z where spd should be computed using
    % expectation and define spd for all Z below boundary
    iJ = find(zGrid <= zBoundary(i),1,'last');
    etaMatrix0 = exp(nu/(2*sigma)*(zMatrix.^2-it)-nu^2/(2*sigma^2)*wMatrix);
    xMatrix0 = X_0*exp((mu-sigma^2/2)*it+sigma*zMatrix);
    spdMatrix0(1:iJ,:) = (lambda+(1-lambda)*etaMatrix0(1:iJ,:))./(xMatrix0(1:iJ,:).*exp(-(mu-sigma^2/2)*(it-T)-...
                 sigma^3/(4*nu)*(exp(2*nu/sigma*(T-it))-1)-sigma*(1-exp(nu/sigma*(T-it)))*zMatrix(1:iJ,:))); % 
    spdMatrix0(nz,:) = spdCompleteMarkets(mu,sigma,nu,lambda,T,xMatrix0(nz,:),...
                           etaMatrix0(nz,:),nu*zMatrix(nz,:),repmat(it,[1 nwDense]));
             
    izStart = iJ+1;
    izEnd = nz-1;
    iWStart = 1;
    iWEnd = ntDense-1;
    p1 = pZDown*(1-pWMatrix(izStart:izEnd,iWStart:iWEnd));
    p2 = pZUp*(1-pWMatrix(izStart:izEnd,iWStart:iWEnd));
    p3 = pZDown*pWMatrix(izStart:izEnd,iWStart:iWEnd);
    p4 = pZUp*pWMatrix(izStart:izEnd,iWStart:iWEnd);
    
    spdMatrix0(izStart:izEnd,iWStart:iWEnd) = p1.*spdMatrix1(izStart-1:izEnd-1,iWStart:iWEnd) +...
                                        p2.*spdMatrix1(izStart+1:izEnd+1,iWStart:iWEnd)+...
                                        p3.*spdMatrix1(izStart-1:izEnd-1,iWStart+1:iWEnd+1)+...
                                        p4.*spdMatrix1(izStart+1:izEnd+1,iWStart+1:iWEnd+1);
    spdMatrix0(izStart:izEnd,ntDense) = pZDown*spdMatrix1(izStart-1:izEnd-1,ntDense) +...
                                      pZUp*spdMatrix1(izStart+1:izEnd+1,ntDense);                            

    if i == k % save spd, calculate price, vol and drifts
        S0 = (lambda+(1-lambda)*etaMatrix0)./spdMatrix0;
        S1 = (lambda+(1-lambda)*etaMatrix1)./spdMatrix1;
        logS0 = log(S0);
        logS1 = log(S1);
        % calculate drift
        driftLogA0(izStart:izEnd,iWStart:iWEnd) = p1.*logS1(izStart-1:izEnd-1,iWStart:iWEnd) +...
                                      p2.*logS1(izStart+1:izEnd+1,iWStart:iWEnd)+...
                                      p3.*logS1(izStart-1:izEnd-1,iWStart+1:iWEnd+1)+...
                                      p4.*logS1(izStart+1:izEnd+1,iWStart+1:iWEnd+1)-logS0(izStart:izEnd,iWStart:iWEnd);
        driftLogA0(izStart:izEnd,ntDense) = pZDown*logS1(izStart-1:izEnd-1,ntDense) +...
                                          pZUp*logS1(izStart+1:izEnd+1,ntDense) - logS0(izStart:izEnd,ntDense); 

        logS01 = logS0(izStart:izEnd,iWStart:iWEnd);
        volatility0(izStart:izEnd,iWStart:iWEnd) = sqrt(p1.*(logS1(izStart-1:izEnd-1,iWStart:iWEnd)-logS01).^2 +...
                                          p2.*(logS1(izStart+1:izEnd+1,iWStart:iWEnd)-logS01).^2+...
                                          p3.*(logS1(izStart-1:izEnd-1,iWStart+1:iWEnd+1)-logS01).^2+...
                                          p4.*(logS1(izStart+1:izEnd+1,iWStart+1:iWEnd+1)-logS01).^2-...
                                          driftLogA0(izStart:izEnd,iWStart:iWEnd).^2);
        logS02 = logS0(izStart:izEnd,ntDense);
        volatility0(izStart:izEnd,ntDense) = sqrt(pZDown*(logS1(izStart-1:izEnd-1,ntDense)-logS02).^2 +...
                                          pZUp*(logS1(izStart+1:izEnd+1,ntDense)-logS02).^2 -...
                                          driftLogA0(izStart:izEnd,ntDense).^2);   
        
        S(:,:,kIndex) = S0(:,subIndexW);
        spd(:,:,kIndex) = spdMatrix0(:,subIndexW);
        xGridArray(:,:,kIndex) = xMatrix0(:,subIndexW);  
        etaGridArray(:,:,kIndex) = etaMatrix0(:,subIndexW);
        zGridArray(:,:,kIndex) = zMatrix(:,subIndexW);
        wGridArray(:,:,kIndex) = wMatrix(:,subIndexW);
        volatility(:,:,kIndex) = volatility0(:,subIndexW)/sqrt(dt);  
        driftA(:,:,kIndex) = (driftLogA0(:,subIndexW)+volatility0(:,subIndexW).^2/2)/dt;  
        
        kIndex = kIndex-1;
        if kIndex>0
            k = subIndexTime(kIndex);
        end
    end
    spdMatrix1 = spdMatrix0;
    etaMatrix1 = etaMatrix0;
end

driftB = driftA+volatility.*zGridArray*nu/sigma;

