function [StCB,diffLogStCB,muRtCB,muRBtCB,thetaAtCB,phiAtCB,phat,tGrid,xBoundary] = uploadFunctionsCircuitBreaker4(params,nt,refinement)
% uploadFunctionsCircuitBreaker3 uploads function in the circuit breaker case
% each function takes 3 arguments: X value, eta value and ! index of time
% in time grid 
%  nt - number of time steps, recommended value: 501 (higher values overload memory)

% Version 4 update: solution is obtained on a denser grid using function  solveViaMarkovChain4
% Version 3 update: volatility is calculated using numerical
% differentiation of spline for stock price (attempt to improve vol evaluation at the boundary)
% Version 2 update: calculating drift and volatility using
% transition probabilities of Markov Chain (rather than numerical differentiation)

[phat,relativeError] = defineThresholdPrice(params,nt);
fprintf('Threshold price is %4.3f with relative error %5.5f \n', phat,relativeError)

params.phat = phat;
%[S,spd,tGridArray,xGridArray,zGridArray,etaGridArray,wGridArray,volatility,driftA,driftB] = solveViaMarkovChain3(params,nt);
[S,spd,tGridArray,xGridArray,zGridArray,etaGridArray,wGridArray,volatility,driftA,driftB] = solveViaMarkovChain4(params,nt,refinement);
%[StCB,diffLogStCB,muRtCB,muRBtCB,thetaAtCB,phiAtCB,xBoundary] = defineInterpolant1(params,spd,zGridArray,etaGridArray,wGridArray,xGridArray,tGridArray);
[StCB,diffLogStCB,muRtCB,muRBtCB,thetaAtCB,phiAtCB,xBoundary] = defineInterpolant2(params,spd,zGridArray,etaGridArray,wGridArray,xGridArray,tGridArray,volatility,driftA,driftB);
tGrid = squeeze(tGridArray(1,1,:));

end

