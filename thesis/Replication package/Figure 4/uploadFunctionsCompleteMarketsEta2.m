function [St,diffLogSt,muRt,muRBt,itoRt,timeRt] = uploadFunctionsCompleteMarketsEta2(X_0,mu,sigma,lambda,nu,T)
% uploadFunctionsCompleteMarketsEta2 uploads stock price, diffusion and
% drift
% Differs from previous version by uploading drift
% muRt - simple return drift
% muRt - simple return drift under agent 2 beliefs

%% this script solves the model with stochastic beliefs and 2 agents
syms X A(t) B(t) W Z t eta

B = sigma/nu*(1-exp(nu/sigma*(T-t)));
A = (t-T)*(mu-sigma^2/2)+sigma^3/(4*nu)*(exp(2*nu/sigma*(T-t))-1);

delta = nu*(log(X)-(mu-sigma^2/2)*t)/sigma;
ksi = exp(nu/(2*sigma)*(Z^2-t)-nu^2/(2*sigma^2)*W); % W = integral_0^t Z_s^2 ds

S = X*(lambda+(1-lambda)*ksi)/(lambda*exp(-(mu-sigma^2)*(T-t))+(1-lambda)*ksi*exp(A+B*delta));
% plot price as function of fundamental value for average realization of
% z^2 path
S1 = subs(S,W,2*sigma^2/nu^2*(nu/(2*sigma)*(Z^2-t)-log(eta)));
S2 = subs(S1,Z,(log(X/X_0)-(mu-sigma^2/2)*t)/sigma);
St = matlabFunction(S2);

S1 = subs(S,X,X_0*exp((mu-sigma^2/2)*t+sigma*Z));
diffusS = diff(S1,Z);
diffLogS = diffusS/S1;
diffLogS = subs(diffLogS,W,2*sigma^2/nu^2*(nu/(2*sigma)*(Z^2-t)-log(eta)));
diffLogS = subs(diffLogS,Z,(log(X/X_0)-(mu-sigma^2/2)*t)/sigma);
diffLogSt = matlabFunction(diffLogS);

% calculate stock price drift
itoTerm = 0.5*diff(diffusS,Z);
S1 = subs(S,X,X_0*exp((mu-sigma^2/2)*t+sigma*Z)); % this is function of W and z only
timeDerivative = diff(S1,t)+diff(S1,W)*Z^2;
muS  = itoTerm+timeDerivative;
muR = muS/S1;
muR  = subs(muR,W,2*sigma^2/nu^2*(nu/(2*sigma)*(Z^2-t)-log(eta)));
muR = subs(muR,Z,(log(X/X_0)-(mu-sigma^2/2)*t)/sigma);
muRt = matlabFunction(muR);

muRB = muR + diffLogS*delta/sigma;
muRBt = matlabFunction(muRB);

itoR = itoTerm/S1;
itoR  = subs(itoR,W,2*sigma^2/nu^2*(nu/(2*sigma)*(Z^2-t)-log(eta)));
itoR = subs(itoR,Z,(log(X/X_0)-(mu-sigma^2/2)*t)/sigma);
itoRt = matlabFunction(itoR);

timeR = timeDerivative/S1;
timeR  = subs(timeR,W,2*sigma^2/nu^2*(nu/(2*sigma)*(Z^2-t)-log(eta)));
timeR = subs(timeR,Z,(log(X/X_0)-(mu-sigma^2/2)*t)/sigma);
timeRt = matlabFunction(timeR);
end

