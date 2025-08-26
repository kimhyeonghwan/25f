function [spd] = spdCompleteMarkets(mu,sigma,nu,lambda,T,X,eta,delta,t)
% spdCompleteMarkets computes SPD in complete markets economy for the case
% when second agent has random walk beliefs


spd = (lambda*exp(-(mu-sigma^2)*(T-t))+(1-lambda)*eta.*exp((t-T)*(mu-sigma^2/2)+(exp(-2*nu/sigma*(t-T))-1)*sigma^3/(4*nu)+delta.*sigma/nu.*(1-exp(nu/sigma*(T-t)))))./...
       X;
end

