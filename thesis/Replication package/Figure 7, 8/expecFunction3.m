function [integr1] = expecFunction3(a,c)
% expecFunction2 calculates integral from minus to plus infty
% exp(-(y+c)^2)/(1+exp(-y*a))

% Version 2: 16-May-2017
% Uses Taylor expansion of order 5 and moments of Normal random variable
% Version 3: 17-May-2017 c can be array


A0 = 1./(1+exp(a*c)); % 0 order term in Taylor expansion
A2 = a^2*exp(a*c).*(-1+exp(a*c))./(1+exp(a*c)).^3;
A4 = a^4*(11*exp(2*a*c)-11*exp(3*a*c)+exp(4*a*c)-exp(a*c))./(1+exp(a*c)).^5;

integr1 = sqrt(pi)*(A0+A2/4+A4/32);
end

