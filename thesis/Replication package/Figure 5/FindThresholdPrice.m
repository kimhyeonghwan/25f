function [phat, rel_err] = FindThresholdPrice(params) 

% parameter values

mu = params.mu; % daily mean
sigma = params.sigma; % daily vol
D_0 = params.D_0; % starting fundamenta value
lambda = params.lambda; % share of stock held by rational agent
T = params.T; % one day
threshold = params.threshold; % price drop
nu = params.nu;

%if threshold < 0.3
    N = 6.5*60*60*4; % number of time step
%else 
    
%end
dt = T/N;   % time step
dz = sqrt(dt*1.1); % space step
t = 0:dt:T; % time grid

a = NaN(1,N+1);
b = NaN(1,N+1);

a = -sigma^2/2.*(T-t) - sigma^3/nu/4 .* (1- exp(2*nu/sigma*(T-t)));
b = sigma/nu * (1-exp(nu/sigma*(T-t)));

% Complete market initial stock price

SA_0 = exp((mu-sigma^2)*T)*D_0;
SB_0 = exp((mu-sigma^2)*T -a(1))*D_0;
S_CM_0 = (lambda/SA_0 + (1-lambda)/SB_0)^(-1);

%disp(['Complete market S_0= ', num2str(S_CM_0)])

% CB threshold

S_0_old = S_CM_0;

rel_err = Inf;


while (abs(rel_err) > 0.05)
S_bar = S_0_old*(1-threshold);

% lower boundary of Z

Zlb_A = 1/sigma * (log(S_bar/D_0) - (mu-sigma^2)*(T-t) -(mu - sigma^2/2)*t);
Zlb_B = 1./(sigma - b*nu) .* (log(S_bar/D_0) - (mu-sigma^2)*(T-t) - (mu-sigma^2/2)*t + a);

Zlb_ma = NaN(2,N+1);
Zlb_ma(1,:) = Zlb_A;
Zlb_ma(2,:) = Zlb_B;

Zlb = max(Zlb_ma); % threshold in Z space for different time
Zlb_ind = NaN(1,N+1);

ZMin = min(Zlb); % lower bound for Z grid
ZMax = 2.5;        % upper bound for Z grid


% solutions matrix

Z = ZMin:dz:ZMax;
M = length(Z);

V = NaN(M, N+1);
VB = NaN(M, N+1);

% boundary values of solution matrix

V(:, N+1) = (D_0 * exp((mu-sigma^2/2)*T + sigma * Z.')).^(-1);
VB(:,N+1) = V(:,N+1);

for n=1:N
for m=1:M
    if Z(m) < Zlb(n)
        V(m,n) = 1/S_bar;
        VB(m,n) = 1/S_bar;
    else
        Zlb_ind(n) = m;
        break;
    end
end
    V(M,n) = (D_0*exp((mu-sigma^2)*(T-t(n)) + (mu-sigma^2/2)*t(n) + sigma*Z(M)))^(-1);
    VB(M,n) = (D_0*exp((mu-sigma^2)*(T-t(n)) + (mu-sigma^2/2)*t(n) - a(n) + (sigma - b(n)*nu)*Z(M)))^(-1);
end

% Solve for V and VB

% Solve for V

c0 = (1-dt/(dz)^2)*ones(1,M);
cu= (dt/(2*dz^2))*ones(1,M);
cd = (dt/(2*dz^2))*ones(1,M);
drift = nu/sigma*Z;
c0B = NaN(1,M);
cuB = NaN(1,M);
cdB = NaN(1,M);

for m=2:(M-1)
    c0B(m) = 1 - dt/(dz^2) - max(drift(m),0)*dt/dz - max(-drift(m),0)*dt/dz;
    cuB(m) = dt/(2*dz^2) + max(drift(m),0)*dt/dz;
    cdB(m) = dt/(2*dz^2) + max(-drift(m),0)*dt/dz;
end

for n=N:-1:1
for m=Zlb_ind(n):(M-1)
    V(m,n) = c0(m)*V(m,n+1) + cu(m)*V(m+1,n+1) + cd(m)*V(m-1,n+1); 
    VB(m,n) = c0B(m)*VB(m,n+1) + cuB(m)*VB(m+1,n+1) + cdB(m)*VB(m-1,n+1);
end
end

% FS_0

FS_0 = (lambda*V(:,1) + (1-lambda)*VB(:,1)).^(-1);
Z_0_ind = find(Z>0,1);
S_CB_0 = FS_0(Z_0_ind);

rel_err = (S_CB_0 - S_0_old)/S_0_old * 100;
S_0_old = S_CB_0;

%disp(['Circuit Breaker S_0= ', num2str(S_CB_0), ' Relative Error (%): ', num2str(rel_err)])

end

phat = S_CB_0 * (1-threshold);

end