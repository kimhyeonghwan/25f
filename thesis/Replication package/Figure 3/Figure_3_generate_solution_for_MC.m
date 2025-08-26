clear
clc

format long

% parameter values

mu = 0.1/250; % daily mean
sigma = 0.03; % daily vol
D_0 = 1; % starting fundamenta value
lambda = 0.9; % share of stock held by rational agent
T = 1; % one day
threshold = 0.05; % price drop
nu = sigma;

SOReps  = 1e-10;
omega   = 1.2;

N = 6.5*60*60*8; % number of time step
dt = T/N;   % time step
dz = sqrt(dt); % space step
t = 0:dt:T; % time grid

tic
%% Find equilibrium price function

% a, b functions

a = NaN(1,N+1);
b = NaN(1,N+1);

a = -sigma^2/2.*(T-t) - sigma^3/nu/4 .* (1- exp(2*nu/sigma*(T-t)));
b = sigma/nu * (1-exp(nu/sigma*(T-t)));

% CB threshold

S_0= 0.9941;
S_bar = S_0*(1-threshold);

% lower boundary of Z

Zlb_A = 1/sigma * (log(S_bar/D_0) - (mu-sigma^2)*(T-t) -(mu - sigma^2/2)*t);
Zlb_B = 1./(sigma - b*nu) .* (log(S_bar/D_0) - (mu-sigma^2)*(T-t) - (mu-sigma^2/2)*t + a);

Zlb_ma = NaN(2,N+1);
Zlb_ma(1,:) = Zlb_A;
Zlb_ma(2,:) = Zlb_B;

Zlb = max(Zlb_ma); % threshold in Z space for different time
Zlb_ind = NaN(1,N+1);

ZMin = min(Zlb); % lower bound for Z grid
ZMax = 3;        % upper bound for Z grid


% solutions matrix

Z = ZMin:dz:ZMax;
M = length(Z);
Z_0_ind = find(Z>0,1);   % index of Z=0

V = NaN(M, N+1); % agent A evaluation
VB = NaN(M, N+1); % agent B evaluation

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

c0 = (1+dt/(dz)^2)*ones(1,M);
cu= (dt/(2*dz^2))*ones(1,M);
cd = (dt/(2*dz^2))*ones(1,M);
drift = nu/sigma*Z;
c0B = NaN(1,M);
cuB = NaN(1,M);
cdB = NaN(1,M);

for m=2:(M-1)
    c0B(m) = 1 + dt/(dz^2) + max(drift(m),0)*dt/dz + max(-drift(m),0)*dt/dz;
    cuB(m) = dt/(2*dz^2) + max(drift(m),0)*dt/dz;
    cdB(m) = dt/(2*dz^2) + max(-drift(m),0)*dt/dz;
end

VB_new = NaN(1,M);
VB_old = NaN(1,M);
VB_tilde = NaN(1,M);

V_new = NaN(1,M);
V_old = NaN(1,M);
V_tilde = NaN(1,M);

for n=N:-1:1

    % Main SOR step
    SORerror = Inf;

    VB_new = VB(:,n);
    VB_new(Zlb_ind(n):M-1) = VB(Zlb_ind(n):M-1,n+1);
        
    while (SORerror>SOReps)
        VB_old = VB_new;
        for m = Zlb_ind(n):M-1
            VB_tilde(m)   = 1/c0B(m) * (VB(m,n+1) + cuB(m)*VB_old(m+1) + cdB(m)*VB_new(m-1));
            VB_new(m)     = VB_old(m) + omega*(VB_tilde(m) - VB_old(m));
        end
        SORerror = norm(VB_new-VB_old, 'fro');
    end
    VB(:,n) = VB_new;
    
    SORerror = Inf;

    V_new = V(:,n);
    V_new(Zlb_ind(n):M-1) = V(Zlb_ind(n):M-1,n+1);
        
    while (SORerror>SOReps)
        V_old = V_new;
        for m = Zlb_ind(n):M-1
            V_tilde(m)   = 1/c0(m) * (V(m,n+1) + cu(m)*V_old(m+1) + cd(m)*V_new(m-1));
            V_new(m)     = V_old(m) + omega*(V_tilde(m) - V_old(m));
        end
        SORerror = norm(V_new-V_old, 'fro');
    end
    V(:,n) = V_new;
    
    if mod(n, 100) == 0
        disp(['time step = ', num2str(n)])
    end
end

% CB Price

eta = 1;

% Price at time 0 with CB 

S_CB_0 = (lambda/(lambda + (1-lambda)*eta) * V(:,1) + (1-lambda)*eta / (lambda + (1-lambda)*eta) * VB(:,1)).^(-1);
S_CB_initial = S_CB_0(Z_0_ind);

disp(['Circuit Breaker S_0= ', num2str(S_CB_initial)])


% Complete market Price

V_CM = NaN(M,N+1);
VB_CM = NaN(M,N+1);
D_Z = NaN(M,N+1);

for n=1:N+1
for m=1:M
    V_CM(m,n) = (exp((mu-sigma^2)*(T-t(n))+(mu-sigma^2/2)*t(n) + sigma*Z(m)))^(-1);
    VB_CM(m,n) = (exp((mu-sigma^2)*(T-t(n))- a(n) - b(n)*nu*Z(m) +(mu-sigma^2/2)*t(n) + sigma*Z(m)))^(-1);
    D_Z(m,n) = D_0 * exp((mu - 1/2*sigma^2)*t(n)+sigma*Z(m));
end
end

S_CM_0 = (lambda/(lambda + (1-lambda)*eta) * V_CM(:,1) + (1-lambda)*eta / (lambda + (1-lambda)*eta) * VB_CM(:,1)).^(-1);
S_CM_initial = S_CM_0(Z_0_ind);

disp(['Complete market S_0= ', num2str(S_CM_initial)])
