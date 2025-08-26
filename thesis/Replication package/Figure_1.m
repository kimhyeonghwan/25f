% instantaneous moments
% stochastic delta
% Figure 1

clear
clc

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

%% Find equilibrium price function

% a, b functions

a = NaN(1,N+1);
b = NaN(1,N+1);

a = -sigma^2/2.*(T-t) - sigma^3/nu/4 .* (1- exp(2*nu/sigma*(T-t)));
b = sigma/nu * (1-exp(nu/sigma*(T-t)));

% CB threshold

S_0= 0.9941;  % obtained from fixed point, agreed with CB S_0 obtained later
S_bar = S_0*(1-threshold);

% lower boundary of Z

Zlb_A = 1/sigma * (log(S_bar/D_0) - (mu-sigma^2)*(T-t) -(mu - sigma^2/2)*t);
Zlb_B = 1./(sigma - b*nu) .* (log(S_bar/D_0) - (mu-sigma^2)*(T-t) - (mu-sigma^2/2)*t + a);

Zlb_ma = NaN(2,N+1);
Zlb_ma(1,:) = Zlb_A;
Zlb_ma(2,:) = Zlb_B;

Zlb = max(Zlb_ma); % threshold in Z space for different time
Zlb_ind = NaN(1,N+1);

ZMin = -3;
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

%% Price related quantities

V_Z = NaN(M,N+1);
VB_Z = NaN(M,N+1);
V_CM_Z = NaN(M,N+1);
VB_CM_Z = NaN(M,N+1);

V_ZZ = NaN(M,N+1);
VB_ZZ = NaN(M,N+1);
V_t = NaN(M,N+1);
VB_t = NaN(M,N+1);
V_CM_ZZ = NaN(M,N+1);
VB_CM_ZZ = NaN(M,N+1);
V_CM_t = NaN(M,N+1);
VB_CM_t = NaN(M,N+1);

for n=1:N
for m=2:M-1
    V_Z(m,n) = (V(m+1,n) - V(m-1,n))/(2*dz);
    VB_Z(m,n) = (VB(m+1,n) - VB(m-1,n))/(2*dz);

    V_CM_Z(m,n) = (V_CM(m+1,n) - V_CM(m-1,n))/(2*dz);
    VB_CM_Z(m,n) = (VB_CM(m+1,n) - VB_CM(m-1,n))/(2*dz);
    
    V_t(m,n) = (V(m,n+1)-V(m,n))/dt;
    VB_t(m,n) = (VB(m,n+1)-VB(m,n))/dt;
    V_ZZ(m,n) = (V(m-1,n)+V(m+1,n)-2*V(m,n))/dz^2;
    VB_ZZ(m,n) = (VB(m-1,n)+VB(m+1,n)-2*VB(m,n))/dz^2;
    V_CM_t(m,n) = (V_CM(m,n+1)-V_CM(m,n))/dt;
    VB_CM_t(m,n) = (VB_CM(m,n+1)-VB_CM(m,n))/dt;
    V_CM_ZZ(m,n) = (V_CM(m-1,n)+V_CM(m+1,n)-2*V_CM(m,n))/dz^2;
    VB_CM_ZZ(m,n) = (VB_CM(m-1,n)+VB_CM(m+1,n)-2*VB_CM(m,n))/dz^2;
   
end
end


%% Monte Carlo simulation

T_MC = 0.1;
t_ind = floor(T_MC/dt+1);
step_MC = 50;
dt_MC = step_MC * dt;
N_MC = floor(T_MC/dt_MC);
Num_sim = 5*1e7;
rng('default')

Z_sim = NaN(1,Num_sim);
W_sim = NaN(1,Num_sim);
Z_sim_ind = NaN(1,Num_sim);
eta_sim = NaN(1,Num_sim);
D_sim = NaN(1,Num_sim);

for sim = 1:Num_sim
    
MC_rand = randn(1, N_MC);
Z_path = 0;  % initial value of Z at time 0
W_path = 0;  % W is time integral of Z^2

for t = 1:N_MC
    Z_path = Z_path + sqrt(dt_MC) * MC_rand(t);
    Z_path = max(min(Z_path, ZMax), ZMin);
    W_path = W_path + dt_MC * Z_path^2;
end

    Z_sim(sim) = Z_path;
    W_sim(sim) = W_path;
    Z_sim_ind(sim) = floor((Z_path-ZMin)/dz)+1;
    eta_sim(sim) = exp(nu/(2*sigma) * (Z_path^2 - T_MC) - nu^2/(2*sigma^2) * W_path);
    D_sim(sim) = exp((mu - sigma^2/2)*T_MC + sigma*Z_path);

end

Z_sim = Z_sim(~isnan(Z_sim));
W_sim = W_sim(~isnan(W_sim));
Z_sim_ind = Z_sim_ind(~isnan(Z_sim_ind));
eta_sim = eta_sim(~isnan(eta_sim));
D_sim = D_sim(~isnan(D_sim));


% Conditioning on D at T_MC

S_CB_path_end = NaN(1,length(Z_sim));
S_CM_path_end = NaN(1,length(Z_sim));
S_CB_repB_path_end = NaN(1,length(Z_sim));
PDiff_CB_path_end = NaN(1,length(Z_sim));
PDiff_CM_path_end = NaN(1,length(Z_sim));
PD_CB_path_end = NaN(1,length(Z_sim));
PD_CB_repB_path_end = NaN(1,length(Z_sim));
PD_CM_path_end = NaN(1,length(Z_sim));
vol_CB_path_end = NaN(1,length(Z_sim));
vol_CM_path_end = NaN(1,length(Z_sim));
mu_CB_path_end = NaN(1,length(Z_sim));
mu_CM_path_end = NaN(1,length(Z_sim));
theta_CB_path_end = NaN(1,length(Z_sim));
theta_CM_path_end = NaN(1,length(Z_sim));
iskew_CB_path_end = NaN(1,length(Z_sim));
iskew_CM_path_end = NaN(1,length(Z_sim));


for n = 1:length(Z_sim)
   S_CB_path_end(n) = PRICE(V(Z_sim_ind(n), t_ind), VB(Z_sim_ind(n), t_ind), lambda, eta_sim(n));
   S_CB_repB_path_end(n) = VB_CM(Z_sim_ind(n), t_ind)^(-1);
   S_CM_path_end(n) = PRICE(V_CM(Z_sim_ind(n), t_ind), VB_CM(Z_sim_ind(n), t_ind), lambda, eta_sim(n));
   PD_CB_path_end(n) = S_CB_path_end(n) / D_sim(n);
   PD_CB_repB_path_end(n) = S_CB_repB_path_end(n) / D_sim(n);
   PD_CM_path_end(n) = S_CM_path_end(n) / D_sim(n);
   PDiff_CB_path_end(n) = (S_CB_path_end(n) - (1-threshold)*S_CB_initial)/S_CB_initial;
   PDiff_CM_path_end(n) = (S_CM_path_end(n) - (1-threshold)*S_CM_initial)/S_CM_initial;
   vol_CB_path_end(n) = VOL(S_CB_path_end(n), eta_sim(n), V_Z(Z_sim_ind(n), t_ind), VB_Z(Z_sim_ind(n), t_ind), V(Z_sim_ind(n), t_ind), VB(Z_sim_ind(n), t_ind), lambda, nu, sigma, Z_sim(n));
   vol_CM_path_end(n) = VOL(S_CM_path_end(n), eta_sim(n), V_CM_Z(Z_sim_ind(n), t_ind), VB_CM_Z(Z_sim_ind(n), t_ind), V_CM(Z_sim_ind(n), t_ind), VB_CM(Z_sim_ind(n), t_ind), lambda, nu, sigma, Z_sim(n));
   theta_CB_path_end(n) = lambda/((1-lambda)*eta_sim(n) + lambda) - lambda*(1-lambda)*eta_sim(n)*Z_sim(n)/((1-lambda)*eta_sim(n) +lambda)^2 / vol_CB_path_end(n);
   theta_CM_path_end(n) = lambda/((1-lambda)*eta_sim(n) + lambda) - lambda*(1-lambda)*eta_sim(n)*Z_sim(n)/((1-lambda)*eta_sim(n) +lambda)^2 / vol_CM_path_end(n);
   mu_CB_path_end(n) = MU(S_CB_path_end(n), eta_sim(n), V_t(Z_sim_ind(n), t_ind), V_ZZ(Z_sim_ind(n), t_ind),...
                        VB_t(Z_sim_ind(n), t_ind), VB_ZZ(Z_sim_ind(n), t_ind), V(Z_sim_ind(n), t_ind), VB(Z_sim_ind(n), t_ind),...
                        V_Z(Z_sim_ind(n), t_ind), VB_Z(Z_sim_ind(n), t_ind), lambda, nu, sigma, Z_sim(n), vol_CB_path_end(n));
   mu_CM_path_end(n) = MU(S_CM_path_end(n), eta_sim(n), V_CM_t(Z_sim_ind(n), t_ind), V_CM_ZZ(Z_sim_ind(n), t_ind),...
                        VB_CM_t(Z_sim_ind(n), t_ind), VB_CM_ZZ(Z_sim_ind(n), t_ind), V_CM(Z_sim_ind(n), t_ind), VB_CM(Z_sim_ind(n), t_ind),...
                        V_CM_Z(Z_sim_ind(n), t_ind), VB_CM_Z(Z_sim_ind(n), t_ind), lambda, nu, sigma, Z_sim(n), vol_CM_path_end(n));
   iskew_CB_path_end(n) = ISKEW(S_CB_path_end(n), vol_CB_path_end(n), V(Z_sim_ind(n), t_ind), VB(Z_sim_ind(n), t_ind),...
       V_Z(Z_sim_ind(n), t_ind), VB_Z(Z_sim_ind(n), t_ind), V_ZZ(Z_sim_ind(n), t_ind), VB_ZZ(Z_sim_ind(n), t_ind), lambda, nu, sigma, eta_sim(n), Z_sim(n));
   
   iskew_CM_path_end(n) = ISKEW(S_CM_path_end(n), vol_CM_path_end(n), V_CM(Z_sim_ind(n), t_ind), VB_CM(Z_sim_ind(n), t_ind),...
       V_CM_Z(Z_sim_ind(n), t_ind), VB_CM_Z(Z_sim_ind(n), t_ind), V_CM_ZZ(Z_sim_ind(n), t_ind), VB_CM_ZZ(Z_sim_ind(n), t_ind), lambda, nu, sigma, eta_sim(n), Z_sim(n));
end


Dl_min = D_Z(Zlb_ind(t_ind), t_ind);  % trigger level in D at time T_MC

dD = 0.005;

Dl_set = Dl_min-0.03:dD:1.056;
Dl_set = [Dl_set(1:7) Dl_min+0.002 Dl_set(8:end)];
S_CB_sim = NaN(1,length(Dl_set));
S_CM_sim = NaN(1,length(Dl_set));
PD_ratio_CB_sim = NaN(1,length(Dl_set));
PD_ratio_CB_repB_sim = NaN(1,length(Dl_set));
PD_ratio_CM_sim = NaN(1,length(Dl_set));
vol_CB_sim = NaN(1,length(Dl_set));
vol_CM_sim = NaN(1,length(Dl_set));
theta_CB_sim = NaN(1,length(Dl_set));
theta_CM_sim = NaN(1,length(Dl_set));
mu_CB_sim = NaN(1,length(Dl_set));
mu_CM_sim = NaN(1,length(Dl_set));
D_mean_sim = NaN(1,length(Dl_set));               

for n = 1:length(Dl_set)

    Dl = Dl_set(n); 
    if n==7 || n==8
        Du = Dl+0.002;
    else
        Du = Dl+dD;
    end

    I = find(D_sim >= Dl & D_sim < Du);
    if ~isempty(I)
    D_mean_sim(n) = mean(D_sim(I));

    % simulated P/D
    PD_ratio_CB_sim(n) = mean(PD_CB_path_end(I), 'omitnan');
    
    PD_ratio_CB_repB_sim(n) = mean(PD_CB_repB_path_end(I), 'omitnan');
   
    PD_ratio_CM_sim(n) = mean(PD_CM_path_end(I), 'omitnan');
    
    % simulated volatility
    vol_CB_sim(n) = mean(vol_CB_path_end(I), 'omitnan');
    
    vol_CM_sim(n) = mean(vol_CM_path_end(I), 'omitnan');
    
    % simulated agent A portfolio holding
    
    theta_CB_sim(n) = mean(theta_CB_path_end(I), 'omitnan');
    
    theta_CM_sim(n) = mean(theta_CM_path_end(I), 'omitnan');
 
    % simulated expected return
    mu_CB_sim(n) = mean(mu_CB_path_end(I), 'omitnan');
   
    mu_CM_sim(n) = mean(mu_CM_path_end(I), 'omitnan');
    
    end
    
    if n==5
       theta_at_CB = mean(lambda./((1-lambda)*eta_sim + lambda), 'omitnan'); 
    end
end


% plot figures

four_std = 4*sqrt(T_MC)*sigma;
D_plot_min = 1 - four_std;
D_plot_max = 1 + four_std;
D_min_ind = 1; %find(D_mean_sim>D_plot_min, 1);
D_max_ind = length(Dl_set); %find(D_mean_sim>D_plot_max, 1);


figure
subplot(1,4,1)
hold on
plot(D_mean_sim(D_min_ind:D_max_ind), PD_ratio_CM_sim(D_min_ind:D_max_ind), ':', 'color',my_color('maroon'),'linewidth',2)
plot(D_mean_sim(7:D_max_ind), PD_ratio_CB_sim(7:D_max_ind), 'color',my_color('pale_blue'),'linewidth',2)
plot(D_mean_sim(D_min_ind:D_max_ind), PD_ratio_CB_repB_sim(D_min_ind:D_max_ind), 'color',[1 1 1]*.85,'linewidth',1.5,'LineStyle','--')
z = axis;
line([Dl_min Dl_min],[z(3) z(4)], 'color',[1 1 1]*.85,'linewidth',1.5,'LineStyle','-');
line([z(1) z(2)],[exp((mu-sigma^2)*(T-T_MC)) exp((mu-sigma^2)*(T-T_MC))], 'color',[1 1 1]*.85,'linewidth',1.5,'LineStyle','--');
xlim([0.95 1.06])
ylim([0.96 1.03])
xlabel('$D_t$', 'interpreter', 'latex')
title('A. $\, S_t/D_t$', 'interpreter', 'latex')
legend({'No CB', 'CB'}, 'fontsize', 10, 'interpreter', 'latex', 'location', 'southeast');
legend boxoff;

subplot(1,4,2)
hold on
plot(D_mean_sim(D_min_ind:D_max_ind), 100*vol_CM_sim(D_min_ind:D_max_ind), ':', 'color',my_color('maroon'),'linewidth',2)
plot(D_mean_sim(7:D_max_ind), 100*vol_CB_sim(7:D_max_ind), 'color',my_color('pale_blue'),'linewidth',2)
vol_repB = -b(t_ind)*nu+sigma;
z = axis;
line([Dl_min Dl_min],[2 z(4)], 'color',[1 1 1]*.85,'linewidth',1.5,'LineStyle','-');
line([z(1) z(2)],[sigma sigma]*100, 'color',[1 1 1]*.85,'linewidth',1.5,'LineStyle','--');
line([z(1) z(2)],[vol_repB vol_repB]*100, 'color',[1 1 1]*.85,'linewidth',1.5,'LineStyle','--');
xlim([0.95 1.06])
ylim([2 z(4)])
xlabel('$D_t$', 'interpreter', 'latex')
title('B. $\, \sigma_{S,t}\, (\%)$', 'interpreter', 'latex')
legend({'No CB', 'CB'}, 'fontsize', 10, 'interpreter', 'latex', 'location', 'northeast');
legend boxoff;

subplot(1,4,3)
hold on
plot(D_mean_sim(D_min_ind:D_max_ind), 100*mu_CM_sim(D_min_ind:D_max_ind), ':', 'color',my_color('maroon'),'linewidth',2)
plot(D_mean_sim(7:D_max_ind), 100*mu_CB_sim(7:D_max_ind), 'color',my_color('pale_blue'),'linewidth',2)
z = axis;
line([Dl_min Dl_min],[z(3) z(4)], 'color',[1 1 1]*.85,'linewidth',1.5,'LineStyle','-');
xlim([0.95 1.06])
ylim([-4 2])
xlabel('$D_t$', 'interpreter', 'latex')
title('C. $\, \mu^A_{S,t}\, (\%)$', 'interpreter', 'latex')
legend({'No CB', 'CB'}, 'fontsize', 10, 'interpreter', 'latex', 'location', 'northeast');
legend boxoff;

subplot(1,4,4)
hold on
plot(D_mean_sim(D_min_ind:D_max_ind), theta_CM_sim(D_min_ind:D_max_ind), ':', 'color',my_color('maroon'),'linewidth',2)
plot(D_mean_sim(7:D_max_ind), theta_CB_sim(7:D_max_ind), 'color',my_color('pale_blue'),'linewidth',2)
plot(Dl_min, theta_at_CB, '.','color',my_color('pale_blue'),'MarkerSize',10)
z = axis;
line([Dl_min Dl_min],[z(3) z(4)], 'color',[1 1 1]*.85,'linewidth',1.5,'LineStyle','-');
xlim([0.95 1.06])
ylim([-5 4])
xlabel('$D_t$', 'interpreter', 'latex')
title('D. $\, \theta^A_{t}$', 'interpreter', 'latex')
legend({'No CB', 'CB'}, 'fontsize', 10, 'interpreter', 'latex', 'location', 'northeast');
legend boxoff;

%{
mainPath = '/Users/haoxing/Documents/Work/My paper/circuit_breaker/stoch_del_codes/important code';
addpath('/Users/haoxing/Documents/MATLAB/TIKZ Figures/src')
destination = [mainPath,'/Fig1'];
print(destination,'-dpdf')  
axoptions ={'scaled ticks=false, x tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, y tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, xticklabel shift={.1cm}, ylabel style={yshift=0.1cm}'}; % first option turns off scientific labeling of axes, 2nd fixes label styles across subplots, third insures that tick labels do not overlap
matlab2tikz([destination,'.tikz'],'height', '\figureheight', 'width', '\figurewidth','extraAxisOptions',axoptions)
close
%}



function S = PRICE(V, VB, lambda, eta)
    S= (lambda/(lambda + (1-lambda)*eta) * V + (1-lambda)*eta / (lambda + (1-lambda)*eta) * VB)^(-1);
end

function vol = VOL(S, eta, V_Z, VB_Z, V, VB, lambda, nu, sigma, Z)
    vol = - S* (lambda/(lambda+(1-lambda)*eta) * V_Z + (1-lambda)*eta/(lambda+(1-lambda)*eta) * VB_Z...
                   -(V-VB)*lambda*(lambda+(1-lambda)*eta)^(-2)*(1-lambda)*eta*nu/sigma*Z);
end

function mu = MU(S, eta, V_t, V_ZZ, VB_t, VB_ZZ, V, VB, V_Z, VB_Z, lambda, nu, sigma, Z, vol)
    mu = - S * (lambda/(lambda+(1-lambda)*eta)*(V_t + 1/2*V_ZZ)...
              +(1-lambda)*eta/(lambda+(1-lambda)*eta)*(VB_t + 1/2*VB_ZZ)...
              +(V - VB)*lambda*(lambda+(1-lambda)*eta)^(-3)*(1-lambda)^2*eta^2*nu^2/sigma^2 * Z^2 ...
              -(V_Z - VB_Z)*lambda*(lambda+(1-lambda)*eta)^(-2)*(1-lambda)*eta*nu/sigma * Z) ...
              + vol^2;  
end

function skew = ISKEW(S, vol, V, VB, V_Z, VB_Z, V_ZZ, VB_ZZ, lambda, nu, sigma, eta, Z)
    vol_mart = vol^2 - S *...
        (lambda/(lambda + (1-lambda)*eta) * V_ZZ + (1-lambda)*eta/(lambda + (1-lambda)*eta) * VB_ZZ...
        -2*(V_Z - VB_Z) * lambda*(1-lambda) * (lambda + (1-lambda)*eta)^(-2) * eta * nu/sigma * Z...
        +2*(V - VB) * lambda*(1-lambda)^2 * (lambda + (1-lambda)*eta)^(-3) * (eta*nu/sigma*Z)^2 ...
        -(V - VB) * lambda*(1-lambda) * (lambda + (1-lambda)*eta)^(-2) * eta * (nu/sigma*Z)^2 ...
        -(V - VB) * lambda*(1-lambda) * (lambda + (1-lambda)*eta)^(-2) * eta * nu/sigma);
    skew = 2* vol^2 * vol_mart;
end