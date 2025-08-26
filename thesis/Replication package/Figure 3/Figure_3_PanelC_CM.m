%% Monte Carlo simulation
% expected return for market without circuit breaker

format long


% parameter values

mu = 0.1/250; % daily mean
sigma = 0.03; % daily vol
D_0 = 1; % starting fundamenta value
lambda = 0.9; % share of stock held by rational agent
T = 1; % one day
threshold = 0.05; % price drop
nu = sigma;

N = 6.5*60*60*8; % number of time step
dt = T/N;   % time step
dz = 0.0023; % space step
t = 0:dt:T; % time grid



parallel_num = 8;
N_sim = floor(2*1e6/8);
T_MC = (15+25/60 - 9.5)/6.5; % fraction of day before 3:25pm
t_end_ind = floor(T_MC/dt+1);
S_CM_initial = 0.9994;
subsample = 8;

dPDiff = 0.005;
PDiff_set = 0.005:dPDiff:0.1;
PDiff_set = [0.001 0.0025 PDiff_set];
obs = subsample*60; % calculate minimal within obs dt
pred = subsample*60; % expected return in next pred dt

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

clear Zlb_A Zlb_B Zlb_ma Zlb
% solutions matrix


tot_exp_ret_MC = zeros(length(PDiff_set),parallel_num);
tot_exp_ret_var = zeros(length(PDiff_set),parallel_num);
count_MC = zeros(length(PDiff_set),parallel_num);

mean_exp_ret = NaN(length(PDiff_set),1);
std_exp_ret = NaN(length(PDiff_set),1);

rng('default')

parfor par = 1:parallel_num
    
temp_tot_exp_ret_MC = zeros(length(PDiff_set),1);
temp_tot_exp_ret_var = zeros(length(PDiff_set),1);
temp_count_MC = zeros(length(PDiff_set),1);

t = 0:dt:T;
Z = ZMin:dz:ZMax;
a = -sigma^2/2.*(T-t) - sigma^3/nu/4 .* (1- exp(2*nu/sigma*(T-t)));
b = sigma/nu * (1-exp(nu/sigma*(T-t)));

for sim = 1:N_sim
% simulate price path
MC_rand = randn(1, t_end_ind-1);

for pair_MC = 1:2
eta_MC = NaN(1, t_end_ind);
weight = zeros(1, t_end_ind);
V_CM_lb = zeros(1, t_end_ind);
V_CM_ub = zeros(1, t_end_ind);
VB_CM_lb = zeros(1, t_end_ind);
VB_CM_ub = zeros(1, t_end_ind);
S_CM_MC = NaN(1, t_end_ind);
eta_MC(1) = 1;



Z_sim = [0 cumsum(MC_rand)*sqrt(dt)];
Z_sim = max(ZMin,min(Z_sim, ZMax));
W_sim = cumsum(Z_sim.^2)*dt;

trigger = t_end_ind;

Z_sim_ind_l = floor((Z_sim(1:trigger) - ZMin)/dz)+1;
Z_sim_ind_u = min(Z_sim_ind_l+1,length(Z));



    eta_MC(2:trigger) = exp(nu/(2*sigma) * (Z_sim(2:trigger).^2 - t(2:trigger)) - nu^2/(2*sigma^2) * W_sim(2:trigger));
    weight(2:trigger) = (Z_sim(2:trigger) - Z(Z_sim_ind_l(2:trigger)))/dz;
    V_CM_lb(2:trigger) = (exp((mu-sigma^2)*(T-t(2:trigger))+(mu-sigma^2/2)*t(2:trigger) + sigma*Z(Z_sim_ind_l(2:trigger)))).^(-1);
    V_CM_ub(2:trigger) = (exp((mu-sigma^2)*(T-t(2:trigger))+(mu-sigma^2/2)*t(2:trigger) + sigma*Z(Z_sim_ind_u(2:trigger)))).^(-1);
    VB_CM_lb(2:trigger) = (exp((mu-sigma^2)*(T-t(2:trigger))- a(2:trigger) - b(2:trigger).*nu.*Z(Z_sim_ind_l(2:trigger)) +(mu-sigma^2/2)*t(2:trigger) + sigma*Z(Z_sim_ind_l(2:trigger)))).^(-1);
    VB_CM_ub(2:trigger) = (exp((mu-sigma^2)*(T-t(2:trigger))- a(2:trigger) - b(2:trigger).*nu.*Z(Z_sim_ind_u(2:trigger)) +(mu-sigma^2/2)*t(2:trigger) + sigma*Z(Z_sim_ind_u(2:trigger)))).^(-1);
    
    S_CM_MC = (1-weight) .* (lambda./(lambda + (1-lambda)*eta_MC) .* V_CM_lb + (1-lambda)*eta_MC ./ (lambda + (1-lambda)*eta_MC) .* VB_CM_lb).^(-1)...
             + weight .* (lambda./(lambda + (1-lambda)*eta_MC) .* V_CM_ub + (1-lambda)*eta_MC ./ (lambda + (1-lambda)*eta_MC) .* VB_CM_ub).^(-1);
    S_CM_MC(1) = S_CM_initial;
    


% find minimal obs sec PDiff
int = obs;
N_min = floor((trigger-1) / int);
S_CM_min = NaN(1, N_min);

for m = 1:N_min
    S_CM_min(m) = min(S_CM_MC(int*(m-1)+1: min(int*m,trigger)));
end

PDiff = (S_CM_min - (1-threshold)*S_CM_initial)./S_CM_initial;

for n = 1:length(PDiff_set)
   PDiff_l = PDiff_set(n)-0.001;
   PDiff_u = PDiff_set(n)+0.001;
   I_min = find(PDiff>PDiff_l & PDiff<=PDiff_u);
   if ~isempty(I_min)
      for i=1:length(I_min)
      
      if (I_min(i) >1) && (int*I_min(i)+1 <= trigger) % exclude the first minute    
        
         temp_count_MC(n) = temp_count_MC(n)+1;
         
            ind_end = min(trigger, int*I_min(i)+pred+1);
          
            Ret_next_1min = log(S_CM_MC(ind_end)) - log(S_CM_MC(int*I_min(i)+1));
            temp_tot_exp_ret_MC(n) = temp_tot_exp_ret_MC(n) + Ret_next_1min;
            temp_tot_exp_ret_var(n) = temp_tot_exp_ret_var(n) + Ret_next_1min^2;
           
      end
      end
   end
end 

MC_rand = - MC_rand;
end

if mod(sim, 1000) == 0
    disp(['sim = ', num2str(sim)])
end

end

count_MC(:, par) = temp_count_MC;
tot_exp_ret_MC(:, par) = temp_tot_exp_ret_MC;
tot_exp_ret_var(:, par) = temp_tot_exp_ret_var;

end

scale = 60*subsample/pred;


    for n = 1:length(PDiff_set)
   
    mean_exp_ret(n) = sum(tot_exp_ret_MC(n,:)) / sum(count_MC(n,:)) * 60*6.5*scale;
    std_exp_ret(n) = sqrt(sum(tot_exp_ret_var(n,:)) / sum(count_MC(n,:)) - (sum(tot_exp_ret_MC(n,:)) / sum(count_MC(n,:)))^2) / sqrt(sum(count_MC(n,:))) *  60*6.5*scale;
    end

   
    figure
    hold on 
    plot(PDiff_set, mean_exp_ret, '-ok', 'Linewidth', 2)
    plot(PDiff_set, mean_exp_ret + 2*std_exp_ret, 'or', 'Linewidth', 1)
    plot(PDiff_set, mean_exp_ret - 2*std_exp_ret, 'or', 'Linewidth', 1)
    title('Expected return')
    xlabel('DTCB')



toc


%save('CM_MC_sim2e6_11042022.mat', 'mean_exp_ret', 'std_exp_ret', 'tot_exp_ret_MC', 'tot_exp_ret_var', 'count_MC');
