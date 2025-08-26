%% Monte Carlo simulation
format long


tic 

parallel_num = 8;
N_sim = floor(2.5*1e6/8);
T_MC = (15+25/60 - 9.5)/6.5; % fraction of day before 3:25pm
t_end_ind = floor(T_MC/dt+1);
subsample = 8;

dPDiff = 0.005;
PDiff_set = 0.005:dPDiff:0.1;
PDiff_set = [0.001 0.0025 PDiff_set];
obs = subsample*60; % calculate minimal within obs dt
pred = subsample*60; % expected return in next pred dt

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

for sim = 1:N_sim
% simulate price path
MC_rand = randn(1, t_end_ind-1);

for pair_MC = 1:2
eta_MC = NaN(1, t_end_ind);
S_CB_MC = NaN(1, t_end_ind);
eta_MC(1) = 1;
S_CB_MC(1) = S_CB_initial;


Z_sim = [0 cumsum(MC_rand)*sqrt(dt)];
Z_sim = min(Z_sim, ZMax);
W_sim = cumsum(Z_sim.^2)*dt;

trigger = find(Z_sim < Zlb(1:t_end_ind), 1);

if isempty(trigger)
    trigger = t_end_ind;
end

Z_sim_ind_l = floor((Z_sim(1:trigger) - ZMin)/dz)+1;
Z_sim_ind_u = min(Z_sim_ind_l+1,M);


for n = 2:trigger
    eta_MC(n) = exp(nu/(2*sigma) * (Z_sim(n)^2 - n*dt) - nu^2/(2*sigma^2) * W_sim(n));
    weight = (Z_sim(n) - Z(Z_sim_ind_l(n)))/dz;
    S_CB_MC(n) = (1-weight) * PRICE(V(Z_sim_ind_l(n),n), VB(Z_sim_ind_l(n),n), lambda, eta_MC(n))...
                 +weight * PRICE(V(Z_sim_ind_u(n),n), VB(Z_sim_ind_u(n),n), lambda, eta_MC(n));
end

% find minimal obs sec PDiff
int = obs;
N_min = floor((trigger-1) / int);
S_CB_min = NaN(1, N_min);

for m = 1:N_min
    S_CB_min(m) = min(S_CB_MC(int*(m-1)+1: min(int*m,trigger)));
end

PDiff = (S_CB_min - (1-threshold)*S_CB_initial)./S_CB_initial;

for n = 1:length(PDiff_set)
   PDiff_l = PDiff_set(n)-0.001;
   PDiff_u = PDiff_set(n)+0.001;
   I_min = find(PDiff>PDiff_l & PDiff<=PDiff_u);
   if ~isempty(I_min)
      for i=1:length(I_min)
      
      if (I_min(i) >1) && (int*I_min(i)+1 <= trigger) % exclude the first minute    
         
         temp_count_MC(n) = temp_count_MC(n)+1;
        
            ind_end = min(trigger, int*I_min(i)+pred+1);
            Ret_next_1min = log(S_CB_MC(ind_end)) - log(S_CB_MC(int*I_min(i)+1));
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





function S = PRICE(V, VB, lambda, eta)
    S= (lambda/(lambda + (1-lambda)*eta) * V + (1-lambda)*eta / (lambda + (1-lambda)*eta) * VB)^(-1);
end