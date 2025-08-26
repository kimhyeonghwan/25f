%% Monte Carlo simulation
% volatility and skewness for market with circuit breaker

format long
 
tic 
 
N_sim = 1*1e5;
T_MC = (15+25/60 - 9.5)/6.5; % fraction of day before 3:25pm
t_end_ind = floor(T_MC/dt+1);
 
dPDiff = 0.005;
PDiff_set = 0.005:dPDiff:0.1;
PDiff_set = [0.001 0.0025 PDiff_set];
 
 
tot_vol_MC = zeros(length(PDiff_set),1);
tot_skew_MC = zeros(length(PDiff_set),1);
tot_third_mom_MC = zeros(length(PDiff_set),1);
tot_vol_var = zeros(length(PDiff_set),1);
tot_skew_var = zeros(length(PDiff_set),1);
tot_third_mom_var = zeros(length(PDiff_set),1);
count_MC = zeros(length(PDiff_set),1);
 
% mean moment (normalized to a day)
mean_vol = NaN(length(PDiff_set),1);
mean_skew = NaN(length(PDiff_set),1);
mean_third_mom = NaN(length(PDiff_set),1);
std_vol = NaN(length(PDiff_set),1);
std_skew = NaN(length(PDiff_set),1);
std_third_mom = NaN(length(PDiff_set),1);
 
 
rng('default')
 
for sim = 1:N_sim
% simulate price path
MC_rand = randn(1, t_end_ind);
 
for pair_MC = 1:2
 
eta_MC = NaN(1, t_end_ind);
S_CB_MC = NaN(1, t_end_ind);
eta_MC(1) = 1;
S_CB_MC(1) = S_CB_initial;
 
Z_sim = cumsum(MC_rand)*sqrt(dt);
Z_sim = max(min(Z_sim, ZMax), ZMin);
 
trigger = find(Z_sim < Zlb(1:t_end_ind), 1);
 
if isempty(trigger)
    trigger = t_end_ind;
end
 
Z_sim_ind_l = floor((Z_sim(1:trigger) - ZMin)/dz)+1;
Z_sim_ind_u = min(Z_sim_ind_l+1,M);
weight = NaN(1,trigger);
 
for n = 2:trigger
    eta_MC(n) = eta_MC(n-1) * (1 + nu/sigma * Z_sim(n-1) * (Z_sim(n)-Z_sim(n-1)));
    weight(n) = (Z_sim(n) - Z(Z_sim_ind_l(n)))/dz;
    S_CB_MC(n) = (1-weight(n)) * PRICE(V(Z_sim_ind_l(n),n), VB(Z_sim_ind_l(n),n), lambda, eta_MC(n))...
                 +weight(n) * PRICE(V(Z_sim_ind_u(n),n), VB(Z_sim_ind_u(n),n), lambda, eta_MC(n));
end
 
 
 
% find minimal minute PDiff
 
N_min = floor(trigger / 60);
S_CB_min = NaN(1, N_min);
 
for m = 1:(N_min-1)
    S_CB_min(m) = min(S_CB_MC(60*(m-1)+1: min(60*m+1,trigger)));
end
 
PDiff = (S_CB_min - (1-threshold)*S_CB_initial)./S_CB_initial;
 
 
for n = 1:length(PDiff_set)
   PDiff_l = PDiff_set(n)-0.001;
   PDiff_u = PDiff_set(n)+0.001;
   I_min = find(PDiff>PDiff_l & PDiff<=PDiff_u);
   if ~isempty(I_min)
      for i=1:length(I_min)
      if I_min(i) >10 % exclude the first minute
         
         Ret_MC = log(S_CB_MC(60*(I_min(i)-1)+2:60*I_min(i)+1)) - log(S_CB_MC(60*(I_min(i)-1)+1:60*I_min(i))); 
         count_MC(n) = count_MC(n)+1;
         sec_momSq = norm(Ret_MC);
         third_mom = sum(Ret_MC.^3);
         tot_vol_MC(n) = tot_vol_MC(n) + sec_momSq;
         tot_vol_var(n) = tot_vol_var(n) + sec_momSq^2;
         tot_third_mom_MC(n) = tot_third_mom_MC(n) + third_mom;
         tot_third_mom_var(n) = tot_third_mom_var(n) + third_mom^2;
         ins_skew = third_mom*sqrt(length(Ret_MC)) / (sec_momSq)^3;
         tot_skew_MC(n) = tot_skew_MC(n) + ins_skew;
         tot_skew_var(n) = tot_skew_var(n) + ins_skew^2;
         
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

    for n = 1:length(PDiff_set)
    mean_vol(n) = tot_vol_MC(n) / count_MC(n) * sqrt(60*6.5);
    std_vol(n) = sqrt(tot_vol_var(n) / count_MC(n) - (tot_vol_MC(n) / count_MC(n))^2) / sqrt(count_MC(n)) * sqrt(60*6.5);
    mean_skew(n) = tot_skew_MC(n) / count_MC(n);
    std_skew(n) = sqrt(tot_skew_var(n) / count_MC(n) - (tot_skew_MC(n) / count_MC(n))^2) / sqrt(count_MC(n));
    mean_third_mom(n) = tot_third_mom_MC(n) / count_MC(n) * (60*6.5)^(3/2);
    std_third_mom(n) = sqrt(tot_third_mom_var(n) / count_MC(n) - (tot_third_mom_MC(n) / count_MC(n))^2) / sqrt(count_MC(n)) * (60*6.5)^(3/2);
    end
 
    figure
    subplot(1,2,1)
    hold on
    plot(PDiff_set, mean_vol, '-ok', 'Linewidth', 2)
    title('Volatility')
 
    subplot(1,2,2)
    hold on
    plot(PDiff_set, mean_skew, '-ok', 'Linewidth', 2)
    title('Skewness')
 
 
function S = PRICE(V, VB, lambda, eta)
    S= (lambda/(lambda + (1-lambda)*eta) * V + (1-lambda)*eta / (lambda + (1-lambda)*eta) * VB)^(-1);
end
