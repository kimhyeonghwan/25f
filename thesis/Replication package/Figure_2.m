% run this script after running Figure_1.m
% function V, VB, and other variables are from there

dS = 0.002;
duration_step = 5000;
S_lb_set = S_bar+dS:dS:1.02;
hitting_prob = NaN(2,length(S_lb_set));

for S_ind = 1:length(S_lb_set)
S_lb = S_lb_set(S_ind);
S_ub = S_lb+dS;
duration = ((S_lb+S_ub)/2 - S_bar)^2 / (0.99 - S_bar)^2 * 30/(60*6.5);
dt_sample = duration/duration_step;


% Circuit Breaker
S_CB_after = NaN(1, duration_step+1);
Z_after_path = NaN(1,duration_step+1);
W_after_path = NaN(1,duration_step+1);
eta_after = NaN(1,duration_step+1);

I = find(S_CB_path_end < S_ub & S_CB_path_end > S_lb);
I_sample = randsample(length(I), 50000, true);
I_select = I(I_sample);
S_CB_select = S_CB_path_end(I_select);
Z_CB_select = Z_sim(I_select);
W_CB_select = W_sim(I_select);
eta_CB_select = eta_sim(I_select);

rng('default')
MC_rand = randn(length(I_sample), duration_step);
hit = NaN(1, length(I_select));

for i=1:length(I_select)
   Z_after_path(1) = Z_CB_select(i);
   W_after_path(1) = W_CB_select(i);
   S_CB_after(1) = S_CB_select(i);
   eta_after(1) = eta_CB_select(i);
   
   for t_step = 1:duration_step
       t_after = T_MC + t_step*dt_sample;
       Z_after_path(t_step+1) = Z_after_path(t_step) + sqrt(dt_sample) * MC_rand(i,t_step);
       Z_after_path(t_step+1) = max(min(Z_after_path(t_step+1), ZMax), ZMin);
       W_after_path(t_step+1) = W_after_path(t_step) + dt_sample * Z_after_path(t_step+1)^2;
       eta_after(t_step+1) = exp(nu/(2*sigma) * (Z_after_path(t_step+1)^2 - t_after) - nu^2/(2*sigma^2) * W_after_path(t_step+1));
       Z_l_ind = floor((Z_after_path(t_step+1) - ZMin)/dz +1);
       Z_u_ind = min(Z_l_ind+1, length(Z));
       t_l_ind = floor(t_after/dt)+1;
       t_u_ind = t_l_ind+1;
       weight_Z = 1 - (Z_after_path(t_step+1) - ((Z_l_ind-1)*dz + ZMin))/dz;
       weight_t = 1 - (t_after - (t_l_ind-1)*dt)/dt;
       S_CB_after(t_step+1) = weight_t*(weight_Z*PRICE(V(Z_l_ind, t_l_ind), VB(Z_l_ind, t_l_ind), lambda, eta_after(t_step+1))...
                            +(1-weight_Z)*PRICE(V(Z_u_ind, t_l_ind), VB(Z_u_ind, t_l_ind), lambda, eta_after(t_step+1)))...
                            +(1-weight_t)*(weight_Z*PRICE(V(Z_l_ind, t_u_ind), VB(Z_l_ind, t_u_ind), lambda, eta_after(t_step+1))...
                            +(1-weight_Z)*PRICE(V(Z_u_ind, t_u_ind), VB(Z_u_ind, t_u_ind), lambda, eta_after(t_step+1)));
   end
   hit(i) = min(S_CB_after) < S_bar;
 
end

hitting_prob(1,S_ind) = mean(hit);

disp(['S_lb: ', num2str(S_lb), ' CB Hitting prob: ', num2str(hitting_prob(1,S_ind))])


% Complete Market
S_CM_after = NaN(1, duration_step+1);
Z_after_path = NaN(1,duration_step+1);
W_after_path = NaN(1,duration_step+1);
eta_after = NaN(1,duration_step+1);

I = find(S_CM_path_end < S_ub & S_CM_path_end > S_lb);
I_sample = randsample(length(I), 50000, true);
I_select = I(I_sample);
S_CM_select = S_CM_path_end(I_select);
Z_CM_select = Z_sim(I_select);
W_CM_select = W_sim(I_select);
eta_CM_select = eta_sim(I_select);

rng('default')
MC_rand = randn(length(I_sample), duration_step);
hit = NaN(1, length(I_select));

for i=1:length(I_select)
   Z_after_path(1) = Z_CM_select(i);
   W_after_path(1) = W_CM_select(i);
   S_CM_after(1) = S_CM_select(i);
   eta_after(1) = eta_CM_select(i);
   
   for t_step = 1:duration_step
       t_after = T_MC + t_step*dt_sample;
       Z_after_path(t_step+1) = Z_after_path(t_step) + sqrt(dt_sample) * MC_rand(i,t_step);
       Z_after_path(t_step+1) = max(min(Z_after_path(t_step+1), ZMax), ZMin);
       W_after_path(t_step+1) = W_after_path(t_step) + dt_sample * Z_after_path(t_step+1)^2;
       eta_after(t_step+1) = exp(nu/(2*sigma) * (Z_after_path(t_step+1)^2 - t_after) - nu^2/(2*sigma^2) * W_after_path(t_step+1));
       Z_l_ind = floor((Z_after_path(t_step+1) - ZMin)/dz +1);
       Z_u_ind = min(Z_l_ind+1, length(Z));
       t_l_ind = floor(t_after/dt)+1;
       t_u_ind = t_l_ind+1;
       weight_Z = 1 - (Z_after_path(t_step+1) - ((Z_l_ind-1)*dz + ZMin))/dz;
       weight_t = 1 - (t_after - (t_l_ind-1)*dt)/dt;
       S_CM_after(t_step+1) = weight_t*(weight_Z*PRICE(V_CM(Z_l_ind, t_l_ind), VB_CM(Z_l_ind, t_l_ind), lambda, eta_after(t_step+1))...
                            +(1-weight_Z)*PRICE(V_CM(Z_u_ind, t_l_ind), VB_CM(Z_u_ind, t_l_ind), lambda, eta_after(t_step+1)))...
                            +(1-weight_t)*(weight_Z*PRICE(V_CM(Z_l_ind, t_u_ind), VB_CM(Z_l_ind, t_u_ind), lambda, eta_after(t_step+1))...
                            +(1-weight_Z)*PRICE(V_CM(Z_u_ind, t_u_ind), VB_CM(Z_u_ind, t_u_ind), lambda, eta_after(t_step+1)));
   end
   hit(i) = min(S_CM_after) < S_bar;
  
end

hitting_prob(2,S_ind) = mean(hit);

disp(['S_lb: ', num2str(S_lb), ' CM Hitting prob: ', num2str(hitting_prob(2,S_ind))])

end

% plot figure

S_plot = 0.944:0.0005:1.02;

p1 = polyfit(S_lb_set, hitting_prob(1,:), 4);
prob_plot1 = polyval(p1, S_plot);

p2 = polyfit(S_lb_set, hitting_prob(2,:), 4);
prob_plot2 = polyval(p2, S_plot);

figure
hold on
h1=plot(S_plot, prob_plot2, ':','color',my_color('maroon'),'linewidth',2);
h2=plot(S_plot, prob_plot1, '-','color',my_color('pale_blue'),'linewidth',2);
plot(0.994, 0, 'ko', 'MarkerSize', 8)
z = axis;
line([0.944 0.944],[z(3) z(4)], 'color',[1 1 1]*.85,'linewidth',1.5,'LineStyle','-');
xlabel('$S_t$', 'interpreter','latex')
ylabel('Probability','interpreter','latex')
ylim([-0.0005 0.12])

legend([h1 h2], {'no CB', 'CB'}, 'fontsize', 10, 'interpreter', 'latex', 'location', 'northeast');
legend boxoff;


%{
mainPath = '/Users/haoxing/Documents/Work/My paper/circuit_breaker/stoch_del_codes/important code';
addpath('/Users/haoxing/Documents/MATLAB/TIKZ Figures/src')
destination = [mainPath,'/hitting_prob'];
print(destination,'-dpdf')  
axoptions ={'scaled ticks=false, x tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, y tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, xticklabel shift={.1cm}, ylabel style={yshift=0.1cm}'}; % first option turns off scientific labeling of axes, 2nd fixes label styles across subplots, third insures that tick labels do not overlap
matlab2tikz([destination,'.tikz'],'height', '\figureheight', 'width', '\figurewidth','extraAxisOptions',axoptions)
close
%}


function S = PRICE(V, VB, lambda, eta)
    S= (lambda/(lambda + (1-lambda)*eta) * V + (1-lambda)*eta / (lambda + (1-lambda)*eta) * VB)^(-1);
end