%% Plot Figure_3

load('CB_vol_skewness.mat')
load('CB_exp_return.mat')

CB_mean_vol = mean_vol;
CB_mean_skew = mean_skew;
CB_mean_exp_ret = mean_exp_ret;

load('CM_vol_skewness.mat')
load('CM_exp_return.mat')

CM_mean_vol = mean_vol;
CM_mean_skew = mean_skew;
CM_mean_exp_ret = mean_exp_ret;

dPDiff = 0.005;
PDiff_set = 0.005:dPDiff:0.1;
PDiff_set = [0.001 0.0025 PDiff_set];
PDiff_plot = 0.001:0.002:0.1;

p = polyfit(PDiff_set, CB_mean_exp_ret', 5);
CB_exp_ret_fit = polyval(p, PDiff_plot);

p = polyfit(PDiff_set, CM_mean_exp_ret', 5);
CM_exp_ret_fit = polyval(p, PDiff_plot);

figure

subplot(1,3,1)
hold on
q=plot(PDiff_set*100, CM_mean_vol*100, ':', 'color',my_color('maroon'),'linewidth',2); 
p=plot(PDiff_set*100, CB_mean_vol*100, '-', 'color', my_color('pale_blue'),'linewidth',2);
z = axis;
line([5 5],[z(3) z(4)], 'color',[0 0 0],'linewidth',1.5,'LineStyle',':');
legend({'No CB', 'CB'}, 'fontsize', 10, 'interpreter', 'latex', 'location', 'northeast');
title('A. Volatility (%)')
xlabel('DTCB (%)')
xticks([0 5 10])
legend boxoff;
 
subplot(1,3,2)
hold on
q=plot(PDiff_set*100, CM_mean_skew, ':', 'color',my_color('maroon'),'linewidth',2);
p=plot(PDiff_set*100, CB_mean_skew, '-', 'color', my_color('pale_blue'),'linewidth',2);
z = axis;
line([5 5],[z(3) z(4)], 'color',[0 0 0],'linewidth',1.5,'LineStyle',':');
legend({'No CB', 'CB'}, 'fontsize', 10, 'interpreter', 'latex', 'location', 'southeast');
title('B. Skewness')
xlabel('DTCB (%)')
xticks([0 5 10])
legend boxoff;

subplot(1,3,3)
hold on
q=plot(PDiff_plot*100, CM_exp_ret_fit*100, ':', 'color',my_color('maroon'),'linewidth',2);
p=plot(PDiff_plot*100, CB_exp_ret_fit*100, '-', 'color', my_color('pale_blue'),'linewidth',2);
z = axis;
line([5 5],[z(3) z(4)], 'color',[0 0 0],'linewidth',1.5,'LineStyle',':');
legend({'No CB', 'CB'}, 'fontsize', 10, 'interpreter', 'latex', 'location', 'northeast');
title('C. Expected return (%)')
xlabel('DTCB (%)')
xticks([0 5 10])
legend boxoff;

%{
mainPath = '/Users/haoxing/Documents/Work/My paper/circuit_breaker/stoch_del_codes/important code';
addpath('/Users/haoxing/Documents/MATLAB/TIKZ Figures/src')
destination = [mainPath,'/ModelPrediction_comparison_new_5e6_curves'];
print(destination,'-dpdf')  
axoptions ={'scaled ticks=false, x tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, y tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, xticklabel shift={.1cm}, ylabel style={yshift=0.1cm}'}; % first option turns off scientific labeling of axes, 2nd fixes label styles across subplots, third insures that tick labels do not overlap
matlab2tikz([destination,'.tikz'],'height', '\figureheight', 'width', '\figurewidth','extraAxisOptions',axoptions)
%close
%}
 