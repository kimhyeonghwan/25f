clear all
clc

Bin = [0.75 1.25 2 3 4 5 6 7 8 9 10 11 12 13 14];
X = [Bin, fliplr(Bin)];
g = [1 1 1]*0.9;


TSTV_Data = readtable('TSTV_09222022.csv');

TSTV = table2array(TSTV_Data(:,2:3));
TSTV_up = TSTV(:,1)-2*sqrt(TSTV(:,2));
TSTV_down = TSTV(:,1)+2*sqrt(TSTV(:,2));

TSTV_Y = [TSTV_up.', fliplr(TSTV_down.')];


Skew_Data = readtable('Skewness_09222022.csv');

Skew = table2array(Skew_Data(:,2:3));
Skew_up = Skew(:,1)-2*sqrt(Skew(:,2));
Skew_down = Skew(:,1)+2*sqrt(Skew(:,2));

Skew_Y = [Skew_up.', fliplr(Skew_down.')];


Volume_Data = readtable('Volume_09222022.csv');

Volume = table2array(Volume_Data(:,2:3));
Volume_up = Volume(:,1)-2*sqrt(Volume(:,2));
Volume_down = Volume(:,1)+2*sqrt(Volume(:,2));

Volume_Y = [Volume_up.', fliplr(Volume_down.')];



exp_ret_Data = readtable('exp_ret_09222022.csv');

exp_ret = table2array(exp_ret_Data(:,2:3));
exp_ret_up = exp_ret(:,1)-2*sqrt(exp_ret(:,2));
exp_ret_down = exp_ret(:,1)+2*sqrt(exp_ret(:,2));

exp_ret_Y = [exp_ret_up.', fliplr(exp_ret_down.')];



figure
subplot(2,2,1)
hold on
fill(X,TSTV_Y, g, 'LineStyle','none');
plot(Bin, TSTV(:,1), 'ok', 'MarkerSize', 4)
z = axis;
line([7 7],[z(3) z(4)], 'color',[0 0 0],'linewidth',1.5,'LineStyle',':');
xlabel('DTCB (%)')
%ylabel('Volatility (%)')
title('\rm Panel A. Volatility (%)')
xticks([0 2 4 6 7 8 10 12 14])

subplot(2,2,2)
hold on
fill(X,Skew_Y, g, 'LineStyle','none');
plot(Bin, Skew(:,1), 'ok', 'MarkerSize', 4)
z = axis;
line([7 7],[z(3) z(4)], 'color',[0 0 0],'linewidth',1.5,'LineStyle',':');
xlabel('DTCB (%)')
%ylabel('Skewness')
title('\rm Panel B. Skewness')
xticks([0 2 4 6 7 8 10 12 14])

subplot(2,2,3)
hold on
fill(X,exp_ret_Y, g, 'LineStyle','none');
plot(Bin, exp_ret(:,1), 'ok', 'MarkerSize', 4)
z = axis;
line([7 7],[z(3) z(4)], 'color',[0 0 0],'linewidth',1.5,'LineStyle',':');
xlabel('DTCB (%)')
%ylabel('Expected Return (%)')
title('\rm Panel C. Expected Return (%)')
xticks([0 2 4 6 7 8 10 12 14])

subplot(2,2,4)
hold on
fill(X,Volume_Y, g, 'LineStyle','none');
plot(Bin, Volume(:,1), 'ok', 'MarkerSize', 4)
z = axis;
line([7 7],[z(3) z(4)], 'color',[0 0 0],'linewidth',1.5,'LineStyle',':');
xlabel('DTCB (%)')
%ylabel('Abnormal Volume')
title('\rm Panel D. Abnormal Volume')
xticks([0 2 4 6 7 8 10 12 14])


%{
mainPath = '/Users/haoxing/Documents/Work/My paper/circuit_breaker/stoch_del_codes/important code';
addpath('/Users/haoxing/Documents/MATLAB/TIKZ Figures/src')
destination = [mainPath,'/Bin_09222022'];
print(destination,'-dpdf')  
axoptions ={'scaled ticks=false, x tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, y tick label style={/pgf/number format/fixed, /pgf/number format/precision=3}, xticklabel shift={.1cm}, ylabel style={yshift=0.1cm}'}; % first option turns off scientific labeling of axes, 2nd fixes label styles across subplots, third insures that tick labels do not overlap
matlab2tikz([destination,'.tikz'],'height', '\figureheight', 'width', '\figurewidth','extraAxisOptions',axoptions)
%close
%}
