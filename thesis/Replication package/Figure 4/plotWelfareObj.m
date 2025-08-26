function [welfareLossA,welfareLossB,totalWelfareLoss] = plotWelfareObj(lambdaGrid,realizedUtilityA,realizedUtilityB,smoothing,xsymbol)
% plotWelfareObj plots fractional certainty equivalent welfare loss due
% to effects non-rationality under objective beliefs

% lambdaGrid - (nGrid,1) array of lambda
% realizedUtilityA/B - (2,1) cell arrays with realised utilities of agents
% A/B, each array of size (nSimulations,nGrid), array 1 -- complete
% markets, array 2 -- circuit breakers
% smoothing -- logical variable, if true cubic smoothing is performed
% before ploting
% confBounds -- logical, if true, add confidence bounds to graph
% xsymbol -- string label for x axis

% calculate fractional equivalent losses
% Loss (gain) of agent A due to irrationality of B in complete markets
welfareLossA1 = 1-exp(mean(realizedUtilityA{2})-mean(realizedUtilityA{1}));
stdErrorA1 = exp(mean(realizedUtilityA{2})-mean(realizedUtilityA{1})).*std(realizedUtilityA{2}-realizedUtilityA{1})/sqrt(size(realizedUtilityA{2},1));
% Loss(gain) of agent A due to irrationality of B with circuit breakers
welfareLossA2 = 1-exp(mean(realizedUtilityA{3})-mean(realizedUtilityA{1}));
stdErrorA2 = exp(mean(realizedUtilityA{3})-mean(realizedUtilityA{1})).*std(realizedUtilityA{3}-realizedUtilityA{1})/sqrt(size(realizedUtilityA{3},1));

welfareLossA = 1-exp(mean(realizedUtilityA{3})-mean(realizedUtilityA{2}));
stdErrorA = exp(mean(realizedUtilityA{3})-mean(realizedUtilityA{2})).*std(realizedUtilityA{3}-realizedUtilityA{2})/sqrt(size(realizedUtilityA{3},1));

% Loss (gain) of agent B due to irrationality of B in complete markets
welfareLossB1 = 1-exp(mean(realizedUtilityB{2})-mean(realizedUtilityB{1}));
stdErrorB1 = exp(mean(realizedUtilityB{2})-mean(realizedUtilityB{1})).*std(realizedUtilityB{2}-realizedUtilityB{1})/sqrt(size(realizedUtilityB{2},1));
% Loss(gain) of agent B due to irrationality of B with circuit breakers
welfareLossB2 = 1-exp(mean(realizedUtilityB{3})-mean(realizedUtilityB{1}));
stdErrorB2 = exp(mean(realizedUtilityB{3})-mean(realizedUtilityB{1})).*std(realizedUtilityB{3}-realizedUtilityB{1})/sqrt(size(realizedUtilityB{3},1));

welfareLossB = 1-exp(mean(realizedUtilityB{3})-mean(realizedUtilityB{2}));
stdErrorB = exp(mean(realizedUtilityB{3})-mean(realizedUtilityB{2})).*std(realizedUtilityB{3}-realizedUtilityB{2})/sqrt(size(realizedUtilityB{3},1));

Nsim = size(realizedUtilityB{2},1);
totalWelfareLoss1 = 1-exp(mean(realizedUtilityB{2}).*(1-lambdaGrid')-mean(realizedUtilityB{1}).*(1-lambdaGrid')+...
                          mean(realizedUtilityA{2}).*lambdaGrid'-mean(realizedUtilityA{1}).*lambdaGrid');
stdError1 = exp(mean(realizedUtilityB{2}).*(1-lambdaGrid')-mean(realizedUtilityB{1}).*(1-lambdaGrid')+...
                          mean(realizedUtilityA{2}).*lambdaGrid'-mean(realizedUtilityA{1}).*lambdaGrid').*...
                          std(realizedUtilityB{2}.*(1-repmat(lambdaGrid',[Nsim,1]))-mean(realizedUtilityB{1}).*(1-repmat(lambdaGrid',[Nsim,1]))+...
                          mean(realizedUtilityA{2}).*repmat(lambdaGrid',[Nsim,1])-mean(realizedUtilityA{1}).*repmat(lambdaGrid',[Nsim,1]))/sqrt(size(realizedUtilityB{2},1));
totalWelfareLoss2 = 1-exp(mean(realizedUtilityB{3}).*(1-lambdaGrid')-mean(realizedUtilityB{1}).*(1-lambdaGrid')+...
                          mean(realizedUtilityA{3}).*lambdaGrid'-mean(realizedUtilityA{1}).*lambdaGrid');
stdError2 = exp(mean(realizedUtilityB{3}).*(1-lambdaGrid')-mean(realizedUtilityB{1}).*(1-lambdaGrid')+...
                          mean(realizedUtilityA{3}).*lambdaGrid'-mean(realizedUtilityA{1}).*lambdaGrid').*...
                          std(realizedUtilityB{3}.*(1-repmat(lambdaGrid',[Nsim,1]))-mean(realizedUtilityB{1}).*(1-repmat(lambdaGrid',[Nsim,1]))+...
                          mean(realizedUtilityA{3}).*repmat(lambdaGrid',[Nsim,1])-mean(realizedUtilityA{1}).*repmat(lambdaGrid',[Nsim,1]))/sqrt(size(realizedUtilityB{3},1));

totalWelfareLoss = 1-exp(mean(realizedUtilityB{3}).*(1-lambdaGrid')-mean(realizedUtilityB{2}).*(1-lambdaGrid')+...
                          mean(realizedUtilityA{3}).*lambdaGrid'-mean(realizedUtilityA{2}).*lambdaGrid');                   
stdError = exp(mean(realizedUtilityB{3}).*(1-lambdaGrid')-mean(realizedUtilityB{2}).*(1-lambdaGrid')+...
                          mean(realizedUtilityA{3}).*lambdaGrid'-mean(realizedUtilityA{2}).*lambdaGrid').*...
                          std(realizedUtilityB{3}.*(1-repmat(lambdaGrid',[Nsim,1]))-mean(realizedUtilityB{2}).*(1-repmat(lambdaGrid',[Nsim,1]))+...
                          mean(realizedUtilityA{3}).*repmat(lambdaGrid',[Nsim,1])-mean(realizedUtilityA{2}).*repmat(lambdaGrid',[Nsim,1]))/sqrt(size(realizedUtilityB{3},1));                      

                      smoothPar = 1/(1+(abs(lambdaGrid(2)-lambdaGrid(1)))^3/0.6);
if smoothing % smoothes the curves and brings them to 0 
     %totalWelfareLoss = csaps(lambdaGrid,totalWelfareLoss,0.999,lambdaGrid);
     %totalWelfareLoss = totalWelfareLoss-min(totalWelfareLoss); % bring minimum to 0, which can be different due to noise
     welfareLossA1 = csaps(lambdaGrid,welfareLossA1,smoothPar,lambdaGrid);
    %welfareLossA1 = welfareLossA1-min(welfareLossA1);
     welfareLossA2 = csaps(lambdaGrid,welfareLossA2,smoothPar,lambdaGrid);
    %welfareLossA2 = welfareLossA2-min(welfareLossA2);
     welfareLossB1 = csaps(lambdaGrid,welfareLossB1,smoothPar,lambdaGrid);
    %welfareLossB1 = welfareLossB1-min(welfareLossB1);
     welfareLossB2 = csaps(lambdaGrid,welfareLossB2,smoothPar,lambdaGrid);
    %welfareLossB2 = welfareLossB2-min(welfareLossB2);
     totalWelfareLoss1 = welfareLossA1.*lambdaGrid+welfareLossB1.*(1-lambdaGrid);
     totalWelfareLoss2 = welfareLossA2.*lambdaGrid+welfareLossB2.*(1-lambdaGrid);
end

figure
hold on
h1=plot(lambdaGrid,-welfareLossA*100,'-.','color',my_color('pale_blue'),'linewidth',2);
h2=plot(lambdaGrid,-welfareLossB*100,':','color',my_color('green'),'linewidth',2);
h3=plot(lambdaGrid,-totalWelfareLoss*100,'-','color',my_color('maroon'),'linewidth',2);
h4=line(xlim, [0,0], 'Color', uint8([17 17 17]), 'LineWidth', 1);
xlabel(xsymbol,'interpreter','latex')
ylabel('Fractional C.E. gain (\%)','interpreter','latex')
legend([h1 h2 h3], {'Agent A', 'Agnet B', 'Total'}, 'location','southwest')


%{
figure
subplot(1,3,1)
hold on
plot(lambdaGrid,welfareLossA1*100,':','color',my_color('maroon'),'linewidth',2)
plot(lambdaGrid,welfareLossA2*100,'-','color',my_color('pale_blue'),'linewidth',2)
xlabel(xsymbol,'interpreter','latex')
ylabel('Fractional C.E. loss (\%)','interpreter','latex')
title('Agent A','interpreter','latex')
legend('No CB','With CB','location','southeast')
subplot(1,3,2)
hold on
plot(lambdaGrid,welfareLossB1*100,':','color',my_color('maroon'),'linewidth',2)
plot(lambdaGrid,welfareLossB2*100,'-','color',my_color('pale_blue'),'linewidth',2)
ylabel('Fractional C.E. loss (\%)','interpreter','latex')
xlabel(xsymbol,'interpreter','latex')
title('Agent B','interpreter','latex')

subplot(1,3,3)
hold on
plot(lambdaGrid, totalWelfareLoss1*100,':','color',my_color('maroon'),'linewidth',2)
plot(lambdaGrid, totalWelfareLoss2*100,'-','color',my_color('pale_blue'),'linewidth',2)
xlabel(xsymbol,'interpreter','latex')
ylabel('Fractional C.E. loss (\%)','interpreter','latex')
title('Total','interpreter','latex')
%}
%legend('Total welfare','Agent A', 'Agent B','location','northwest')

