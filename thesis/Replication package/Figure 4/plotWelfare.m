function [welfareLossA,welfareLossB,totalWelfareLoss] = plotWelfare(lambdaGrid,realizedUtilityA,realizedUtilityB,smoothing,confBounds,xsymbol)
% plotWelfare plots fractional certainty equivalent welfare loss due
% effects of the circuit breakers

% lambdaGrid - (nGrid,1) array of lambda
% realizedUtilityA/B - (2,1) cell arrays with realised utilities of agents
% A/B, each array of size (nSimulations,nGrid), array 1 -- complete
% markets, array 2 -- circuit breakers
% smoothing -- logical variable, if true cubic smoothing is performed
% before ploting
% confBounds -- logical, if true, add confidence bounds to graph
% xsymbol -- string label for x axis

% calculate fractional equivalent losses
welfareLossA = 1-exp(mean(realizedUtilityA{2})-mean(realizedUtilityA{1}));
stdErrorA = exp(mean(realizedUtilityA{2})-mean(realizedUtilityA{1})).*std(realizedUtilityA{2}-realizedUtilityA{1})/sqrt(size(realizedUtilityA{2},1));

welfareLossB = 1-exp(mean(realizedUtilityB{2})-mean(realizedUtilityB{1}));
stdErrorB = exp(mean(realizedUtilityB{2})-mean(realizedUtilityB{1})).*std(realizedUtilityB{2}-realizedUtilityB{1})/sqrt(size(realizedUtilityB{2},1));

Nsim = size(realizedUtilityB{2},1);
totalWelfareLoss = 1-exp(mean(realizedUtilityB{2}).*(1-lambdaGrid')-mean(realizedUtilityB{1}).*(1-lambdaGrid')+...
                          mean(realizedUtilityA{2}).*lambdaGrid'-mean(realizedUtilityA{1}).*lambdaGrid');
stdError = exp(mean(realizedUtilityB{2}).*(1-lambdaGrid')-mean(realizedUtilityB{1}).*(1-lambdaGrid')+...
                          mean(realizedUtilityA{2}).*lambdaGrid'-mean(realizedUtilityA{1}).*lambdaGrid').*...
                          std(realizedUtilityB{2}.*(1-repmat(lambdaGrid',[Nsim,1]))-mean(realizedUtilityB{1}).*(1-repmat(lambdaGrid',[Nsim,1]))+...
                          mean(realizedUtilityA{2}).*repmat(lambdaGrid',[Nsim,1])-mean(realizedUtilityA{1}).*repmat(lambdaGrid',[Nsim,1]))/sqrt(size(realizedUtilityB{2},1));
smoothPar = 1/(1+(abs(lambdaGrid(2)-lambdaGrid(1)))^3/0.6);
if smoothing % smoothes the curves and brings them to 0 
     %totalWelfareLoss = csaps(lambdaGrid,totalWelfareLoss,0.999,lambdaGrid);
     %totalWelfareLoss = totalWelfareLoss-min(totalWelfareLoss); % bring minimum to 0, which can be different due to noise
     welfareLossA = csaps(lambdaGrid,welfareLossA,smoothPar,lambdaGrid);
     welfareLossA = welfareLossA-min(welfareLossA);
     welfareLossB = csaps(lambdaGrid,welfareLossB,smoothPar,lambdaGrid);
     welfareLossB = welfareLossB-min(welfareLossB);
     totalWelfareLoss = welfareLossA.*lambdaGrid+welfareLossB.*(1-lambdaGrid);
end

figure
hold on
plot(lambdaGrid,-welfareLossA*100,'-.','color',my_color('pale_blue'),'linewidth',2)
plot(lambdaGrid,-welfareLossB*100,':','color',my_color('green'),'linewidth',2)
plot(lambdaGrid, -totalWelfareLoss*100,'-','color',my_color('maroon'),'linewidth',2)
if confBounds
    plot(lambdaGrid, (totalWelfareLoss+1.96*stdError)*100,':','color',my_color('pale_blue'),'linewidth',1)
    plot(lambdaGrid, (totalWelfareLoss-1.96*stdError)*100,':','color',my_color('pale_blue'),'linewidth',1)
    plot(lambdaGrid,(welfareLossA+1.96*stdErrorA)*100,':','color',my_color('pale_blue'),'linewidth',1)
    plot(lambdaGrid,(welfareLossA-1.96*stdErrorA)*100,':','color',my_color('pale_blue'),'linewidth',1)
    plot(lambdaGrid,(welfareLossB+1.96*stdErrorB)*100,':','color',my_color('pale_blue'),'linewidth',1)
    plot(lambdaGrid,(welfareLossB-1.96*stdErrorB)*100,':','color',my_color('pale_blue'),'linewidth',1)
end
xlabel(xsymbol,'interpreter','latex')
ylabel('Fractional C.E. gain (\%)','interpreter','latex')
box on
legend('Agent A', 'Agent B','Total','location','southwest')

