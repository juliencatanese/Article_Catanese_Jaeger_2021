function pub_fig3_Distrib_Bootstrap_JCfun(dzSMA, group_ID, group_bool, parfig)
% function pub_fig2_Distrib_TaskMod_Bootstrap_JCfun(Tcombo_z,iter,parfig)
% OUTPUT: % figures: 3 curves with confidence interval (95%)  
% INPUTS: % Tcombo_z = table obtained using : pub_table6_Tcombo_JCscript;
% iter is number of iteration for the bootstrap. 
% parfig is a structure that need: Ztrh, pre, post(e.g. parfig.Zthr=3)  
% written by Julien Catanese 01/01/2019
% last updates: 03/31/2019

figure,
ncellmax = sum(sum(group_bool))
sdzSMA = dzSMA(logical(sum(group_bool,2)),:);
groupsize = parfig.bootstrap.groupsize; % 30-40 cells per group (randomly selected among the all pool (allbool)
iter= parfig.bootstrap.iter % number of iteration for Boostrap
Yall = [];
for ii=1:iter
    idx_rand=[]; BarSig=[]; Yconv=[];
    
    idx_rand = randi(ncellmax, 1, groupsize);
    GrpsdzSMA = sdzSMA(idx_rand,:);
    BarSig= zeros(size(GrpsdzSMA)); nsc=size(GrpsdzSMA,1);
    idxB=find(GrpsdzSMA>parfig.Zthr);
    BarSig(idxB)= 1;
    parfig.plot = 0;
    parfig.plotEpochLines = 0;
    hold on;
    [Yconv ] = pub_fig2_distrib_percZthr_JCfun(BarSig, parfig);
    Yall(ii,:) = Yconv;
    
end
X=[-parfig.pre:1:parfig.post];
figure,
plotshaded(X, [mean(Yall)-(2*std(Yall)) ; mean(Yall)+(2*std(Yall))],  'k');
xlim(parfig.xlim);
% ylim(parfig.ylim);
ylabel('% cell modulated Above thr');
xlabel('time from GoCue (ms)');

sdzSMA=[]; 
for Group = 1:size(group_bool,2)
    sdzSMA= dzSMA(group_bool(:,Group),:);
    BarSig= zeros(size(sdzSMA)); nsc=size(sdzSMA,1);
    idxB=find(sdzSMA>parfig.Zthr); BarSig(idxB)= 1;
    
    parfig.plot=1;
    hold on,
    [Yconv ]=pub_fig2_distrib_percZthr_JCfun(BarSig, parfig);
    Yall = Yconv
    
end

parfig.xlim = [0-2250 0+750];
parfig.ylim = [0  100]
parfig.title = ['Distrib of Modulation (Zthr >' num2str(parfig.Zthr) ') trials-' parfig.trial_type{1} ];

pre=0;
if parfig.center_evt == 'GoCue'
    hold on, plot([pre-1500 pre-1500], ylim,'k--','LineWidth',2.5);
    hold on, plot([pre-750 pre-750], ylim,'k--', 'LineWidth',2.5);
    hold on, plot([pre pre], ylim,'k--','LineWidth',2.5);
elseif parfig.center_evt == 'Delay'
    hold on, plot([pre-750 pre-750], ylim,'k--', 'LineWidth',2.5);
    hold on, plot([pre pre], ylim,'k--','LineWidth',2.5);
    hold on, plot([pre+750 pre+750], ylim,'k--','LineWidth',2.5);
end
legend(group_ID, 'Location', 'northwest') ;
