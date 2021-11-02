function pub_fig2_Distrib_TaskMod_Bootstrap_JCfun(zSMA, Tcombo_z, iter, parfig)
% function pub_fig2_Distrib_TaskMod_Bootstrap_JCfun(Tcombo_z,iter,parfig)
% OUTPUT: % figures: 3 curves with confidence interval (95%)  
% INPUTS: % Tcombo_z = table obtained using : pub_table6_Tcombo_JCscript;
% iter is number of iteration for the bootstrap. 
% parfig is a structure that need: Ztrh, pre, post(e.g. parfig.Zthr=3)  
% written by Julien Catanese 01/01/2019
% last updates: 03/31/2019

close all, figure,
legend_all = [];
% parfig.Zthr = 3 ;

allbool=  Tcombo_z.VMVL & Tcombo_z.Opto_post_sess & Tcombo_z.z_exct;
ncellmax = sum(allbool)
szSMA = zSMA(allbool,:);
groupsize = 30; % 30 cells per group (randomly selected among the all pool (allbool)
Yall = [];
for ii=1:iter
    idx_rand=[]; BarSig=[]; Yconv=[];
    
    idx_rand = randi(ncellmax, 1, groupsize);
    GrpszSMA = szSMA(idx_rand,:);
    BarSig= zeros(size(GrpszSMA)); nsc=size(GrpszSMA,1);
    idxB=find(GrpszSMA>parfig.Zthr);
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


%
bool_ori = [Tcombo_z.VMVL & Tcombo_z.Opto_post_sess & Tcombo_z.z_exct & ~(Tcombo_z.Opto_inib | Tcombo_z.Opto_exct), ...
            Tcombo_z.VMVL & Tcombo_z.Opto_post_sess & Tcombo_z.z_exct & Tcombo_z.Opto_inib, ...
            Tcombo_z.VMVL & Tcombo_z.Opto_post_sess & Tcombo_z.z_exct & Tcombo_z.Opto_exct ]; % nonopto=237 % opto-=64 % opto+=37
GroupID = {'conf int limit (95%)';...
    ['Thal NOL cells (' num2str(sum(Tcombo_z.VMVL & Tcombo_z.Opto_post_sess & Tcombo_z.z_exct & ~(Tcombo_z.Opto_inib | Tcombo_z.Opto_exct))) ')'];...
    ['SvTh NEG cells (0' num2str(num2str(sum(Tcombo_z.VMVL & Tcombo_z.Opto_post_sess & Tcombo_z.z_exct & Tcombo_z.Opto_inib))) ')' ];...
    ['SvTh POS cells (0' num2str(num2str(sum(Tcombo_z.VMVL & Tcombo_z.Opto_post_sess & Tcombo_z.z_exct & Tcombo_z.Opto_exct))) ')' ]};

szSMA=[]; 
for Group = 1:size(bool_ori,2)
    szSMA= zSMA(bool_ori(:,Group),:);
    BarSig= zeros(size(szSMA)); nsc=size(szSMA,1);
    idxB=find(szSMA>parfig.Zthr); BarSig(idxB)= 1;
    
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
legend(GroupID) ;
