% Get_Tables_JCscript
% written by Julien Catanese 2019 in JaegerLAB 
% last updates 5/10/2020 by Juien Catanese

%% Table1 : list of recorded cells + anatomy + opto Y/N
% pub_table1_listcell_JCscript
% load('listcell')
% idx_list=[]; idx_list(:,1)=1:size(Tfig1_VMopto,1);
% Tfig1_VMopto = addvars(Tfig1_VMopto, idx_list, 'Before', 1);Tfig1_VMopto(1,:)
% save('Tfig1_VMopto.mat','Tfig1_VMopto')
% save('listcell', 'Tfig1_VMopto' ,'-append')


%% Table2 : Ephys (correct trials) Firing rate, resp types based on zscore(Tephys_z) or ttest(Tephys_tt)
% clearvars -except mypath parfig
% load('listcell.mat'); load('Tcombo.mat');
parfig.plot=0;
parfig.saveSMA=1;
parfig.saveTABLE =1;
parfig.center_evt = 'GoCue';    % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.BaselineEpoch= [150:2150];
parfig.pre = parfig.BaselineEpoch(end) + 1500 % define how much time before zero
parfig.post = 1500 % define how much time after zero (zero will be defined by trigtimes.cor)
parfig.trial_type = {'cor'};
parfig.col = {'k'};
% pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig)
%
% load('SMA_cor_GoCue545.mat');
% pub_table2_Tfig2_cor_JCscript; close all;
%
% idx_list=[]; idx_list(:,1)=1:size(Tfig2_cor,1);
% Tfig2_cor = addvars(Tfig2_cor, idx_list, 'Before', 1);Tfig2_cor(1,:)
% save('Tfig2_cor.mat','Tfig2_cor')
% save('listcell', 'Tfig2_cor','-append')

%% Table3 : Ephys (Zipsi-Zcontra) resp types based on zscore(Tephys_z)
% disp('Start Table3: ipsi-contra')
% clearvars -except mypath parfig
% load('listcell.mat'); load('Tcombo.mat');
% parfig.trial_type = {'cCL'};
% parfig.col = {'r'};
% [SMA_CL FRepoch_CL ttestEpoch_CL] = pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig);
% parfig.trial_type = {'cCR'};
% parfig.col = {'b'};
% [SMA_CR FRepoch_CR ttestEpoch_CR] = pub_comp_ephySMA_task_ttest_JCfun(listcell, parfig);

% pub_table3_Tfig3_RvL_JCscript; close all;
%
% idx_list=[]; idx_list(:,1)=1:size(Tfig3_RvL,1);
% Tfig3_RvL = addvars(Tfig3_RvL, idx_list, 'Before', 1);Tfig3_RvL(1,:)
% save('Tfig3_RvL.mat','Tfig3_RvL')
% save('listcell', 'Tfig3_RvL' , '-append')

%% Table4: Tcombo
% disp('Start Table4: Tcombo')
% load('listcell2')
% pub_table4_Tcombo_JCscript

%% Table5: Ephys (impulse-omissions trials) to define dzic=zImp-zCor and dzoc=zOmi-zCor
disp('Start Table5: impulse omission')
clearvars -except mypath parfig
load('listcell.mat'); load('Tcombo.mat');
parfig.trial_type = {'imp'};
parfig.col = {'m'};
% sub_listcell = listcell(logical(sum(Tcombo.nSess == [4,7,10,13,15],2)),:)
sub_listcell = listcell(logical(Tfig2_cor.Opto_post_sess),:);
[SMA_imp FRepoch_imp ttestEpoch_imp] = pub_comp_ephySMA_task_ttest_JCfun2(sub_listcell, parfig);


parfig.trial_type =  {'omi'};
parfig.col = {'c'};
% sub_listcell = listcell(logical(sum(Tcombo.nSess == [4,7,10,13,15],2)),:)
sub_listcell = listcell(logical(Tcombo.Opto_post_sess),:);
[SMA_omi FRepoch_omi ttestEpoch_omi] = pub_comp_ephySMA_task_ttest_JCfun(sub_listcell, parfig);
%%
load('listcell.mat', 'listcell', 'Tfig1_VMopto'); load('Tcombo.mat');
%%
sub_listcell = listcell(logical(Tcombo.Opto_post_sess),:);
pub_table5_Tfig5_imp_JCscript;

idx_list=[]; idx_list(:,1)=1:size(Tfig5_imp,1);
Tfig5_imp = addvars(Tfig5_imp, idx_list, 'Before', 1);Tfig5_imp(1,:)
% save('Tfig5_imp.mat','Tfig5_imp','Zthr');
% save('listcell2', 'Tfig5_imp' , '-append')

%% ALIGN TO LICK and GET:  SMA_imp_Licks360.mat
% ATTENTION THE EPOCH ARE NOT USABLE FOR LICK ALIGN so ignore all the TTest
% and other epoch related analysis in the SMA_lick
parfig.saveSMA=1;  parfig.plot=0;
parfig.BaselineEpoch= [150:2150];
parfig.pre = parfig.BaselineEpoch(end) + 1500 % define how much time before zero
parfig.post = 1500 % define how much time after zero (zero will be defined by trigtimes.cor)

sub_listcell = listcell(logical(Tcombo.Opto_post_sess),:);

parfig.center_evt = 'Licks';    % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.trial_type = {'imp'};
pub_comp_ephySMA_task_ttest_JCfun(sub_listcell, parfig)

parfig.center_evt = 'Licks';    % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.trial_type = {'cor'};
pub_comp_ephySMA_task_ttest_JCfun(sub_listcell, parfig)
% ATTENTION THE EPOCH ARE NOT USABLE FOR LICK ALIGN so ignore all the TTest
% and other epoch related analysis in The SMA_lick
