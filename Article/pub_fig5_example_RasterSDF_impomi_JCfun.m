function pub_fig_example_RasterSDF_ipsicontra_JCfun(Tfig5_imp,dzSMA_imp, dzSMA_omi, listcell, parfig)
% function pub_fig_example_RasterSDF_ipsicontra_JCfun(Tcombo, SDF, listcell, parfig)
% plot Raster and Spike density function for single units (SUA). 
% Screen using Zscore for example of each type. 
% Written by Julien Catanese 4/15/2019

%% during testing (then delete)
load('Tcombo.mat')
SDF=dzSMA_imp;

imp_all =  Tfig5_imp.dzic_exctAll  & Tfig5_imp.VMVL ;  Nimp_all = sum(imp_all)
omi_all = Tfig5_imp.dzoc_inibAll & Tfig5_imp.VMVL ;  Nomi_all = sum(omi_all)

imp = Tfig5_imp.dzic_exct & Tfig5_imp.VMVL; Nimp= sum( imp)
omi = Tfig5_imp.dzoc_inib & Tfig5_imp.VMVL ;  Nomi= sum( omi)

imp_puf = Tfig5_imp.dzic_exct_pufAll & Tfig5_imp.VMVL; Nimp_puf = sum(imp_puf )
omi_puf = Tfig5_imp.dzoc_inib_pufAll & Tfig5_imp.VMVL; Nomi_puf = sum(omi_puf)

imp_del = Tfig5_imp.dzic_exct_delAll & Tfig5_imp.VMVL; Nimp_del = sum( imp_del )
omi_del = Tfig5_imp.dzoc_inib_delAll & Tfig5_imp.VMVL;Nomi_del= sum(omi_del)

imp_res = Tfig5_imp.dzic_exct_resAll & Tfig5_imp.VMVL;  Nimp_res = sum(imp_res )
omi_res = Tfig5_imp.dzoc_inib_resAll & Tfig5_imp.VMVL; Nomi_res = sum(omi_res )

%% TEST  

LT=Tfig5_imp(imp & omi,:); 
nn= size(LT,1)
for nc = 1:nn;
RowName = LT.Row(nc);
listc1= listcell(RowName{1}, :)  
figure, hold on, 
% listc1= listcell('vgat12w11d5S3Ch6clu#01', :)  
% listc1= listcell('vgat14w14d2S1Ch5clu#01', :)  
% listc1= listcell('vgat12w11d5S4Ch1clu#02', :)  
% listc1= listcell('vgat17w10d8S4Ch3clu#03', :)  

% the 3 final example 
listc1= listcell('vgat12w11d5S3Ch5clu#02', :)  
listc1= listcell('vgat14w14d8S3Ch2clu#01', :)  
listc1= listcell('vgat17w10d4S2Ch3clu#01', :)  
%Imp
parfig.col = parfig.parfigimp.col 
parfig.trial_type = parfig.parfigimp.trial_type
parfig.currentpanelNb = 1; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
% omi
parfig.col = parfig.parfigomi.col 
parfig.trial_type = parfig.parfigomi.trial_type
parfig.currentpanelNb = 2; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
disp('done Raster')
% cor
parfig.col = {'k'} 
parfig.trial_type = {'cor'}
parfig.currentpanelNb = 3; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
disp('done Raster')

%% IMP ALL (Screening dzicSMA)
SDF = dzSMA_imp; parfig.title=['IMP all (#' num2str(Nimp_all) ')' ]; 
pub_fig_SMA_SUA_JCfun2(Tfig5_imp(imp_all,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(imp_all, SDF, parfig)

%% OMI ALL (Screening dzocSMA)
SDF = dzSMA_omi; parfig.title=['OMI all (#' num2str(Nomi_all) ')' ]; 
pub_fig_SMA_SUA_JCfun2(Tfig5_imp(omi_all,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(omi_all, SDF, parfig)

%% IMP Type1 PUFF 
% Screening IMP
SDF = dzSMA_imp; parfig.title=['IMP PUF (#' num2str(Nimp_puf) ')' ]; 
pub_fig_SMA_SUA_JCfun2(Tcombo(imp_puf,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(imp_puf, SDF, parfig)

%% OMI Type1 PUFF 
SDF = dzSMA_omi; parfig.title=['OMI PUF (#' num2str(Nomi_puf) ')' ]; 
pub_fig_SMA_SUA_JCfun2(Tcombo(omi_puf,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(omi_puf, SDF, parfig)

%% IMP Type2 DEL 
% Screening IMP
SDF = dzSMA_imp; parfig.title=['IMP DEL (#' num2str(Nimp_del) ')' ]; 
pub_fig_SMA_SUA_JCfun2(Tcombo(imp_del,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(imp_del, SDF, parfig)

%% OMI Type2 DEL 
SDF = dzSMA_omi; parfig.title=['OMI DEL (#' num2str(Nomi_del) ')' ]; 
pub_fig_SMA_SUA_JCfun2(Tcombo(omi_del,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(omi_del, SDF, parfig)

%% IMP Type3 RES
% Screening IMP
SDF = dzSMA_imp; parfig.title=['IMP RES (#' num2str(Nimp_res) ')' ]; 
pub_fig_SMA_SUA_JCfun2(Tcombo(imp_res,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(imp_res, SDF, parfig)

%% OMI Type3 RES 
SDF = dzSMA_omi; parfig.title=['OMI RES (#' num2str(Nomi_res) ')' ]; 
pub_fig_SMA_SUA_JCfun2(Tcombo(omi_res,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(omi_res, SDF, parfig)

%% Example from Fig2 
%Type1 PUFF 
listc1 = listcell('vgat15w10d3S1Ch4clu#01', :)
figure,
%contra
parfig.col = parfig.parfigcont.col 
parfig.trial_type = parfig.parfigcont.trial_type
parfig.currentpanelNb = 1; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%ipsi
parfig.col = parfig.parfigipsi.col 
parfig.trial_type = parfig.parfigipsi.trial_type
parfig.currentpanelNb = 2; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%Type2 DELAY
listc1= listcell('vgat17w10d4S3Ch4clu#01', :);
figure,
%contra
parfig.col = parfig.parfigcont.col 
parfig.trial_type = parfig.parfigcont.trial_type
parfig.currentpanelNb = 1; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%ipsi
parfig.col = parfig.parfigipsi.col 
parfig.trial_type = parfig.parfigipsi.trial_type
parfig.currentpanelNb = 2; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%Type 3 RESP 
listc1= listcell('vgat11w10d4S4Ch6clu#01', :)
figure,
%contra
parfig.col = parfig.parfigcont.col 
parfig.trial_type = parfig.parfigcont.trial_type
parfig.currentpanelNb = 1; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%ipsi
parfig.col = parfig.parfigipsi.col 
parfig.trial_type = parfig.parfigipsi.trial_type
parfig.currentpanelNb = 2; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%Type 4 (puff+delay)
listc1= listcell('vgat15w10d7S1Ch2clu#01', :)
figure,
%contra
parfig.col = parfig.parfigcont.col 
parfig.trial_type = parfig.parfigcont.trial_type
parfig.currentpanelNb = 1; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%ipsi
parfig.col = parfig.parfigipsi.col 
parfig.trial_type = parfig.parfigipsi.trial_type
parfig.currentpanelNb = 2; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%Type 5  (RAMP)
listc1= listcell('vgat12w11d5S3Ch6clu#02', :)
figure,
%contra
parfig.col = parfig.parfigcont.col 
parfig.trial_type = parfig.parfigcont.trial_type
parfig.currentpanelNb = 1; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%ipsi
parfig.col = parfig.parfigipsi.col 
parfig.trial_type = parfig.parfigipsi.trial_type
parfig.currentpanelNb = 2; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

%  Type 6  (BI-MOD)
listc1= listcell('vgat15w10d3S4Ch5clu#01', :)
figure,
%contra
parfig.col = parfig.parfigcont.col 
parfig.trial_type = parfig.parfigcont.trial_type
parfig.currentpanelNb = 1; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%ipsi
parfig.col = parfig.parfigipsi.col 
parfig.trial_type = parfig.parfigipsi.trial_type
parfig.currentpanelNb = 2; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

% Type 7 (TRI-MOD)
listc1= listcell('vgat15w10d7S1Ch1clu#02', :)
figure,
%contra
parfig.col = parfig.parfigcont.col 
parfig.trial_type = parfig.parfigcont.trial_type
parfig.currentpanelNb = 1; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%ipsi
parfig.col = parfig.parfigipsi.col 
parfig.trial_type = parfig.parfigipsi.trial_type
parfig.currentpanelNb = 2; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

listc1= listcell('vgat17w10d7S2Ch6clu#01', :)
figure,
%contra
parfig.col = parfig.parfigcont.col 
parfig.trial_type = parfig.parfigcont.trial_type
parfig.currentpanelNb = 1; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%ipsi
parfig.col = parfig.parfigipsi.col 
parfig.trial_type = parfig.parfigipsi.trial_type
parfig.currentpanelNb = 2; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

% Type 0  
listc1= listcell('vgat11w10d4S3Ch5clu#03',:)
figure,
%contra
parfig.col = parfig.parfigcont.col 
parfig.trial_type = parfig.parfigcont.trial_type
parfig.currentpanelNb = 1; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
%ipsi
parfig.col = parfig.parfigipsi.col 
parfig.trial_type = parfig.parfigipsi.trial_type
parfig.currentpanelNb = 2; 
pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);

