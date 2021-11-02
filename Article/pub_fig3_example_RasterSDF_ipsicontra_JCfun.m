function pub_fig_example_RasterSDF_ipsicontra_JCfun(Tfig3_RvL, SDF, listcell, parfig)
% function pub_fig_example_RasterSDF_ipsicontra_JCfun(Tcombo, SDF, listcell, parfig)
% plot Raster and Spike density function for single units (SUA). 
% Screen using Zscore for example of each type. 
% Written by Julien Catanese 3/15/2019
% last updated 4/2/19 JC

%% during testing (then delete)
load('Tcombo.mat')
% SDF=dzSMA;


%% Identify each Type (using Ttest (FrBLE vs FrEpoch)) 
type0      = Tcombo.VMVL & Tfig3_RvL.dzRL_nosig; 
type1_exct = Tcombo.VMVL & Tfig3_RvL.dzRL_exct_UniMod_1puf;
type2_exct = Tcombo.VMVL & Tfig3_RvL.dzRL_exct_UniMod_2del;
type3_exct = Tcombo.VMVL & Tfig3_RvL.dzRL_exct_UniMod_3res;   
type4_exct = Tcombo.VMVL & Tfig3_RvL.dzRL_exct_biMod_12;
type5_exct = Tcombo.VMVL & Tfig3_RvL.dzRL_exct_biMod_23;
type6_exct = Tcombo.VMVL & Tfig3_RvL.dzRL_exct_biMod_13;
type7_exct = Tcombo.VMVL & Tfig3_RvL.dzRL_exct_triMod_123; 


Ntype0 = sum(type0)
Ntype1_exct = sum(type1_exct)
Ntype2_exct = sum(type2_exct)
Ntype3_exct = sum(type3_exct)
Ntype4_exct = sum(type4_exct)
Ntype5_exct = sum(type5_exct)
Ntype6_exct = sum(type6_exct)
Ntype7_exct = sum(type7_exct)


type0     = Tcombo.VMVL & Tfig3_RvL.dzRL_nosig; 
type1_inib = Tcombo.VMVL & Tfig3_RvL.dzRL_inib_UniMod_1puf;
type2_inib = Tcombo.VMVL & Tfig3_RvL.dzRL_inib_UniMod_2del;
type3_inib = Tcombo.VMVL & Tfig3_RvL.dzRL_inib_UniMod_3res;   
type4_inib = Tcombo.VMVL & Tfig3_RvL.dzRL_inib_biMod_12;
type5_inib = Tcombo.VMVL & Tfig3_RvL.dzRL_inib_biMod_23;
type6_inib = Tcombo.VMVL & Tfig3_RvL.dzRL_inib_biMod_13;
type7_inib = Tcombo.VMVL & Tfig3_RvL.dzRL_inib_triMod_123; 


Ntype0 = sum(type0)
Ntype1_inib = sum(type1_inib)
Ntype2_inib = sum(type2_inib)
Ntype3_inib = sum(type3_inib)
Ntype4_inib = sum(type4_inib)
Ntype5_inib = sum(type5_inib)
Ntype6_inib = sum(type6_inib)
Ntype7_inib = sum(type7_inib)

figure,
bar([Ntype1_exct Ntype1_inib; Ntype2_exct Ntype2_inib; Ntype3_exct Ntype3_inib;...
    Ntype4_exct Ntype4_inib; Ntype5_exct Ntype5_inib; Ntype6_exct Ntype6_inib; ...
    Ntype7_exct Ntype7_inib])

%% TEST  

figure,
listc1= listcell('vgat14w14d7S1Ch3clu#03', :)  
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

%% Type1 PUFF
% Screening CONTRA
parfig.title=['Type1 Cont (#' num2str(Ntype1_exct) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type1_exct,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type1_exct, SDF, parfig)
% Screening IPSI
parfig.title=['Type1 ipsi (#' num2str(Ntype1_inib) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type1_inib,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type1_inib, SDF, parfig)

pause(2)
% example cells to plot
% listc1= listcell('vgat17w10d8S4Ch6clu#02', :)
% see example from Figure2 BiMod


%% Type2 DELAY
% Screening CONTRA
parfig.title=['Type2 Cont (#' num2str(Ntype1_exct) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type2_exct,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type2_exct, SDF, parfig)
% Screening IPSI
parfig.title=['Type2 ipsi (#' num2str(Ntype1_inib) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type2_inib,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type2_inib, SDF, parfig)

% example cells to plot
listc1= listcell('vgat17w10d7S1Ch8clu#01', :)
listc1= listcell('vgat14w14d7S1Ch3clu#03', :) 
pause(2)

%% Type 3 RESP 
% Screening
parfig.title=['Type3 Cont (#' num2str(Ntype3_exct) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type3_exct,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type3_exct, SDF, parfig)
%  Screening IPSI
parfig.title=['Type3 ipsi (#' num2str(Ntype3_inib) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type3_inib,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type3_inib, SDF, parfig)
% example cells to plot
listc1= listcell('vgat17w10d8S2Ch5clu#01', :);
listc1= listcell('vgat17w10d4S4Ch8clu#01', :);
pause(2)

%% Type 4 (puff+delay)
% Screening
parfig.title=['Type4 Cont (#' num2str(Ntype4_exct) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type4_exct,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type4_exct, SDF, parfig)
% example cells to plot
listc1= listcell('vgat17w10d3S3Ch5clu#04', :); 
listc1= listcell('vgat15w10d8S2Ch2clu#02', :);
listc1= listcell('vgat15w10d7S1Ch2clu#01', :); 
pause(2)

%% Type 5  (RAMP)
% Screening
parfig.title=['Type5 Cont(#' num2str(Ntype5_exct) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type5_exct,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type5_exct, SDF, parfig)
% example cells to plot
listc1= listcell('vgat17w10d3S2Ch5clu#01', :) 
listc1= listcell('vgat17w10d3S2Ch5clu#01', :) 
pause(2)

%% Type 6  (BI-MOD)
% Screening
parfig.title=['Type6 Cont(#' num2str(Ntype6_exct) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type6_exct,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type6_exct, SDF, parfig)
% example cells to plot
listc1= listcell('vgat15w10d3S3Ch4clu#01', :); 
listc1= listcell('vgat17w10d3S4Ch1clu#03', :) 
pause(2)

%% Type 7 (TRI-MOD)
% Screening
parfig.title=['Type7 Cont(#' num2str(Ntype7_exct) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type7_exct,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type7_exct, SDF, parfig)
% example cells to plot


%% Type 0  
% Screening
parfig.title=['Type0 (#' num2str(Ntype0) ')' ]
pub_fig_SMA_SUA_JCfun2(Tcombo(type0,:), SDF, parfig)
pub_fig_SMA_MUA_JCfun(type0, SDF, parfig)
% example cells to plot
listc1= listcell('vgat17w10d3S2Ch5clu#01', :) 

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

