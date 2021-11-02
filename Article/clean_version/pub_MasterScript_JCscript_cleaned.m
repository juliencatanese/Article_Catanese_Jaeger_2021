% pub_MasterScript_JCscript_cleaned
% by JC 05/10/2020

mypath = 'D:\DATA EMORY\JC_Analysis'
cd(mypath)

%% Get all events and task trials
pub_Behavior_MasterScript_JCscript
% cd(mypath);

%% Preprocess data: Create Tables and SMA    
Get_Tables_JCscript

%% Figure1: BEHAVIOR
pub_fig1_Behavior_JCscript;

%% Figure2: Ephys during corr trials (SUA and pop)
close all, clear all,
load('listcell2.mat'); load('Tcombo');
load('SMA_cor_GoCue545.mat');
parfig=parfigUsed;
parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.trial_type = {'cor'};
parfig.col = {'k'};
parfig.xlim = [-2750 1250];

zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); zSMA(find(zSMA==inf))=NaN;
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSA,2); zSMA(find(zSMA==inf))=NaN;
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSA,2); zSMA(find(zSMA==inf))=NaN;


%% Fig 2A Venn Diagrams (2 types: +/-/complex/nosig and 7celltypes)
%colormap ( %'jet'   %  'flag' % 'copper'  % 'cool' % 'hot' % 'colorcube'   % 'bone'  % 'autumn' )
colmap = 'bone'
resolution = 0.05
pub_fig2_venndiagrams_zcor_JCfun(Tfig2_cor, colmap, resolution)

%% Fig 2G Matrice Populations
zSMA_VMexct= zSMA(logical(Tfig2_cor.VMVL) & Tfig2_cor.z_exct,:);
ZorTT= 'z_thr'
parfig.colormap = 'jet';
parfig.caxis = [-1 1];
parfig.sort_Xepoch = [parfig.pre-1500 : parfig.pre+1000];
parfig.trial_type = parfig.trial_type;
parfig.sort_variable = 'peak';
pub_fig2_MatricePop_JCfun(ZorTT, zSMA_VMexct, Tcombo, parfig)
% saveas()
parfig.sort_variable = 'Z-tr'
parfig.caxis = [0 9];
parfig.sort_Xepoch = [parfig.pre-750 : parfig.pre+1000];
pub_fig2_MatricePop_JCfun(ZorTT, zSMA_VMexct, Tcombo, parfig)
% saveas()

%% Fig 3R : Population firing rate
close all,
VMVL_all = Tcombo.VMVL;
Exct_only = Tcombo.z_exct & Tcombo.VMVL;
Inib_only = Tcombo.z_inib  & Tcombo.VMVL;
complex_only = Tcombo.z_complex  & Tcombo.VMVL;
nosig_only = Tcombo.z_nosig  & Tcombo.VMVL;
idx_list = logical([VMVL_all, Exct_only, Inib_only, complex_only, nosig_only ]);
parfig.ylabel = 'Zscore'
parfig.title = 'grand average'
parfig.col = {'k','r','g','c','y'}
parfig.plotshaded='sem'
parfig.k=2;
% parfig.plotshaded='std'
parfig.xlim = [0-2250 0+750];
parfig.ylim= [-5 7];
pub_fig_SMA_MUA_1by1_JCfun(idx_list, zSMA, parfig)
pub_fig_SMA_MUA_JCfun(logical(idx_list), zSMA, parfig)
legend ('VMVL all', ['(n=' num2str(sum(VMVL_all)) ')'],...
    'Exct only', ['(n=' num2str(sum(Exct_only)) ')'],...
    'inib only',  ['(n=' num2str(sum(Inib_only)) ')'], ...
    'complex', ['(n=' num2str(sum(complex_only)) ')'],...
    'nosig', ['(n=' num2str(sum(nosig_only)) ')'],...
    'Location', 'Best')

%% SUPP FIG S3
close all
pub_figS3_Supp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4 EPYS IPSI CONTRA (opto-Post)
close all, clearvars -except mypath
load('listcell2.mat'); load('Tcombo.mat')
parfig.WorkFolder = mypath;
parfig.saveSMA=0;
parfig.saveTABLE=0;
parfig.center_evt =  'GoCue';
parfig.BaselineEpoch= [150:2150];
parfig.pre = parfig.BaselineEpoch(end) + 1500 % define how much time before zero
parfig.post = 1500 % define how much time after zero (zero will be defined by trigtimes.cor)
parfig.xlim = [-2750 1250];
parfig.merge=1;
for iside = 1:2
    % IPSI
    if iside == 1
        load('SMA_cCL_GoCue545.mat');
        parfigipsi.trial_type = {'cCL'};
        parfigipsi.col = {'r'};
        zSMAipsi = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
        SMAipsi= SMA; SSemAipsi=SSemA; clear SMA SSemA SSA;
        % CONTRA
    elseif iside == 2
        load('SMA_cCR_GoCue545.mat');
        parfigcont.trial_type = {'cCR'};
        parfigcont.col = {'b'};
        zSMAcont = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
        SMAcont= SMA; SSemAcont=SSemA; clear SMA SSemA SSA;
    end
end
%
Zthr=3
dzSMA = zSMAcont-zSMAipsi;
contra_dz = Tcombo.z_exct & Tcombo.VMVL & (Tfig3_RvL.dzRL_exctAll) % all Contra
ipsi_dz = Tcombo.z_exct & Tcombo.VMVL & (Tfig3_RvL.dzRL_inibAll); % all Ipsi

%% Fig 4A Venn Diagrams (2 types: +/-/complex/nosig and 7celltypes)
close all, colmap = 'jet';  resolution = 0.05;
Ncomplex_dz = sum(contra_dz & ipsi_dz)
data= [sum(contra_dz)-Ncomplex_dz , Ncomplex_dz, sum(ipsi_dz) - Ncomplex_dz]
vennX(data, resolution, colmap); colormapeditor

%% Fig 4B Bar Graphs + stat
% contra_dz include Ncomplex
pub_fig4_Bar_stat_ipsicontra_Epoch_JCscript

%% Fig 4C-H : Examples SUA => SDF-RASTER for each Type (from I to VII)
close all
parfig.saveSMA=0;
parfig.saveTABLE=0;
parfig.plotmerge=1;
parfig.nRaster = 2;
parfig.xlim= [-1750 1000];
parfig.plot=1;
parfig.k=2;
parfig.ylabel = 'z';
parfig.plotshaded= 'sem';
parfig.parfigcont=parfigcont;
parfig.parfigipsi = parfigipsi;
dzSMA=  zSMAcont-zSMAipsi;
pub_fig3_example_RasterSDF_ipsicontra_JCfun(Tfig3_RvL, dzSMA, listcell, parfig)
% saveas(gcf, ['D:\JC_Figures\SDF_RASTER_cor\SDF_RASTER_cor_listcell_ncell#' num2str(ncell) ],'jpg')

%% Fig 4I-J Mean DZ
Xt=[1:1:size(dzSMA,2)]-parfig.pre;
parfig.xlim=[-2000 750]
K=2;

figure,
hold on, plot(Xt, mean(dzSMA(contra_dz,:)))
hold on, plot(Xt, mean(dzSMA(ipsi_dz,:)))

pub_fig3_LvR_shaded_dzSMA_JCscript

%% Fig 4K-L Distribution Bootstrap Opto
% contra (Fig 3K)
parfig.Zthr = 3;
parfig.bootstrap.iter = 500
parfig.bootstrap.groupsize=30
parfig.title = 'CONTRA Dist bootstrap'; parfig.ylim=[0 70];
parfig.trial_type = {'cCR'};
group_bool= [contra_dz & Tcombo.Opto_post_sess &  ~(Tcombo.Opto_inib | Tcombo.Opto_exct), ...
    contra_dz & Tcombo.Opto_post_sess & Tcombo.Opto_inib, ...
    contra_dz & Tcombo.Opto_post_sess & Tcombo.Opto_exct ];
group_ID = {'conf int limit (3std)';...
    ['Thal NOL cells (' num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & ~(Tcombo.Opto_inib | Tcombo.Opto_exct))) ')'];...
    ['SvTh NEG cells (0' num2str(num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & Tcombo.Opto_inib))) ')' ];...
    ['SvTh POS cells (0' num2str(num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & Tcombo.Opto_exct))) ')' ]};
pub_fig3_Distrib_Bootstrap_JCfun(dzSMA, group_ID, group_bool, parfig)
xlim(parfig.xlim);

% ipsi (Fig 3L)
parfig.Zthr = 3;
parfig.bootstrap.iter = 500
parfig.bootstrap.groupsize=30
parfig.title = ' IPSI Dist bootstrap'; parfig.ylim=[0 70];
parfig.trial_type = {'cCL'};
group_bool= [ipsi_dz & Tcombo.Opto_post_sess &  ~(Tcombo.Opto_inib | Tcombo.Opto_exct), ...
    ipsi_dz & Tcombo.Opto_post_sess & Tcombo.Opto_inib, ...
    ipsi_dz & Tcombo.Opto_post_sess & Tcombo.Opto_exct ];
GroupID = {'conf int limit (3std)';...
    ['Thal NOL cells (' num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & ~(Tcombo.Opto_inib | Tcombo.Opto_exct))) ')'];...
    ['SvTh NEG cells (0' num2str(num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & Tcombo.Opto_inib))) ')' ];...
    ['SvTh POS cells (0' num2str(num2str(sum(Tcombo.VMVL & Tcombo.Opto_post_sess & Tcombo.z_exct & Tcombo.Opto_exct))) ')' ]};
pub_fig3_Distrib_Bootstrap_JCfun(-dzSMA, group_ID, group_bool, parfig)
xlim(parfig.xlim);

%% Fig 4M-N Matrices ipsi vs Contra
close all,
dzSMA = zSMAcont-zSMAipsi;
dzSMA_contra = abs(dzSMA(contra_dz,:));
dzSMA_ipsi = abs(dzSMA(ipsi_dz,:));

ZorTT= 'z_thr'
parfig.sort_variable = 'peak';
parfig.colormap = 'jet';
parfig.caxis = [-1 1];
parfig.sort_Xepoch = [parfig.pre-1500 : parfig.pre+750];
parfig.trial_type = parfigcont.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, -dzSMA_contra, Tcombo, parfig)
parfig.trial_type = parfigipsi.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, dzSMA_ipsi, Tcombo, parfig)

parfig.sort_variable = 'Z-tr'
parfig.colormap = 'jet';
parfig.caxis = [-5 5];
parfig.sort_Xepoch = [parfig.pre-1500 : parfig.pre+1000];
parfig.trial_type = parfigcont.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, -dzSMA_contra, Tcombo, parfig)
parfig.trial_type = parfigipsi.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, dzSMA_ipsi, Tcombo, parfig)
pub_fig2_MatricePop_JCfun(ZorTT, -dzSMA(Tcombo.z_exct & Tcombo.VMVL,:), Tcombo, parfig)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUPP FIG S4
close all
pub_figS4_Supp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fig5: OPTO Behavior
% need to clean up and consolidate all the scripts
%% Fig 5A
% Task Schematic made in illustrator
%% Fig 5B Behavior opto task
pub_fig_BehavOpto_BarMeanStat_JCscript
pub_fig5_optoBehavior_JCscript
pub_fig5_optoBehaviorRatio_JCscript3
%% Fig 5C Ephys opto task
close all 
clearvars -except mypath parfig
load('listcell3.mat'); 
load('Tfig2_cor.mat');
parfig.plot=0;
parfig.saveSMA=1;
parfig.saveTABLE =1;
parfig.BaselineEpoch= [150:2150]; BLE = parfig.BaselineEpoch;
parfig.pre = parfig.BaselineEpoch(end) + 1500 % define how much time before zero
parfig.post = 1500 % define how much time after zero (zero will be defined by trigtimes.cor)
parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)

pub_fig5_optoEphys_JCscript







%% PlotShaded Average Pop Z-score of the Firing rate during TASK : OPTO+ vs OPTO-
figure,
K=2;
load SMA_opt_GoCue360.mat  % to compare different cells types during opto trial
zSMA=[]; zSMA1=[]; SMA2=[];
% zSMA1 = (SMA-FRepoch.BLE.mean)./FRepoch.BLE.Sem; zSMA1(find(zSMA1==inf))=NaN;
zSMA1 = (SMA-FRepoch.BLE.mean)./mean(SSemA,2); 
Tfig2_cor=Tfig2_cor(logical(Tfig2_cor.Opto_post_sess),:)

idx=logical(Tfig2_cor.VMVL & Tfig2_cor.z_exct & ~Tfig2_cor.Opto_inib & ~Tfig2_cor.Opto_exct & Tfig2_cor.Opto_post_sess); 
SMA2=zSMA1(idx,:); col='k'; fstr= col; idx1=idx;
Xt=[1:max(size(nanmean(SMA2)))]-parfig.pre;
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;nanmean(SMA2) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;  ], fstr)

idx=logical(Tfig2_cor.VMVL & Tfig2_cor.Opto_inib & Tfig2_cor.z_exct); 
SMA2=zSMA1(idx,:); col='c'; fstr= col; idx2=idx;
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;nanmean(SMA2) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;  ], fstr)

idx=logical(Tfig2_cor.VMVL & Tfig2_cor.Opto_exct & Tfig2_cor.z_exct); 
SMA2=zSMA1(idx,:); col='m'; fstr= col; idx3=idx;
hold on, plotshaded(Xt,[nanmean(SMA2)+(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;nanmean(SMA2) ; nanmean(SMA2)-(K*(nanstd(SMA2)/sqrt(sum(idx)))) ;  ], fstr)

XoptoStim = Xt(parfig.pre-500: parfig.pre+500);
hold on, plot(XoptoStim , zeros(size(XoptoStim))+4, 'c','LineWidth' , 3)

ylabel('Zscore'); xlabel('time');
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);

xlim([-2500 1000]);
ylim([-2 7])
legend(['sem n=' num2str(sum(idx1)) ],' Thal NOL',['sem n=' num2str(sum(idx2))],' SvTh- ', ['sem n=' num2str(sum(idx3))],' SvTh+ ', 'Location','northwest')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fig6: Ephys Impusle Omission SUA and MUA
% close all, clearvars -except mypath sublistcell
% load('listcell.mat'); load('Tcombo.mat')
parfig.saveSMA=0;
parfig.saveTABLE=0;
parfig.center_evt =  'GoCue';
parfig.BaselineEpoch= [150:2150];
parfig.pre = parfig.BaselineEpoch(end) + 1500 % define how much time before zero
parfig.post = 1500 % define how much time after zero (zero will be defined by trigtimes.cor)
parfig.xlim = [-2750 1250];
parfig.merge=1;

% IMPULSE
load('SMA_imp_GoCue360.mat');
parfigimp.trial_type = {'imp'};
parfigimp.col = {'m'};
zSMAimp = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
SMAimp= SMA; SSemAimp=SSemA; clear SMA SSemA SSA;
% OMISSION
load('SMA_omi_GoCue360.mat');
parfigomi.trial_type = {'omi'};
parfigomi.col = {'c'};
zSMAomi = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
SMAomi= SMA; SSemAomi=SSemA; clear SMA SSemA SSA;
% CORRECT
load('SMA_cor_GoCue545.mat');
zSMA = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
zSMAcor= zSMA(logical(Tfig1_VMopto.Opto_post_sess),:);

%% For future figure
Zthr=3
dzSMA_imp = zSMAimp-zSMAcor;
dzSMA_omi = zSMAomi-zSMAcor;

%% TFig6_imp contain only cells tested with opto (360 cells)
Select_Sess_trialTrsh_JCscript
sub_T = listcell(logical(bool_sel),:);
% sub_T = listcell(logical(Tfig1_VMopto.Opto_post_sess),:);
% imp_dz = Tcombo.z_exct(sub_T.ncell) & Tfig5_imp.VMVL & Tfig5_imp.dzic_exctAll ; % dzic ==> dzSMA_imp > 3
% omi_dz = Tcombo.z_exct(sub_T.ncell) & Tfig5_imp.VMVL & Tfig5_imp.dzoc_inibAll; % dzoc ==> dzSMA_omi < 3
imp_dz = Tfig2_cor.z_exct(sub_T.ncell) & Tfig5_imp.VMVL & Tfig5_imp.dzic_exct_delAll ; % dzic ==> dzSMA_imp > 3
omi_dz = Tfig2_cor.z_exct(sub_T.ncell) & Tfig5_imp.VMVL & Tfig5_imp.dzoc_inib_delAll; % dzoc ==> dzSMA_omi < 3
disp('done')

%% Fig 6A Venn Diagrams (2 types: +/-/complex/nosig and 7celltypes)
close all, colmap = 'jet';  resolution = 0.01;
Ncomplex_dz = sum(imp_dz & omi_dz)
data= [sum(imp_dz)-Ncomplex_dz , Ncomplex_dz, sum(omi_dz) - Ncomplex_dz]
vennX(data, resolution, colmap); colormapeditor
Nnosig=sum(~imp_dz & ~omi_dz & Tcombo.z_exct(sub_T.ncell) & Tfig5_imp.VMVL)

%% Fig6B BAR Identify each Type (using Ttest (FrBLE vs FrEpoch))
figure,
Nimp_puf = sum(imp_dz & Tfig5_imp.dzic_exct_pufAll) ;
Nimp_del = sum(imp_dz & Tfig5_imp.dzic_exct_delAll) ;
Nimp_res = sum(imp_dz & Tfig5_imp.dzic_exct_resAll) ;
Nomi_puf = sum(omi_dz & Tfig5_imp.dzoc_inib_pufAll) ;
Nomi_del = sum(omi_dz & Tfig5_imp.dzoc_inib_delAll) ;
Nomi_res = sum(omi_dz & Tfig5_imp.dzoc_inib_resAll);

data = [Nimp_puf Nomi_puf; Nimp_del Nomi_del ;  Nimp_res Nomi_res ]
bar(data)

%% Fig 6C-H : Examples SUA => SDF-RASTER for each Type (from I to VII)
close all
parfig.saveSMA=0;
parfig.saveTABLE=0;
parfig.parfigomi=parfigomi;
parfig.parfigimp = parfigimp;
parfig.xlim= [-1750 1000];
parfig.ylabel = 'z';
parfig.plotshaded= 'sem';
parfig.k=2;
parfig.plot=1;
parfig.plotmerge=1;
parfig.nRaster = 3;

%% PLOT 3 examples for Illustrator
parfig.k=2;
parfig.nRaster = 1;
parfig.currentpanelNb = 1;
pub_fig5_example_impomi_Raster1by1_JCscript

% Screening based on dz
% pub_fig5_example_RasterSDF_impomi_JCfun(Tfig5_imp, dzSMA_imp, dzSMA_omi,  listcell, parfig)

%% PLOT ALL RASTERS
close all;
SaveFig_Fodler = 'D:\JC_Figures\SDF_RASTER\'

LT=Tfig5_imp(imp_dz & omi_dz,:)
nn= size(LT,1)
for nc = 1:nn;
    RowName = LT.Row(nc);
    listc1= listcell(RowName{1}, :)
    figure, hold on,
    %Imp
    parfig.col = parfig.parfigimp.col;
    parfig.trial_type = parfig.parfigimp.trial_type;
    parfig.currentpanelNb = 1;
    pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
    % omi
    parfig.col = parfig.parfigomi.col;
    parfig.trial_type = parfig.parfigomi.trial_type;
    parfig.currentpanelNb = 2;
    pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
    disp('done Raster')
    % cor
    parfig.col = {'k'};
    parfig.trial_type = {'cor'};
    parfig.currentpanelNb = 3;
    pub_comp_ephySMA_task_ttest_JCfun(listc1, parfig);
    disp('done Raster')
    saveas(gcf, [SaveFig_Fodler 'impomi_' RowName{1}],'png')
end


%% Fig 6I-J Mean DZ
Xt=[1:1:size(dzSMA_imp,2)]-parfig.pre;
parfig.xlim=[-2000 750]
k=2;

figure,
hold on, plot(Xt, mean(dzSMA_imp(imp_dz,:)), 'm')
hold on, plot(Xt, mean(dzSMA_omi(omi_dz,:)), 'c')
SDF1 = dzSMA_imp(imp_dz,:);
SDF2 = dzSMA_omi(omi_dz,:);
hold on, plotshaded(Xt, [mean(SDF1)+ (k*(std(SDF1 )/sqrt(sum(imp_dz)))); mean(SDF1)- (k*(std(SDF1 )/sqrt(sum(imp_dz))))] ,'m');
hold on, plotshaded(Xt, [mean(SDF2)+ (k*(std(SDF2 )/sqrt(sum(omi_dz)))); mean(SDF2)- (k*(std(SDF2 )/sqrt(sum(omi_dz))))] ,'c');

ylabel('mean(dz vs cor))'), xlabel('time ms')
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
xlim([parfig.xlim(1)-150 parfig.xlim(2)])
legend(['imp (#' num2str(sum(imp_dz))  ')' ], ['omi (#' num2str(sum(omi_dz))  ')'], 'Location','Best')

%% Fig 6K-L Distribution Bootstrap Opto

%  DIST IMPULSES opto Bootstrap (Fig 5L)
parfig.Zthr = 3;
parfig.bootstrap.iter = 1000
parfig.bootstrap.groupsize=10
parfig.trial_type = parfig.parfigimp.trial_type
parfig.title = 'IMP Dist bootstrap'; parfig.ylim=[0 70];
group_bool= [imp_dz & Tfig5_imp.Opto_post_sess &  ~(Tfig5_imp.Opto_inib | Tfig5_imp.Opto_exct), ...
    imp_dz & Tfig5_imp.Opto_post_sess & Tfig5_imp.Opto_inib, ...
    imp_dz & Tfig5_imp.Opto_post_sess & Tfig5_imp.Opto_exct ];
group_ID = {'conf int limit (3std)';...
    ['Thal NOL cells (' num2str(sum( imp_dz & Tfig5_imp.VMVL & Tfig5_imp.Opto_post_sess & ~(Tfig5_imp.Opto_inib | Tfig5_imp.Opto_exct))) ')'];...
    ['SvTh NEG cells (' num2str(num2str(sum( imp_dz & Tfig5_imp.VMVL & Tfig5_imp.Opto_post_sess  & Tfig5_imp.Opto_inib))) ')' ];...
    ['SvTh POS cells (0' num2str(num2str(sum( imp_dz & Tfig5_imp.VMVL & Tfig5_imp.Opto_post_sess & Tfig5_imp.Opto_exct))) ')' ]};
pub_fig3_Distrib_Bootstrap_JCfun(dzSMA_imp, group_ID, group_bool, parfig)
xlim(parfig.xlim);
ylim([0 100])

% DIST OMISSION opto Bootstrap (Fig 5L)
parfig.Zthr = 3;
parfig.bootstrap.iter = 1000
parfig.bootstrap.groupsize=10
parfig.trial_type = parfig.parfigomi.trial_type
parfig.title = ' OMI Dist bootstrap'; parfig.ylim=[0 70];
group_bool= [omi_dz & Tfig5_imp.Opto_post_sess &  ~(Tfig5_imp.Opto_inib | Tfig5_imp.Opto_exct), ...
    omi_dz & Tfig5_imp.Opto_post_sess & Tfig5_imp.Opto_inib, ...
    omi_dz & Tfig5_imp.Opto_post_sess & Tfig5_imp.Opto_exct ];
group_ID = {'conf int limit (3std)';...
    ['Thal NOL cells (' num2str(sum(omi_dz & Tfig5_imp.Opto_post_sess &  ~(Tfig5_imp.Opto_inib | Tfig5_imp.Opto_exct))) ')'];...
    ['SvTh NEG cells (' num2str(sum(omi_dz & Tfig5_imp.Opto_post_sess & Tfig5_imp.Opto_inib)) ')' ];...
    ['SvTh POS cells (0' num2str(num2str(sum(omi_dz & Tfig5_imp.Opto_post_sess & Tfig5_imp.Opto_exct))) ')' ]};
pub_fig3_Distrib_Bootstrap_JCfun(-dzSMA_omi, group_ID, group_bool, parfig)
xlim(parfig.xlim);
ylim([0 100])

%% Fig 6M-N Matrices 
close all,
% dzSMA_imp =
ZorTT= 'z_thr'
parfig.sort_variable = 'peak';
parfig.colormap = 'jet';
parfig.caxis = [-1 1];
parfig.sort_Xepoch = [parfig.pre-750 : parfig.pre+750];
parfig.trial_type = parfigimp.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, dzSMA_imp(imp_dz,:), Tfig5_imp, parfig)
parfig.trial_type = parfigomi.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, dzSMA_omi(omi_dz,:), Tfig5_imp, parfig)

parfig.sort_variable = 'Z-tr'
parfig.colormap = 'jet';
parfig.caxis = [-3.5 3.5];
parfig.sort_Xepoch = [parfig.pre-1500 : parfig.pre+1000];
parfig.trial_type = parfigimp.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT,  dzSMA_imp(imp_dz,:), Tfig5_imp, parfig)
parfig.trial_type = parfigomi.trial_type;
pub_fig2_MatricePop_JCfun(ZorTT, dzSMA_omi(omi_dz & ~(imp_dz & ~(imp_dz & omi_dz)) ,:), Tfig5_imp, parfig)


%% SUPP Fig6S : LICK ALIGNS Ephys Impusle Omission SUA and MUA
close all, clearvars -except mypath sublistcell
load('listcell3.mat'); load('Tfig1_VMopto.mat'); load('Tfig5_imp.mat'); load('Tfig2_cor.mat'); 
parfig.saveSMA=0;
parfig.saveTABLE=0;
parfig.center_evt =  'Licks';
parfig.BaselineEpoch= [150:2150];
parfig.pre = parfig.BaselineEpoch(end) + 1500 % define how much time before zero
parfig.post = 1500 % define how much time after zero (zero will be defined by trigtimes.cor)
parfig.xlim = [-2750 1250];
parfig.merge=1;
sub_T = listcell(logical(Tfig1_VMopto. Opto_post_sess),:);

% IMPULSE
load('SMA_imp_Licks360.mat');
parfigimp.trial_type = {'imp'};
parfigimp.col = {'m'};
zSMAimp = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
SMAimp= SMA; SSemAimp=SSemA; clear SMA SSemA SSA;
% OMISSION
load('SMA_omi_GoCue360.mat');
parfigomi.trial_type = {'omi'};
parfigomi.col = {'c'};
zSMAomi = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
SMAomi= SMA; SSemAomi=SSemA; clear SMA SSemA SSA;
% CORRECT
load('SMA_cor_Licks360.mat');
zSMAcor = (SMA-FRepoch.BLE.mean)./mean(SSemA,2);
SMAcor= SMA; SSemAcor=SSemA; clear SMA SSemA SSA;

Zthr=3
dzSMA_imp = zSMAimp-zSMAcor;


% imp_dz = Tcombo.z_exct(sub_T.ncell) & Tfig5_imp.VMVL & Tfig5_imp.dzic_exctAll; % all Contra
% omi_dz = Tcombo.z_exct(sub_T.ncell) & Tfig5_imp.VMVL & Tfig5_imp.dzoc_inibAll; % all Ipsi
imp_dz = Tfig2_cor.z_exct(sub_T.ncell) & Tfig5_imp.VMVL & Tfig5_imp.dzic_exct_delAll; % all Contra
omi_dz = Tfig2_cor.z_exct(sub_T.ncell) & Tfig5_imp.VMVL & Tfig5_imp.dzoc_inib_delAll; % all Ipsi
disp('done')

%% Fig 6I-J Mean DZ
Xt=[1:1:size(dzSMA_imp,2)]-parfig.pre;
parfig.xlim=[-2000 750]
k=2;

figure,
% hold on, plot(Xt, mean(dzSMA_imp(imp_dz,:)), 'r')
% hold on, plot(Xt, mean(dzSMA_omi(omi_dz,:)), 'c')

hold on, plot(Xt, mean(zSMAimp(imp_dz,:)), 'm')
hold on, plot(Xt, mean(zSMAcor(imp_dz,:)), 'k')

SDF1 = zSMAimp(imp_dz,:);
SDF2 = zSMAcor(imp_dz,:);
hold on, plotshaded(Xt, [mean(SDF1)+ (k*(std(SDF1 )/sqrt(sum(imp_dz)))); mean(SDF1)- (k*(std(SDF1 )/sqrt(sum(imp_dz))))] ,'m');
hold on, plotshaded(Xt, [mean(SDF2)+ (k*(std(SDF2 )/sqrt(sum(omi_dz)))); mean(SDF2)- (k*(std(SDF2 )/sqrt(sum(omi_dz))))] ,'k');

hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);

xlim([parfig.xlim(1)-150 parfig.xlim(2)])
legend(['imp (#' num2str(sum(imp_dz)) ')' ],['cor (#' num2str(sum(imp_dz)) ')'], 'Location','Best')
ylabel('mean zscore'), xlabel('time ms')


%% Find SDF Peak time relative to Lick onset
pre=parfig.pre; 
MAXimp = max(SDF1'); % peak Value for 60 cells 
MAXcor = max(SDF2'); 
for cc=1:60
  AA(cc)=find(SDF1(cc,:)==MAXimp(cc))  ;
  BB(cc)=find(SDF2(cc,:)==MAXcor(cc));
end

AAA= AA-parfig.pre
BBB= BB-parfig.pre
mean(AAA)
mean(BBB)
std(AAA)/sqrt(60)
std(BBB)/sqrt(60)


disp(['Time of Peak before imp Lick = ' num2str(mean(AAA)) '+/- ' num2str(std(AAA)/sqrt(60)) 'ms  (ncell = ' num2str(cc) ') windows around Lick = +/-' num2str(w) 'ms'])
disp(['Time of Peak before cor Lick = ' num2str(mean(BBB)) '+/- ' num2str(std(BBB)/sqrt(60)) 'ms  (ncell = ' num2str(cc) ') windows around Lick = +/-' num2str(w) 'ms'])
% FIND THE START of THE RAMP 
SDF1PNorm = (SDF1'./MAXimp);
SDF2PNorm = (SDF2'./MAXcor);
w=750
SDF1w = SDF1PNorm(pre-w:pre+w,:)';
SDF2w = SDF2PNorm(pre-w:pre+w,:)';

% SDF1w = SDF1(:,pre-750:pre);
% SDF2w = SDF2(:,pre-750:pre);
RAMPimp=[]; RAMPcor=[];
for cc=1:60
    cc;
  RAMPimp(cc) = min(find(SDF1w(cc,:)>0.20));
  RAMPcor(cc) = min(find(SDF2w(cc,:)>0.20));
%   THR2 = (SDF2(cc,:)/MAXcor(cc))>0.25
%   SDF1w = SDF1(:,pre-750:pre);
%   SDF2w = SDF2(:,pre-750:pre);
%   RAMPimp(cc) = min(find()
%   RAMPcor(cc) = min(find((SDF2(cc,:)/MAXcor(cc))>0.25))
end

CCC= RAMPimp-w;
DDD= RAMPcor-w;
mean(CCC);
mean(DDD);
std(CCC)/sqrt(60);
std(DDD)/sqrt(60);

disp(['Time Ramp start before imp Lick = ' num2str(mean(CCC)) '+/- ' num2str(std(CCC)/sqrt(60)) 'ms  (ncell = ' num2str(cc) ') windows around Lick = +/-' num2str(w) 'ms'])
disp(['Time Ramp start before cor Lick = ' num2str(mean(DDD)) '+/- ' num2str(std(DDD)/sqrt(60)) 'ms  (ncell = ' num2str(cc) ') windows around Lick = +/-' num2str(w) 'ms'])
%% 


A1 = mean(SDF1);
AA1=A1/max(A1);
AAA1 = min(find(AA1(pre-750:pre+750)>0.20));
AAAA1 = AAA1-750

B1 = mean(SDF2);
BB1=B1/max(B1);
BBB1 = min(find(BB1(pre-750:pre+750)>0.20));
BBBB1 = BBB1-750

A1 = mean(SDF1);
AA1=A1;
AAA1 = min(find(AA1(pre-750:pre+750)>1));
AAAA1 = AAA1-750

B1 = mean(SDF2);
BB1=B1;
BBB1 = min(find(BB1(pre-750:pre+750)>1));
BBBB1 = BBB1-750






















%% Fig7  plot Classifier Left vs Right during Delay
% close all;
clearvars -except mypath parfig
rng('shuffle');
load('listcell.mat'); load('Tcombo.mat');
parfig.plot=0; parfig.saveSMA=0; parfig.saveTABLE =0;

parfig.center_evt = 'GoCue';    % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.epoch  = 'delay'
if parfig.epoch  =='delay'
    parfig.pre = 750 ; % define how much time before zero (in ms)
    parfig.post= 0;
end

Mall=[];
SelectSess = Tcombo.nSess

figure(1), close(1), figure (1),
Nb_sess= max(Tcombo.nSess )
%
for ii=1:Nb_sess
    listc=listcell(Tcombo.nSess==ii & Tcombo.VMVL,:);
    ncell = size(listc,1);
    if ncell >= 20;
        
        parfig.trial_type = {'cCL'}; trialtype1= parfig.trial_type{1};
        [nspx_A] = pub_get_nspx_trialtype_JCfun(listc, parfig);
        parfig.trial_type = {'cCR'}; trialtype2= parfig.trial_type{1};
        [nspx_B] = pub_get_nspx_trialtype_JCfun(listc, parfig);
        %
        parfig.Nrepeat = 10;
        parfig.typeClass='FoldXVal' ; %     typeClass='%HoldOut'
        parfig.learner= 'logistic';%, 'svm'}
        parfig.Nfold = 10;
        
        Lnall=[];
        for Nr = 1:parfig.Nrepeat;
            %     [Ln] = play the function n times (ncell, )
            parfig.ControlShuffle=1;
            [Ln] = pub_Classifier_Ln_JCfun(ncell, nspx_A, nspx_B, parfig); % nspx_A = ntrials x ncells
            Lnall=[Lnall;Ln];
        end
        
        % PLOT
        MouseID= listc.MouseID(1,:)
        Day=listc.Day(1,:)
        M=mean(Lnall);
        S=std(Lnall);
        
        figure(1), hold on,
        plot( [1:1:ncell] , M, 'LineWidth', 2)
        
        % M1=M; save(['ClassM1_' psth_center_evtID '_' trial_type_list{1} '_' trial_type_list{2}],'M1','trial_type_list')
        ylim([0 1])
        xlabel('#Neurons', 'FontSize', 11)
        ylabel('Proba errors', 'FontSize', 11)
        title(parfig.center_evt, 'FontWeight','bold' ,'FontSize', 12)
        
        Mall = [Mall; M(1:20)];
    else
        disp('not enough cells in this session')
    end
end
%%
hold on,
figure,
plot( [1:1:ncell] , ones(1,ncell)/2, 'k--')
%%
% save Mall
if parfig.ControlShuffle==0;
    save(['Class10Fold_Mall_'  trialtype1 'v' trialtype2  'delay_Orig.mat'], 'Mall')
elseif parfig.ControlShuffle==1;
    save(['Class10Fold_Mall_' trialtype1 'v' trialtype2  '_delay_ShuffleControl.mat'], 'Mall')
end

%% Fig7  Classifier
% DELAY epoch, GOCue Centered, cor vs omi, cor vs imp
% Computing the Classifier10FOld _Mall.mat (saved in JC_Analysis)

Classifier_2020 % need to be renamed 


close all;
SaveFigFolder = 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig7_classifier'
clearvars -except mypath parfig SaveFigFolder
parfig.ControlShuffle=0;
rng('shuffle');
load('listcell.mat'); load('Tfig1_VMopto.mat');

parfig.Nrepeat = 100;
parfig.typeClass='FoldXVal' ; %     typeClass='%HoldOut'
parfig.learner= 'logistic';%, 'svm'}
parfig.Nfold = 10;

parfig.center_evt = 'GoCue';% center ('Delay' ; 'GoCue' ;'APuff'; 'Licks')
parfig.epoch = 'delay'

if parfig.epoch  =='delay'
    parfig.pre = 750 ; % define how much time before zero (in ms)
    parfig.post= 0;
elseif parfig.epoch == 'Apuff'
    parfig.pre = 1500 ; % define how much time before zero (in ms)
    parfig.post= -750;
end

%% Th- CELLS
celltype = 'SvTh-'
idxcell = Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib
parfig.minNbcell = 9;
parfig.ControlShuffle=0;

% trial_type_1 = 'cCL'
% trial_type_2 = 'cCR'
% pub_fig7_Classifier_JCscript

trial_type_1 = 'cor'
trial_type_2 = 'omi'
pub_fig7_Classifier_JCscript

trial_type_1 = 'cor'
trial_type_2 = 'imp'
pub_fig7_Classifier_JCscript

%% NOL CELLS
celltype = 'NOL'
idxcell = Tfig1_VMopto.VMVL & ~(Tfig1_VMopto.Opto_exct | Tfig1_VMopto.Opto_inib)
parfig.minNbcell = 15;
parfig.ControlShuffle=0;
%
% trial_type_1 = 'cCL'
% trial_type_2 = 'cCR'
% pub_fig7_Classifier_JCscript

trial_type_1 = 'cor'
trial_type_2 = 'omi'
pub_fig7_Classifier_JCscript

trial_type_1 = 'cor'
trial_type_2 = 'imp'
pub_fig7_Classifier_JCscript

%% Shuffle control
celltype = 'ShuffleControl'
idxcell = Tfig1_VMopto.VMVL
parfig.minNbcell = 20;
parfig.ControlShuffle=1;

% trial_type_1 = 'cCL'
% trial_type_2 = 'cCR'
% pub_fig7_Classifier_JCscript

trial_type_1 = 'cor'
trial_type_2 = 'omi'
pub_fig7_Classifier_JCscript

trial_type_1 = 'cor'
trial_type_2 = 'imp'
pub_fig7_Classifier_JCscript

% All raw traces per sessions
%% VMVL all
celltype = 'VMVL'
idxcell = Tfig1_VMopto.VMVL
parfig.minNbcell = 15;
parfig.ControlShuffle=0;

% trial_type_1 = 'cCL'
% trial_type_2 = 'cCR'
% pub_fig7_Classifier_JCscript

% trial_type_1 = 'cor'
% trial_type_2 = 'omi'
% pub_fig7_Classifier_JCscript
% 
% trial_type_1 = 'cor'
% trial_type_2 = 'imp'
% pub_fig7_Classifier_JCscript

trial_type_1 = 'imp'
trial_type_2 = 'omi'
pub_fig7_Classifier_JCscript

%% Average for distinct population
% DELAY epoch, GOCue Centered, cor vs omi, cor vs imp
% PLOTING Mean +/- 1SEM (shaded)
SaveFigFolder = 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig7_classifier'
for jj = 1:2 % 1=omi ; 2=imp
clearvars -except mypath parfig Tfig1_VMopto listcell SaveFigFolder jj
parfig.center_evt = 'GoCue';% center ('Delay' ; 'GoCue' ;'APuff'; 'Licks')
parfig.epoch = 'delay'
figure(10) ; close (10); figure(10)
    for ii=1:3 %1=NOL; 2=Th-; 3=Shuffle
        if ii==1
            celltype = 'NOL'
            parfig.minNbcell = 15;
            col = 'k'
        elseif ii==2
            celltype = 'SvTh-'
            parfig.minNbcell = 9;
            col = 'c'
        elseif ii==3
            celltype = 'ShuffleControl'
            parfig.minNbcell = 20;
            col = 'r'
        end
        
        trial_type_1 = 'cor'
        if jj==1; trial_type_2 = 'omi'; end
        if jj==2  trial_type_2 = 'imp'; end
        
        
        load(['Class10Fold_Mall_' trial_type_1 'v' trial_type_2  '_' parfig.epoch '_' parfig.center_evt  '_' celltype '.mat'], 'Mall')
        figure (10), hold on,
        K = 1;
        
        Mm=mean(Mall)
        NSessFinal = size(Mall,1)
        Mstd = std(Mall)/sqrt(NSessFinal)
        XX=[1:1:parfig.minNbcell];
        plot( XX , Mm, col, 'LineWidth', 3)
        plotshaded(XX, [Mm-K*Mstd ; Mm+K*Mstd], col)
    end
    
    hold on,
    plot( XX , ones(1,parfig.minNbcell)/2, 'k--')
    ylim([0.2 0.8]);
    title(['Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt  ' ' celltype ])
    legend('NOL', '', 'SvTh-', '','Shuffle','')
    
    saveas(gcf, [SaveFigFolder '\Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt  ' ' celltype ],'png')
    saveas(gcf, [SaveFigFolder '\Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt  ' ' celltype ],'emf')
    saveas(gcf, [SaveFigFolder '\Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt  ' ' celltype ],'fig')
    
end

%% LICKS centered for IMP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 750ms PRE-LICK epoch, Licks Centered, cor vs imp
% Computing the Classifier10FOld _Mall.mat (saved in JC_Analysis)
clearvars -except mypath parfig SaveFigFolder
parfig.ControlShuffle=0;
rng('shuffle');
load('listcell.mat'); load('Tfig1_VMopto.mat');

parfig.Nrepeat = 10;
parfig.typeClass='FoldXVal' ; %     typeClass='%HoldOut'
parfig.learner= 'logistic';%, 'svm'}
parfig.Nfold = 10;


parfig.center_evt = 'Licks';    % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.epoch = 'delay'

if parfig.epoch  =='delay'
    parfig.pre = 750 ; % define how much time before zero (in ms)
    parfig.post= 0;
elseif parfig.epoch == 'Apuff'
    parfig.pre = 1500 ; % define how much time before zero (in ms)
    parfig.post= -750;
end

% Th- CELLS
celltype = 'SvTh-'
idxcell = Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib
parfig.minNbcell = 9;
parfig.ControlShuffle=0;

% trial_type_1 = 'cCL'
% trial_type_2 = 'cCR'
% pub_fig7_Classifier_JCscript

trial_type_1 = 'cor'
trial_type_2 = 'imp'
pub_fig7_Classifier_JCscript

% NOL CELLS
celltype = 'NOL'
idxcell = Tfig1_VMopto.VMVL & ~(Tfig1_VMopto.Opto_exct | Tfig1_VMopto.Opto_inib)
parfig.minNbcell = 15;
parfig.ControlShuffle=0;
%
% trial_type_1 = 'cCL'
% trial_type_2 = 'cCR'
% pub_fig7_Classifier_JCscript

trial_type_1 = 'cor'
trial_type_2 = 'imp'
pub_fig7_Classifier_JCscript


% Shuffle control
celltype = 'ShuffleControl'
idxcell = Tfig1_VMopto.VMVL
parfig.minNbcell = 20;
parfig.ControlShuffle=1;

% trial_type_1 = 'cCL'
% trial_type_2 = 'cCR'
% pub_fig7_Classifier_JCscript

trial_type_1 = 'cor'
trial_type_2 = 'imp'
pub_fig7_Classifier_JCscript

% VMVL all raw plots for each sessions
celltype = 'VMVL'
idxcell = Tfig1_VMopto.VMVL
parfig.minNbcell = 15;
parfig.ControlShuffle=0;

trial_type_1 = 'cor'
trial_type_2 = 'imp'
pub_fig7_Classifier_JCscript

%% Average for distinct population
% PRE-LICK epoch, Lick Centered, cor vs imp
% PLOTING Mean +/- 1SEM (shaded)
SaveFigFolder = 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig7_classifier'
clearvars -except mypath parfig Tfig1_VMopto listcell SaveFigFolder jj
figure(10) ; close (10); figure(10)
parfig.center_evt = 'Licks';    % center for SDF ('Delay' ; 'GOcue' ; Licks)
parfig.epoch = 'delay'

for ii=1:3
    if ii==1
        celltype = 'NOL'
        parfig.minNbcell = 15;
        col = 'k'
    elseif ii==2
        celltype = 'SvTh-'
        parfig.minNbcell = 9;
        col = 'c'
    elseif ii==3
        celltype = 'ShuffleControl'
        parfig.minNbcell = 20;
        col = 'r'
    end
    
    trial_type_1 = 'cor'
    trial_type_2 = 'imp';
    
    load(['Class10Fold_Mall_' trial_type_1 'v' trial_type_2  '_' parfig.epoch '_' parfig.center_evt '_' celltype '.mat'], 'Mall')
    figure (10), hold on,
    K = 1;
    dddd
    Mm=mean(Mall)
    NSessFinal = size(Mall,1)
    Mstd = std(Mall)/sqrt(NSessFinal)
    XX=[1:1:parfig.minNbcell];
    plot( XX , Mm, col, 'LineWidth', 3)
    plotshaded(XX, [Mm-K*Mstd ; Mm+K*Mstd], col)
end

hold on,
plot( XX , ones(1,parfig.minNbcell)/2, 'k--')
ylim([0.2 0.8]);
title(['Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt ' '  celltype ])
legend('NOL', '', 'SvTh-', '','Shuffle','')

saveas(gcf, [SaveFigFolder '\Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt  ' ' celltype ],'png')
saveas(gcf, [SaveFigFolder '\Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt  ' ' celltype ],'emf')
saveas(gcf, [SaveFigFolder '\Mean Class10Fold ' trial_type_1 ' vs ' trial_type_2  ' ' parfig.epoch  ' ' parfig.center_evt  ' ' celltype ],'fig')



