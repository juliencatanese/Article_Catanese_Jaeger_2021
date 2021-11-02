% Plot_DLC_opto_AverageMultiSess_JCScript
% function plot_Abs_Average_DLCvar_JCfun(VAR1, VAR2, VARname_str1, VARname_str2)
% Julien Catanese 9/23/2020

clc, close all, clear all,

ncor  = 0; ncoropto = 0;
corW = []; corWopto = []; corWsem = []; corWoptosem = []; 
corT = []; corTopto = []; corTsem = []; corToptosem = [];   
corN = []; corNopto = []; corNsem = []; corNoptosem = [];  


nimp  = 0; nimpopto = 0;
impW = []; impWopto = []; impWsem = []; impWoptosem = []; 
impT = []; impTopto = []; impTsem = []; impToptosem = [];   
impN = []; impNopto = []; impNsem = []; impNoptosem = [];  

nomi  = 0; nomiopto = 0;
omiW = []; omiWopto = []; omiWsem = []; omiWoptosem = []; 
omiT = []; omiTopto = []; omiTsem = []; omiToptosem = [];   
omiN = []; omiNopto = []; omiNsem = []; omiNoptosem = [];  

for nm=1:4
    nm
    if nm==1
        FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW';
        MouseID =  'vgat15w10d7';
    elseif nm==2
        FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d8_z4900_VM_taskopto_optopost_G912_180715_vidY_150tr_20cel_13mW_tooDeep';
        MouseID =  'vgat14w14d8';
    elseif nm==3
        FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW'
        MouseID =  'vgat17w10d7'
    elseif nm==4
        FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d5_z4300_VM_taskopto_optopost_EOBD_180728_vidY_200tr_20cel_12mW_Q4';
        MouseID =  'vgat17w10d5';
    end
    cd(FileLocation);
    load([ MouseID '_DLCresults.mat']);
    load('Ntrial_type.mat');
    load ('evt.mat','evt_trial', 'evt_lick_L', 'evt_lick_R');
    load ('time.mat');
    
    % Get_DLC_VAR_Trial_Mat_script
    %% average over trials (100 frames)
    NoseY_tr = allTab.NoseY(1:100); TongueY_tr = allTab.TongueY(1:100); WhiskerX_tr = allTab.WhiskerX(1:100);
    for itr = 1:(size(allTab,1)/100)-1;
        NoseY_tr = [NoseY_tr allTab.NoseY(itr*100:((itr+1)*100)-1)];
        TongueY_tr = [TongueY_tr allTab.TongueY(itr*100:((itr+1)*100)-1)];
        WhiskerX_tr = [WhiskerX_tr allTab.WhiskerX(itr*100:((itr+1)*100)-1)];
    end
    
    %% cor vs opto trials selection
    idx1 = [trial.idx_correct_L trial.idx_correct_R];
    idx2 = [trial.idx_correct_L_opto trial.idx_correct_R_opto];
    idx_str1 = 'cor';
    idx_str2 = 'cor-opto';
    
    VAR_DLC = WhiskerX_tr; VAR_str = 'WhiskerX'
    [VAR1norm VAR2norm VAR1sem VAR2sem n1 n2]= plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2);
    corW = [corW; VAR1norm ];
    corWopto = [corWopto; VAR2norm ];
    corWsem = [corWsem; VAR1sem ];
    corWoptosem = [corWoptosem; VAR2sem ];
    ncor  = ncor + n1;
    ncoropto = ncoropto + n2;
    
    VAR_DLC=TongueY_tr; VAR_str = 'TongueY'
    [VAR1norm VAR2norm VAR1sem VAR2sem n1 n2]= plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2);
    corT = [corT; VAR1norm ];
    corTopto = [corTopto; VAR2norm ];
    corTsem = [corTsem; VAR1sem ];
    corToptosem = [corToptosem; VAR2sem ];
    
    VAR_DLC=NoseY_tr; VAR_str = 'NoseY'
    [VAR1norm VAR2norm VAR1sem VAR2sem n1 n2]= plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2);
    corN = [corN; VAR1norm ];
    corNopto = [corNopto; VAR2norm ];
    corNsem = [corNsem; VAR1sem ];
    corNoptosem = [corNoptosem; VAR2sem ];
    
    
    %% imp vs omiopto trials selection
    idx1 = sort([trial.idx_errorDelay_PL_CL, trial.idx_errorDelay_PR_CR, trial.idx_errorDelay_PL_CR,...
        trial.idx_errorDelay_PR_CL]);
    idx2 = sort([trial.idx_errorDelay_PL_CL_opto, trial.idx_errorDelay_PR_CR_opto,...
        trial.idx_errorDelay_PL_CR_opto, trial.idx_errorDelay_PR_CL_opto]);
    idx_str1 = 'imp';
    idx_str2 = 'imp-opto';
    
    VAR_DLC = WhiskerX_tr; VAR_str = 'WhiskerX'
    [VAR1norm VAR2norm VAR1sem VAR2sem ntr1 ntr2]= plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2);
    impW = [impW; VAR1norm ];
    impWopto = [impWopto; VAR2norm ];
    impWsem = [impWsem; VAR1sem ];
    impWoptosem = [impWoptosem; VAR2sem ];
    nimp  = nimp + n1; 
    nimpopto = nimpopto + n2; Ni(nm)=n2; 
    
    VAR_DLC=TongueY_tr; VAR_str = 'TongueY'
    [VAR1norm VAR2norm VAR1sem VAR2sem ntr1 ntr2]= plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2);
    impT = [impT; VAR1norm ];
    impTopto = [impTopto; VAR2norm ];
    impTsem = [impTsem; VAR1sem ];
    impToptosem = [impToptosem; VAR2sem ];
    
    VAR_DLC=NoseY_tr; VAR_str = 'NoseY'
    [VAR1norm VAR2norm VAR1sem VAR2sem ntr1 ntr2]= plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2);
    impN = [impN; VAR1norm ];
    impNopto = [impNopto; VAR2norm ];
    impNsem = [impNsem; VAR1sem ];
    impNoptosem = [impNoptosem; VAR2sem ];
    
    %% omi vs omiopto trials selection
    idx1 = trial.idx_NoLick;
    idx2 = trial.idx_NoLick_opto;
    idx_str1 = 'omi';
    idx_str2 = 'omi-opto';
    
    VAR_DLC = WhiskerX_tr; VAR_str = 'WhiskerX'
    [VAR1norm VAR2norm VAR1sem VAR2sem n1 n2]= plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2);
    omiW = [omiW; VAR1norm ];
    omiWopto = [omiWopto; VAR2norm ];
    omiWsem = [omiWsem; VAR1sem ];
    omiWoptosem = [omiWoptosem; VAR2sem ];
    nomi  = nomi + n1;
    nomiopto = nomiopto + n2;
    
    VAR_DLC=TongueY_tr; VAR_str = 'TongueY'
    [VAR1norm VAR2norm VAR1sem VAR2sem n1 n2]= plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2);
    omiT = [omiT; VAR1norm ];
    omiTopto = [omiTopto; VAR2norm ];
    omiTsem = [omiTsem; VAR1sem ];
    omiToptosem = [omiToptosem; VAR2sem ];
    
    VAR_DLC=NoseY_tr; VAR_str = 'NoseY'
    [VAR1norm VAR2norm VAR1sem VAR2sem n1 n2]= plot_Abs_Average_DLCvar_JCfun(VAR_DLC, idx1, idx2, VAR_str, idx_str1, idx_str2);
    omiN = [omiN; VAR1norm ];
    omiNopto = [omiNopto; VAR2norm ];
    omiNsem = [omiNsem; VAR1sem ];
    omiNoptosem = [omiNoptosem; VAR2sem ];
    
    
end

%% PLOT AVERAGE
close all 
Ylimitas=[-2 25]
Xlimitas=[0 75]
% PLOT COR
% SNOUT cor opto
figure, hold on,
title(['SNOUT cor (n=' num2str(ncor) ')'  'opto (n=' num2str(ncoropto)  ')'])
K=2
X=[1:1:100]';
nanmean1 = nanmean((corN'),2);
nanmean2 = nanmean((corNopto'),2);
std1 = nanmean((corNsem'),2);
std2 = nanmean((corNoptosem'),2);
plot(abs(nanmean1),'k','LineWidth',2);
plot(abs(nanmean2),'c','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
plotshaded(X', abs([nanmean1-(K*std1) nanmean1+(K*std1)]),'k');
plotshaded(X', abs([nanmean2-(K*std2) nanmean2+(K*std2)]),'c');
ylim([Ylimitas])
xlim([Xlimitas])
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_corvopto_Snout','png')
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_corvopto_Snout','emf')

% Tongue cor opto
figure, hold on,
title(['Tongue cor (n=' num2str(ncor) ')'  'opto (n=' num2str(ncoropto)  ')'])
K=2
X=[1:1:100]';
nanmean1 = nanmean((corT'),2);
nanmean2 = nanmean((corTopto'),2);
std1 = nanmean((corTsem'),2);
std2 = nanmean((corToptosem'),2);
plot(abs(nanmean1),'k','LineWidth',2);
plot(abs(nanmean2),'c','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
plotshaded(X', abs([nanmean1-(K*std1) nanmean1+(K*std1)]),'k');
plotshaded(X', abs([nanmean2-(K*std2) nanmean2+(K*std2)]),'c');
ylim([Ylimitas])
xlim([Xlimitas])
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_corvopto_Tongue','png')
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_corvopto_Tongue','emf')

% Whisker cor opto
figure, hold on,
title(['Whisker cor (n=' num2str(ncor) ')'  'opto (n=' num2str(ncoropto)  ')'])
K=2
X=[1:1:100]';
nanmean1 = nanmean((corW'),2);
nanmean2 = nanmean((corWopto'),2);
std1 = nanmean((corWsem'),2);
std2 = nanmean((corWoptosem'),2);
plot(abs(nanmean1),'k','LineWidth',2);
plot(abs(nanmean2),'c','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
plotshaded(X', abs([nanmean1-(K*std1) nanmean1+(K*std1)]),'k');
plotshaded(X', abs([nanmean2-(K*std2) nanmean2+(K*std2)]),'c');
ylim([Ylimitas])
xlim([Xlimitas])
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_corvopto_Whisker','png')
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_corvopto_Whisker','emf')

% PLOT IMP
% SNOUT imp opto
figure, hold on,
title(['SNOUT imp (n=' num2str(nimp) ')'  'opto (n=' num2str(nimpopto)  ')'])
K=2
X=[1:1:100]';
nanmean1 = nanmean((impN'),2);
nanmean2 = nanmean((impNopto'),2);
std1 = nanmean((impNsem'),2);
std2 = nanmean((impNoptosem'),2);
plot(abs(nanmean1),'m','LineWidth',2);
plot(abs(nanmean2),'c','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
plotshaded(X', abs([nanmean1-(K*std1) nanmean1+(K*std1)]),'m');
plotshaded(X', abs([nanmean2-(K*std2) nanmean2+(K*std2)]),'c');
ylim([Ylimitas])
xlim([Xlimitas])
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_impvopto_Snout','png')
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_impvopto_Snout','emf')

% Tongue imp opto
figure, hold on,
title(['Tongue imp (n=' num2str(nimp) ')'  'opto (n=' num2str(nimpopto)  ')'])
K=2
X=[1:1:100]';
nanmean1 = nanmean((impT'),2);
nanmean2 = nanmean((impTopto'),2);
std1 = nanmean((impTsem'),2);
std2 = nanmean((impToptosem'),2);
plot(abs(nanmean1),'m','LineWidth',2);
plot(abs(nanmean2),'c','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
plotshaded(X', abs([nanmean1-(K*std1) nanmean1+(K*std1)]),'m');
plotshaded(X', abs([nanmean2-(K*std2) nanmean2+(K*std2)]),'c');
ylim([Ylimitas])
xlim([Xlimitas])
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_impvopto_Tongue','png')
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_impvopto_Tongue','emf')

% Whisker imp opto
figure, hold on,
title(['Whisker imp (n=' num2str(nimp) ')'  'opto (n=' num2str(nimpopto)  ')'])
K=2
X=[1:1:100]';
nanmean1 = nanmean((impW'),2);
nanmean2 = nanmean((impWopto'),2);
std1 = nanmean((impWsem'),2);
std2 = nanmean((impWoptosem'),2);
plot(abs(nanmean1),'m','LineWidth',2);
plot(abs(nanmean2),'c','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
plotshaded(X', abs([nanmean1-(K*std1) nanmean1+(K*std1)]),'m');
plotshaded(X', abs([nanmean2-(K*std2) nanmean2+(K*std2)]),'c');
ylim([Ylimitas])
xlim([Xlimitas])
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_impvopto_Whisker','png')
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_impvopto_Whisker','emf')

% PLOT OMI
% SNOUT omi opto
figure, hold on,
title(['SNOUT omi (n=' num2str(nomi) ')'  'opto (n=' num2str(nomiopto)  ')'])
K=2
X=[1:1:100]';
nanmean1 = nanmean((omiN'),2);
nanmean2 = nanmean((omiNopto'),2);
std1 = nanmean((omiNsem'),2);
std2 = nanmean((omiNoptosem'),2);
plot(abs(nanmean1),'g','LineWidth',2);
plot(abs(nanmean2),'c','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
plotshaded(X', abs([nanmean1-(K*std1) nanmean1+(K*std1)]),'g');
plotshaded(X', abs([nanmean2-(K*std2) nanmean2+(K*std2)]),'c');
ylim([Ylimitas])
xlim([Xlimitas])
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_omivopto_Snout','png')
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_omivopto_Snout','emf')

% Tongue omi opto
figure, hold on,
title(['Tongue omi (n=' num2str(nomi) ')'  'opto (n=' num2str(nomiopto)  ')'])
K=2
X=[1:1:100]';
nanmean1 = nanmean((omiT'),2);
nanmean2 = nanmean((omiTopto'),2);
std1 = nanmean((omiTsem'),2);
std2 = nanmean((omiToptosem'),2);
plot(abs(nanmean1),'g','LineWidth',2);
plot(abs(nanmean2),'c','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
plotshaded(X', abs([nanmean1-(K*std1) nanmean1+(K*std1)]),'g');
plotshaded(X', abs([nanmean2-(K*std2) nanmean2+(K*std2)]),'c');
ylim([Ylimitas])
xlim([Xlimitas])
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_omivopto_Tongue','png')
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_omivopto_Tongue','emf')

% Whisker omi opto
figure, hold on,
title(['Whisker omi (n=' num2str(nomi) ')'  'opto (n=' num2str(nomiopto)  ')'])
K=2
X=[1:1:100]';
nanmean1 = nanmean((omiW'),2);
nanmean2 = nanmean((omiWopto'),2);
std1 = nanmean((omiWsem'),2);
std2 = nanmean((omiWoptosem'),2);
plot(abs(nanmean1),'g','LineWidth',2);
plot(abs(nanmean2),'c','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);
plotshaded(X', abs([nanmean1-(K*std1) nanmean1+(K*std1)]),'g');
plotshaded(X', abs([nanmean2-(K*std2) nanmean2+(K*std2)]),'c');
ylim([Ylimitas])
xlim([Xlimitas])
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_omivopto_Whisker','png')
saveas(gcf, 'D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Fig6OPTO_omivopto_Whisker','emf')




