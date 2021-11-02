% pub_fig_BehavOpto_BarMeanStat_JCscript
% Written by Julien Catanese 10/26/2018
% last update: by JC 10/27/2018

%% Define Session list
clear all
close all
cd('D:\JC_Analysis');
% SessList = dir(['**/*' '_' 'taskopto' '_' '*']);
SessList = dir('**/*task*_optopost_*mW*');
NSess= max(size(SessList)) % Number of Sessions
Nol=[]; Nolo=[]; cCR=[]; cCRo=[]; cCL=[]; cCLo=[];
imp=[]; impo=[]; impR=[]; impRo=[]; impL=[]; impLo=[];
erP=[]; erPo=[]; erPR=[]; erPRo=[]; erPL=[]; erPLo=[];
cor=[]; coro=[]; TOTAL_opt_trial = []; TOTAL_imp_trial = []; TOTAL_omi_trial = []; TOTAL_cor_trial = [];TOTAL_eSi_trial=[]; 

%% loop trhough all Sessions named "taskopto"
for of=1:NSess
    SessID= SessList(of).name
    SessPath=SessList(of).folder;
    cd([SessPath '\' SessID])
    
    load ('Ntrial_type.mat')
    
    Ncorrect = trial.Nb_correct_L + trial.Nb_correct_R; % only trial without opto stim
    Ncorrect_opto = trial.Nb_correct_L_opto + trial.Nb_correct_R_opto; % only trial without opto stim
    Ncorrect_tot = Ncorrect + Ncorrect_opto;
    
    Nnolick = trial.Nb_NoLick; % Nolick non-opto
    Nnolick_opto = trial.Nb_NoLick_opto; %Nolick opto
    Nnolick_tot = Nnolick + Nnolick_opto;
    
    Nimpulse =   trial.Nb_errorDelay_PL_CL + trial.Nb_errorDelay_PR_CR + trial.Nb_errorDelay_PL_CR + trial.Nb_errorDelay_PR_CL;
    Nimpulse_opto = trial.Nb_errorDelay_PL_CL_opto + trial.Nb_errorDelay_PR_CR_opto + trial.Nb_errorDelay_PL_CR_opto + trial.Nb_errorDelay_PR_CL_opto;
    Nimpulse_tot = Nimpulse + Nimpulse_opto;
    
    Nimpulse_CR =   trial.Nb_errorDelay_PR_CR + trial.Nb_errorDelay_PL_CR ;
    Nimpulse_CR_opto =  trial.Nb_errorDelay_PR_CR_opto + trial.Nb_errorDelay_PL_CR_opto ;
    Nimpulse_CL =   trial.Nb_errorDelay_PL_CL + trial.Nb_errorDelay_PR_CL;
    Nimpulse_CL_opto = trial.Nb_errorDelay_PL_CL_opto + trial.Nb_errorDelay_PR_CL_opto;
    
    NerrSide =  trial.Nb_errorResp_PR +  trial.Nb_errorResp_PL;
    NerrSide_opto = trial.Nb_errorResp_PR_opto +  trial.Nb_errorResp_PL_opto;
    NerrSide_tot = NerrSide + NerrSide_opto;
    
    NerrSide_contra =   trial.Nb_errorResp_PR;
    NerrSide_ipsi = trial.Nb_errorResp_PL ;
    NerrSide_contra_opto = trial.Nb_errorResp_PR_opto;
    NerrSide_ipsi_opto = trial.Nb_errorResp_PL_opto;
    
    Ntrial_tot = trial.Ntrial;
    Ntrial_opto = trial.Nb_all_opto ;
    Ntrial = Ntrial_tot - Ntrial_opto;
    
    TOTAL_opt_trial =  [TOTAL_opt_trial  Ntrial_opto];
    TOTAL_imp_trial = [TOTAL_imp_trial  Nimpulse ];
    TOTAL_omi_trial = [TOTAL_omi_trial  Nnolick ];
    TOTAL_cor_trial = [TOTAL_cor_trial  Ncorrect ];
    TOTAL_eSi_trial = [TOTAL_eSi_trial  NerrSide ];

    %% Convert in percent
    cor = [cor 100*Ncorrect/Ntrial];
    coro=[coro 100*Ncorrect_opto/Ntrial_opto];
    
    Nol=[Nol  100*Nnolick/Ntrial];
    Nolo=[Nolo 100*Nnolick_opto/Ntrial_opto];
    
    imp=[imp  100*Nimpulse/Ntrial];
    impo = [impo 100*Nimpulse_opto/Ntrial_opto];
    
    erP = [erP 100*NerrSide/Ntrial] ;
    erPo = [erPo 100*NerrSide_opto/Ntrial_opto];
    
    
    cCR=[cCR 100*trial.Nb_correct_R/Ncorrect_tot];
    cCL=[cCL 100*trial.Nb_correct_L/Ncorrect_tot];
    cCRo=[cCRo 100*trial.Nb_correct_R_opto/Ncorrect_tot];
    cCLo=[cCLo 100*trial.Nb_correct_L_opto/Ncorrect_tot];
    
    impR=[impR  100*Nimpulse_CR/Nimpulse_tot];
    impRo=[impRo  100*Nimpulse_CR_opto/Nimpulse_tot];
    impL=[impL  100*Nimpulse_CL/Nimpulse_tot];
    impLo=[impLo  100*Nimpulse_CL_opto/Nimpulse_tot];
    
    erPR= [erPR 100*NerrSide_contra/NerrSide_tot] ;
    erPRo= [erPRo 100*NerrSide_contra_opto/NerrSide_tot]  ;
    erPL= [erPL 100*NerrSide_ipsi/NerrSide_tot] ;
    erPLo= [erPLo 100*NerrSide_ipsi_opto/NerrSide_tot]  ;
end

%% STATS : unpaired ttest (trials stim and no stim are independent)

[Pcor,Hcor] = ranksum(cor,coro, 'alpha' ,0.05)
[PerP,HerP] = ranksum(erP,erPo, 'alpha' ,0.05)
[Pnol,Hnol] = ranksum(Nol,Nolo, 'alpha' ,0.05)
[Pimp,Himp] = ranksum(imp,impo, 'alpha' ,0.05)

[P_cor_eSi,H_cor_eSi] = ranksum(cor,erP, 'alpha' ,0.05)

%% Display Count
tot_opt = sum(TOTAL_opt_trial)
tot_imp = sum(TOTAL_imp_trial)
tot_omi = sum(TOTAL_omi_trial)
tot_cor = sum(TOTAL_cor_trial)
tot_eSi = sum(TOTAL_eSi_trial) 

mean_opt = mean(TOTAL_opt_trial)
mean_imp = mean(TOTAL_imp_trial)
mean_omi = mean(TOTAL_omi_trial)
mean_cor = mean(TOTAL_cor_trial)
mean_eSi = mean(TOTAL_eSi_trial) 


%% Calculate mean and std
cor_mean = mean(cor);
cor_std =std(cor,1);
coro_mean= mean(coro);;
coro_std =std(coro,1);
;
Nol_mean=mean(Nol);
Nol_std=std(Nol,1);
Nolo_mean=mean(Nolo);
Nolo_std=std(Nolo,1);

imp_mean=mean(imp);
imp_std=std(imp,1);
impo_mean=mean(impo);
impo_std=std(impo,1);

erP_mean = mean(erP);
erPo_mean = mean(erPo);
erP_std = std(erP);
erPo_std = std(erPo);
;
cCR_mean=mean(cCR);
cCR_std=std(cCR,1);
cCRo_mean=mean(cCRo);
cCRo_std=std(cCRo,1);
cCL_mean=mean(cCL);
cCL_std=std(cCL,1);
cCLo_mean=mean(cCLo);
cCLo_std=std(cCLo,1);

impR_mean=mean(impR);
impR_std=std(impR,1);
impRo_mean=mean(impRo);
impRo_std=std(impRo,1);
impL_mean=mean(impL);
impL_std=std(impL,1);
impLo_mean=mean(impLo);
impLo_std=std(impLo,1);

erPR_mean = mean(erPR);
erPRo_mean = mean(erPRo) ;
erPL_mean = mean(erPL) ;
erPLo_mean = mean(erPLo);
erPR_std = std(erPR,1);
erPRo_std = std(erPRo,1) ;
erPL_std = std(erPL,1) ;
erPLo_std = std(erPLo,1);



%% PLOT1 bar plot with Std
close all
%Figure
figH = figure;
axes1 = axes;
hold on,
y =         [cor_mean,  coro_mean;  erP_mean,   erPo_mean;  imp_mean,   impo_mean;  Nol_mean,   Nolo_mean];
std_dev =   [cor_std,   coro_std;   erP_std,    erPo_std;   imp_std,    impo_std;   Nol_std,    Nolo_std ];
num = 4; %number of different subcategories
c = 1:num;
%Bar(s) %You can not color differently the same bar.
for i = 1:num
    bar(c(i)-0.15,y(i,1),0.2, 'FaceColor','k','EdgeColor','k','LineWidth',1.5);
    bar(c(i)+0.15,y(i,2),0.2, 'FaceColor','y','EdgeColor','c','LineWidth',1.5);
end
%Errorbar
errH1 = errorbar(c-0.15,y(:,1),std_dev(:,1),'.','Color','k');
errH2 = errorbar(c+0.15,y(:,2),std_dev(:,2),'.','Color','k');
errH1.LineWidth = 1;
errH2.LineWidth = 1;
% errH1.Color = [1 0.5 0];
% errH2.Color = [1 0.3 1];
% Set x-ticks
set(axes1,'Ylim',[0 100]);
set(axes1,'Xlim',[0.5 4.5]);
set(axes1,'XTick',[1 1.5 2 2.5 3 3.5 4 4.5],'XTickLabel',...
    {'Correct',' ','Error Side',' ','Error Impulse ',' ','Error Omission', ' '});
% plot stat
SigX = [];
if Hcor == 1; SigX=[SigX 1]; end;
if HerP == 1; SigX=[SigX 2]; end;
if Himp == 1; SigX=[SigX 3]; end;
if Hnol == 1; SigX=[SigX 4]; end;

if ~isempty(SigX);
    SigY = ones(size(SigX))*70;
    plot(SigX,SigY, '*k', 'LineWidth',1);
end

% Label Legend Title
xlabel('trials types');
ylabel('Percent (% among trial type)');
title (['Average over ' num2str(NSess) ' sessions (n=4)' ]) ;
legend('No-Stim trials','Opto-Stim trials' ,'Location','northwest');
%% SAVING
saveas(gcf,['BarPlot_4TrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'jpeg');
saveas(gcf,['D:\JC_Figures\BarPlot_4TrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'tif');
saveas(gcf,['D:\JC_Figures\BarPlot_4TrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'eps');
saveas(gcf,['D:\JC_Figures\BarPlot_4TrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'emf');






