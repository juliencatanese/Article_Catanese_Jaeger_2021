% BarPlot_4TrialType_MeanSTAT_BehavOpto_JCscript
% Written by Julien Catanese 10/26/2018
% last update: by JC 10/27/2018

%% Define Session list
clear all
close all
cd('C:\Users\catan\Documents\EMORY\JC_Analysis');
% SessList = dir(['**/*' '_' 'taskopto' '_' '*']);
SessList = dir('**/*VM_task*');
NSess= max(size(SessList)) % Number of Sessions
cor=[]; coro=[]; corL=[];corLo=[]; corR=[]; corRo=[];
omi=[]; omio=[]; omiL=[]; omiLo=[]; omiR=[]; omiRo=[];
imp=[]; impo=[]; impR=[]; impRo=[]; impL=[]; impLo=[];
eSi=[]; eSio=[]; eSiR=[]; eSiRo=[]; eSiL=[]; eSiLo=[];

cor_RLratio=[]; coro_RLratio =[]; imp_RLratio = [];impo_RLratio = []; 
eSi_RLratio = []; eSio_RLratio = []; omi_RLratio = []; omio_RLratio = []
cor_perc=[]; coro_perc =[]; imp_perc = [];impo_perc = []; 
eSi_perc = []; eSio_perc = []; omi_perc = []; omio_perc = []

%% loop trhough all Sessions named "taskopto"
for of=1:NSess
    SessID= SessList(of).name
    SessPath=SessList(of).folder;
    cd([SessPath '\' SessID])
    
    load ('Ntrial_type.mat')
    
    Ncorrect = trial.Nb_correct_L + trial.Nb_correct_R; % only trial without opto stim
    Ncorrect_opto = trial.Nb_correct_L_opto + trial.Nb_correct_R_opto; % only trial without opto stim
    Ncorrect_tot = Ncorrect + Ncorrect_opto;
    
    Ncorrect_L = trial.Nb_correct_L; 
    Ncorrect_Lopto = trial.Nb_correct_L_opto
    Ncorrect_R = trial.Nb_correct_R;
    Ncorrect_Ropto = trial.Nb_correct_R_opto
    
    Nnolick = trial.Nb_NoLick; % Nolick non-opto
    Nnolick_opto = trial.Nb_NoLick_opto; %Nolick opto
    Nnolick_tot = Nnolick + Nnolick_opto;
    
    Nnolick_R =   trial.Nb_NoLick_R;
    Nnolick_R_opto =  trial.Nb_NoLick_R_opto;
    Nnolick_L =   trial.Nb_NoLick_L ;
    Nnolick_L_opto = trial.Nb_NoLick_L_opto;
    
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
    
    NerrSide_PR =   trial.Nb_errorResp_PR;
    NerrSide_PL = trial.Nb_errorResp_PL ;
    NerrSide_PR_opto = trial.Nb_errorResp_PR_opto;
    NerrSide_PL_opto = trial.Nb_errorResp_PL_opto;
    
    Ntrial = trial.Ntrial - trial.Nb_all_opto;
    Ntrial_opto = trial.Nb_all_opto ;
    Ntrial_tot = trial.Ntrial;
    
    %%
    cor = [cor Ncorrect];
    coro=[coro Ncorrect_opto];
    corL = [corL Ncorrect_L]  
    corR = [corR Ncorrect_R]
    corLo = [corLo Ncorrect_Lopto]  
    corRo = [corRo Ncorrect_Ropto]
    
    eSi = [eSi NerrSide] ;
    eSio = [eSio NerrSide_opto]; 
    eSiR = [eSiR NerrSide_PR] ;
    eSiRo = [eSiRo NerrSide_PR_opto]; 
    eSiL = [eSiL NerrSide_PL] ;
    eSiLo = [eSiLo NerrSide_PL_opto]; 
    
    imp=[imp  Nimpulse];
    impo = [impo Nimpulse_opto];
    impL=[impL  Nimpulse_CL];
    impLo = [impLo Nimpulse_CL_opto];
    impR=[impR  Nimpulse_CR];
    impRo = [impRo Nimpulse_CR_opto];
    
    omi=[omi  Nnolick];
    omio=[omio Nnolick_opto];
    omiL=[omiL  Nnolick_L];
    omiLo=[omiLo Nnolick_L_opto];
    omiR=[omiR  Nnolick_R];
    omiRo=[omiRo Nnolick_R_opto];
    
    %% Convert in percent
    
    cor_perc  = [cor_perc 100*Ncorrect/Ntrial];
    coro_perc =[coro_perc 100*Ncorrect_opto/Ntrial_opto];
    
    omi_perc =[omi_perc  100*Nnolick/Ntrial];
    omio_perc =[omio_perc 100*Nnolick_opto/Ntrial_opto];
    
    imp_perc=[imp  100*Nimpulse/Ntrial];
    impo_perc = [impo 100*Nimpulse_opto/Ntrial_opto];
    
    eSi_perc = [eSi_perc  100*NerrSide/Ntrial] ;
    eSio_perc = [eSio_perc 100*NerrSide_opto/Ntrial_opto];
    
    cor_RLratio=[cor_RLratio ((trial.Nb_correct_R+1)/(trial.Nb_correct_L+1))-1];
    coro_RLratio =[coro_RLratio (trial.Nb_correct_R_opto+1)/(trial.Nb_correct_L_opto+1)];
    
    imp_RLratio = [imp_RLratio ((Nimpulse_CR+1)/(Nimpulse_CL+1))-1]
    impo_RLratio = [impo_RLratio (Nimpulse_CR_opto+1)/(Nimpulse_CL_opto+1)]
    
    eSi_RLratio = [eSi_RLratio ((NerrSide_PR+1)/(NerrSide_PL+1))-1]
    eSio_RLratio = [eSio_RLratio (NerrSide_PR_opto+1)/(NerrSide_PL_opto+1)]
    
    omi_RLratio = [omi_RLratio ((Nnolick_R+1)/(Nnolick_L+1))-1]
    omio_RLratio = [omio_RLratio (Nnolick_R_opto+1)/(Nnolick_L_opto+1)]
    
end


%% STATS : 2 sample ttest 

A= corL-corR;
B= eSiL-eSiR;
C= impL-impR;
D= omiL-omiR;

[HA,PA] = ttest(A) % (alpha 5%)
[HB,PB] = ttest(B)  % (alpha 5%)
[HC,PC] = ttest(C) % (alpha 5%)
[HD,PD] = ttest(D) % (alpha 5%)


%% PLOT2 bar plot: Contra-ipsi Ratio 
figH2 = figure; axes2 = axes; hold on,
y= [ mean(A);  mean(B); mean(C); mean(D) ];
std_dev = [ std(A);  std(B); std(C); std(D) ];
num = 4; %number of different subcategories
c = 1:num;
for i = 1:num
    bar(c(i)-0.15,y(i,1),0.2, 'FaceColor','k','EdgeColor','k','LineWidth',1.5);
%     bar(c(i)+0.15,y(i,2),0.2, 'FaceColor','y','EdgeColor','c','LineWidth',1.5);
end
errH1 = errorbar(c-0.15,y(:,1),std_dev(:,1),'.','Color','k');
% errH2 = errorbar(c+0.15,y(:,2),std_dev(:,2),'.','Color','k');
errH1.LineWidth = 1;
% errH2.LineWidth = 1;
set(axes2,'XTick',[1 1.5 2 2.5 3 3.5 4],'XTickLabel',...
    {'Correct', '', 'errorSide', '', 'impulse','','omission'});
SigX = [];
if HA == 1; SigX=[SigX 1-0.15]; end;
if HB == 1; SigX=[SigX 2-0.15]; end;
if HC == 1; SigX=[SigX 3-0.15]; end;
if HD == 1; SigX=[SigX 4-0.15]; end;

if ~isempty(SigX)
    SigY = ones(size(SigX))*10
    plot(SigX,SigY, '*k', 'LineWidth',1)
end
% Label Legend Title
xlabel('trials types');
ylabel('ratio (contra/ipsi)');
title (['Average over ' num2str(NSess) ' sessions (n=4)' ]) ;


