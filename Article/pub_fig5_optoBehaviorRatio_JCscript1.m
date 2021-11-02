% BarPlot_4TrialType_MeanSTAT_BehavOpto_JCscript
% Written by Julien Catanese 10/26/2018
% last update: by JC 10/27/2018

%% Define Session list
clear all
close all
cd('C:\Users\catan\Documents\EMORY\JC_Analysis');
% SessList = dir(['**/*' '_' 'taskopto' '_' '*']);
SessList = dir('**/*taskopto*');
NSess= max(size(SessList)) % Number of Sessions
Nol=[]; Nolo=[]; cCR=[]; cCRo=[]; cCL=[]; cCLo=[];
imp=[]; impo=[]; impR=[]; impRo=[]; impL=[]; impLo=[];
erP=[]; erPo=[]; erPR=[]; erPRo=[]; erPL=[]; erPLo=[];
cor=[]; coro=[]; 
CRvsCL=[]; CRovsCLo =[]; impRvsimpL = [];impRovsimpLo = []; 
erPRvserPL = []; erPRovserPLo = []; omiRvsomiL = []; omiRovsomiLo = []
%% loop trhough all Sessions named "taskopto"
for of=1:NSess-1
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
    
    
    %% Convert in percent
    
    cor = [cor 100*Ncorrect/Ntrial];
    coro=[coro 100*Ncorrect_opto/Ntrial_opto];
    
    Nol=[Nol  100*Nnolick/Ntrial];
    Nolo=[Nolo 100*Nnolick_opto/Ntrial_opto];
    
    imp=[imp  100*Nimpulse/Ntrial];
    impo = [impo 100*Nimpulse_opto/Ntrial_opto];
    
    erP = [erP 100*NerrSide/Ntrial] ;
    erPo = [erPo 100*NerrSide_opto/Ntrial_opto];
    
    CRvsCL=[CRvsCL (trial.Nb_correct_R+1)/(trial.Nb_correct_L+1)];
    CRovsCLo =[CRovsCLo (trial.Nb_correct_R_opto+1)/(trial.Nb_correct_L_opto+1)];
    
    impRvsimpL = [impRvsimpL (Nimpulse_CR+1)/(Nimpulse_CL+1)]
    impRovsimpLo = [impRovsimpLo (Nimpulse_CR_opto+1)/(Nimpulse_CL_opto+1)]
    
    erPRvserPL = [erPRvserPL (NerrSide_PR+1)/(NerrSide_PL+1)]
    erPRovserPLo = [erPRovserPLo (NerrSide_PR_opto+1)/(NerrSide_PL_opto+1)]
    
    omiRvsomiL = [omiRvsomiL (Nnolick_R+1)/(Nnolick_L+1)]
    omiRovsomiLo = [omiRovsomiLo (Nnolick_R_opto+1)/(Nnolick_L_opto+1)]
    
end


%% Calculate mean and std
cor_mean = mean(cor)
cor_std =std(cor,1)
coro_mean= mean(coro)
coro_std =std(coro,1)

Nol_mean=mean(Nol);
Nol_std=std(Nol,1);
Nolo_mean=mean(Nolo);
Nolo_std=std(Nolo,1);

imp_mean=mean(imp);
imp_std=std(imp,1);
impo_mean=mean(impo);
impo_std=std(impo,1);

erP_mean = mean(erP)
erPo_mean = mean(erPo)
erP_std = std(erP)
erPo_std = std(erPo)

CRvsCL_mean= mean(CRvsCL) 
CRovsCLo_mean = mean(CRovsCLo) 

impRvsimpL_mean = mean(impRvsimpL) 
impRovsimpLo_mean = mean(impRovsimpLo) 

erPRvserPL_mean = mean(erPRvserPL)
erPRovserPLo_mean = mean(erPRovserPLo)

omiRvsomiL_mean = mean(omiRvsomiL)
omiRovsomiLo_mean = mean(omiRovsomiLo)

CRvsCL_std= std(CRvsCL,1) 
CRovsCLo_std = std(CRovsCLo,1) 

impRvsimpL_std = std(impRvsimpL,1) 
impRovsimpLo_std = std(impRovsimpLo,1) 

erPRvserPL_std = std(erPRvserPL,1)
erPRovserPLo_std = std(erPRovserPLo,1)

omiRvsomiL_std = std(omiRvsomiL,1)
omiRovsomiLo_std = std(omiRovsomiLo,1)

%% STATS : 2 sample ttest 

[Hcor,Pcor] = ttest2(cor,coro, 'alpha' ,0.025)
[HerP,PerP] = ttest2(erP,erPo, 'alpha' ,0.025)
[Hnol,Pnol] = ttest2(Nol,Nolo, 'alpha' ,0.025)
[Himp,Pimp] = ttest2(imp,impo, 'alpha' ,0.025)

[Hcor_RvLratio,Pcor_RvLratio] = ttest2(CRvsCL,CRovsCLo, 'alpha' ,0.1)
[HerP_RvLratio,PerP_RvLratio] = ttest2(erPRvserPL,erPRovserPLo, 'alpha' ,0.1)
[Himp_RvLratio,Pimp_RvLratio] = ttest2(impRvsimpL,impRovsimpLo, 'alpha' ,0.1)
[Himp_RvLratio,Pimp_RvLratio] = ttest2(omiRvsomiL,omiRovsomiLo, 'alpha' ,0.1)


A= CRovsCLo-CRvsCL;
B= erPRovserPLo-erPRvserPL;
C= impRovsimpLo-impRvsimpL;
D= omiRovsomiLo-omiRvsomiL;

[HA0,PA0] = ttest(A) % (alpha 5%)
[HB0,PB0] = ttest(B)  % (alpha 5%)
[HC0,PC0] = ttest(C) % (alpha 5%)
[HD0,PD0] = ttest(D) % (alpha 5%)

[HAB,PAB] = ttest(A,B) % (alpha 5%)
[HBC,PBC] = ttest(B,C) % (alpha 5%)
[HAC,PAC] = ttest(A,C) % (alpha 5%)
[HDC,PDC] = ttest(D,C) % (alpha 5%)

%% PLOT1 bar plot : Mean opto vs Non opto with Std
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

if ~isempty(SigX)
    SigY = ones(size(SigX))*70
    plot(SigX,SigY, '*k', 'LineWidth',1)
end

% Label Legend Title
xlabel('trials types');
ylabel('Percent (% among trial type)');
title (['Average over ' num2str(NSess) ' sessions (n=4)' ]) ;
legend('No-Stim trials','Opto-Stim trials' ,'Location','northwest')
%% SAVING
saveas(gcf,['BarPlot_4TrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'jpeg');
saveas(gcf,['D:\JC_Figures\BarPlot_4TrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'tif');
saveas(gcf,['D:\JC_Figures\BarPlot_4TrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'eps');
saveas(gcf,['D:\JC_Figures\BarPlot_4TrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'emf');




%% PLOT2 bar plot: Contra-ipsi Ratio 
figH2 = figure; axes2 = axes; hold on,
y= [ CRvsCL_mean, CRovsCLo_mean ;  ; erPRvserPL_mean, erPRovserPLo_mean; impRvsimpL_mean, impRovsimpLo_mean ];
std_dev =   [CRvsCL_std, CRovsCLo_std ; erPRvserPL_std, erPRovserPLo_std; impRvsimpL_std, impRovsimpLo_std ];
num = 3; %number of different subcategories
c = 1:num;
for i = 1:num
    bar(c(i)-0.15,y(i,1),0.2, 'FaceColor','k','EdgeColor','k','LineWidth',1.5);
    bar(c(i)+0.15,y(i,2),0.2, 'FaceColor','y','EdgeColor','c','LineWidth',1.5);
end
errH1 = errorbar(c-0.15,y(:,1),std_dev(:,1),'.','Color','k');
errH2 = errorbar(c+0.15,y(:,2),std_dev(:,2),'.','Color','k');
errH1.LineWidth = 1;
errH2.LineWidth = 1;
set(axes2,'Ylim',[0 10]);
set(axes2,'Xlim',[0.5 3.5]);
set(axes2,'XTick',[1 1.5 2 2.5 3 3.5],'XTickLabel',...
    {'Correct', ' ', 'errorSide', ' ', 'impulse',''});
SigX = [];
if Hcor_RvLratio == 1; SigX=[SigX 1]; end;
if HerP_RvLratio == 1; SigX=[SigX 2]; end;
if Himp_RvLratio == 1; SigX=[SigX 3]; end;

if ~isempty(SigX)
    SigY = ones(size(SigX))*7
    plot(SigX,SigY, '*k', 'LineWidth',1)
end
% Label Legend Title
xlabel('trials types');
ylabel('ratio (contra/ipsi)');
title (['Average over ' num2str(NSess) ' sessions (n=4)' ]) ;
legend('No-Stim trials','Opto-Stim trials' ,'Location','northwest')
% Saving
saveas(gcf,['BarPlot_ImpTrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'jpeg');
saveas(gcf,['D:\JC_Figures\BarPlot_ratioLRTrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'tif');
saveas(gcf,['D:\JC_Figures\BarPlot_ratioLRTrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'eps');
saveas(gcf,['D:\JC_Figures\BarPlot_ratioLRTrialType_MeanSTAT_' num2str(NSess) 'Sess_BehavOpto'],'emf');




%% PLOT3 bar plot:  Diff Contra-ipsi ratio 

figure(3), close(3), figure(3), hold on

pos=[0.5 1 1.5 2]; 
y= [ mean(A) ; mean(B); mean(C); mean(D)];
std_dev =   [ std(A) ; std(B); std(C); std(D)];

hold on, bar(pos, y ,0.4, 'FaceColor','k','EdgeColor','k','LineWidth',1.5);
hold on, errorbar(pos,y,std_dev,'.','Color','k','LineWidth',1);

set(gca,...
    'Ylim',[-10 10], ... 
    'Xlim',[pos(1)-0.25 pos(4)+0.25],...
    'Xtick',[pos(1) pos(2) pos(3) pos(4)],...
    'XTickLabel',{'Correct', 'errorSide', 'impulse', 'omission'} );

SigX = [];SigX2=[];
if HA0 == 1; SigX=[SigX pos(1)]; end;
if HB0 == 1; SigX=[SigX pos(2)]; end;
if HC0 == 1; SigX=[SigX pos(3)]; end;
if HD0 == 1; SigX=[SigX pos(4)]; end;

if HAB == 1; SigX2=[SigX2; pos(1) pos(2)]; end;
if HBC == 1; SigX2=[SigX2; pos(2) pos(3)]; end;
if HAC == 1; SigX2=[SigX2; pos(1) pos(3)]; end;
if HDC == 1; SigX2=[SigX2; pos(3) pos(4)]; end;

if ~isempty(SigX)
    SigY = ones(size(SigX))*6
    hold on, plot(SigX,SigY, '*k', 'LineWidth',1)
end
if ~isempty(SigX2)
    for bb=1:size(SigX2,1) 
    SigY2 = ones(size(SigX2))*5+bb
    hold on, plot(SigX2(bb,:),SigY2, 'k', 'LineWidth',1)
    hold on, plot(mean((SigX2(bb,:))),mean(SigY2)+0.25, '*k', 'LineWidth',1)

    end
end

% Label Legend Title
ylabel('contra/ipsi ratio diff (stim - nostim)');
title (['Average over ' num2str(NSess) ' sessions (n=4)' ]) ;


