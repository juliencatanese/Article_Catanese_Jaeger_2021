clear all
close all
cd('C:\Users\catan\Documents\EMORY\JC_Analysis')
% cd('C:\Users\Julien\Documents\WORK\JC_Analysis')
% cd('D:\JC_Analysis');
% SessList = dir(['**/*' '_' 'taskopto' '_' '*']);
% SessList = dir('**/*taskopto*');
SessList = dir('*/*task*');
NSess= max(size(SessList)) % Number of Sessions
    Ntot_all = [] ;
    Ncor_all = [] ;
    NeSi_all =[];
    Nomi_all = [] ;
    Nimp_all = [] ;

% loop trhough all Sessions named "taskopto"
for of=1:NSess
    
    SessID= SessList(of).name
    SessPath=SessList(of).folder;
    cd([SessPath '\' SessID])
    
    load ('info.mat')
    load ('Ntrial_type.mat')
    
    Ncorrect = trial.Nb_correct_L + trial.Nb_correct_R; % only trial without opto stim
    
    NerSide = trial.Nb_errorResp_PL + trial.Nb_errorResp_PR
    
    NerSide2 =  NerSide + trial.Nb_errorDelay_PL_CR + trial.Nb_errorDelay_PR_CL;
    
    Nnolick = trial.Nb_NoLick
    
    Nimpulse =   trial.Nb_errorDelay_PL_CL + trial.Nb_errorDelay_PR_CR 
    
    Nerror_tot = NerSide + Nnolick + Nimpulse
    
    Ntot = Ncorrect + NerSide + Nnolick + Nimpulse
    %% PLOT
%     figure;
%     c = categorical({'Correct Side';'Error Side'});
%     y = [      100*Ncorrect /(NerSide2 + Ncorrect) ; 100*NerSide2 /(NerSide2 + Ncorrect) ];
%     b = bar(c,y, 'FaceColor','flat'); ylabel('% N/(NerSide + Ncorrect)')
%     saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day  '_'  info.info_notes.task '_Bar_2PercSide_cor_eSi'],'jpeg');
%     
%     figure;
%     c = categorical({'Error1 Side'; 'Error2 Impulse' ; 'Error3 Omission' });
%     y = [       100*NerSide /Ntot ; 100*Nimpulse/Ntot ; 100*Nnolick/Ntot ];
%     b = bar(c,y, 'FaceColor','flat'); ylabel('% N/(Ntot)')
%     saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day  '_'  info.info_notes.task '_Bar_3PercErrors_imp_omi_eSi'],'jpeg');
    
    %%
    Ntot_all = [Ntot_all Ntot] ;
    Ncor_all = [Ncor_all Ncorrect] ;
    NeSi_all =[NeSi_all NerSide2] ;
    Nomi_all = [Nomi_all Nnolick] ;
    Nimp_all = [Nimp_all Nimpulse] ;
    
end
%%

[H1cor P1cor] = ttest(((Ncor_all./(Ncor_all+NeSi_all))/0.5)*100, 100)
[H1eSi P1eSi] = ttest(((NeSi_all./(Ncor_all+NeSi_all))/0.5)*100, 100)
fffff
%%

    Ntot_all = [Ntot_all(1) Ntot_all(2) mean(Ntot_all(3:6)) mean(Ntot_all(7:10)) mean(Ntot_all(11:15))] ;
    Ncor_all = [Ncor_all(1) Ncor_all(2) mean(Ncor_all(3:6)) mean(Ncor_all(7:10)) mean(Ncor_all(11:15))] ; 
    NeSi_all = [NeSi_all(1) NeSi_all(2) mean(NeSi_all(3:6)) mean(NeSi_all(7:10)) mean(NeSi_all(11:15))] ; 
    Nomi_all = [Nomi_all(1) Nomi_all(2) mean(Nomi_all(3:6)) mean(Nomi_all(7:10)) mean(Nomi_all(11:15))] ; 
    Nimp_all = [Nimp_all(1) Nimp_all(2) mean(Nimp_all(3:6)) mean(Nimp_all(7:10)) mean(Nimp_all(11:15))] ; 


%% PLOT1 bar plot with Std
Mcor = mean(Ncor_all./ (Ncor_all+NeSi_all))
std_cor= std(Ncor_all./ (Ncor_all+NeSi_all))
MeSi = mean(NeSi_all ./ (Ncor_all+NeSi_all))
std_eSi= std(NeSi_all./ (Ncor_all+NeSi_all))

% [Pcor Hcor] = ranksum(Ncor_all, NeSi_all)

[H1cor P1cor] = ttest(((Ncor_all./(Ncor_all+NeSi_all))/0.5)*100, 100, 'tail', 'right')
[H1eSi P1eSi] = ttest(((NeSi_all./(Ncor_all+NeSi_all))/0.5)*100, 100, 'tail', 'left')
[H2 P2] = ttest(((Ncor_all./(Ncor_all+NeSi_all))/0.5)*100, ((NeSi_all./(Ncor_all+NeSi_all))/0.5)*100)
%%

[H1cor P1cor] = ttest(((Ncor_all./(Ncor_all+NeSi_all))/0.5)*100, 100)
[H1eSi P1eSi] = ttest(((NeSi_all./(Ncor_all+NeSi_all))/0.5)*100, 100)


%%
[H1cor P1cor] = ttest(Ncor_all, mean((Ncor_all+NeSi_all)./2), 'tail', 'right')
[H1eSi P1eSi] = ttest(NeSi_all, mean((Ncor_all+NeSi_all)./2), 'tail', 'left')


[H1cor P1cor] = ttest(Ncor_all./(Ncor_all+NeSi_all), 0.5, 'tail', 'right')
[H1eSi P1eSi] = ttest(NeSi_all./(Ncor_all+NeSi_all), 0.5, 'tail', 'left')

[H1cor P1cor] = ttest(Ncor_all, mean((Ncor_all+NeSi_all)./2), 'tail', 'right')
[H1eSi P1eSi] = ttest(NeSi_all, mean((Ncor_all+NeSi_all)./2), 'tail', 'left')



%%



y =       [Mcor*100,  MeSi*100];
std_dev = [std_cor*100,   std_eSi*100];
X = [ 0.25, 0.50] ; %number of different subcategories
figure, 
axes1=gca, 
hold on, B1= bar(X(1), y(1) , 0.2, 'FaceColor','k','EdgeColor','k','LineWidth',1.5);
hold on, E1 = errorbar(X(1),y(1),std_dev(1),'.','Color','r');
hold on, B2 = bar(X(2), y(2) , 0.2, 'FaceColor','k','EdgeColor','k','LineWidth',1.5);
hold on, E2 = errorbar(X(2), y(2),std_dev(2),'.','Color','r');


E1.LineWidth = 1;
E2.LineWidth = 1;
set(axes1,'Xlim',[0.1 0.65]);
set(axes1,'XTick',[0.25 0.5],'XTickLabel',{ 'Correct','Error Side'});
% plot stat
SigX = [];
if Hcor; SigX=[SigX 0.38]; end;

if ~isempty(SigX);
    SigY = ones(size(SigX))*95;
    plot(SigX,SigY, '*k', 'LineWidth',1);
end

% SAVING
mkdir D:\JC_Figures\behavior;
saveas(gcf,['D:\JC_Figures\Average_Bar_2PercSide_cor_eSi'],'png');
saveas(gcf,['D:\JC_Figures\Average_Bar_2PercSide_cor_eSi'],'eps');
saveas(gcf,['D:\JC_Figures\Average_Bar_2PercSide_cor_eSi'],'emf');

%% PLOT1 bar plot with Std
Mimp = mean(Nimp_all./ (Ntot_all))
MeSi = mean(NeSi_all./ (Ntot_all))
Momi = mean(Nomi_all./ (Ntot_all))

Simp = std(Nimp_all./ (Ntot_all))
SeSi = std(NeSi_all./ (Ntot_all))
Somi = std(Nomi_all./ (Ntot_all))

[Himp Pimp] = ttest(Nimp_all)
[HeSi PeSi] = ttest(NeSi_all)
[Homi Pomi] = ttest(Nomi_all)

y =       [MeSi*100,  Mimp*100, Momi*100];
std_dev = [SeSi*100,  Simp*100, Somi*100];
X = [ 0.25, 0.50, 0.75] ; %number of different subcategories
figure, 
axes1=gca, 
hold on, B1= bar(X(1), y(1) , 0.2, 'FaceColor','k','EdgeColor','k','LineWidth',1.5);
hold on, E1 = errorbar(X(1),y(1),std_dev(1),'.','Color','r');
hold on, B2 = bar(X(2), y(2) , 0.2, 'FaceColor','k','EdgeColor','k','LineWidth',1.5);
hold on, E2 = errorbar(X(2), y(2),std_dev(2),'.','Color','r');
hold on, B3 = bar(X(3), y(3) , 0.2, 'FaceColor','k','EdgeColor','k','LineWidth',1.5);
hold on, E3 = errorbar(X(3), y(3),std_dev(3),'.','Color','r');

E1.LineWidth = 1;
E2.LineWidth = 1;
E3.LineWidth = 1;
set(axes1,'Xlim',[0.1 0.9]);
set(axes1,'XTick',[0.25 0.5 0.75],'XTickLabel',{ 'eSi', 'imp', 'omi'});
% plot stat
SigX = [];
if HeSi; SigX=[SigX 0.25]; end;
if Himp; SigX=[SigX 0.5]; end;
if Homi; SigX=[SigX 0.75]; end;

if ~isempty(SigX);
    SigY = ones(size(SigX))*100;
    plot(SigX,SigY, '*k', 'LineWidth',1);
end

% SAVING
mkdir D:\JC_Figures\behavior;
saveas(gcf,['D:\JC_Figures\Average_Bar_3PercErors_eSi_imp_omi'],'png');
saveas(gcf,['D:\JC_Figures\Average_Bar_3PercErrors_eSi_imp_omi'],'eps');
saveas(gcf,['D:\JC_Figures\Average_Bar_3PercErrors_eSi_imp_omi'],'emf');

%% Additional data training 
TS = [0 3 4 4 5 2 2 0; 2 6 0 3 3 3 4 6; 0 0 3 3 5 6 4 4; 3 5 4 2 4 5 4 0; 4 4 4 4 4 4 4 0]
NsessTOT= sum(TS') 
AvPerMouse = mean(NsessTOT)
StdPerMouse = std(NsessTOT)


average= mean(TS') 
standev = std(TS')