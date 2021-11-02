% BarPlot_4TrialType_1Sess_BehavOpto_JCscript
% Written by Julien Catanese 10/26/2018
% last update: by JC 10/27/2018

%% Calculate trial types
load ('info.mat')
load ('Ntrial_type.mat')

Ncorrect = trial.Nb_correct_L + trial.Nb_correct_R; % only trial without opto stim 
Ncorrect_opto = trial.Nb_correct_L_opto + trial.Nb_correct_R_opto; % only trial without opto stim 
Ncorrect_tot = Ncorrect + Ncorrect_opto;

Nnolick = trial.Nb_NoLick; % Nolick non-opto 
Nnolick_opto = trial.Nb_NoLick_opto; %Nolick opto 
Nnolick_tot = Nnolick + Nnolick_opto

Nimpulse =   trial.Nb_errorDelay_PL_CL + trial.Nb_errorDelay_PR_CR + trial.Nb_errorDelay_PL_CR + trial.Nb_errorDelay_PR_CL;
Nimpulse_opto = trial.Nb_errorDelay_PL_CL_opto + trial.Nb_errorDelay_PR_CR_opto + trial.Nb_errorDelay_PL_CR_opto + trial.Nb_errorDelay_PR_CL_opto;
Nimpulse_tot = Nimpulse + Nimpulse_opto 

Nerror_tot = Nnolick + Nnolick_opto + Nimpulse + Nimpulse_opto 

Nopto_tot = trial.Nb_all_opto ; % all opto trial
Ntrial_tot = trial.Ntrial;

%% PLOT 
figure;
c = categorical({'Error Omission';'Correct Ipsi';'Correct Contra'; 'Error Impulse'});
y = [       100*Nnolick/(Ntrial_tot-Nopto_tot)            ,   100*Nnolick_opto/Nopto_tot                ;...
            100*trial.Nb_correct_L/(Ntrial_tot-Nopto_tot) ,   100*trial.Nb_correct_L_opto/Nopto_tot     ;...
            100*trial.Nb_correct_R/(Ntrial_tot-Nopto_tot) ,   100*trial.Nb_correct_R_opto/Nopto_tot     ;...
            100*Nimpulse/(Ntrial_tot-Nopto_tot)           ,   100*Nimpulse_opto/Nopto_tot              ];
b = bar(c, y, 'FaceColor','flat');

% color
colormap('copper')
for k = 1:size(y,2)
    b(k).CData = k;
end

b(2).LineWidth = 2;
b(2).EdgeColor = 'cyan';

% labels legends titles 
ylabel('Percent (N/Ntot)')
ylim([0 90])

legend('%No Stim', '%Opto Stim', 'Location', 'best')
title ([ info.info_notes.MouseID ''  info.info_notes.Day  ' #correct=' num2str(Ncorrect_tot) '  #omission=' num2str(Nnolick_tot) '  #impulse=' num2str(Nimpulse_tot) '  #opto=' num2str(Nopto_tot)]) 

%% SAVING 
saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day  '_'  info.info_notes.task '_Bar_trial_type4'],'jpeg');
mkdir D:\JC_Figures\behavior;
saveas(gcf,['D:\JC_Figures\behavior\' info.info_notes.MouseID '_'  info.info_notes.Day '_'  info.info_notes.task '_Bar_Count_trial_type4'],'tif');


