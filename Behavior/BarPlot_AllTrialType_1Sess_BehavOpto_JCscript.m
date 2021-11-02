% BarPlot_AllTrialType_1Sess_BehavOpto_JCscript
% by Julien Catanese; 9/25/2018
% last update: by JC 10/27/2018

load ('info.mat')
load ('Ntrial_type.mat')

Ncorrect = trial.Nb_correct_L + trial.Nb_correct_R; % only trial without opto stim 
Nopto = trial.Nb_all_opto ;      
Nerror = trial.Nb_errorDelay_PL_CL + trial.Nb_errorDelay_PL_CR + trial.Nb_errorDelay_PR_CL + trial.Nb_errorDelay_PR_CR + trial.Nb_errorResp_PR + trial.Nb_errorResp_PL; 
Nnolick = trial.Nb_NoLick; 

figure;
subplot(5,5,[1 2 6 7]);
c = categorical({'not licked';'Correct PL';'Correct PR';});
bar(c, [    trial.Nb_NoLick , trial.Nb_NoLick_opto              ;...
            trial.Nb_correct_L, trial.Nb_correct_L_opto     ;...
            trial.Nb_correct_R, trial.Nb_correct_R_opto])   ;
ylabel('Count')
if info.info_notes.MouseID=='vgat07'
    title ([info.info_notes.MouseID ' ' info.info_notes.Day ' '  info.info_notes.task  '   (ARCH in right VM)'])
else
    title ([info.info_notes.MouseID ' ' info.info_notes.Day ' '  info.info_notes.task '   (ChR2 in left VM)']) 
end

subplot(5,5,[3 4 5 8 9 10 ])
c = categorical({'ErrorSide PL';'ErrorSide PR';'ErrorDelay PL CL';'ErrorDelay PR CR';'ErrorDelay PL CR';'ErrorDelay PR CL'});
bar(c, [    trial.Nb_errorResp_PL , trial.Nb_errorResp_PL_opto             ;...
            trial.Nb_errorResp_PR, trial.Nb_errorResp_PR_opto              ;...
            trial.Nb_errorDelay_PL_CL , trial.Nb_errorDelay_PL_CL_opto ;...
            trial.Nb_errorDelay_PR_CR, trial.Nb_errorDelay_PR_CR_opto  ;...
            trial.Nb_errorDelay_PL_CR, trial.Nb_errorDelay_PL_CR_opto  ;...
            trial.Nb_errorDelay_PR_CL, trial.Nb_errorDelay_PR_CL_opto]);
%     ylabel('Count')
legend('non-Opto', 'Opto', 'Location', 'best');


subplot(5,5, [16 17 21 22]);
c = categorical({'not licked';'Correct PL';'Correct PR';});
bar(c, [    100*(trial.Nb_NoLick)/(trial.Ntrial-trial.Nb_all_opto),       100*trial.Nb_NoLick_opto/trial.Nb_all_opto        ;...
            100*(trial.Nb_correct_L)/(trial.Ntrial-trial.Nb_all_opto) ,   100*trial.Nb_correct_L_opto/trial.Nb_all_opto     ;...
            100*(trial.Nb_correct_R)/(trial.Ntrial-trial.Nb_all_opto) ,   100*trial.Nb_correct_R_opto/trial.Nb_all_opto  ]) ;
ylabel('Percent (N/Ntot)')
ylim([0 90])
%     legend('%non-Opto', '%Opto', 'Location', 'eastoutside')
title ([ 'Ncorrect=' num2str(Ncorrect) '    Nnotlicked=' num2str(Nnolick) ]) 

subplot(5,5,[18 19 20 23 24 25])
c = categorical({'ErrorSide PL';'ErrorSide PR';'ErrorDelay PL CL';'ErrorDelay PR CR';'ErrorDelay PL CR';'ErrorDelay PR CL'});
bar(c, [    100*(trial.Nb_errorResp_PL )/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorResp_PL_opto/trial.Nb_all_opto             ;...
            100*(trial.Nb_errorResp_PR)/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorResp_PR_opto/trial.Nb_all_opto              ;...
            100*(trial.Nb_errorDelay_PL_CL)/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorDelay_PL_CL_opto/trial.Nb_all_opto ;...
            100*(trial.Nb_errorDelay_PR_CR)/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorDelay_PR_CR_opto/trial.Nb_all_opto  ;...
            100*(trial.Nb_errorDelay_PL_CR)/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorDelay_PL_CR_opto/trial.Nb_all_opto  ;...
           100* (trial.Nb_errorDelay_PR_CL)/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorDelay_PR_CL_opto/trial.Nb_all_opto  ]);
%     ylabel('Percent (N/Ntot)')
legend('%non-Opto', '%Opto', 'Location', 'best');
ylim([0 90])
title ([' Nerror=' num2str(Nerror) '     Nopto=' num2str(Nopto)]) ;

saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day  '_'  info.info_notes.task '_Bar_Count_trial_type'],'jpeg');

mkdir D:\JC_Figures\behavior;
saveas(gcf,['D:\JC_Figures\behavior\' info.info_notes.MouseID '_'  info.info_notes.Day '_'  info.info_notes.task '_Bar_Count_trial_type'],'png');

load ('info.mat')
load ('Ntrial_type.mat')
trial

Ncorrect = trial.Nb_correct_L + trial.Nb_correct_R; % only trial without opto stim 
Nopto = trial.Nb_all_opto ;      
Nerror = trial.Nb_errorDelay_PL_CL + trial.Nb_errorDelay_PL_CR + trial.Nb_errorDelay_PR_CL + trial.Nb_errorDelay_PR_CR + trial.Nb_errorResp_PR + trial.Nb_errorResp_PL; 
Nnolick = trial.Nb_NoLick; 

figure;
subplot(5,5,[1 2 6 7]);
c = categorical({'not licked';'Correct PL';'Correct PR';});
bar(c, [    trial.Nb_NoLick , trial.Nb_NoLick_opto              ;...
            trial.Nb_correct_L, trial.Nb_correct_L_opto     ;...
            trial.Nb_correct_R, trial.Nb_correct_R_opto])   ;
ylabel('Count')
if info.info_notes.MouseID=='vgat07'
    title ([info.info_notes.MouseID ' ' info.info_notes.Day ' '  info.info_notes.task  '   (ARCH in right VM)'])
else
    title ([info.info_notes.MouseID ' ' info.info_notes.Day ' '  info.info_notes.task '   (ChR2 in left VM)']) 
end

subplot(5,5,[3 4 5 8 9 10 ])
c = categorical({'ErrorSide PL';'ErrorSide PR';'ErrorDelay PL CL';'ErrorDelay PR CR';'ErrorDelay PL CR';'ErrorDelay PR CL'});
bar(c, [    trial.Nb_errorResp_PL , trial.Nb_errorResp_PL_opto             ;...
            trial.Nb_errorResp_PR, trial.Nb_errorResp_PR_opto              ;...
            trial.Nb_errorDelay_PL_CL , trial.Nb_errorDelay_PL_CL_opto ;...
            trial.Nb_errorDelay_PR_CR, trial.Nb_errorDelay_PR_CR_opto  ;...
            trial.Nb_errorDelay_PL_CR, trial.Nb_errorDelay_PL_CR_opto  ;...
            trial.Nb_errorDelay_PR_CL, trial.Nb_errorDelay_PR_CL_opto]);
%     ylabel('Count')
legend('non-Opto', 'Opto', 'Location', 'best');


subplot(5,5, [16 17 21 22]);
c = categorical({'not licked';'Correct PL';'Correct PR';});
bar(c, [    100*(trial.Nb_NoLick)/(trial.Ntrial-trial.Nb_all_opto),       100*trial.Nb_NoLick_opto/trial.Nb_all_opto        ;...
            100*(trial.Nb_correct_L)/(trial.Ntrial-trial.Nb_all_opto) ,   100*trial.Nb_correct_L_opto/trial.Nb_all_opto     ;...
            100*(trial.Nb_correct_R)/(trial.Ntrial-trial.Nb_all_opto) ,   100*trial.Nb_correct_R_opto/trial.Nb_all_opto  ]) ;
ylabel('Percent (N/Ntot)')
ylim([0 90])
%     legend('%non-Opto', '%Opto', 'Location', 'eastoutside')
title ([ 'Ncorrect=' num2str(Ncorrect) '    Nnotlicked=' num2str(Nnolick) ]) 

subplot(5,5,[18 19 20 23 24 25])
c = categorical({'ErrorSide PL';'ErrorSide PR';'ErrorDelay PL CL';'ErrorDelay PR CR';'ErrorDelay PL CR';'ErrorDelay PR CL'});
bar(c, [    100*(trial.Nb_errorResp_PL )/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorResp_PL_opto/trial.Nb_all_opto             ;...
            100*(trial.Nb_errorResp_PR)/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorResp_PR_opto/trial.Nb_all_opto              ;...
            100*(trial.Nb_errorDelay_PL_CL)/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorDelay_PL_CL_opto/trial.Nb_all_opto ;...
            100*(trial.Nb_errorDelay_PR_CR)/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorDelay_PR_CR_opto/trial.Nb_all_opto  ;...
            100*(trial.Nb_errorDelay_PL_CR)/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorDelay_PL_CR_opto/trial.Nb_all_opto  ;...
           100* (trial.Nb_errorDelay_PR_CL)/(trial.Ntrial-trial.Nb_all_opto) , 100*trial.Nb_errorDelay_PR_CL_opto/trial.Nb_all_opto  ]);
%     ylabel('Percent (N/Ntot)')
legend('%non-Opto', '%Opto', 'Location', 'best');
ylim([0 90])
title ([' Nerror=' num2str(Nerror) '     Nopto=' num2str(Nopto)]) ;

saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day  '_'  info.info_notes.task '_Bar_Count_trial_typeAll'],'jpeg');

mkdir D:\JC_Figures\behavior;
saveas(gcf,['D:\JC_Figures\behavior\' info.info_notes.MouseID '_'  info.info_notes.Day '_'  info.info_notes.task '_Bar_Count_trial_typeAll'],'png');


% %% SECOND FIGURE 
% figure;
% 
% c = categorical({'not licked';'Correct PL';'Correct PR';});
% y = [    100*(trial.Nb_NoLick)/(trial.Ntrial-trial.Nb_all_opto),       100*trial.Nb_NoLick_opto/trial.Nb_all_opto        ;...
%             100*(trial.Nb_correct_L)/(trial.Ntrial-trial.Nb_all_opto) ,   100*trial.Nb_correct_L_opto/trial.Nb_all_opto     ;...
%             100*(trial.Nb_correct_R)/(trial.Ntrial-trial.Nb_all_opto) ,   100*trial.Nb_correct_R_opto/trial.Nb_all_opto  ];
% b = bar(c, y, 'FaceColor','flat');
% colormap('copper')
% for k = 1:size(y,2)
%     b(k).CData = k;
% end
% 
% b(2).LineWidth = 2;
% b(2).EdgeColor = 'cyan';
% 
% ylabel('Percent (N/Ntot)')
% ylim([0 90])
% %     legend('%non-Opto', '%Opto', 'Location', 'eastoutside')
% title ([ 'Ncorrect=' num2str(Ncorrect) '  Nnotlicked=' num2str(Nnolick) ' Nerror=' num2str(Nerror) '  Nopto=' num2str(Nopto)]) 
% 
% saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day  '_'  info.info_notes.task '_Bar_trial_type2'],'jpeg');
% 
% mkdir D:\JC_Figures\behavior;
% saveas(gcf,['D:\JC_Figures\behavior\' info.info_notes.MouseID '_'  info.info_notes.Day '_'  info.info_notes.task '_Bar_Count_trial_type2'],'png');
