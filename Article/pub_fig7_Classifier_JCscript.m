% pub_fig7_Classifier_JCscript
% JC 4/23/2019

figure(1), close(1), 
Mall=[];
Nb_sess= max(Tfig1_VMopto.nSess )
for ii=1:Nb_sess
    listc=listcell(Tfig1_VMopto.nSess==ii & idxcell,:);
    ncell = size(listc,1);
    if ncell >= parfig.minNbcell
        parfig.trial_type = {trial_type_1};
        [nspx_A] = pub_get_nspx_trialtype_JCfun(listc, parfig);
        parfig.trial_type =  {trial_type_2};
        [nspx_B] = pub_get_nspx_trialtype_JCfun(listc, parfig);
        
        if ~sum(sum(isnan(nspx_A))) & ~sum(sum(isnan(nspx_B)))
            
            Lnall=[];
            for Nr = 1:parfig.Nrepeat;
                [Ln] = pub_Classifier_Ln_JCfun(ncell, nspx_A, nspx_B, parfig);  %parfig.ControlShuffle = 1 or 0
                Lnall=[Lnall;Ln];
            end
            
            % PLOT
            MouseID= listc.MouseID(1,:)
            Day=listc.Day(1,:)
            M=mean(Lnall);
            S=std(Lnall);
            
            figure(1), hold on,
            plot( [1:1:ncell] , M, 'LineWidth', 2)
            ylim([0 1])
            xlabel('#Neurons', 'FontSize', 11)
            ylabel('Proba errors', 'FontSize', 11)
            title(['LogClass-' trial_type_1 'v' trial_type_2  '-' parfig.epoch '-' parfig.center_evt '-' celltype], 'FontWeight','bold' ,'FontSize', 12)
            
            
            Mall = [Mall; M(1:parfig.minNbcell)];
            
        else
            disp('not enough trials in this session')
        end
    else
        disp('not enough cells in this session')
    end
end

figure(1), hold on,
plot( [1:1:ncell] , ones(1,ncell)/2, 'k--')
saveas(gcf, ['LogClass-' trial_type_1 'v' trial_type_2  '-' parfig.epoch '-' parfig.center_evt '-' celltype], 'fig')
saveas(gcf, ['LogClass-' trial_type_1 'v' trial_type_2  '-' parfig.epoch '-' parfig.center_evt '-' celltype], 'png')
%% save Mall
save([mypath '\Class10Fold_Mall_'  trial_type_1 'v' trial_type_2 '_' parfig.epoch  '_' parfig.center_evt '_' celltype '_minCell' num2str(parfig.minNbcell) '.mat'], 'Mall')
disp('done all')