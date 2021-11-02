% Opto_Script
% Written by Julien Catanese 10/27/2018

%% Define Session list
clear all
close all
cd('D:\JC_Analysis');
% SessList = dir(['**/*' '_' 'taskopto' '_' '*']);
% SessList = dir(['**/*taskopto*']);
SessList = dir(['*/*taskopto*']);

NSess= max(size(SessList)); % Number of Sessions
Nol=[]; Nolo=[]; cCR=[]; cCRo=[]; cCL=[]; cCLo=[];
ALL_SDFALL = []; 
%% loop trhough all Sessions named "taskopto"
for of=1:NSess
    SessID= SessList(of).name;
    SessPath=SessList(of).folder;
    cd([SessPath '\' SessID]);
 
    %% CLASSIFIER 
    % linClass_JC_test1
    
    %% SEQUENCES Matrix of Pop activity 
     Pop_Seq_DistMatrix_1session_GOcent_corr_mlibJCscript;
     ALL_SDFALL = [ALL_SDFALL; sort_sdfall] ;
    
end
%%
ALL_SDFALL_sorted = sortrows(ALL_SDFALL,'ascend');
figure, imagesc(ALL_SDFALL_sorted);
colorbar, 
caxis([0.5 1])
title('ALL SESSIONS')

hold on, plot([750 750], ylim,'w--', 'LineWidth',2.5)
hold on, plot([1500 1500], ylim,'w--','LineWidth',2.5)
hold on, plot([2250 2250], ylim,'w--','LineWidth',2.5)

ylabel ('# neurons (sorted)')
xlabel('time (ms)')

