% Pop_Seq_DistMatrix_ALLsession_JCscript
% Plot population Sequences of peak normalized Activity sorted by time in task
% see Fig2B
% dependency:     Pop_Seq_DistMatrix_1session_mlibJCscript;
% Written by Julien Catanese 10/27/2018

%% Define Session list
cd('D:\JC_Analysis');
SessList = dir(['*/*taskopto*']);

Ntt = max(size(psth_trial_type))
NSess= max(size(SessList)); % Number of Sessions
Nol=[]; Nolo=[]; cCR=[]; cCRo=[]; cCL=[]; cCLo=[];
ALL_SDFALL = [];
%% loop trhough all Sessions named "taskopto"
for ns=1:NSess
    SessID= SessList(ns).name;
    SessPath=SessList(ns).folder;
    cd([SessPath '\' SessID]);
    if Ntt == 1
        Pop_Seq_DistMatrix_1session_mlibJCscript;
    else
        Pop_Seq_Matrice_1session_JCfun;
    end
    ALL_SDFALL = [ALL_SDFALL; sort_sdfall] ;
    
end
%%
ALL_SDFALL_sorted = sortrows(ALL_SDFALL,'ascend');
ALL_SDFALL_sort = ALL_SDFALL_sorted(:,2:end)
figure, imagesc(ALL_SDFALL_sort);
if Ntt == 1
    colormap('default')
    colorbar, caxis([0.5 1])
else
    colormap('jet')
colorbar, caxis([-1.2 1.2])
end

title('ALL SESSIONS')

hold on, plot([750 750], ylim,'w--', 'LineWidth',2.5)
hold on, plot([1500 1500], ylim,'w--','LineWidth',2.5)
hold on, plot([2250 2250], ylim,'w--','LineWidth',2.5)

ylabel ('# neurons (sorted)')
xlabel('time (ms)')
