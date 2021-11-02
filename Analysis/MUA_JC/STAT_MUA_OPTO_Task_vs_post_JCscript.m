% Find All BINs of Opto ON and compare to the equivalent number of BINs of opto OFF.
% Nbin and Binsize should be the same
% Ttest firing rate ON vs OFF
% Written by JC 11/13/2018

%% Define Session list
clear all
close all
cd('D:\JC_Analysis');
SessList = dir(['**/*taskopto*']);
NSess= max(size(SessList)) % Number of Sessions
recap = 1 

%% loop trhough all Sessions named "taskopto"
for ns=recap:NSess;
    SessID= SessList(ns).name
    SessPath=SessList(ns).folder;
    cd([SessPath '\' SessID]);
    
    load('info.mat');
    load('time.mat');
    load('evt.mat','evt_opto', 'evt_trial', 'evt_lick_L', 'evt_lick_R');
    evt_lick = evt_lick_L + evt_lick_R;
    clear evt_lick_R evt_lick_L;
    
    %% Manually Define and plot task epoch
        
    try
        Manual_Epochs_St_End_JCspecScript;
    catch
        disp(['stopped at session Nb ' num2str(ns) ])
        recap = ns
        open Manual_Epochs_St_End_JCScript;
        Manual_Epochs_St_End_JCscript;
    end
    
    
% %     ChID_list = dir([path '\S*Ch*_raw.mat']);
% %     for nc= 1:max(size(ChID_list))
% %         ChanID=ChID_list(nc)
% %         
% %         Plot_MUA_OPTO_JC_Script3
% %         
% %     end
end



