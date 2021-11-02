% Find All BINs of Opto ON and compare to the equivalent number of BINs of opto OFF.
% Nbin and Binsize should be the same
% Ttest firing rate ON vs OFF
% Written by JC 11/13/2018

%% Define Session list

clear all
close all

cd('D:\JC_Analysis');
SessList = dir(['*JCVGAT12/*VM*post*']);
NSess= max(size(SessList)) % Number of Sessions
recap = 1
% post_pre_task=3;
for post_pre_task=1:1 % post=1 pre=2 task=3
    
    %% loop trhough all Sessions named "taskopto"
    for ns=recap:NSess;
        clearvars -except SessList NSess  ns recap post_pre_task
        close all,
        SessID= SessList(ns).name
        SessPath=SessList(ns).folder;
        cd([SessPath '\' SessID]);
        
        %% Manually Define and plot task epoch
        %                      pause(0.001)
        %         StimSize = median(diff(find(diff(find(evo))>1)))
        %                      pause(0.001)
        %         if ~isnan(StimSize)
        ChID_list = dir(['S*Ch*_sub.mat']);
        for nc= 1:max(size(ChID_list))
            ChanID=ChID_list(nc).name(1:5)
            MUA_OPTO_1Chan_Stats_JCfun(ChanID, post_pre_task)
            pause(0.001)
        end
        %         else
        %             display('ATTENTION: No STIM pulses during this Epoch')
        %             display('ATTENTION: No Figure being made Jumoing to next ')
        %         end
    end
end