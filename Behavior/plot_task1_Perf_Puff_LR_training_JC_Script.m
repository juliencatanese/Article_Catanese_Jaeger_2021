% plot_task1_Perf_Puff_LR_training_JC_Script
% Require the use of 'dat2mat_DIN_JC_Script.m' to generate .mat files
% By Julien Catanese in JaegerLab 
% last updated 3/25/2018


%%
% dat2mat_DIN_JC_Script
%%
% clear all
% close all

load('info.mat');
load('evt.mat');
load('time.mat');

%% define trial idx
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0)+1;

idx_trial_start = trig_end(1:end-1);
idx_trial_end = trig_st(2:end);

idx_trial_resp_start = find(diff(evt_delay)<0); % the end of the delay 
idx_trial_resp_end = find(diff(evt_delay)<0) + 4*20000;  % + 4sec  

Ntrial = max(size(idx_trial_start));
TrialperBlock=20;

%%
Nblock = ceil(Ntrial/TrialperBlock);

perc_Lfail = [];
perc_Rfail = [];
perc_Lfail_no_lick = [];
perc_Rfail_no_lick = [];
perc_Lfail_wrong_side = [];
perc_Rfail_wrong_side = [];

for Nb = 1:Nblock;
    if Ntrial < TrialperBlock*Nb;
        endblock= Ntrial;
    else
        endblock = TrialperBlock*Nb;
    end
    
    N_Lpuff = 0;
    N_Lfail =0;
    N_Lfail_no_lick =0;
    N_Lfail_wrong_side =0;
    
    N_Rpuff = 0;
    N_Rfail =0;
    N_Rfail_no_lick =0;
    N_Rfail_wrong_side =0;
    
    
    for tr=1+(50*(Nb-1)):endblock
        % for LEft Trials (only trial with a Puff)
        if sum(evt_puff_L(idx_trial_start(tr):idx_trial_end(tr)))>0;
            N_Lpuff = N_Lpuff+1;
            if sum(evt_rwd_L(idx_trial_start(tr):idx_trial_end(tr)))==0;
                N_Lfail = N_Lfail +1;
                if sum(evt_lick_L(idx_trial_resp_start(tr):idx_trial_resp_end(tr)))==0 & sum(evt_lick_R(idx_trial_resp_start(tr):idx_trial_resp_end(tr)))==0;
                    N_Lfail_no_lick = N_Lfail_no_lick+1;
                else
                    N_Lfail_wrong_side = N_Lfail_wrong_side +1;
                end
            end
        end
        % for Right Trials (only trial with a Puff)
        if sum(evt_puff_R(idx_trial_start(tr):idx_trial_end(tr)))>0;
            N_Rpuff = N_Rpuff+1;
            if sum(evt_rwd_R(idx_trial_start(tr):idx_trial_end(tr)))==0;
                N_Rfail = N_Rfail +1;
                if sum(evt_lick_L(idx_trial_resp_start(tr):idx_trial_resp_end(tr)))==0 & sum(evt_lick_R(idx_trial_resp_start(tr):idx_trial_resp_end(tr)))==0;
                    N_Rfail_no_lick = N_Rfail_no_lick +1;
                else
                    N_Rfail_wrong_side = N_Rfail_wrong_side +1;
                end
            end
        end
        
        
        
        
    end
    perc_Lfail = [perc_Lfail; (N_Lfail/(N_Lpuff+0.1))*100];
    perc_Rfail = [perc_Rfail; (N_Rfail/(N_Rpuff+0.1))*100];
    
    perc_Lfail_no_lick = [perc_Lfail_no_lick; (N_Lfail_no_lick/(N_Lpuff+0.1))*100];
    perc_Rfail_no_lick = [perc_Rfail_no_lick; (N_Rfail_no_lick/(N_Rpuff+0.1))*100];
    
    perc_Lfail_wrong_side = [perc_Lfail_wrong_side; (N_Lfail_wrong_side/(N_Lpuff+0.1))*100];
    perc_Rfail_wrong_side = [perc_Rfail_wrong_side; (N_Rfail_wrong_side/(N_Rpuff+0.1))*100];
    
end

%%
figure,

plot( perc_Lfail , 'r'), hold on;
plot(perc_Rfail, 'b'), hold on;

plot(perc_Lfail_no_lick, 'm--'), hold on;
plot(perc_Rfail_no_lick, 'c--'), hold on;

plot(perc_Lfail_wrong_side, 'm.-'), hold on;
plot(perc_Rfail_wrong_side, 'c.-'), hold on;

legend('Lfail','Rfail','L no lick', 'R no lick' , 'L wrong side','R wrong side','Location','NorthEast')
xlabel(['#blocks (' num2str(TrialperBlock) 'trials/block)']);
ylabel('errors perc (%)');
ylim([0 100]);

saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day '_PERF'],'tif')




