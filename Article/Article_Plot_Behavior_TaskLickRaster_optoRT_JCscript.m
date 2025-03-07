% plot_Task_LickRASTER_JC_Script
% Require the use of 'dat2mat_DIN_JC_Script.m' to generate .mat files
% By Julien Catanese in JaegerLab 2020

clear all
close all
load('info.mat');
load('evt.mat');
load('time.mat');

%% define trial idx
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);

idx_trial_start = trig_end(1:end-1);
idx_trial_end = trig_st(2:end);

Ntrial= length(find(diff(evt_trial)<0))-1
NtrialR= length(find(diff(evt_puff_R)<0))
NtrialL= length(find(diff(evt_puff_L)<0))
NtrialOpto = length(find(diff(evt_opto)<0))
opto_tr=0;
% find long gap between trials.
% A=diff(idx_trial_start)>mean(diff(idx_trial_start));
% AA=find(A);

%% transform 0 into NaN
evt_delay(~evt_delay)=NaN;
evt_opto(~evt_opto)=NaN;

evt_lick_L(~evt_lick_L)=NaN;
evt_lick_R(~evt_lick_R)=NaN;

evt_puff_L(~evt_puff_L)=NaN;
evt_puff_R(~evt_puff_R)=NaN;

evt_rwd_L(~evt_rwd_L)=NaN;
evt_rwd_R(~evt_rwd_R)=NaN;

%%
close all

Nblock = ceil(Ntrial/50)

for Nb = 1:Nblock
    if Ntrial < 50*Nb
        endblock= Ntrial;
    else
        endblock = 50*Nb;
    end
    
    if Nb<10 
        zz = '0'; 
    else
        zz='';
    end
    
    figure,
    pause(0.01);  
    opto_tr;
    for tr=1+(50*(Nb-1)):endblock
        pause(0.01)
        time_tr = time(idx_trial_start(tr):idx_trial_end(tr))-time(idx_trial_start(tr));
        
        plot(time_tr,tr*evt_delay(idx_trial_start(tr):idx_trial_end(tr)),'color',[0.5 0.5 0.5],'LineWidth',2);
        xlim([0,6]);
        
        hold on,
        plot(time_tr,tr*evt_puff_L(idx_trial_start(tr):idx_trial_end(tr)),'r'); %m
        hold on,
        plot(time_tr,tr*evt_puff_R(idx_trial_start(tr):idx_trial_end(tr)),'b'); %c
        hold on,
        
        plot(time_tr,tr*evt_rwd_L(idx_trial_start(tr):idx_trial_end(tr)),'*k');
        hold on,
        plot(time_tr,tr*evt_rwd_R(idx_trial_start(tr):idx_trial_end(tr)),'*k');
        hold on,
        
        plot(time_tr,tr*evt_lick_L(idx_trial_start(tr):idx_trial_end(tr)),'.r') %color',[1 0.45 0.65])
        hold on,
        plot(time_tr,tr*evt_lick_R(idx_trial_start(tr):idx_trial_end(tr)),'.b')
        hold on,
        
        plot(time_tr,tr*evt_opto(idx_trial_start(tr):idx_trial_end(tr)),'c','LineWidth',2)
        hold on,
        
        if nansum(evt_opto(idx_trial_start(tr):idx_trial_end(tr)))>1
            opto_tr=opto_tr+1; 
        end
       
        xlabel('time (s)');
        ylabel('# trials');
        title( [info.info_notes.MouseID '  '  info.info_notes.Day ]);
%         xlim([0 4.2])        
    end
%     xlim([0,13])
%     saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day '_BEHAVIOR_block#' zz num2str(Nb) ],'tif');
%     saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day '_BEHAVIOR_block#' zz num2str(Nb) ],'png');
%     saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day '_BEHAVIOR_block#' zz num2str(Nb) ],'emf');
end

