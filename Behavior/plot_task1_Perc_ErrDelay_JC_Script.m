
% clear all
% close all

load('info.mat')
load('evt.mat')
load('time.mat')
load task_param

param_delaytime = Delay_sec*info.info_freq_parameters.board_dig_in_sample_rate


%% define trial idx
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);

idx_trial_start = trig_end(1:end-1);
idx_trial_end = trig_st(2:end);

Ntrial = max(size(idx_trial_start))

%%
Nblock = ceil(Ntrial/50)
DelayLicks_perc =[]; 

for Nb = 1:Nblock
    if Ntrial < 50*Nb
        endblock= Ntrial
    else
        endblock = 50*Nb
    end
   pause(0.01)
   ErrDelay = 0;

    
    for tr=1+(50*(Nb-1)):endblock
        
   ErrDelay = ErrDelay + double(max(size(find(evt_delay(idx_trial_start(tr):idx_trial_end(tr)))))<param_delaytime) 
       
    end
        DelayLicks_perc = [DelayLicks_perc; ((ErrDelay)/50)*100];
end

%%
figure,
plot(DelayLicks_perc, 'k'), hold on
legend('ErrorDelay','Location','NorthWest')
xlabel('#blocks (x50trials)')
ylabel('perc errors (%)')
ylim([0 100])

saveas(gcf,[info.info_notes.MouseID '_'  info.info_notes.Day '_PERF_Delay'],'tif')




