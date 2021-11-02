function [trigtimes2]=Get_trigtimes(psth_trig_evt, psth_trial_type, FileLocation)
% function [trigtime]=Get_trigtimes(psth_trig_evt, psth_trial_type, FileLocation)
% for each trial of a given type (e.g correct Left trials)
% return the time center on a task event (e.g. delay)
% INPUT :
%       psth_trig_evt : 'Delay' or 'GoCue' or 'APuff' or 'Licks'
%       psth_trial_type : 'cor', 'imp', 'cCR', 'cCL', 'eNo' ...
%       FileLOcation :  'D\vgat12_w11d5_z4300_12cells_100tr_taskopto_CCAE3Sh_VM_180609_231954_vid_optopost'
% OUTPUT:
%       trigtimes = vector(1,ntrial) of times (sec) needed as input to compute psth
% dependency: Count_tr_types_JC_Script.m
% written by Julien Catanese 11/15/2018
% last updated 12/05/2018

try
    psth_trig_evt ;
    psth_trial_type = psth_trial_type{1};
catch
    display ('Wrong INPUT ARG');
    stop;
end

load([FileLocation '\evt.mat']);
load([FileLocation '\Ntrial_type.mat']);
load([FileLocation '\time.mat']);
load([FileLocation '\Epochs_pre_post_task_st_end.mat'])
%% GET trigtimes : times of the psth_trig_evt is sec   (1 vector of times)
% define trial idx
% evt_trial2 = evt_trial(idx_task_st_end(1):idx_task_st_end(end));

trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);
idx_trial_start = trig_end(1:end-1);
idx_trial_end = trig_st(2:end) ;
Ntrials= trial.Ntrial;

% define common events
evt_valve = evt_rwd_L + evt_rwd_R;
evt_lick = evt_lick_L + evt_lick_R;
evt_puff =evt_puff_L + evt_puff_R;
evt_prelick = evt_delay + evt_puff;

%% define trigtimes for each trials
trigtimes=[];
for tr=1:Ntrials
    time_tr = time(idx_trial_start(tr):idx_trial_end(tr));
    if psth_trig_evt=='APuff'
        puff_tr=  evt_puff(idx_trial_start(tr):idx_trial_end(tr));
        idx_puff_st = find(diff(puff_tr)>0); % start of the puff
        time_puff_st = time_tr(idx_puff_st);
        trigtimes = [trigtimes time_puff_st]; % in sec
        
    elseif psth_trig_evt=='Delay'
        if FileLocation(25:30)=='vgat11' | FileLocation(25:30)=='vgat12'
            puff_tr=  evt_puff(idx_trial_start(tr):idx_trial_end(tr));
            idx_puff_st = find(diff(puff_tr)>0); % start of the puff
            time_delay_st = time_tr(idx_puff_st)+0.75;
            if isempty(time_delay_st); time_delay_st=NaN; end
            trigtimes = [trigtimes time_delay_st]; % in sec
        else
            delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
            idx_delay_st = find(diff(delay_tr)>0);  % start of the delay
            time_delay_st = time_tr(idx_delay_st);
            trigtimes = [trigtimes time_delay_st]; % in sec
        end
        
    elseif psth_trig_evt=='GoCue'
        
        if FileLocation(25:30)=='vgat11' | FileLocation(25:30)=='vgat12'
            puff_tr=  evt_puff(idx_trial_start(tr):idx_trial_end(tr));
            idx_puff_st = find(diff(puff_tr)>0); % start of the puff
            time_GO_st = time_tr(idx_puff_st)+1.5;
            trigtimes = [trigtimes time_GO_st]; % in sec
        else
            delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
            idx_GO_st = find(diff(delay_tr)>0);  % start of the delay
            time_GO_st = time_tr(idx_GO_st)+0.75; % 750ms after start of the delay
            trigtimes = [trigtimes time_GO_st]; % in sec
        end
        
    elseif psth_trig_evt=='Licks'
        lick_tr=  evt_lick(idx_trial_start(tr):idx_trial_end(tr));
        idx_lick_st = min(find(diff(lick_tr)>0));% start of the first lick
        if isempty(idx_lick_st) % replace by Lick evt by GOcue evt for centering during NOlick trial
            delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
            idx_GO_st = find(diff(delay_tr)<0);
            idx_lick_st=idx_GO_st;
        end
        time_lick_st = time_tr(idx_lick_st);
        trigtimes = [trigtimes time_lick_st]; % in sec
        
    else
        disp(' WRONG INPUT ARG pick one: psth_trig_evt = GoCue, APuff, Delay, Licks or Valve' )
    end
end

%%
trigtimes_all = trigtimes';
trigtimes_cCL = trigtimes(trial.idx_correct_L);
trigtimes_cCR = trigtimes(trial.idx_correct_R);
trigtimes_iCL = trigtimes([trial.idx_errorDelay_PL_CL trial.idx_errorDelay_PR_CL]);
trigtimes_iCR = trigtimes([trial.idx_errorDelay_PL_CR trial.idx_errorDelay_PR_CR]);
trigtimes_oCR = trigtimes([trial.idx_correct_R_opto]);
trigtimes_oCL = trigtimes([trial.idx_correct_L_opto]);
trigtimes_ocC = trigtimes([trial.idx_correct_L_opto trial.idx_correct_R_opto]);
trigtimes_oNO = trigtimes([trial.idx_NoLick_opto]);
trigtimes_eNO = trigtimes([trial.idx_NoLick]);
trigtimes_ePL = trigtimes([trial.idx_errorResp_PL]);
trigtimes_ePR = trigtimes([trial.idx_errorResp_PR]);
trigtimes_opI = trigtimes([trial.idx_errorDelay_PL_CL_opto trial.idx_errorDelay_PR_CR_opto trial.idx_errorDelay_PR_CL_opto trial.idx_errorDelay_PL_CR_opto])


trigtimes_cor = sort([trigtimes_cCR trigtimes_cCL]);
trigtimes_imp = sort([trigtimes_iCL trigtimes_iCR]);
trigtimes_oOC = sort([trigtimes_oNO trigtimes_ocC]);

trigtimes_opt = trigtimes([trial.idx_all_opto]);

if psth_trial_type == 'all' % All Trials
    trigtimes2 = sort(trigtimes_all) ;
elseif psth_trial_type == 'opt' % All Opto
    trigtimes2 = sort(trigtimes_opt) ;
elseif psth_trial_type == 'cor' % All Correct 
    trigtimes2 = sort(trigtimes_cor) ;
elseif psth_trial_type == 'imp' % All impulse Error
    trigtimes2 = sort(trigtimes_imp) ;
elseif psth_trial_type == 'omi' % All impulse Error
    trigtimes2 = sort(trigtimes_eNO) ;
elseif psth_trial_type == 'eNO' % All Omission Error
    trigtimes2 = sort(trigtimes_eNO) ;
elseif psth_trial_type == 'cCL' % correct Choice Left
    trigtimes2 = sort(trigtimes_cCL) ;
elseif psth_trial_type == 'cCR'
    trigtimes2 = sort(trigtimes_cCR) ;
elseif psth_trial_type == 'iCL' % impulse Choice Left
    trigtimes2 = sort(trigtimes_iCL) ;
elseif psth_trial_type == 'iCR'
    trigtimes2 = sort(trigtimes_iCR) ;
elseif psth_trial_type == 'ocC' % opto correct Choices (left+right)
    trigtimes2 = sort(trigtimes_ocC) ;
elseif psth_trial_type == 'oNO' % opto error Omission 
    trigtimes2 = sort(trigtimes_oNO) ;
elseif psth_trial_type == 'oOC' % Opto [Omission + Correct]
    trigtimes2 = sort(trigtimes_oOC) ;
elseif psth_trial_type == 'oCR' % opto correct Choice right
    trigtimes2 = sort(trigtimes_oCR) ;
elseif psth_trial_type == 'oCL' % opto correct Choice right
    trigtimes2 = sort(trigtimes_oCL) ;
elseif psth_trial_type == 'ePL' 
    trigtimes2 = sort(trigtimes_ePL) ;
elseif psth_trial_type == 'ePR' % error side Puff Right
    trigtimes2 = sort(trigtimes_ePR) ;
elseif psth_trial_type == 'opC'
    trigtimes2 = sort(trigtimes_ocC)
elseif psth_trial_type == 'opO'
    trigtimes2 = sort(trigtimes_oNO)
elseif psth_trial_type == 'opI'    
    trigtimes2 = sort(trigtimes_opI)
    
end




