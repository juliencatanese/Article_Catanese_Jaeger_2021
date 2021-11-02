% Get_TrialsBinVec_TASKrestrict_JCscript
% Generate TrialsBinVec structure with 8 field containing booolean vector for each type of trials,
% eg: TrialBinVec.optoY1N0  
% each row correspond to 1 trial in order of appearence during the task. 
% Written by Julien Catanese 
% last updated 2/22/2019 by JC   

load('info.mat');
load('evt_Restrict_TaskEpoch.mat');
load('time_Restrict_TaskEpoch.mat');
load('Epochs_pre_post_task_st_end.mat');

evt_lick = evt_lick_L + evt_lick_R;
evt_rwd = evt_rwd_L+ evt_rwd_R;

trig_end = find(diff(evt_trial)<0); % RESTRICT TO TASK EPOCH 
trig_st = find(diff(evt_trial)>0); % RESTRICT TO TASK EPOCH 

delay_end = find(diff(evt_delay)<0);
delay_start = find(diff(evt_delay)>0);

sr = info.info_freq_parameters.board_dig_in_sample_rate;

idx_trial_start = trig_end(1:end-1); % from the end of the pulse to the start pulse of next trial
idx_trial_end = trig_st(2:end);
idx_delay_start = delay_start(1:end-1);
idx_delay_end = delay_end(1:end-1);

%% TEMPORARY FIX to be corrected, pblem with the trg_trial correction  ----------------------------------------------
% idx_trial_start= delay_start(1:end-1);
% idx_trial_end= delay_start(2:end)-0.01*sr;% sec
%%--------------------------------------------------------------

idx_resp_start = delay_end(1:end-1);
idx_resp_end = idx_resp_start + 1.5*sr;% + 1.5 sec

Ntrial = max(size(idx_trial_start))
if strcmp(FolderID(nf).name(1:12), 'vgat14_w14d8');
    disp('Ntrial corrected for vgat14_w14d8: mouse stopped task at 170')
    Ntrial=150
end

% Ntrial_opto = max(size(find(diff(evt_opto)<0)))-1
trNb =[];
PuffL1R0= [];
OptoY1N0=[];
ErrDelayY1N0=[];
ErrDelayL1R0=[];
ErrRespY1N0=[];
NorespY1N0=[];
FirstRespL1R0=[];

for tr = 1:Ntrial
    pause(0.01)
    trNb =  [trNb tr];
    PuffL1R0= [PuffL1R0 sum(evt_puff_L((idx_trial_start(tr):idx_trial_end(tr))))~=0 ];
    OptoY1N0=[OptoY1N0 sum(evt_opto((idx_trial_start(tr):idx_trial_end(tr))))~=0];
    NorespY1N0=[NorespY1N0 sum(evt_lick((idx_delay_start(tr):idx_resp_end(tr))))==0];
    %     Rwd1YN0 = [NorespY1N0  sum(evt_rwd(idx_resp_start(tr):idx_resp_end(tr)))~=0];
    
    %Error during delay  (impuls lick)
    impulsR = sum(evt_lick_R(idx_delay_start(tr):idx_delay_end(tr)))~=0;
    impulsL = sum(evt_lick_L(idx_delay_start(tr):idx_delay_end(tr)))~=0;
    
    % Find the First lick Resp L1 or R0  from the start of the delay to the end of the response period.
    FirstlickRespL1R0 = min(find(evt_lick_L(idx_delay_start(tr):idx_resp_end(tr))))<min(find(evt_lick_R(idx_delay_start(tr):idx_resp_end(tr))));
    
    if impulsR && impulsL
        FirstlickRespL1R0 = min(find(evt_lick_L(idx_delay_start(tr):idx_resp_end(tr))))<min(find(evt_lick_R(idx_delay_start(tr):idx_resp_end(tr))));
        FirstRespL1R0=[FirstRespL1R0 FirstlickRespL1R0];
        ErrDelayY1N0=[ErrDelayY1N0 1];
        if PuffL1R0(tr)==1 %left trials
            ErrRespY1N0=[ErrRespY1N0 FirstlickRespL1R0==0];
        else  %Right trials
            ErrRespY1N0=[ErrRespY1N0 FirstlickRespL1R0==1];
        end
        
    elseif impulsR && ~impulsL
        FirstlickRespL1R0 = 0;
        FirstRespL1R0=[FirstRespL1R0 FirstlickRespL1R0];
        ErrDelayL1R0=[ErrDelayL1R0 0];
        ErrDelayY1N0=[ErrDelayY1N0 1];
        if PuffL1R0(tr)==1 %left trials
            ErrRespY1N0=[ErrRespY1N0 1];
        else  %Right trials
            ErrRespY1N0=[ErrRespY1N0 0];
        end
        
    elseif ~impulsR && impulsL
        FirstlickRespL1R0 = 1;
        FirstRespL1R0=[FirstRespL1R0 FirstlickRespL1R0];
        ErrDelayL1R0=[ErrDelayL1R0 1];
        ErrDelayY1N0=[ErrDelayY1N0 1];
        if PuffL1R0(tr)==1 %left trials
            ErrRespY1N0=[ErrRespY1N0 0];
        else  %Right trials
            ErrRespY1N0=[ErrRespY1N0 1];
        end
        
    elseif ~impulsR && ~impulsL   % if no error delay
        ErrDelayL1R0=[ErrDelayL1R0 NaN];
        ErrDelayY1N0=[ErrDelayY1N0 0];
        
        RespR = sum(evt_lick_R(idx_resp_start(tr):idx_resp_end(tr)))~=0; % is there a lick R during resp period
        RespL = sum(evt_lick_L(idx_resp_start(tr):idx_resp_end(tr)))~=0; % is there a lick L during resp period
        
        if RespR && RespL
            FirstlickRespL1R0 = min(find(evt_lick_L(idx_delay_start(tr):idx_resp_end(tr))))<min(find(evt_lick_R(idx_delay_start(tr):idx_resp_end(tr))));
            FirstRespL1R0=[FirstRespL1R0 FirstlickRespL1R0];
            if PuffL1R0(tr)==1 %left trials
                ErrRespY1N0=[ErrRespY1N0 FirstlickRespL1R0==0];
            else  %Right trials
                ErrRespY1N0=[ErrRespY1N0 FirstlickRespL1R0==1];
            end
            
        elseif  RespR && ~RespL
            FirstlickRespL1R0 = 0;
            FirstRespL1R0=[FirstRespL1R0 FirstlickRespL1R0];
            if PuffL1R0(tr)==1 %left trials
                ErrRespY1N0=[ErrRespY1N0 1];
            else  %Right trials
                ErrRespY1N0=[ErrRespY1N0 0];
            end
            
        elseif ~RespR && RespL
            FirstlickRespL1R0 = 1;
            FirstRespL1R0=[FirstRespL1R0 FirstlickRespL1R0];
            if PuffL1R0(tr)==1 %left trials
                ErrRespY1N0=[ErrRespY1N0 0];
            else  %Right trials
                ErrRespY1N0=[ErrRespY1N0 1];
            end
            
        elseif ~RespR && ~RespL
            FirstlickRespL1R0 = NaN;
            FirstRespL1R0=[FirstRespL1R0 FirstlickRespL1R0];
            ErrRespY1N0=[ErrRespY1N0 NaN];                   
        end
        
    end
    
end
%% make a Table Trial structure and save it:
table.puffL1R0 = PuffL1R0;
table.firstrespL1R0 = FirstRespL1R0;
table.optoY1N0 = OptoY1N0;
table.errrespY1N0 = ErrRespY1N0;
table.norespY1N0 = NorespY1N0;
table.errdelayL1R0=ErrDelayL1R0;
table.errdelayY1N0=ErrDelayY1N0;
table.trNb=trNb; 
table
save('table_trials', 'table')
disp('saved as table_trials.mat')
% 
% disp(['trNb  trL1/R0  Resp  erSide erDelay NoResp opto'])
% [table.trNb' table.puffL1R0' table.firstrespL1R0'  table.errrespY1N0' table.errdelayY1N0' table.norespY1N0' table.optoY1N0' ]
% 
