% Count_tr_types_JC_Script
% Count and find idx of each trial type for the JC task. 
% 9 trials type (18 types with opto)
% Dependency: Get_TrialsBinVec_TASKrestrict_JCscript.m or Get_Table_trials_behav_JC_SCript.m
% Written by Julien Catanese 
% last updated 9/23/2018 by JC   

load('table_trials.mat');

% disp(['trNb  trL1/R0  Resp  erSide erDelay NoResp opto']);
% [table.trNb' table.puffL1R0' table.firstrespL1R0'  table.errrespY1N0' table.errdelayY1N0' table.norespY1N0' table.optoY1N0' ];
% disp(['trNb  trL1/R0  Resp  erSide erDelay NoResp opto']);
%%
Ntrial = table.trNb(end)
trial.Ntrial = Ntrial;

trial.NPL= sum(table.puffL1R0);
trial.NPR= sum(~table.puffL1R0);

%  case1 = Correct PL (PL=Puff Left P                                                                    
L=1)
idx_correct_L = find(table.puffL1R0==1  & table.firstrespL1R0==1 & table.errdelayY1N0==0  & table.errrespY1N0==0 & table.optoY1N0==0);
trial.Nb_correct_L = size(idx_correct_L,2);
trial.idx_correct_L= idx_correct_L;

% case2 = Correct PR (PR=Puff Right  PR=0)
idx_correct_R = find(table.puffL1R0==0  & table.firstrespL1R0==0 & table.errdelayY1N0==0  & table.errrespY1N0==0 & table.optoY1N0==0);
trial.Nb_correct_R = size(idx_correct_R,2);
trial.idx_correct_R = idx_correct_R;

% case3 = wrong side PL
idx_errorResp_PL = find(table.puffL1R0==1 & table.errrespY1N0==1  & table.errdelayY1N0==0 & table.optoY1N0==0);
trial.Nb_errorResp_PL = size(idx_errorResp_PL,2) ;
trial.idx_errorResp_PL = idx_errorResp_PL;

% case4 = wrong side PR
idx_errorResp_PR = find(table.puffL1R0==0 & table.errrespY1N0==1  & table.errdelayY1N0==0 & table.optoY1N0==0);
trial.Nb_errorResp_PR = size(idx_errorResp_PR,2) ;
trial.idx_errorResp_PR = idx_errorResp_PR;

% case5 = errorDelay PL_CL (CL=Choice Left;  PL=1)
idx_errorDelay_PL_CL = find(table.errdelayY1N0==1 & table.puffL1R0==1  & table.firstrespL1R0==1 & table.optoY1N0==0);
trial.Nb_errorDelay_PL_CL = size(idx_errorDelay_PL_CL,2);
trial.idx_errorDelay_PL_CL = idx_errorDelay_PL_CL; 

% case6 = errorDelay PR_CR (CR=Choice Right;  PR=0)
idx_errorDelay_PR_CR = find(table.errdelayY1N0==1 & table.puffL1R0==0  & table.firstrespL1R0==0 & table.optoY1N0==0);
trial.Nb_errorDelay_PR_CR = size(idx_errorDelay_PR_CR,2);
trial.idx_errorDelay_PR_CR =  idx_errorDelay_PR_CR ; 

% case7 = errorDelay PL_CR (CR=Choice Right;  PL=1)
idx_errorDelay_PL_CR = find(table.errdelayY1N0==1 & table.puffL1R0==1  & table.firstrespL1R0==0 & table.optoY1N0==0);
trial.Nb_errorDelay_PL_CR = size(idx_errorDelay_PL_CR,2);
trial.idx_errorDelay_PL_CR =idx_errorDelay_PL_CR ; 

% case8 = errorDelay PR_CL (CL=Choice Left;  PR=0)
idx_errorDelay_PR_CL = find(table.errdelayY1N0==1 & table.puffL1R0==0  & table.firstrespL1R0==1 & table.optoY1N0==0);
trial.Nb_errorDelay_PR_CL = size(idx_errorDelay_PR_CL,2);
trial.idx_errorDelay_PR_CL =idx_errorDelay_PR_CL;

% case9 = NoLick 
idx_NoLick = find(table.norespY1N0==1 & table.optoY1N0==0);
trial.Nb_NoLick = size(idx_NoLick,2);
trial.idx_NoLick = idx_NoLick; 

%  case10 = omission PL (PL=Puff Left PL=1)
idx_NoLick_L = find(table.puffL1R0==1  & table.norespY1N0==1 & table.optoY1N0==0);
trial.Nb_NoLick_L = size(idx_NoLick_L,2);
trial.idx_NoLick_L= idx_NoLick_L;

% case11 = Omission PR (PR=Puff Right  PR=0)
idx_NoLick_R = find(table.puffL1R0==0  & table.norespY1N0==1 & table.optoY1N0==0);
trial.Nb_NoLick_R = size(idx_NoLick_R,2);
trial.idx_NoLick_R = idx_NoLick_R;

%% OPTO ;
Nopto = sum(table.optoY1N0);
trial.Nb_all_opto= Nopto;
trial.idx_all_opto= find(table.optoY1N0==1);

%  case1 = Correct PL (PL=Puff Left PL=1)
idx_correct_L_opto = find(table.puffL1R0==1  & table.firstrespL1R0==1 & table.optoY1N0==1 & table.errdelayY1N0==0  & table.errrespY1N0==0);
trial.Nb_correct_L_opto = size(idx_correct_L_opto,2);
trial.idx_correct_L_opto = idx_correct_L_opto; 

% case2 = Correct PR (PR=Puff Right  PR=0);
idx_correct_R_opto = find(table.puffL1R0==0  & table.firstrespL1R0==0 & table.errdelayY1N0==0  & table.errrespY1N0==0 & table.optoY1N0==1);
trial.Nb_correct_R_opto = size(idx_correct_R_opto,2);
trial.idx_correct_R_opto = idx_correct_R_opto;

% case3 = wrong side PL
idx_errorResp_PL_opto = find(table.puffL1R0==1 & table.errrespY1N0==1  & table.errdelayY1N0==0 & table.optoY1N0==1);
trial.Nb_errorResp_PL_opto = size(idx_errorResp_PL_opto,2) ;
trial.idx_errorResp_PL_opto = idx_errorResp_PL_opto; 

% case4 = wrong side PR
idx_errorResp_PR_opto = find(table.puffL1R0==0 & table.errrespY1N0==1   & table.errdelayY1N0==0 & table.optoY1N0==1);
trial.Nb_errorResp_PR_opto = size(idx_errorResp_PR_opto,2) ;
trial.idx_errorResp_PR_opto = idx_errorResp_PR_opto; 

% case5 = errorDelay PL_CL (CL=Choice Left;  PL=1)
idx_errorDelay_PL_CL_opto = find(table.errdelayY1N0==1 & table.puffL1R0==1  & table.firstrespL1R0==1 & table.optoY1N0==1);
trial.Nb_errorDelay_PL_CL_opto = size(idx_errorDelay_PL_CL_opto,2);
trial.idx_errorDelay_PL_CL_opto = idx_errorDelay_PL_CL_opto; 

% case6 = errorDelay PR_CR (CR=Choice Right;  PR=0)
idx_errorDelay_PR_CR_opto = find(table.errdelayY1N0==1 & table.puffL1R0==0  & table.firstrespL1R0==0 & table.optoY1N0==1);
trial.Nb_errorDelay_PR_CR_opto = size(idx_errorDelay_PR_CR_opto,2);
trial.idx_errorDelay_PR_CR_opto = idx_errorDelay_PR_CR_opto; 

% case7 = errorDelay PL_CR (CL=Choice Right;  PL=1)
idx_errorDelay_PL_CR_opto = find(table.errdelayY1N0==1 & table.puffL1R0==1  & table.firstrespL1R0==0 & table.optoY1N0==1);
trial.Nb_errorDelay_PL_CR_opto = size(idx_errorDelay_PL_CR_opto,2);
trial.idx_errorDelay_PL_CR_opto = idx_errorDelay_PL_CR_opto; 

% case8 = errorDelay PR_CL (CL=Choice Left;  PL=1)
idx_errorDelay_PR_CL_opto = find(table.errdelayY1N0==1 & table.puffL1R0==0  & table.firstrespL1R0==1 & table.optoY1N0==1);
trial.Nb_errorDelay_PR_CL_opto = size(idx_errorDelay_PR_CL_opto,2);
trial.idx_errorDelay_PR_CL_opto = idx_errorDelay_PR_CL_opto; 

% case9 = NoLick 
idx_NoLick_opto = find(table.norespY1N0==1 & table.optoY1N0==1);
trial.Nb_NoLick_opto = size(idx_NoLick_opto,2);;
trial.idx_NoLick_opto = idx_NoLick_opto; 

%  case10 = omission PL (PL=Puff Left PL=1)
idx_NoLick_L_opto = find(table.puffL1R0==1  & table.norespY1N0==1 & table.optoY1N0==1);
trial.Nb_NoLick_L_opto = size(idx_NoLick_L_opto,2);
trial.idx_NoLick_L_opto= idx_NoLick_L_opto;

% case11 = Omission PR (PR=Puff Right  PR=0)
idx_NoLick_R_opto = find(table.puffL1R0==0  & table.norespY1N0==1 & table.optoY1N0==1);
trial.Nb_NoLick_R_opto = size(idx_NoLick_R_opto,2);
trial.idx_NoLick_R_opto = idx_NoLick_R_opto;

save('Ntrial_type', 'trial')

% % plot
% ErrRespY1N0(isnan(ErrRespY1N0))=0;
% PuffL1R0(find(ErrRespY1N0));
% firstrespL1R0(find(ErrRespY1N0));
%
% Nb_no_resp=sum(NorespY1N0);
% Nb_errorResp=  size(PuffL1R0(find(ErrRespY1N0)),2);
% Nb_errREsp_left = sum(PuffL1R0(find(ErrRespY1N0)));
% Nb_errREsp_right = sum(~PuffL1R0(find(ErrRespY1N0)));
%
% figure,
% c = categorical({'total_ErrResp','Left_errResp','Right_errResp'});
% bar(c, [Nb_errorResp, Nb_errREsp_left, Nb_errREsp_right]);
%
% trial_idx_err_resp = find(ErrRespY1N0)