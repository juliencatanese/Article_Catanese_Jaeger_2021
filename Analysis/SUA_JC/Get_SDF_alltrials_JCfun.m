function [sdf_alltr BLE_trspx puf_trspx del_trspx res_trspx] = Get_SDF_alltrials_JCfun(spxtimes, trigtimes, parfig, ntrial)
% function [sdf_alltrial BLE_trspx puf_trspx del_trspx res_trspx] = Get_SDF_alltrials_JCfun(spxtimes, trigtimes, parfig, ntrial)
% by JC 1/16/2019

BLE = parfig.BaselineEpoch;
pre=parfig.pre;
post=parfig.post;

if parfig.center_evt == 'Delay'
    GO=750; 
elseif parfig.center_evt == 'GoCue'
    GO=0; 
elseif parfig.center_evt == 'Licks'
    GO=0;
end

sdf=[]; sdf_tt=[];  BLE_trspx= [];   puf_trspx= [];   del_trspx= [];   res_trspx= [];
for tt=1:ntrial  %  loop across trials
    ntrial;
    fr=1; %if '1', normalizes to firing rate (Hz); if '0', does nothing (default)
    binsz=1; %ms  bin size of psth (default: 1 ms)
    [psth_tt trialspx] = mpsth(spxtimes, trigtimes(tt), 'pre', pre, 'post', post, 'fr',fr, 'binsz', binsz, 'chart',0, 'tb',0);
    %for example if DELAY START CENTERED (delay = 0 to 750) trialspx{1}% from -3000 to +1500
    
    BLE_trspx= [BLE_trspx sum(trialspx{1}> (BLE(1)-pre) & trialspx{1}< (BLE(end)-pre))];  % from -2850 to 850
    puf_trspx= [puf_trspx sum(trialspx{1}> GO-1500 & trialspx{1}< GO-750)];
    del_trspx= [del_trspx sum(trialspx{1}> GO-750 & trialspx{1}< GO)];
    res_trspx= [res_trspx sum(trialspx{1}> GO & trialspx{1}< GO+750)];
    
    ftype = 'Gauss' ;% boxcar, Gauss, exp, exGauss
    w = 75; %75ms
    [sdf_tt kernel] = msdf(psth_tt,ftype,w);
    sdf=[sdf sdf_tt];
    
end
sdf_alltr = sdf';
% Baseline_fr_alltr= (trspx_tt)
% Baseline_fr_alltr= (trspx_tt)./((BLE(end)-BLE(1))/1000); % in Hz
