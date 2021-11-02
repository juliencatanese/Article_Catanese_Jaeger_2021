function [all_nspx ] = pub_get_nspx_trialtype_JCfun(listc, parfig) 
% function [all_nspx ] = pub_get_nspx_trialtype_JCfun(listc, parfig)  
% Julien Catanese 10/31/2018
% last modified 4/23/2019

MouseID= listc.MouseID(1,:)
Day=listc.Day(1,:)
SessPath = listc.SessPath(1,:)
SessID = listc.SessID(1,:);
FileLocation = [SessPath '\' SessID{1}];

%% GET trigtimes : times of the psth_trig_evt is sec   (1 vector of times)
[trigtimes] = Get_trigtimes(parfig.center_evt, parfig.trial_type, FileLocation); AA=trigtimes;  % trigtimes in sec
ntrial = size(trigtimes,2)
if ntrial > 5
%% loop trough each cells of the tabel/session
 all_nspx =[]; ncell = size(listc,1)
for nc=1:ncell;     disp([num2str(nc) '/' num2str(ncell) ' ' listc.ChanID(nc,:) ' ' SessID{1} ]);   
    CLUST = listc.CLUST(nc,:);  ;
    CluFile2load= listc.CluFile2load(nc,:); 
    % Load cluster file containing cluster_class
    load([FileLocation '\' CluFile2load]); 
    idx_spk = find(cluster_class(:,1)==CLUST);
    spxtimes = cluster_class(idx_spk,2)/10^3; spxtimes = spxtimes'; % convert time to sec
    
    nspx = mnspx(spxtimes,trigtimes, parfig.pre, parfig.post); % 1 sec before center event (GO cue) 
    all_nspx = [all_nspx nspx]; 
end

disp('done')
else 
    disp('not enough trials')
    all_nspx = nan; 
end


