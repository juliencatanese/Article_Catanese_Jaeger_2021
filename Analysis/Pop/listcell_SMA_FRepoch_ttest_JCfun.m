function [SMA SSA SSemA FRepoch ttestEpoch trigtimes spxtimes sdf_alltr] = listcell_SMA_FRepoch_ttest_JCfun(listcell, parfig)
% function [SMA SSA SSemA Sig_idx BarSig] = listcell_SMA_SdfMeanAll_JCfun(pre, post, center_evt, trial_type)
% written by Julien Catanese 12/05/2018
cd('D:\JC_Analysis');
disp('starting: listcell_SMA_SdfMeanAll_JCfun')

SMA = []; SSA = [] ; SSemA = [];  M_BLE = []; S_BLE = []; Sem_BLE= [];
Hinf_puf_BLE_all = [];   Pinf_puf_BLE_all = [];  Hinf_del_BLE_all = [];  Pinf_del_BLE_all = [];  Hinf_res_BLE_all = [];   Pinf_res_BLE_all = [];
Hsup_puf_BLE_all = [];   Psup_puf_BLE_all = [];  Hsup_del_BLE_all = [];  Psup_del_BLE_all = [];  Hsup_res_BLE_all = [];   Psup_res_BLE_all = [];
M_puf = [];    S_puf = [];  Sem_puf= [];  M_del = [];  S_del = [];   Sem_del= [];  M_res = [];  S_res = [];  Sem_res= [];

pre = parfig.pre;
post = parfig.post;
center_evt = parfig.center_evt;
trial_type = parfig.trial_type;
BLE= parfig.BaselineEpoch;

ncell = size(listcell,1);
for nc=1:ncell
    disp([num2str(nc) '/' num2str(ncell)])
    MouseID= listcell{nc,2}; Day=listcell{nc,3}; ChanID = listcell{nc,4}; CLUST = listcell{nc,5}; 
    CluFile2load= listcell{nc,6}; SessPath = listcell{nc,7}; SessID = listcell{nc,8}; 
    FileLocation = [SessPath '\' SessID{1}];
    
    load([FileLocation '\' CluFile2load]);
    
    disp([num2str(nc) '/' num2str(ncell) ':' num2str(nc/ncell*100) '%'] )
    if nc>1 && sum(MouseID~=listcell{nc-1,2})==0 && sum(Day ~= listcell{nc-1,3})==0;
        [ MouseID ' ' Day ' Chan' ChanID ' clust#' num2str(CLUST) ]
    else
        [trigtimes] = Get_trigtimes(center_evt, trial_type, FileLocation); AA=trigtimes;  % trigtimes in sec
    end
    
    ntrial = size(trigtimes,2);
    idx_spk = find(cluster_class(:,1)==CLUST);
    spxtimes = cluster_class(idx_spk,2)/10^3; % convert time to sec
    spxtimes = spxtimes';
    
    [sdf_alltr BLE_trspx puf_trspx del_trspx res_trspx] = Get_SDF_alltrials_JCfun(spxtimes, trigtimes, parfig, ntrial);
    
    
    sdf_mean = mean(sdf_alltr);         sdf_std = std(sdf_alltr);           sdf_sem =  sdf_std/sqrt(ntrial) ;
    SMA = [SMA; sdf_mean];             SSA = [SSA; sdf_std];               SSemA = [SSemA; sdf_sem];
    
    BLE_fr = BLE_trspx./((BLE(end)-BLE(1))/1000);
    puf_fr = puf_trspx./(750/1000);
    del_fr = del_trspx./(750/1000);
    res_fr = res_trspx./(750/1000);

    M_BLE = [M_BLE; mean(BLE_fr)];       S_BLE = [S_BLE; std(BLE_fr)];      Sem_BLE= [Sem_BLE; std(BLE_fr)/sqrt(ntrial)];
    M_puf = [M_puf; mean(puf_fr)];       S_puf = [S_puf; std(puf_fr)];      Sem_puf= [Sem_puf; std(puf_fr)/sqrt(ntrial)];
    M_del = [M_del; mean(del_fr)];       S_del = [S_del; std(del_fr)];      Sem_del= [Sem_del; std(del_fr)/sqrt(ntrial)];
    M_res = [M_res; mean(res_fr)];       S_res = [S_res; std(res_fr)];      Sem_res= [Sem_res; std(res_fr)/sqrt(ntrial)];
    
% TTEST STATISTIC 
    [Hsup_puf_BLE Psup_puf_BLE] = ttest(puf_fr, BLE_fr,'tail', 'right');
    [Hsup_del_BLE Psup_del_BLE] = ttest(del_fr, BLE_fr,'tail', 'right');
    [Hsup_res_BLE Psup_res_BLE] = ttest(res_fr, BLE_fr,'tail', 'right');  
    
    [Hinf_puf_BLE Pinf_puf_BLE] = ttest(puf_fr, BLE_fr,'tail', 'left');
    [Hinf_del_BLE Pinf_del_BLE] = ttest(del_fr, BLE_fr,'tail', 'left');
    [Hinf_res_BLE Pinf_res_BLE] = ttest(res_fr, BLE_fr,'tail', 'left');      
    
    Hsup_puf_BLE_all = [Hsup_puf_BLE_all;  Hsup_puf_BLE];
    Hsup_del_BLE_all = [Hsup_del_BLE_all;  Hsup_del_BLE];
    Hsup_res_BLE_all = [Hsup_res_BLE_all;  Hsup_res_BLE];
    Psup_puf_BLE_all = [Psup_puf_BLE_all;  Psup_puf_BLE];
    Psup_del_BLE_all = [Psup_del_BLE_all;  Psup_del_BLE];
    Psup_res_BLE_all = [Psup_res_BLE_all;  Psup_res_BLE];
    
    Hinf_puf_BLE_all = [Hinf_puf_BLE_all;  Hinf_puf_BLE];
    Hinf_del_BLE_all = [Hinf_del_BLE_all;  Hinf_del_BLE];
    Hinf_res_BLE_all = [Hinf_res_BLE_all;  Hinf_res_BLE];
    Pinf_puf_BLE_all = [Pinf_puf_BLE_all;  Pinf_puf_BLE];
    Pinf_del_BLE_all = [Pinf_del_BLE_all;  Pinf_del_BLE];
    Pinf_res_BLE_all = [Pinf_res_BLE_all;  Pinf_res_BLE];
    
    if parfig.plot==1
        cell2plot=listcell(nc,:);
        Plot_Raster_SDFshaded_JCfun(SMA, SSA, spxtimes, trigtimes, cell2plot, parfig)
    end
    %     [Sig_idx, BarSig] = Select_sigCell_1std_JCfun_v2(sdf_mean, sdf_std, Sig_idx, BarSig, pre, post, ncell) ;
end

ttestEpoch.Hsup_puf_BLE= Hsup_puf_BLE_all;
ttestEpoch.Psup_puf_BLE= Psup_puf_BLE_all;
ttestEpoch.Hsup_del_BLE= Hsup_del_BLE_all;
ttestEpoch.Psup_del_BLE= Psup_del_BLE_all;
ttestEpoch.Hsup_res_BLE= Hsup_res_BLE_all;
ttestEpoch.Psup_res_BLE= Psup_res_BLE_all;

ttestEpoch.Hinf_puf_BLE= Hinf_puf_BLE_all;
ttestEpoch.Pinf_puf_BLE= Pinf_puf_BLE_all;
ttestEpoch.Hinf_del_BLE= Hinf_del_BLE_all;
ttestEpoch.Pinf_del_BLE= Pinf_del_BLE_all;
ttestEpoch.Hinf_res_BLE= Hinf_res_BLE_all;
ttestEpoch.Pinf_res_BLE= Pinf_res_BLE_all;


FRepoch.BLE.mean=M_BLE;
FRepoch.BLE.std=S_BLE;
FRepoch.BLE.Sem=Sem_BLE;

FRepoch.res.mean=M_res;
FRepoch.res.std=S_res;;
FRepoch.res.Sem=Sem_res;

FRepoch.puf.mean=M_puf;
FRepoch.puf.std=S_puf;
FRepoch.puf.Sem=Sem_puf;

FRepoch.del.mean=M_del;
FRepoch.del.std=S_del;
FRepoch.del.Sem=Sem_del;

save('SMA_BLE_lastsaved.mat', 'SMA', 'SSA', 'SSemA', 'FRepoch','ttestEpoch', 'parfig')

end

