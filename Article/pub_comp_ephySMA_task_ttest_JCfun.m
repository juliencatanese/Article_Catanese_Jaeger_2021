function [SMA FRepoch ttestEpoch trigtimes spxtimes sdf_alltr] = pub_comp_ephySMA_task_ttest_JCfun(listc, parfig)
% function [SMA FRepoch ttestEpoch trigtimes spxtimes sdf_alltr] = pub_comp_ephySMA_task_ttest_JCfun(listc, parfig)
% written by Julien Catanese 12/05/2018
% last update : JC 03/12/2019

disp('starting: pub_comp_ephySMA_task_ttest_JCfun')

try
    if parfig.saveSMA==1;
        disp('ATTENTION WILL SAVE A NEW SMA')
    end
catch
    parfig.saveSMA=0;
end

trigtimes = []; Ntrials = []; 
SMA = []; SSA = [] ; SSemA = [];  M_BLE = []; S_BLE = []; Sem_BLE= [];
H_puf_BLE_all = [];   P_puf_BLE_all = [];  H_del_BLE_all = [];  P_del_BLE_all = [];  H_res_BLE_all = [];   P_res_BLE_all = [];
M_puf = [];    S_puf = [];  Sem_puf= [];  M_del = [];  S_del = [];   Sem_del= [];  M_res = [];  S_res = [];  Sem_res= [];

trial_type = parfig.trial_type;
BLE= parfig.BaselineEpoch;
center_evt = parfig.center_evt;
pre=parfig.pre;
post=parfig.post;

ncell = size(listc,1);
for nc=1:ncell;
    
    if trial_type{1} == 'opt' & listc.Opto_post_sess(nc)==0;
        SMA = [SMA;  NaN(1,post + pre +1)];
        SSA = [SSA; NaN(1,post + pre +1)];
        SSemA = [SSemA; NaN(1,post + pre +1)];
        
        M_BLE = [M_BLE; NaN];       S_BLE = [S_BLE; NaN];      Sem_BLE= [Sem_BLE; NaN];
        M_puf = [M_puf; NaN];       S_puf = [S_puf; NaN];      Sem_puf= [Sem_puf; NaN];
        M_del = [M_del; NaN];       S_del = [S_del; NaN];      Sem_del= [Sem_del; NaN];
        M_res = [M_res; NaN];       S_res = [S_res; NaN];      Sem_res= [Sem_res; NaN];
        
        H_puf_BLE_all = [H_puf_BLE_all;  0];
        H_del_BLE_all = [H_del_BLE_all;  0];
        H_res_BLE_all = [H_res_BLE_all;  0];
        
        P_puf_BLE_all = [P_puf_BLE_all;  NaN];
        P_del_BLE_all = [P_del_BLE_all;  NaN];
        P_res_BLE_all = [P_res_BLE_all;  NaN];
        
    else
        disp([num2str(nc) '/' num2str(ncell)]);
        MouseID= listc.MouseID(nc,:); Day=listc.Day(nc,:);
        ChanID = listc.ChanID(nc,:); 
        SessPath = [parfig.WorkFolder listc.SessPath(nc,end-8:end)];       
        SessID = listc.SessID(nc,:)        
        CLUST = listc.CLUST(nc,:);
        CluFile2load= listc.CluFile2load(nc,:);
        FileLocation = [SessPath '\' SessID{1}];
        load([FileLocation '\' CluFile2load]);
       
        disp([num2str(nc) '/' num2str(ncell) ':' num2str(nc/ncell*100) '%'] );
        if nc>1 && sum(MouseID~=listc{nc-1,2})==0 && sum(Day ~= listc{nc-1,3})==0;
            [ MouseID ' ' Day ' Chan' ChanID ' clust#' num2str(CLUST) ]
        else
            [trigtimes] = Get_trigtimes(center_evt, trial_type, FileLocation); AA=trigtimes;  % trigtimes in sec
        end
        
        ntrial = size(trigtimes,2)
        idx_spk = find(cluster_class(:,1)==CLUST);
        spxtimes = cluster_class(idx_spk,2)/10^3; % convert time to sec
        spxtimes = spxtimes';
        
        [sdf_alltr BLE_trspx puf_trspx del_trspx res_trspx] = Get_SDF_alltrials_JCfun(spxtimes, trigtimes, parfig, ntrial);
        
        
        sdf_mean = mean(sdf_alltr);        sdf_std = std(sdf_alltr);           sdf_sem =  sdf_std/sqrt(ntrial);
        SMA = [SMA; sdf_mean];             SSA = [SSA; sdf_std];               SSemA = [SSemA; sdf_sem];
        
        BLE_fr = BLE_trspx./((BLE(end)-BLE(1))/1000);
        puf_fr = puf_trspx./(750/1000);
        del_fr = del_trspx./(750/1000);
        res_fr = res_trspx./(750/1000);
        
        M_BLE = [M_BLE; mean(BLE_fr)];       S_BLE = [S_BLE; std(BLE_fr)];      Sem_BLE= [Sem_BLE; std(BLE_fr)/sqrt(ntrial)];
        M_puf = [M_puf; mean(puf_fr)];       S_puf = [S_puf; std(puf_fr)];      Sem_puf= [Sem_puf; std(puf_fr)/sqrt(ntrial)];
        M_del = [M_del; mean(del_fr)];       S_del = [S_del; std(del_fr)];      Sem_del= [Sem_del; std(del_fr)/sqrt(ntrial)];
        M_res = [M_res; mean(res_fr)];       S_res = [S_res; std(res_fr)];      Sem_res= [Sem_res; std(res_fr)/sqrt(ntrial)];
        
        % TEST STATISTIC (non param)
        
        [P_puf_BLE, H_puf_BLE] = ranksum(puf_fr, BLE_fr);
        [P_del_BLE, H_del_BLE] = ranksum(del_fr, BLE_fr);
        [P_res_BLE, H_res_BLE] = ranksum(res_fr, BLE_fr);
        
        H_puf_BLE_all = [H_puf_BLE_all;  H_puf_BLE];
        H_del_BLE_all = [H_del_BLE_all;  H_del_BLE];
        H_res_BLE_all = [H_res_BLE_all;  H_res_BLE];
        
        P_puf_BLE_all = [P_puf_BLE_all;  P_puf_BLE];
        P_del_BLE_all = [P_del_BLE_all;  P_del_BLE];
        P_res_BLE_all = [P_res_BLE_all;  P_res_BLE];
        
        Ntrials = [Ntrials ; ntrial]; 
        
        if parfig.plot==1 & ntrial>5
            pause(0.5);
            cell2plot=listc(nc,:)
%             Plot_Raster_SDFshaded_JCfun1(SMA(nc,:), SSemA(nc,:), spxtimes, trigtimes, cell2plot, parfig);
            pub_fig_Raster_SDFshaded_JCfun(SMA(nc,:), SSemA(nc,:), spxtimes, trigtimes, cell2plot, parfig);
            % saveas(gcf, ['D:\JC_Figures\SDF_RASTER_cor\SDF_RASTER_cor_listcell_ncell#' num2str(ncell) ],'jpg')
        end
    end
    %     [Sig_idx, BarSig] = Select_sigCell_1std_JCfun_v2(sdf_mean, sdf_std, Sig_idx, BarSig, pre, post, ncell) ;
end

ttestEpoch.H_puf_BLE= H_puf_BLE_all;
ttestEpoch.H_del_BLE= H_del_BLE_all;
ttestEpoch.H_res_BLE= H_res_BLE_all;

ttestEpoch.P_puf_BLE= P_puf_BLE_all;
ttestEpoch.P_del_BLE= P_del_BLE_all;
ttestEpoch.P_res_BLE= P_res_BLE_all;

FRepoch.BLE.mean=M_BLE;
FRepoch.BLE.std=S_BLE;
FRepoch.BLE.Sem=Sem_BLE;

FRepoch.res.mean=M_res;
FRepoch.res.std=S_res;
FRepoch.res.Sem=Sem_res;

FRepoch.puf.mean=M_puf;
FRepoch.puf.std=S_puf;
FRepoch.puf.Sem=Sem_puf;

FRepoch.del.mean=M_del;
FRepoch.del.std=S_del;
FRepoch.del.Sem=Sem_del;


%%
if parfig.saveSMA==1
    parfigUsed=parfig; 
    save('SMA_BLE_lastsaved.mat', 'SMA', 'SSA', 'SSemA', 'FRepoch','ttestEpoch', 'parfigUsed')
    save(['SMA_' trial_type{1} '_' center_evt num2str(size(SMA,1)) '.mat'], 'SMA', 'SSA', 'SSemA', 'FRepoch','ttestEpoch', 'Ntrials', 'parfigUsed')
    disp('SMA SAVED')
end

