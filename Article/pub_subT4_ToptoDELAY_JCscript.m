% pub_table4_optoSUA_post_stat_JCscript
% create and save a table ('Topto') that contain new cells fields
% Fields: Hopto_OFFon,  Popto_OFFon,  Hopto_ONoff, Popto_ONoff, OPTO_POST
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Topto')
% written by Julien Catanese 2/5/2019
% last updated JC 3/12/2019.
% last updated JC 4/12/2019. Add Tcombo, reduce Topto.

%%
close all
clearvars -except mypath parfig
load('Tfig1_VMopto.mat', 'Tfig1_VMopto');
load listcell.mat;
CList=listcell(Tfig1_VMopto.Opto_inib,:)
ncell = max(size(CList))
%% loop over cells
H_OFFon_all = []; P_OFFon_all = []; P_ONoff_all = []; H_ONoff_all = [];
MeanNspk_ON_all = []; MeanNspk_OFF_all = []; Hopto_all = []; Popto_all = [];
OPTO_POST = []; OPTO_POST_Nstim = [];

for nc=1:ncell
    disp([num2str(nc) '/' num2str(ncell) '  table4' ]);
    MouseID= CList.MouseID(nc,:); Day=CList.Day(nc,:);
    ChanID = CList.ChanID(nc,:); CLUST = CList.CLUST(nc,:);
    CluFile2load= CList.CluFile2load(nc,:);
    SessPath = CList.SessPath(nc,:);
    SessID = CList.SessID(nc,:);
    FileLocation = [SessPath '\' SessID{1}];
    load([FileLocation '\' CluFile2load]);
    disp(['loading evt for ' MouseID ' ' Day ' Chan' ChanID ' clust#' num2str(CLUST)]);
    load([FileLocation '\evt.mat'],'evt_opto');
    load([FileLocation '\Epochs_pre_post_task_st_end.mat']);
    load([FileLocation '\time.mat']);
    
    %% Restrict to POST epoch
    %  opto stim time
    evt_opto_post = evt_opto(idx_postStim_st_end(1):idx_postStim_st_end(2));
    time_opto_post = time(idx_postStim_st_end(1):idx_postStim_st_end(2));
    % Spk times
    spk_time_all = cluster_class(find(cluster_class(:,1)==CLUST),2)/10^3; % ATTENTION in converted to sec from ms (/10^3)
    spk_time_post = spk_time_all(spk_time_all>time_opto_post(1) & spk_time_all<time_opto_post(end));
    
    %% Define stimtrial: stimON stimOFF start and end   
    idx_stim_st = find(diff(evt_opto_post)>0); % idx relative to the start of Epoch: POST opto (from 0 to 80sec*sr)
    idx_stim_end = find(diff(evt_opto_post)<0);
    
    time_trON_st=   time_opto_post(idx_stim_st(1:end));  % time absolute in sec
    time_trON_end=  time_opto_post(idx_stim_end(1:end));
    time_trOFF_st=  time_opto_post(idx_stim_end(1:end));
    time_trOFF_end= time_opto_post(idx_stim_end(1:end))+median(time_trON_end-time_trON_st);
    
    Ntrial = max(size(idx_stim_st));
    
    %% Compute minimum delay between stimOn and spk 
    Aon=[]; Aspk=[]; Bidx=[]; Bspk = []; Delos = []; 
    Aon= time_trON_st;
    Aspk=spk_time_post;
    Bidx=zeros(max(size(Aon)),1); 
    for ii=1:max(size(Aon)); 
        if sum(Aspk>Aon(ii))>1
            Bidx(ii)=min(find(Aspk>Aon(ii)));
        else 
           Bidx(ii)=1 
        end
    end
        
    Bspk= Aspk(Bidx);
    Delos = Bspk-Aon; 
    
    MINd(nc)  = min(Delos)
    MEANd(nc) = mean(Delos)
    STDd(nc) = std(Delos)
    MEDd(nc)  = median(Delos)
    
    % distribution visu
    figure, 
    histogram(Delos,'BinWidth',0.01, 'BinLimits', [0.01 0.3]) 
    
    
    %% STATISTIC ttest
    Nspk_trON=[];
    Nspk_trOFF=[];
    for nt=1:Ntrial
        Nspk_trON =   [Nspk_trON  ; sum(spk_time_post > time_trON_st(nt) & spk_time_post < time_trON_end(nt))  ];
        Nspk_trOFF =  [Nspk_trOFF ; sum(spk_time_post > time_trOFF_st(nt) & spk_time_post < time_trOFF_end(nt))];
    end
    
    [Popto,Hopto]=ttest2(Nspk_trON, Nspk_trOFF); % Wilcoxon or U-Mann-Wihtney
    
    MeanNspk_ON = mean(Nspk_trON);
    MeanNspk_OFF = mean(Nspk_trOFF);
    
    Hopto_all = [Hopto_all Hopto];
    Popto_all = [Popto_all Popto];
    MeanNspk_ON_all = [MeanNspk_ON_all MeanNspk_ON];
    MeanNspk_OFF_all = [MeanNspk_OFF_all MeanNspk_OFF];
    
    OPTO_POST = [OPTO_POST; 1];
    OPTO_POST_Nstim = [OPTO_POST_Nstim; Ntrial] ;
    
    disp([ 'Popto=' num2str(Hopto) ' ' MouseID ' ' Day ' ' ChanID ' CLUST' num2str(CLUST) '( nstim=' num2str(Ntrial) ' ;   Hopto=' num2str(Popto) ')'])
    
end


