% pub_table4_optoSUA_post_stat_JCscript
% create and save a table ('Topto') that contain new cells fields
% Fields: Hopto_OFFon,  Popto_OFFon,  Hopto_ONoff, Popto_ONoff, OPTO_POST
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Topto')
% written by Julien Catanese 2/5/2019
% last updated JC 3/12/2019.

%%
clearvars -except parfig nn listcell  st av tv

try
    if parfig.saveTABLE ==1
        disp('ATTENTION WILL SAVE A NEW TABLE Tephys')
    end
catch
    parfig.saveTABLE=0
end

% cd('D:\JC_Analysis');
% load('listcell.mat');
T1=listcell;
% T1=T0(~strcmp(cellstr(T0.MouseID),'vgat11') & strcmp(cellstr(T0.Day),'w10d7'),:)
ncell=size(T1,1)
%%
H_OFFon_all = [];
P_OFFon_all = [];
P_ONoff_all = [];
H_ONoff_all = [];
MeanNspk_ON_all = [];
MeanNspk_OFF_all = [];
Hopto_all = [];
Popto_all = [];
OPTO_POST = [];
OPTO_POST_Nstim = [];

for nc=1:ncell
    disp([num2str(nc) '/' num2str(ncell) '  table4' ]);
    
    MouseID= T1.MouseID(nc,:); Day=T1.Day(nc,:)
    ChanID = T1.ChanID(nc,:); CLUST = T1.CLUST(nc,:)
    CluFile2load= T1.CluFile2load(nc,:)
    SessPath = T1.SessPath(nc,:)
    SessID = T1.SessID(nc,:)
    FileLocation = [SessPath '\' SessID{1}];
    
    load([FileLocation '\' CluFile2load]);
    
    if nc>1 && sum(MouseID~=T1{nc-1,2})==0 && sum(Day ~=T1{nc-1,3})==0;
        %         disp([ MouseID ' ' Day ' Chan' ChanID ' clust#' num2str(CLUST) ])
    else
        disp(['loading evt for ' MouseID ' ' Day ' Chan' ChanID ' clust#' num2str(CLUST)])
        load([FileLocation '\evt.mat'])
        load([FileLocation '\Epochs_pre_post_task_st_end.mat'])
        load([FileLocation '\time.mat'])
        figure, plot(time(idx_postStim_st_end(1):idx_postStim_st_end(2)), evt_opto(idx_postStim_st_end(1):idx_postStim_st_end(2))), xlabel('sec')
    end
    
    %% Restrict to POST epoch
    evt_opto_post = evt_opto(idx_postStim_st_end(1):idx_postStim_st_end(2));
    time_opto_post = time(idx_postStim_st_end(1):idx_postStim_st_end(2));
    
    if sum(evt_opto_post)==0
        OPTO_POST = [OPTO_POST; 0];
        OPTO_POST_Nstim = [OPTO_POST_Nstim; 0];
        
        Hopto_all = [Hopto_all 0];
        Popto_all = [Popto_all 1];
        MeanNspk_ON_all = [MeanNspk_ON_all 0];
        MeanNspk_OFF_all = [MeanNspk_OFF_all 0];
        
        disp([ 'no_opto_post '  MouseID ' ' Day ])
        
    else
        spk_time_all = cluster_class(find(cluster_class(:,1)==CLUST),2)/10^3; % ATTENTION in converted to sec from ms (/10^3)
        spk_time_post = spk_time_all(spk_time_all>time_opto_post(1) & spk_time_all<time_opto_post(end));
        
        %% Define trials
        idx_stim_st = find(diff(evt_opto_post)>0); % idx relative to the start of Epoch: POST opto (from 0 to 80sec*sr)
        idx_stim_end = find(diff(evt_opto_post)<0);
        
        time_trON_st=   time_opto_post(idx_stim_st(1:end));  % time absolute in sec
        time_trON_end=  time_opto_post(idx_stim_end(1:end));
        time_trOFF_st=  time_opto_post(idx_stim_end(1:end));
        time_trOFF_end= time_opto_post(idx_stim_end(1:end))+median(time_trON_end-time_trON_st);
        
        Ntrial = max(size(idx_stim_st));
        %% STAT
        Nspk_trON=[];
        Nspk_trOFF=[];
        for nt=1:Ntrial
            Nspk_trON =   [Nspk_trON  ; sum(spk_time_post > time_trON_st(nt) & spk_time_post < time_trON_end(nt))  ];
            Nspk_trOFF =  [Nspk_trOFF ; sum(spk_time_post > time_trOFF_st(nt) & spk_time_post < time_trOFF_end(nt))];
        end
        
        
        %         ZspkON = (Nspk_trON-mean(Nspk_trON))/std(Nspk_trON);
        %         ZspkOFF = (Nspk_trOFF-mean(Nspk_trOFF))/std(Nspk_trOFF);
        % %         %         figure, subplot(1,2,1); hist(ZspkON) ; subplot(1,2,2); hist(ZspkOFF);
        %         try
        %             [HKSoff]=kstest(ZspkOFF);
        %             [HKSon]=kstest(ZspkON);
        %         catch
        %             HKSoff=1
        %             HKSon=1
        %         end
        %
        %         if ~HKSoff & ~HKSon
        %             [Hopto,Popto]=ttest(Nspk_trON, Nspk_trOFF);
        %             disp('ttest')
        %         else
        [Popto,Hopto]=ranksum(Nspk_trON, Nspk_trOFF); % Wilcoxon or U-Mann-Wihtney
        %             disp('Wilcoxon')
        %         end
        
        MeanNspk_ON = mean(Nspk_trON);
        MeanNspk_OFF = mean(Nspk_trOFF);
        
        Hopto_all = [Hopto_all Hopto];
        Popto_all = [Popto_all Popto];
        MeanNspk_ON_all = [MeanNspk_ON_all MeanNspk_ON];
        MeanNspk_OFF_all = [MeanNspk_OFF_all MeanNspk_OFF];
        
        OPTO_POST = [OPTO_POST; 1];
        OPTO_POST_Nstim = [OPTO_POST_Nstim; Ntrial] ;
        
        disp([ 'Hopto=' num2str(Hopto) ' ' MouseID ' ' Day ' ' ChanID ' CLUST' num2str(CLUST) '( nstim=' num2str(Ntrial) ' ;   Popto=' num2str(Popto) ')'])
        
    end
    
end
%%
H_opto_post =  Hopto_all';
P_opto_post =  Popto_all';
MeanNspk_optoON_post = MeanNspk_ON_all';
MeanNspk_optoOFF_post = MeanNspk_OFF_all';

T2=[]; T2=removevars(listcell,[6, 7, 8]);
T2=addvars(T2, H_opto_post,  P_opto_post, MeanNspk_optoON_post, MeanNspk_optoOFF_post, OPTO_POST, OPTO_POST_Nstim ) ;
Topto=T2;

if parfig.saveTABLE ==1
    save('D:\JC_Analysis\listcell.mat','Topto', '-append')
    disp('Topto saved')
end

%%
Topto.H_opto_post(isnan(Topto.H_opto_post))=0;
Ncell_opto =  sum(Topto.H_opto_post)
Ncell_opto_inib = sum(Topto.H_opto_post & (Topto.MeanNspk_optoON_post<Topto.MeanNspk_optoOFF_post))
Ncell_opto_exct = sum(Topto.H_opto_post & (Topto.MeanNspk_optoON_post>Topto.MeanNspk_optoOFF_post))

disp('done Topto')
close all;

% %%
% T2(1:5,:)
% close all;
% nansum(T2.Hopto_OFFon)
% nansum(T2.Hopto_ONoff)
%
% nansum(T2.Hopto_OFFon((strcmp(cellstr(T2.Day),'w11d5'))))
% nansum(T2.Hopto_OFFon((strcmp(cellstr(T2.Day),'w10d4'))))
% size(T2.Hopto_OFFon(strcmp(cellstr(T2.MouseID),'vgat17') & strcmp(cellstr(T2.Day),'w10d4')))
% nansum(T2.Hopto_OFFon(strcmp(cellstr(T2.MouseID),'vgat17') & strcmp(cellstr(T2.Day),'w10d4')))
% %%
% sum(T2.OPTO_POST)
%
% sum(T2.OPTO_POST & strcmp(cellstr(T2.MouseID),'vgat12'))
% sum(T2.OPTO_POST & strcmp(cellstr(T2.MouseID),'vgat14'))
% sum(T2.OPTO_POST & strcmp(cellstr(T2.MouseID),'vgat15'))
% sum(T2.OPTO_POST & strcmp(cellstr(T2.MouseID),'vgat17'))

