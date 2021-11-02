% pub_table4_optoSUA_post_stat_JCscript
% create and save a table ('Topto') that contain new cells fields
% Fields: Hopto_OFFon,  Popto_OFFon,  Hopto_ONoff, Popto_ONoff, OPTO_POST
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Topto')
% written by Julien Catanese 2/5/2019
% last updated JC 3/12/2019.
% last updated JC 4/12/2019. Add Tcombo, reduce Topto. 

%%
clearvars -except mypath parfig

load('listcell.mat');
T1=[]; T1=listcell;
ncell=size(T1,1)

%% loop over cells
H_OFFon_all = []; P_OFFon_all = []; P_ONoff_all = []; H_ONoff_all = [];
MeanNspk_ON_all = []; MeanNspk_OFF_all = []; Hopto_all = []; Popto_all = [];
OPTO_POST = []; OPTO_POST_Nstim = [];

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
        
        [Popto,Hopto]=ranksum(Nspk_trON, Nspk_trOFF); % Wilcoxon or U-Mann-Wihtney
        
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
H_opto_post =  Hopto_all'; H_opto_post(isnan(H_opto_post))=0;
P_opto_post =  Popto_all';
MeanNspk_optoON_post = MeanNspk_ON_all';
MeanNspk_optoOFF_post = MeanNspk_OFF_all';

Opto_inib = H_opto_post & (MeanNspk_optoON_post<MeanNspk_optoOFF_post);
Opto_exct = H_opto_post & (MeanNspk_optoON_post>MeanNspk_optoOFF_post);
Opto_post_sess = OPTO_POST;
Opto_post_Nstim = OPTO_POST_Nstim;

Topto=addvars(Tcoord, Opto_inib, Opto_exct, Opto_post_sess, Opto_post_Nstim,...
    H_opto_post,  P_opto_post, MeanNspk_optoON_post, MeanNspk_optoOFF_post) ;
%%
save([mypath '\listcell.mat'],'Topto', '-append'); Topto(1,:)
disp('Topto saved');

%%

Tcombo=[]; 
Tcombo=listcell(:,1:5); AreaID=Tcoord.AreaID; VMVL=Tcoord.VMVL; 
ncell =  Tcoord.ncell; nSess =  Tcoord.nSess; nMouse = Tcoord.nMouse;
Tcombo = addvars(Tcombo, nMouse, nSess ,VMVL, AreaID);
Opto_inib = Topto.Opto_inib; Opto_exct = Topto.Opto_exct; Topto.Opto_post_sess; 
Tcombo= addvars(Tcombo, Opto_inib, Opto_exct, Opto_post_sess);
save([mypath '\Tcombo.mat'],'Tcombo'); Tcombo(1,:)
disp('Tcombo SAVED')


%% count 
Topto.H_opto_post(isnan(Topto.H_opto_post))=0;
Ncell_opto =  sum(Topto.H_opto_post)
Ncell_opto_inib = sum(Topto.H_opto_post & (Topto.MeanNspk_optoON_post<Topto.MeanNspk_optoOFF_post))
Ncell_opto_exct = sum(Topto.H_opto_post & (Topto.MeanNspk_optoON_post>Topto.MeanNspk_optoOFF_post))

disp('done Topto')
close all;

