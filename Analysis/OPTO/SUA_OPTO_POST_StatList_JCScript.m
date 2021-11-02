
close all, clear all,
cd('D:\JC_Analysis');
load('listcell.mat');
T1=listcell;
% T1=T0(~strcmp(cellstr(T0.MouseID),'vgat11') & strcmp(cellstr(T0.Day),'w10d7'),:)
ncell=size(T1,1)
%%
H_OFFon_all = [];
P_OFFon_all = [];
P_ONoff_all = [];
H_ONoff_all = [];
OPTO_POST = [];

for nc=1:ncell
    SessID =T1{nc,9};
    SessPath =T1{nc,8};
    CluFile2load=T1{nc,7};
    MouseID=T1{nc,2};
    Day=T1{nc,3};
    ChanID =T1{nc,5};
    CLUST =T1{nc,6};
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
        
        H_OFFon_all = [H_OFFon_all NaN];
        H_ONoff_all = [H_ONoff_all NaN];
        P_OFFon_all = [P_OFFon_all NaN];
        P_ONoff_all = [P_ONoff_all NaN];
        
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
        %%
        Nspk_trON=[]; 
        Nspk_trOFF=[]; 
        for nt=1:Ntrial
            Nspk_trON =   [Nspk_trON  ; sum(spk_time_post > time_trON_st(nt) & spk_time_post < time_trON_end(nt))  ];
            Nspk_trOFF =  [Nspk_trOFF ; sum(spk_time_post > time_trOFF_st(nt) & spk_time_post < time_trOFF_end(nt))];
        end
        
        [H_OFFon,P_OFFon]=ttest(Nspk_trON, Nspk_trOFF, 'tail', 'left');
        [H_ONoff,P_ONoff]=ttest(Nspk_trON, Nspk_trOFF, 'tail', 'right');
        
        %     [Pr_OFFon,Hr_OFFon]=ranksum(Nspk_trON', Nspk_trOFF', 'tail', 'left');
        %     [Pr_ONoff,Hr_ONoff]=ranksum(Nspk_trON', Nspk_trOFF', 'tail', 'right');
        
        H_OFFon_all = [H_OFFon_all H_OFFon];
        P_OFFon_all = [P_OFFon_all P_OFFon];
        P_ONoff_all = [P_ONoff_all P_ONoff];
        H_ONoff_all = [H_ONoff_all H_ONoff];
        OPTO_POST = [OPTO_POST; 1];
        
        disp([ 'H_OFFon=' num2str(H_OFFon) '  H_ONoff=' num2str(H_ONoff) '  ' MouseID ' ' Day ' ' ChanID ' CLUST' num2str(CLUST) '( nstim=' num2str(Ntrial) ' ;   P_OFFon=' num2str(P_OFFon) ')'])
        
    end
    
end
%%
 Hopto_OFFon =  H_OFFon_all';  
 Popto_OFFon =  P_OFFon_all';  
 Hopto_ONoff =  H_ONoff_all'; 
 Popto_ONoff =  P_ONoff_all'; 

T2=[]; 
T2=addvars(T1, Hopto_OFFon,  Popto_OFFon,  Hopto_ONoff, Popto_ONoff, OPTO_POST) ; 
T2(1:5,:)
%%
close all; 
nansum(T2.Hopto_OFFon)
nansum(T2.Hopto_ONoff)

nansum(T2.Hopto_OFFon((strcmp(cellstr(T2.Day),'w11d5'))))
nansum(T2.Hopto_OFFon((strcmp(cellstr(T2.Day),'w10d4'))))
size(T2.Hopto_OFFon(strcmp(cellstr(T2.MouseID),'vgat17') & strcmp(cellstr(T2.Day),'w10d4')))
nansum(T2.Hopto_OFFon(strcmp(cellstr(T2.MouseID),'vgat17') & strcmp(cellstr(T2.Day),'w10d4')))
%%
sum(T2.OPTO_POST)

sum(T2.OPTO_POST & strcmp(cellstr(T2.MouseID),'vgat12'))
sum(T2.OPTO_POST & strcmp(cellstr(T2.MouseID),'vgat14'))
sum(T2.OPTO_POST & strcmp(cellstr(T2.MouseID),'vgat15'))
sum(T2.OPTO_POST & strcmp(cellstr(T2.MouseID),'vgat17'))

