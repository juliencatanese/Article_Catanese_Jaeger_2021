function [psth trialspx] = SDF_Raster_PSTH_LvRvO_mlibJC_fun(psth_center_evtID, Other_trial_type)
% function [psth trialspx] = SDF_Raster_PSTH_LvRvO_mlibJC_fun(psth_center_evtID, Other_trial_type)
% Plot a SDF and Raster blue ('cCR' = correct Choice Right), red ('cCL' = correct Choice Left)
% and grey (Other: Define by INPUT argument.  
% INPUT1 : evt_center = 
% 'Delay'; 'APuff'; 'GoCue'; 'Licks'; 'Valve'
% INPUT2 : Other_trial_type =    
% 'iCR' = impulse Choice Right (or Left: 'iCL')
% 'oCR' = opto correct Choice Right (or Left: 'oCL')
% 'ocC' = opto correct Choices (L+R)
% 'oNO' = opto Not Licked. 
% 'ePL' = error Puff Left (lick R instead of L)
% This function use the mlib toolbox for mraster, mpsth and msdf.  
% Written by Julien Catanese, 10/03/2018, in JaegerLab. 
% Last Updated: 10/08/2018

% psth_center_evtID = 'Delay'
% Other_trial_type = 

clearvars -except  psth_center_evtID  Other_trial_type  MouseID FolderID
close all,

%% Define PSTH parameter

if nargin==0;
        display('NOT ENOUGH INPUT ARGUMENTS: min=1')
elseif nargin==1; 
    disp(['psth Center = ' psth_center_evtID]);
    disp('NO Other trial type, plotting LvR only');
    psth_trial_type = {'cCL', 'cCR'};
    Special_trial= 'x';
elseif nargin==2;   
     display(['Other trial type = ' Other_trial_type]);
     psth_trial_type = {'cCL', 'cCR', Other_trial_type};
     Special_trial= Other_trial_type;
elseif nargin>2; 
    display('TOO MANY INPUT ARGUMENTS: max=2');
end

psth_trig_evt =  psth_center_evtID; % 'Delay', % 'APuff' % 'GoCue' % 'Valve' %'Licks'
NbTtrialType = max(size(psth_trial_type));
%% load data
load('info.mat'); MouseID = info.info_notes.MouseID;  Day=info.info_notes.Day; 
load('evt.mat');
load('time.mat');
load('Ntrial_type.mat');
sr=info.info_freq_parameters.board_dig_in_sample_rate;
FileLocation=pwd; 

%% GET trigtimes : times of the psth_trig_evt is sec   (1 vector of times)
% define trial idx
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);

idx_trial_start = trig_end(1:end-1); %-(1*sr);
idx_trial_end = trig_st(2:end) ;     %-(1*sr);

Ntrials= trial.Ntrial; %if Ntrials==trial.Ntrial; disp('good Ntrial'); else forcestop; end ;

% define common events
evt_valve = evt_rwd_L + evt_rwd_R;
evt_lick = evt_lick_L + evt_lick_R;
evt_puff =evt_puff_L + evt_puff_R;
evt_prelick = evt_delay + evt_puff;
%%

trigtimes=[];
for tr=1:Ntrials
    time_tr = time(idx_trial_start(tr):idx_trial_end(tr));
    if psth_trig_evt=='APuff'
        puff_tr=  evt_puff(idx_trial_start(tr):idx_trial_end(tr));
        idx_puff_st = find(diff(puff_tr)>0); % start of the puff
        time_puff_st = time_tr(idx_puff_st);
        trigtimes = [trigtimes time_puff_st]; % in sec
        
    elseif psth_trig_evt=='Delay'
        if FileLocation(25:30)=='vgat11' | FileLocation(25:30)=='vgat12'
            puff_tr=  evt_puff(idx_trial_start(tr):idx_trial_end(tr));
            idx_puff_st = find(diff(puff_tr)>0); % start of the puff
            time_delay_st = time_tr(idx_puff_st)+0.75;
            trigtimes = [trigtimes time_delay_st]; % in sec
        else
            delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
            idx_delay_st = find(diff(delay_tr)>0)  % start of the delay
            time_delay_st = time_tr(idx_delay_st);
            trigtimes = [trigtimes time_delay_st]; % in sec
        end
        
    elseif psth_trig_evt=='GoCue'
        
        if FileLocation(25:30)=='vgat11' | FileLocation(25:30)=='vgat12'
            puff_tr=  evt_puff(idx_trial_start(tr):idx_trial_end(tr));
            idx_puff_st = find(diff(puff_tr)>0); % start of the puff
            time_GO_st = time_tr(idx_puff_st)+1.5;
            trigtimes = [trigtimes time_GO_st]; % in sec         
        else
            delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
            idx_GO_st = find(diff(delay_tr)<0); % end of the delay
            time_GO_st = time_tr(idx_GO_st);
            trigtimes = [trigtimes time_GO_st]; % in sec
        end
        
    elseif psth_trig_evt=='Licks'
        lick_tr=  evt_lick(idx_trial_start(tr):idx_trial_end(tr));
        idx_lick_st = min(find(diff(lick_tr)>0));% start of the first lick
        if isempty(idx_lick_st) % replace by Lick evt by GOcue evt for centering during NOlick trial
            delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
            idx_GO_st = find(diff(delay_tr)<0)
            idx_lick_st=idx_GO_st;
        end
        time_lick_st = time_tr(idx_lick_st);
        trigtimes = [trigtimes time_lick_st]; % in sec
        
    else
        disp(' WRONG INPUT ARG pick one: psth_trig_evt = GoCue, APuff, Delay, Licks or Valve' )
    end
end

%%
trigtimes_all=trigtimes';
trigtimes_cCL = trigtimes(trial.idx_correct_L);
trigtimes_cCR = trigtimes(trial.idx_correct_R);
trigtimes_iCL = trigtimes([trial.idx_errorDelay_PL_CL trial.idx_errorDelay_PR_CL]);
trigtimes_iCR = trigtimes([trial.idx_errorDelay_PL_CR trial.idx_errorDelay_PR_CR]);
trigtimes_oCR = trigtimes([trial.idx_correct_R_opto]);
trigtimes_oCL = trigtimes([trial.idx_correct_L_opto]);
trigtimes_ocC = trigtimes([trial.idx_correct_L_opto trial.idx_correct_R_opto]);
trigtimes_oNO = trigtimes([trial.idx_NoLick_opto]);
trigtimes_eNO = trigtimes([trial.idx_NoLick]);
trigtimes_ePL = trigtimes([trial.idx_errorResp_PL]);
trigtimes_ePR = trigtimes([trial.idx_errorResp_PR]);

 if MouseID == 'vgat11'; trigtimes=trigtimes(1:20); end 


%% GET spxtimes = SPIKE times in Sec (1 vector of times)
close all,
clust_file = dir('times_*S*Ch*_sub.mat');
for ncluf = 1:max(size(clust_file)) %  channel loop
    
    load(clust_file(ncluf).name); %  e.g.: load('times_S2Ch6_man.mat')
    Nclust = max(cluster_class(:,1));
    
    for CLUST= 1:Nclust %Cluster loop within each channel
        idx_spk = find(cluster_class(:,1)==CLUST);
        spxtimes_ms = cluster_class(idx_spk,2);
        % convert time to sec
        spxtimes = spxtimes_ms/10^3;
        Nspk=max(size(spxtimes));
        
        figure,
        for ii=1:NbTtrialType
            
            if psth_trial_type{ii} == 'all'
                trigtimes = trigtimes_all ;
            elseif psth_trial_type{ii} == 'cCL'
                trigtimes = trigtimes_cCL ;
            elseif psth_trial_type{ii} == 'cCR'
                trigtimes = trigtimes_cCR ;
            elseif psth_trial_type{ii} == 'iCL'
                trigtimes = trigtimes_iCL ;
            elseif psth_trial_type{ii} == 'iCR'
                trigtimes = trigtimes_iCR ;
            elseif psth_trial_type{ii} == 'oCR' % opto correct Choice right
                trigtimes = trigtimes_oCR ;
            elseif psth_trial_type{ii} == 'oCL' 
                trigtimes = trigtimes_oCL ;
            elseif psth_trial_type{ii} == 'ocC' % opto correct Choices l+r 
                trigtimes = trigtimes_ocC ;
            elseif psth_trial_type{ii} == 'oNO' % opto NOlick
                trigtimes = trigtimes_oNO ;
            elseif psth_trial_type{ii} == 'ePL' % error Side Puff left 
                trigtimes = trigtimes_ePL ;
            elseif psth_trial_type{ii} == 'ePR'
                trigtimes = trigtimes_ePR ;
            elseif psth_trial_type{ii} == 'eNO' % error NOlick
                trigtimes = trigtimes_eNO ; 
            end
            
            if MouseID == 'vgat11'; trigtimes=trigtimes(1:25); end
            
            pre=2500; %ms
            post=1750;%ms
            fr=1; %if '1', normalizes to firing rate (Hz); if '0', does nothing (default)
            binsz=1; %ms  bin size of psth (default: 1 ms)
            [psth trialspx] = mpsth(spxtimes, trigtimes, 'pre', pre, 'post', post, 'fr',fr, 'binsz', binsz, 'chart',0);
            
            subplot(NbTtrialType+2,1,ii)
            % preallocate
            rastmat = zeros(numel(trialspx),pre+1+post);
            timevec = -pre:1:post;
            % generate and plot raster
            for i = 1:numel(trialspx)
                rastmat(i,trialspx{i}+pre+1) = 1;
                if ii==1
                    plot(timevec(rastmat(i,:)~=0),rastmat(i,rastmat(i,:)~=0)*i,'r.','MarkerSize',5),hold on
%                     title(['raster center on ' psth_trig_evt ' ' psth_trial_type{ii} 'trials'])
                    legend({['Correct ' psth_trial_type{ii}]}, 'FontSize',10,'FontWeight', 'normal') 
                elseif ii==2
                    plot(timevec(rastmat(i,:)~=0),rastmat(i,rastmat(i,:)~=0)*i,'b.','MarkerSize',5),hold on
%                     title(['raster center on ' psth_trig_evt ' ' psth_trial_type{ii} 'trials'])      
                    legend({['Correct ' psth_trial_type{ii}]}, 'FontSize',10,'FontWeight', 'normal') 
                else
                    col = [0.3 0.3 0.3]+((ii/30)*(ii-2))
                    plot(timevec(rastmat(i,:)~=0),rastmat(i,rastmat(i,:)~=0)*i,'.','Color',col,'MarkerSize',5),hold on
                    title([ 'raster center on ' psth_trig_evt ' ' psth_trial_type{ii} 'trials'])
                    if Special_trial(1)=='i' 
                        legend({['Impulse ' psth_trial_type{ii}]}, 'FontSize',10,'FontWeight', 'normal')
                    elseif Special_trial(1)=='o' 
                        legend({['Opto ' psth_trial_type{ii}]}, 'FontSize',10,'FontWeight', 'normal') 
                    elseif Special_trial(1)=='e'
                        legend({['Error ' psth_trial_type{ii}]}, 'FontSize',10,'FontWeight', 'normal')
                    end
                    
                end
            end
            xlim([-2250 1500])
            gca, axis([-pre+10 post+10 0.5 numel(trialspx)+0.5]);
            xlabel(['time from ' psth_trig_evt '(ms)']);
            ylabel('trials');
            % Set the 'visible' property 'off'
            ax = gca;
            ax.Visible = 'off';
            xlim([-2250 1500])
             
            subplot(NbTtrialType+2,1,[NbTtrialType+1 NbTtrialType+2]),
            ftype = 'Gauss' % boxcar, Gauss, exp, exGauss
            w = 75 %ms
            [sdf kernel] = msdf(psth,ftype,w);
            if ii==1
                hold on, plot(sdf(:,1), sdf(:,2),'r','LineWidth', 2);
            elseif ii==2
                hold on, plot(sdf(:,1), sdf(:,2),'b','LineWidth', 2);
            else
                col = [0.3 0.3 0.3]+((ii/30)*(ii-2));
                hold on, plot(sdf(:,1), sdf(:,2),'color',col,'LineWidth', 2)
            end
            
            if  fr==0
                ylabel(['counts per ' num2str(binsz) 'ms bin'],'FontSize',11,'FontWeight', 'normal')
            elseif fr==1
                ylabel(['fr (Hz)'],'FontSize',11,'FontWeight', 'normal')
            end
            xlabel(['time from ' psth_trig_evt ' start (ms)'], 'FontSize',11,'FontWeight', 'normal'); 
            
        end
        tagsave = [];
        for ptt=1:NbTtrialType
            tagsave = [tagsave psth_trial_type{ptt} '_'];
        end
        xlim([-2250 1500])
%         legend('cCLeft', 'cCRight' , ['error' psth_trial_type{ptt}],   'Location', 'northeast')
        title([info.info_notes.MouseID ' ' info.info_notes.Day ' '  clust_file(ncluf).name(7:11) ' clust#' num2str(CLUST)],'FontSize',11,'FontWeight', 'bold')
        saveas(gcf, ['fig2print_'  clust_file(ncluf).name(7:11) '_clust#' num2str(CLUST) '_SDF_' tagsave 'center_' psth_trig_evt ], 'jpg')
    end
end

