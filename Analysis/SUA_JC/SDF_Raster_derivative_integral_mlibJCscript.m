%  SDF_Raster_Correct_mlibJCscript
%  by JC 11/22/2018

%% Define NbTtrialType
try
    NbTtrialType = max(size(psth_trial_type))
catch
    psth_trig_evt =  'GoCue'; % center for SDF
    psth_trial_type = {'cor'};
    NbTtrialType = max(size(psth_trial_type))
end

%% load data
load('info.mat'); Day=info.info_notes.Day;  sr=info.info_freq_parameters.board_dig_in_sample_rate;
load('evt.mat');
load('time.mat');
load('Ntrial_type.mat');

%% GET trigtimes : times of the psth_trig_evt is sec   (1 vector of times)
% define trial idx
trig_end = find(diff(evt_trial)<0);
trig_st = find(diff(evt_trial)>0);

idx_trial_start = trig_end(1:end-1); %-(1*sr);
idx_trial_end = trig_st(2:end) ;     %-(1*sr);

Ntrials= trial.Ntrial;

% define common events
evt_valve = evt_rwd_L + evt_rwd_R;
evt_lick = evt_lick_L + evt_lick_R;
evt_puff =evt_puff_L + evt_puff_R;
evt_prelick = evt_delay + evt_puff;

%% define trigtimes for each trials
trigtimes=[];
for tr=1:Ntrials
    time_tr = time(idx_trial_start(tr):idx_trial_end(tr));
    if  psth_trig_evt=='Delay'
        delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
        idx_delay_st = find(diff(delay_tr)>0)  % start of the delay
        time_delay_st = time_tr(idx_delay_st);
        trigtimes = [trigtimes time_delay_st]; % in sec
        
    elseif psth_trig_evt=='GoCue'
        delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
        idx_GO_st = find(diff(delay_tr)<0); % end of the delay
        time_GO_st = time_tr(idx_GO_st);
        trigtimes = [trigtimes time_GO_st]; % in sec
        
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

%% define trigtime

trigtimes_all=trigtimes';
trigtimes_cor = sort([trigtimes(trial.idx_correct_L) trigtimes(trial.idx_correct_R)]);

%% GET spxtimes = SPIKE times in Sec (1 vector of times)
clust_file = dir('times_*S*Ch*_sub.mat');
for ncluf = 1:max(size(clust_file)) %  channel loop
    load(clust_file(ncluf).name); %  e.g.: load('times_S2Ch6_man.mat')
    Nclust = max(cluster_class(:,1));
    ChanID = clust_file(ncluf).name(7:11)
    
    for CLUST= 1:Nclust %Cluster loop within each channel
        idx_spk = find(cluster_class(:,1)==CLUST);
        spxtimes_ms = cluster_class(idx_spk,2);
        spxtimes = spxtimes_ms/10^3;        % convert time to sec
        
        figure,
        for ii=1:NbTtrialType
            
            if psth_trial_type{ii} == 'all'
                trigtimes = trigtimes_all ;
            elseif psth_trial_type{ii} == 'cor' % error NOlick
                trigtimes = trigtimes_cor ;
            end
            
            %% SET PARAMETERS
            pre=1500; %ms
            post=1500;%ms
            fr=1; %if '1', normalizes to firing rate (Hz); if '0', does nothing (default)
            binsz=1; %ms  bin size of psth (default: 1 ms)
            [psth trialspx] = mpsth(spxtimes, trigtimes, 'pre', pre, 'post', post, 'fr',fr, 'binsz', binsz, 'chart',0);
            
            %% PLOT RASTER (top pannel)
            subplot(4,1,ii)
            rastmat = zeros(numel(trialspx),pre+1+post);
            timevec = -pre:1:post;
            for i = 1:numel(trialspx)
                rastmat(i,trialspx{i}+pre+1) = 1;
                col = [0.1 0.1 0.3]+((ii/30)*(ii-2))
                plot(timevec(rastmat(i,:)~=0),rastmat(i,rastmat(i,:)~=0)*i,'b.','MarkerSize',5),hold on
                legend({['Correct ' psth_trial_type{ii}]}, 'FontSize',11,'FontWeight', 'normal')
            end
            gca, axis([-pre+10 post+10 0.5 numel(trialspx)+0.5]);
            xlabel(['time from ' psth_trig_evt '(ms)']);
            ylabel('trials');
            ax = gca; ax.Visible = 'off'; % Set the 'visible' property 'off'
            
            %% Calculate SDF and derivative and Intergral
            ftype = 'Gauss' % boxcar, Gauss, exp, exGauss
            w = 100 %ms
            [sdf kernel] = msdf(psth,ftype,w);
            sdf_norm = sdf(:,2)/max(sdf(:,2));
            sdf_der = diff(sdf(:,2))/max(diff(sdf(:,2)));
%             sdf_int = int(sdf(:,2))/max(diff(sdf(:,2)));



            %% PLOT SDF derivative and Intergral
            subplot(4,1,[4]),
            hold on, plot(sdf(:,1), sdf_norm,'color','b','LineWidth',1 );
            hold on, plot(sdf(1:end-1,1), abs(sdf_der),'--m','LineWidth',1 );
%             hold on, plot(sdf(1:end-1,1), sdf_int,'color','--y','LineWidth',1 );
            
            pos= [-1500 -750 0 750 1500]
            YY = [0 1]
            for ipos=1:max(size(pos))
                hold on, line([pos(ipos) pos(ipos)], YY, 'Color','k','LineStyle','--')
            end
            xlabel(['time from ' psth_trig_evt ' start (ms)'], 'FontSize',11,'FontWeight', 'normal');
            ylabel(['norm fr'],'FontSize',11,'FontWeight', 'normal')
                         %% PLOT SDF 
            subplot(4,1,[2 3]),
            hold on, plot(sdf(:,1), sdf(:,2),'color','b','LineWidth',2 );
            
            pos= [-1500 -750 0 750 1500]
            ym = max(sdf(:,2))
            YY = [0 ym]
            for ipos=1:max(size(pos))
                hold on, line([pos(ipos) pos(ipos)], YY, 'Color','k','LineStyle','--')
            end
            ylabel(['fr (Hz)'],'FontSize',11,'FontWeight', 'normal')
            
        end
        title([MouseID ' ' Day ' ' ChanID ' clust#' num2str(CLUST)],'FontSize',11,'FontWeight', 'bold')
        
        
        tag_trtype = [];
        for ptt=1:NbTtrialType
            tag_trtype = [tag_trtype psth_trial_type{ptt}];
        end
        
%         saveas(gcf, ['fig2print_'  clust_file(ncluf).name(7:11) '_sub_clust#' num2str(CLUST) '_SDF_' tag_trtype '_' psth_trig_evt ], 'jpg')
mkdir('D:\JC_Figures\SUA\SDF\derivative')
saveas(gcf, ['D:\JC_Figures\SUA\SDF\derivative\' MouseID '_' Day  '_' ChanID '_clust#' num2str(CLUST) '_SDFder_' tag_trtype '_' psth_trig_evt ], 'png')
        
    end
end

