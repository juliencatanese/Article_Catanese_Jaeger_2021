function Plot_Raster_SDFshaded_JCfun(SMA, SSA, spxtimes, trigtimes, cell2plot, parfig)
% function Plot_Raster_SDFshaded_JCfun(SMA, SSA, spxtimes, trigtimes)

% param 
pre = parfig.pre;
post = parfig.post;
center_evt = parfig.center_evt;
trial_type = parfig.trial_type;
col= parfig.col;

MouseID = cell2plot{1,2};
Day = cell2plot{1,3};
ChanID = cell2plot{1,4};
CLUST = cell2plot{1,5};

NbTtrialType = size(SMA,1)
fr=1; %if '1', normalizes to firing rate (Hz); if '0', does nothing (default)
binsz=1; %ms  bin size of psth (default: 1 ms)



%% PLOT RASTER

figure, 
for ii=1:NbTtrialType
    
    % Compute Raster/PSTH
    [psth trialspx] = mpsth(spxtimes(ii,:), trigtimes(ii,:), 'pre', pre, 'post', post, 'fr',fr, 'binsz', binsz, 'chart',0, 'tb',0);

    subplot(NbTtrialType+2,1,ii)
    rastmat = zeros(numel(trialspx),pre+1+post); size(rastmat)
    timevec = -pre:1:post;
    for i = 1:numel(trialspx)
        rastmat(i,trialspx{i}+pre+1) = 1;
        plot(timevec(rastmat(i,:)~=0),rastmat(i,rastmat(i,:)~=0)*i,'.','Color',col{ii},'MarkerSize',5),hold on
%         legend({['trial ' trial_type{ii}]}, 'FontSize',11,'FontWeight', 'normal')
    end
    gca, axis([-pre+10 post+10 0.5 numel(trialspx)+0.5]);
    xlabel(['time from ' center_evt '(ms)']);
    ylabel('trials');
    xlim([parfig.xlim])
    ax = gca; ax.Visible = 'off';  % Set the 'visible' property 'off'
    
    
    %% PLOT SDF
    subplot(NbTtrialType+2,1,[NbTtrialType+1 NbTtrialType+2]),
    
    [sdf kernel] = msdf(psth,'Gauss',100);
    hold on, plot([-pre:1:post], sdf,'color',col{ii},'LineWidth',2);
    hold on, plotshaded([-pre:1:post], [sdf'+ (SSA(ii,:)); sdf'- (SSA(ii,:))] ,col{ii});

    pos=[-1500; -750; 0; +750; +1500]; 
    ym = max(ylim);
    YY = [0 ym];
    for ipos=1:max(size(pos))
        hold on, line([pos(ipos) pos(ipos)], YY, 'Color','k','LineStyle','--','LineWidth',1 )
    end
    
end
%%
ylabel(['fr (Hz)'],'FontSize',11,'FontWeight', 'normal')
xlabel(['time from ' center_evt ' start (ms)'], 'FontSize',11,'FontWeight', 'normal');

xlim([parfig.xlim])
ylim([0 ym])

title([MouseID ' ' Day ' ' ChanID ' clust#' num2str(CLUST)],'FontSize',11,'FontWeight', 'bold')
