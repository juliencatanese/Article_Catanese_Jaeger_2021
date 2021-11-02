mypath = 'D:\DATA EMORY\JC_Analysis'
cd(mypath)
mkdir('ISI')
close all,
load ('listcell3')
load('Tfig1_VMopto.mat')
ncell= size(listcell,1)
figure(1),
percentRFS_all =[]; all_CLUperChan = []
for ii=1:ncell
    load ([listcell.SessPath(ii,:) '\' listcell.SessID{ii} '\' listcell.CluFile2load(ii,:) ])
    
    % read cluster_class (clu#, time in ms)
    CLU=listcell.CLUST(ii);
    clu_idx = find(cluster_class(:,1)==CLU);
    clu_times = cluster_class(clu_idx,2)/1000; % convert ms to sec
    
    % Calculate percent of spike in refractory period (3ms)
    Refractory = 0.003; % 2.5ms
    
    nSPKtot = size(clu_idx,1);
    ISI = diff(clu_times);
    
    nRFS = sum(ISI<Refractory); % nb spikes in refractory period
    percentRFS = (nRFS/nSPKtot)*100; % nb spikes in refractory period
    percentRFS_all =[percentRFS_all; percentRFS];
    
    
    if ii~=ncell && sum(listcell.ChanID(ii,:)~= listcell.ChanID(ii+1,:))~=0
        all_CLUperChan = [all_CLUperChan; CLU];
    elseif ii==ncell
        all_CLUperChan = [all_CLUperChan; CLU];
    end
    
    BW = 1;%ms
    BMIN = 1
    BMAX = 100
    close(1),figure(1),h=histogram(ISI*1000,'BinWidth',BW,'BinLimits',[BMIN, BMAX])
    h.FaceColor = 'k'; h.EdgeColor = 'k'; h.BinEdges([1:100])
    legend([ 'rpv =' num2str(percentRFS) '%  (with rp=' num2str(Refractory*1000) 'ms)'])
    title([ Tfig1_VMopto.AreaID{ii} ' ' listcell.Properties.RowNames{ii} '  binW=' num2str(BW)])
    %     saveas (gcf,[listcell.SessPath(ii,:) '\' listcell.SessID{ii} '\ISI_cell#' num2str(ii) '_bmax' num2str(BMAX) 'ms'],'png')
%     saveas (gcf,[ mypath '\ISI\ISI_cell#' num2str(ii) '_bmax' num2str(BMAX) 'ms'] ,'png');
        saveas (gcf,[ mypath '\ISI\ISI_cell#' num2str(ii) '_bmax' num2str(BMAX) 'ms'] ,'emf');

end
%%
AA=sum(percentRFS_all>1.5)
BB= max(percentRFS_all);
CC=find(percentRFS_all>1.5);
DD=percentRFS_all(CC)
listcell(CC,:)

figure,
BW = 0.1;  %
BMIN = -0.1;
BMAX = 7;
h2=histogram(percentRFS_all,'BinWidth',BW,'BinLimits',[BMIN, BMAX])
h2.FaceColor = 'y'; h2.EdgeColor = 'k'; 
ylabel('cell count')
xlabel('% rpv')
MedianRFS = median(percentRFS_all)
MeanRFS = mean(percentRFS_all)
StdRFS = std(percentRFS_all)
title(['refractory period = ' num2str(Refractory*1000) 'ms'])
xlim([0 BMAX ])
saveas (gcf,[ mypath '\ISI\DISTRIBUTION_rpv'] ,'emf');
saveas (gcf,[ mypath '\ISI\DISTRIBUTION_rpv'] ,'png');
saveas (gcf,[ mypath '\ISI\DISTRIBUTION_rpv'] ,'fig');

%% Nb of cluster per channel
nclu= sum(all_CLUperChan)
meannclu= mean(all_CLUperChan)
stdnclu = std(all_CLUperChan)
minnclu= min(all_CLUperChan)
maxclu = max(all_CLUperChan)
figure, histogram(all_CLUperChan)
ylabel('count')
xlabel('#cluster per channel')