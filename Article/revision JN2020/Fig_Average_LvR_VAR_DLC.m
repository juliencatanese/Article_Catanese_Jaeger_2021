function [idxL, idxR] = Fig_Average_LvR_VAR_DLC(MouthY_tr, TongueProba_tr, VAR_DLC, VARname_str)
% function [idxL, idxR] = Fig_Average_LvR_VAR_DLC(MouthY_tr, TongueProba_tr, VAR_DLC, VARname_str);
% Julien Catanese 9/23/2020

VAR=MouthY_tr;
TITLE_EVT = 'MouthY';
Thr = median(mean(VAR));

% figure, plot(VAR), hold on, plot(mean(VAR,2),'r','LineWidth',2) 

AV=[]; 
for itr=1:125;
idxT = find(TongueProba_tr(:,itr)>0.95);
AV = [AV mean(VAR(idxT,itr))];
end

idxL = find(AV>Thr);
idxR = find(AV<Thr);

figure, hold on, 
plot(mean(VAR(:,idxL),2),'r','LineWidth',2);
plot(mean(VAR(:,idxR),2),'b','LineWidth',2);

xline(1.5*25,'k--','LineWidth',2);

plot(ones(1,100)*Thr,'k--');
for ii=5:25
plot(VAR(:,idxL(ii)),'m');
plot(VAR(:,idxR(ii)),'c');
end
plot(mean(VAR(:,idxL),2),'r','LineWidth',2);
plot(mean(VAR(:,idxR),2),'b','LineWidth',2);


maxR=max(abs(mean(VAR(:,idxR)-Thr,2))) ;
maxL=max(abs(mean(VAR(:,idxL)-Thr,2)));
ylim([Thr-max(maxL,maxR)-5 Thr+max(maxL,maxR)+5]);

legend('Left','Right','GoCue');
title(TITLE_EVT);
saveas(gcf,['D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Trial_Average_DLC_' TITLE_EVT  '.png']);

%% Plot VAR using idxL and idxR

VAR=VAR_DLC;
TITLE_EVT = VARname_str;
Thr = median(mean(VAR));

figure, hold on, 
plot(mean(VAR(:,idxL),2),'r','LineWidth',2);
plot(mean(VAR(:,idxR),2),'b','LineWidth',2);
xline(1.5*25,'k--','LineWidth',2);

plot(ones(1,100)*Thr,'k--');
for ii=5:25
plot(VAR(:,idxL(ii)),'m');
plot(VAR(:,idxR(ii)),'c');
end
plot(mean(VAR(:,idxL),2),'r','LineWidth',2);
plot(mean(VAR(:,idxR),2),'b','LineWidth',2);

maxR=max(abs(mean(VAR(:,idxR)-Thr,2))) ;
maxL=max(abs(mean(VAR(:,idxL)-Thr,2)));
ylim([Thr-max(maxL,maxR)-5 Thr+max(maxL,maxR)+5]);

legend('Left','Right','GoCue');
title(TITLE_EVT);

saveas(gcf,['D:\DATA EMORY\JC_Analysis\deeplabcut_fig\Trial_Average_DLC_' TITLE_EVT  '.png']);

