% pub_fig3_LvR_shaded_dzSMA_JCscript
% plot and compare pair of shaded curve (mean +/- k*sem)
% written by Julien Catanese 04/09/2019


%% Contra + ipsi (ABS MEAN without complex)
figure,
idx = [ (contra_dz & ~ipsi_dz), (~contra_dz & ipsi_dz)]
col = {'b','r'}
for ii= 1:2
    SDF2=dzSMA(idx(:,ii),:); fstr= col{ii};  idx1=idx(:,ii); hold on,
    plotshaded(Xt,[abs(mean(SDF2))+(K*(nanstd(SDF2)/sqrt(sum(idx1)))) ;abs(mean(SDF2)); abs(mean(SDF2))-(K*(nanstd(SDF2)/sqrt(sum(idx1))));], fstr);
end
xlim(parfig.xlim); ylim([0 2.5])
ylabel('abs(mean(dz))'), xlabel('time ms')
legend('cont 100 cells','', 'ipsi 99 cells', '', 'Location','northwest')
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);

%% contra + ipsi (ABS MEAN including complex)  
figure,
idx = [ contra_dz, ipsi_dz];
col = {'b','r'};
for ii= 1:size(idx,2); ii
    SDF2=dzSMA(idx(:,ii),:); idx1=idx(:,ii); fstr= col{ii}; hold on;
    plotshaded(Xt,[abs(mean(SDF2))+(K*(nanstd(SDF2)/sqrt(sum(idx1)))) ;abs(mean(SDF2)); abs(mean(SDF2))-(K*(nanstd(SDF2)/sqrt(sum(idx1))));], fstr);
end
legend('cont 210 cells', 'ipsi 209 cells', 'Location', 'northwest');
xlim(parfig.xlim);  ylim([0 2.5])
ylabel('abs(mean(dz))'), xlabel('time ms')
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);


%% contra + ipsi (MEAN without complex)
figure,
idx = [ (contra_dz & ~ipsi_dz), (~contra_dz & ipsi_dz), ]
Xt=[1:1:size(dzSMA,2)]-parfig.pre;
col = {'b','r'}
for ii= 1:2; %size(idx,2); ii
    SDF2=dzSMA(idx(:,ii),:); fstr= col{ii};  idx1=idx(:,ii); hold on,
    plotshaded(Xt,[mean(SDF2)+(K*(nanstd(SDF2)/sqrt(sum(idx1)))) ;mean(SDF2); mean(SDF2)-(K*(nanstd(SDF2)/sqrt(sum(idx1))));], fstr);
end
xlim(parfig.xlim)
ylim([-2.5 2.5])
ylabel('mean(dz)'), xlabel('time ms')
legend('cont 100 cells','', 'ipsi 99 cells', '', 'Location','northwest')
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);


%% Grand Average + Complex only (MEAN)
figure,
idx = [(contra_dz & ipsi_dz), Tcombo.VMVL & Tcombo.z_exct] % (~contra_dz & ~ipsi_dz & Tcombo.VMVL)]
Xt=[1:1:size(dzSMA,2)]-parfig.pre;
col = {'m','k'}
for ii= 1:2; %size(idx,2); ii
    SDF2=dzSMA(idx(:,ii),:); fstr= col{ii};  idx1=idx(:,ii); hold on,
    plotshaded(Xt,[mean(SDF2)+(K*(nanstd(SDF2)/sqrt(sum(idx1)))) ;mean(SDF2); mean(SDF2)-(K*(nanstd(SDF2)/sqrt(sum(idx1))));], fstr);
end
xlim(parfig.xlim)
ylim([-2.5 2.5])
ylabel('mean(dz)'), xlabel('time ms')
legend('ALL 338 cells','', 'Complex 110 cells', '', 'Location','northwest')
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);

%% No sig only (MEAN) 
figure,
idx = [ ~contra_dz & ~ipsi_dz & Tcombo.VMVL & Tcombo.z_exct ]
Xt=[1:1:size(dzSMA,2)]-parfig.pre;
col = {'y'}
for ii= 1 %size(idx,2); ii
    SDF2=dzSMA(idx(:,ii),:); fstr= col{ii};  idx1=idx(:,ii); hold on,
    plotshaded(Xt,[mean(SDF2)+(K*(nanstd(SDF2)/sqrt(sum(idx1)))) ;mean(SDF2); mean(SDF2)-(K*(nanstd(SDF2)/sqrt(sum(idx1))));], fstr);
end
xlim(parfig.xlim)
ylim([-2.5 2.5])
ylabel('mean(dz)'), xlabel('time ms')
legend('nosig 28 cells','', 'Location','northwest')
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);


%% contra + ipsi (MEAN with complex + all)
figure,
idx = [ (contra_dz & ~ipsi_dz), (~contra_dz & ipsi_dz), ]
Xt=[1:1:size(dzSMA,2)]-parfig.pre;
col = {'b','r'}
for ii= 1:2; %size(idx,2); ii
    SDF2=dzSMA(idx(:,ii),:); fstr= col{ii};  idx1=idx(:,ii); hold on,
    plotshaded(Xt,[mean(SDF2)+(K*(nanstd(SDF2)/sqrt(sum(idx1)))) ;mean(SDF2); mean(SDF2)-(K*(nanstd(SDF2)/sqrt(sum(idx1))));], fstr);
end
xlim(parfig.xlim)
ylim([-2.5 2.5])
ylabel('mean(dz)'), xlabel('time ms')
% legend('ALL 338 cells','', 'Complex 110 cells', '', 'Location','northwest')
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);

figure

idx = [(contra_dz & ipsi_dz), Tcombo.VMVL & Tcombo.z_exct] % (~contra_dz & ~ipsi_dz & Tcombo.VMVL)]
Xt=[1:1:size(dzSMA,2)]-parfig.pre;
col = {'m','k'}
for ii= 1:2; %size(idx,2); ii
    SDF2=dzSMA(idx(:,ii),:); fstr= col{ii};  idx1=idx(:,ii); hold on,
    plotshaded(Xt,[mean(SDF2)+(K*(nanstd(SDF2)/sqrt(sum(idx1)))) ;mean(SDF2); mean(SDF2)-(K*(nanstd(SDF2)/sqrt(sum(idx1))));], fstr);
end
xlim(parfig.xlim)
ylim([-2.5 2.5])
ylabel('mean(dz)'), xlabel('time ms')
% legend('ALL 338 cells','', 'Complex 110 cells', '', 'Location','northwest')
hold on, plot([-1500 -1500], ylim,'k--','LineWidth',1.5);
hold on, plot([-750 -750], ylim,'k--', 'LineWidth',1.5);
hold on, plot([0 0], ylim,'k--','LineWidth',1.5);
hold on, plot(Xt, zeros(size(Xt)),'k--','LineWidth',1.5);
