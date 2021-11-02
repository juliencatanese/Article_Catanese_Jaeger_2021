% pub_fig_FRpop_JCscript
% Fig 2A : Population firing rate
% written by Julien Catanese 3/15/2019

% parfig.BaselineEpoch= [150:2150]; BLE = parfig.BaselineEpoch;
% parfig.center_evt =  'GoCue'; % center for SDF ('Delay' ; 'GOcue' ; Licks)
% parfig.trial_type = {'cor'};
% parfig.col = {'k'};
% parfig.xlim = [-2750 1250];
% if parfig.center_evt == 'Delay'
%     parfig.pre = BLE(end) + 750  ;
%     parfig.post = 1500 + 750 ;
% elseif parfig.center_evt == 'GoCue'
%     parfig.pre = BLE(end) + 1500  ;
%     parfig.post = 1500 ;
% end
% listcell2= []; listcell2=listcell;
% 
% % Fig 2A : Population firing rate
% % load('listcell.mat')
% % load('SMA_cor_GoCue.mat', 'SMA')

Hpuf_pos = (Tephys.H_puf & (Tephys.Fr_puf_Hz > Tephys.Fr_BLE_Hz)) &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz)));

Hdel_pos = (Tephys.H_del & (Tephys.Fr_del_Hz > Tephys.Fr_BLE_Hz)) &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz)));

Hres_pos = (Tephys.H_res & (Tephys.Fr_res_Hz > Tephys.Fr_BLE_Hz)) &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz)));

idx = [ logical(Tcombo.VMVL), ...
    logical(Tcombo.VMVL & (Hpuf_pos | Hdel_pos | Hdel_pos)), ...
    logical(Tcombo.VMVL & (~Hpuf_pos | ~Hdel_pos | ~Hdel_pos))];

figure,
Xt=[1:1:size(SMA,2)]-parfig.pre;
col = {'k','m','g'}
LegText={'ALL VMVL', 'Exc VMVL' , 'Inh VMVL'}
legend_all = []; idx1=[]; SMA2=[];
for ii= 1:size(idx,2); ii
    SMA2=SMA(idx(:,ii),:);
    fstr= col{ii};
    idx1=idx(:,ii);
    hold on,
    plotshaded(Xt,[nanmean(SMA2)+(nanstd(SMA2)/sqrt(sum(idx1))) ;nanmean(SMA2); nanmean(SMA2)-(nanstd(SMA2)/sqrt(sum(idx1)));], fstr);
    legend_all = [legend_all ;  LegText{ii}];
end

xlim(parfig.xlim);
ylabel('Fr (Hz)') ;
xlabel('time (ms)');
title('Average Firing rate over all cell');

pre = 0; %parfig.pre; 
if parfig.center_evt == 'GoCue'
    hold on, plot([pre-1500 pre-1500], ylim,'k--','LineWidth',2.5);
    hold on, plot([pre-750 pre-750], ylim,'k--', 'LineWidth',2.5);
    hold on, plot([pre pre], ylim,'k--','LineWidth',2.5);
elseif parfig.center_evt == 'Delay'
    hold on, plot([pre-750 pre-750], ylim,'k--', 'LineWidth',2.5);
    hold on, plot([pre pre], ylim,'k--','LineWidth',2.5);
    hold on, plot([pre+750 pre+750], ylim,'k--','LineWidth',2.5);
end

legend('All VMVL', '' , 'Exc VMVL' , '' , 'Inh VMVL' , '' )
