function pub_fig_SMA_SUA_JCfun2(listc, SDF, parfig)
% function pub_fig_SMA_JCfun(listcell, parfig)
% fast plot SMA only for vizualisation of a specific list of cells
% plot SDF all
% written by Julien Catanese the 12 March 2019

try 
    close 100 101 102 103 104 105 106 120
    close 101 
    close 102 
    close 103 
    close 104
catch 
end

legendtext= [];
Xpos = [-parfig.pre:1:parfig.post];
nctot= size(listc,1)
for nc = 1:nctot;
    %     listc(nc,:)
    if nc<8,
        figure(100);hold on; 
    elseif nc >= 8 & nc < 15;
        gcf; xlim([parfig.xlim]); legend(legendtext(1:7,:)); title(parfig.title);
         figure(101);hold on;
    elseif nc >= 15 & nc < 21
        gcf; xlim([parfig.xlim]); legend(legendtext(8:14,:)); title(parfig.title);
        figure(102);hold on; 
    elseif nc >= 21 & nc < 28;
        gcf; xlim([parfig.xlim]); legend(legendtext(15:20,:)); title(parfig.title);
         figure(103),hold on;
    elseif nc >= 100 & nc < 107;
        gcf; xlim([parfig.xlim]); legend(legendtext(15:20,:)); title(parfig.title);
         figure(104),hold on;
         
    else
        gcf; xlim([parfig.xlim]); legend(legendtext(21:27,:)); title(parfig.title);
        figure(120); hold on;
    end
    hold on;
    plot(Xpos, SDF(listc.idx_list(nc),:));
    
    RowName = listc.Row(nc);
    legendtext= [legendtext ; RowName{1} ];
   
end
nrest = ceil(nctot/7);  

title(parfig.title);
legend(legendtext((nc-nrest):nc,:));
ylabel(parfig.ylabel,'FontSize',11,'FontWeight', 'normal');
xlabel(['time from ' parfig.center_evt ' start (ms)'], 'FontSize',11,'FontWeight', 'normal');



