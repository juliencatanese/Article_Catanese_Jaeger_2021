%% Curr Biol Lick article
% question are all lick the sames or goal directed differ from regular?
% compare Lick1 vs Other licks
%

% Choose 1 Session
cd C:\Users\catan\Documents\EMORY\JC_Analysis\JCVGAT17\vgat17_w10d8_z4500_VM_taskopto_optopost_G912_180730_vidY_250tr_43cel_13mW

% Choose 1 cell
Sh=1, Ch=1; CLUST=3;

% Choose which Lick
for Lii = 1:5
    LickNb = Lii;
    LickID = ['Lick' num2str(LickNb)];
    
    % Define time (ms) around Lick
    pre=250 %ms
    post=500 %ms
    
    % plot raster and PSTH center on lick
    Raster_1cell_lick12345(Sh, Ch, CLUST, LickID , {'cor'}, {'k'}, pre, post, LickID)
    
    % save
    ChanID = ['S' num2str(Sh) 'Ch' num2str(Ch) ];
    
    saveas(gcf, ['./Licks/fig_LickSpk_'  ChanID '_clust#' num2str(CLUST) LickID '_' DirSess(ii).name(1:12) ], 'jpg')
end
end
