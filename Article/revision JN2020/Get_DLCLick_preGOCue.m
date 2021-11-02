% Plot_DLCLick_vs_TTLSensor_align
% JC Script 10/02/2020
% cd (DATA_folder)
% load('Ntrial_type.mat')
% load ('evt.mat','evt_trial', 'evt_lick_L', 'evt_lick_R'); load ('time.mat');

clear all
close all
NTRtot = 0; 
PERC_ttl=[]; PERC_go=[]; TTLminusDLC = [];  TTLminusGOcue =[]; meanFrame_fromLick1=[]; meanFrame_fromGoCue=[]; TTLminusGOcue2 =[] 
for nm=1:5
    if nm==1
        FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d8_z4900_VM_taskopto_optopost_G912_180715_vidY_150tr_20cel_13mW_tooDeep';
        MouseID =  'vgat14w14d8';
    elseif nm==7
        FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW';
        MouseID =  'vgat17w10d7';
    elseif nm==3
        FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW';
        MouseID =  'vgat15w10d7';
    elseif nm==4
        FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w08d5_z4700_VM_taskopto_optopost_G912_180706_vidY_300tr_18cel_12mW';
        MouseID =  'vgat15w08d5';
    elseif nm==6
        FileLocation = 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d2_z4670_VM_taskopto_optopost_G912_180709_vidY_100tr_29cel_10mW';
        MouseID =  'vgat14w14d2';
    elseif nm==2
        FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d5_z4300_VM_taskopto_optopost_EOBD_180728_vidY_200tr_20cel_12mW_Q4';
        MouseID =  'vgat17w10d5';
    elseif nm==5
        FileLocation= 'D:\DATA EMORY\JC_Analysis\JCVGAT11\vgat11_w10d4_z4300_VM_taskonly_nonepost_CC4F_180505_vidY_150tr_40cel_00mW';
        MouseID =  'vgat11w10d4';
    end

    cd(FileLocation);
    load([ MouseID '_DLCresults.mat'])
    load('Ntrial_type.mat')
    load ('evt.mat','evt_trial', 'evt_lick_L', 'evt_lick_R', 'evt_delay')
    load ('time.mat');
    
    %% define trial start and end
    trig_end = find(diff(evt_trial)<0);
    trig_st = find(diff(evt_trial)>0);
    TTL_lenght_sec=trig_end(2)/20000 - trig_st(2)/20000
    %
    idx_trial_start = trig_st(1:end-1);
    idx_trial_end = idx_trial_start + 4*20000;
    evt_lick = evt_lick_L + evt_lick_R;
    
    %% average impulsive alignement
    plotON = 0;
    K=2
    trtype = 'cor'
    % for tt=1:2
    %     if tt==1
    %         trtype = 'cor'
    %     else
    %         trtype = 'imp'
    %     end
    
 idxtr =[]; 
    if trtype=='imp'
        idxtr = sort([trial.idx_errorDelay_PL_CL,...
            trial.idx_errorDelay_PR_CR,...
            trial.idx_errorDelay_PL_CR,...
            trial.idx_errorDelay_PR_CL]);
    elseif trtype=='cor'
        idxtr = sort([trial.idx_correct_R trial.idx_correct_L]);
    elseif trtype=='all'
        idxtr =[trial.idx_correct_R trial.idx_correct_L trial.idx_errorDelay_PL_CL,...
            trial.idx_errorDelay_PR_CR trial.idx_errorDelay_PL_CR trial.idx_errorDelay_PR_CL];
    end
    
    if MouseID(5:6)=='11'; idxtr(idxtr>150)=[]; end 
    ntrial = max(size(idxtr)); 
    NTRtot = NTRtot + ntrial; 
    
    TTL =[]; DLC = []; YDLC =[];BaseLine =[]; THR=[]; Y=[]; DLC_lick1frame = []; 
    TTLminusDLC =  []; TTLminusGOcue =  [];  TTLminusGOcue2 =  []; 
    TTL_lick1frame=[];  TTL_GOframe= [];  TTLGOframe_all =[]; TTLLick1frame_all = [];
    
    for ii= 1 : ntrial;
        tr=idxtr(ii);
        idx_tr = [idx_trial_start(tr)+1:1:idx_trial_end(tr)];
        Xt =time(idx_tr); Xt(1);
        XtLick1st(ii)= Xt(min(find(evt_lick(idx_tr)==1)));
        XtGoCue(ii)= Xt(max(find(evt_delay(idx_tr)==1)));

        TTL_lick1frame= (XtLick1st(ii)-Xt(1))*25
        TTL_GOframe= (XtGoCue(ii)-Xt(1))*25
        
        YDLC = allTab.TongueY(((tr-1)*100)+1:tr*100);
        BaseLine = median(abs(YDLC(1:30)));
        Y  = abs(YDLC-BaseLine);
        THR = median(abs(Y(1:30)))+K*std(abs(Y)',1);
        DLC_lick1frame = 29+min(find(Y(30:end)>THR));
        
        if plotON==1 & ii<10
            figure, hold on,
            plot((Xt-Xt(1))*25, evt_lick(idx_tr)*10,'k');
            plot(Y');
            yline(THR);
            
            legend('TTL LickSensor','DLC TongueY','DLC TongueX'  );
            title(['SENSOR vs DLC trial' num2str(tr)]);
            xlabel('frames (40ms/frame)');
            ylabel('Y (pixels)');
            %         ylim([80 180])
            %         xlim([1 70])
            %         saveas(gcf,[SAVEFIG_folder 'LickSensor vs DLC ' trtype ' trial#' num2str(tr) '_' MouseID],'png')
        end
        
        TTL = [TTL evt_lick(idx_tr)];
        DLC = [DLC allTab.TongueY(((tr-1)*100)+1:tr*100)];

        TTLminusDLC(ii) =  DLC_lick1frame - TTL_lick1frame;
        TTLminusGOcue(ii) =  DLC_lick1frame - 37.5;
        TTLminusGOcue2(ii) =  DLC_lick1frame - TTL_GOframe;

        TTLGOframe_all(ii) =TTL_GOframe;
        TTLLick1frame_all(ii) = TTL_lick1frame ;
    end
    
    figure,
    histogram(TTLminusDLC,-5:1:5)
    ntrial = max(size(TTLminusDLC))
    percEarly_ttl = 100*sum(TTLminusDLC<0)/ntrial
    title(['TTL Lick1 ' trtype ' ntrial=' num2str(ntrial) '  percEarly=' num2str(percEarly_ttl) '% std=' num2str(K)  ' ' MouseID  ])
    saveas(gcf, ['D:\DATA EMORY\JC_Analysis\deeplabcut_fig\ttl-dlc dist\TTLminusDLC' MouseID '_' trtype '-K' num2str(K)],'png')
   
    figure,
    histogram(TTLminusGOcue	,-5:1:5)
    ntrial = max(size(TTLminusGOcue	))
    percEarly_go = 100*sum(TTLminusGOcue<0)/ntrial
    title(['GOcue ' trtype ' ntrial=' num2str(ntrial) '  percEarly=' num2str(percEarly_go) '% std=' num2str(K)  ' ' MouseID  ])
    saveas(gcf, ['D:\DATA EMORY\JC_Analysis\deeplabcut_fig\ttl-dlc dist\TTLminusGOCUE' MouseID '_' trtype '-K' num2str(K)],'png')
    
    figure,
    histogram(TTLminusGOcue2	,-5:1:5)
    ntrial = max(size(TTLminusGOcue2	))
    percEarly_go2 = 100*sum(TTLminusGOcue2<0)/ntrial
    title(['GOcue2 ' trtype ' ntrial=' num2str(ntrial) '  percEarly=' num2str(percEarly_go) '% std=' num2str(K)  ' ' MouseID  ])
    saveas(gcf, ['D:\DATA EMORY\JC_Analysis\deeplabcut_fig\ttl-dlc dist\TTLminusGOCUE' MouseID '_' trtype '-K' num2str(K)],'png')
    
    
    
    PERC_ttl(nm)=   percEarly_ttl
    PERC_go(nm)=   percEarly_go
    PERC_go2(nm)=   percEarly_go2
    meanFrame_fromLick1(nm)= mean(TTLminusDLC(TTLminusDLC<0)) 
    meanFrame_fromGoCue(nm)=mean(TTLminusGOcue(TTLminusGOcue<0)) 
    meanFrame_fromGoCue2(nm)=mean(TTLminusGOcue2(TTLminusGOcue2<0))     
    
end

%%
clc
STD_PERC_TTL = std( PERC_ttl)
STD_PERC_GO = std( PERC_go)
STD_FRAME_2TTL = std(meanFrame_fromLick1)
STD_FRAME_2GO = std(meanFrame_fromGoCue)

MEAN_PERC_TTL = mean( PERC_ttl)
MEAN_PERC_GO = mean( PERC_go)
MEAN_PERC_GO2 = mean( PERC_go2)
MEAN_FRAME_2TTL = mean(meanFrame_fromLick1)
MEAN_FRAME_2GO = mean(meanFrame_fromGoCue)
MEAN_FRAME_2GO2 = mean(meanFrame_fromGoCue2)

disp(['#' num2str(nm) 'sessions ' 'ntrial=' num2str(NTRtot) ' trialtype ' trtype '  K=' num2str(K) ])
disp(['Nframe before TTL ' num2str(40*MEAN_FRAME_2TTL) 'ms +/-' num2str(40*STD_FRAME_2TTL) ' occuring in ' num2str(MEAN_PERC_TTL) '+/-' num2str(STD_PERC_TTL) '% of the trials'])
disp(['Nframe before Go-cue ' num2str(40*MEAN_FRAME_2GO) 'ms +/-' num2str(40*STD_FRAME_2GO) ' occuring in ' num2str(MEAN_PERC_GO) '+/-' num2str(STD_PERC_GO) '% of the trials'])
disp(['Nframe before Go-cue2 ' num2str(40*MEAN_FRAME_2GO2) 'ms +/-' num2str(40*STD_FRAME_2GO) ' occuring in ' num2str(MEAN_PERC_GO2) '+/-' num2str(STD_PERC_GO) '% of the trials'])

%% plot average from start trial
figure, hold on,
plot((Xt-Xt(1))*25,  mean(TTL,2)*-100,'k')
plot(mean(DLC,2)- mean(median(DLC,2)),'r')
title(['SENSOR vs DLC trial Average n=' num2str(max(size(idxtr)))])
xlabel('frames (40ms/frame)')
ylabel('Y (pixels)')


