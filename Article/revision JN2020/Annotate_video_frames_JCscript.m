% Annotate_video_frames_JCscript
close all, clear all, clc
cd ('D:\DATA EMORY\JC_Analysis')
% load a full raw video

%% load evts and trials

% VGAT14-w14d2
% DATA_folder = 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d2_z4670_VM_taskopto_optopost_G912_180709_vidY_100tr_29cel_10mW'
% MouseIDfull =  'vgat14w14d2'
% MouseID = 'Mouse#3'
% cd (DATA_folder)
% load('vid14raw.mat'); vid=vid14raw; clear vid14raw;


% VGAT14-w14d8
% DATA_folder = 'D:\DATA EMORY\JC_Analysis\JCVGAT14\vgat14_w14d8_z4900_VM_taskopto_optopost_G912_180715_vidY_150tr_20cel_13mW_tooDeep'
% MouseIDfull =  'vgat14w14d8'
% MouseID = 'Mouse#3'
% cd (DATA_folder)
% load('vid14raw.mat'); vid=vid14raw; clear vid14raw;


% VGAT15-w10d7
% DATA_folder = 'D:\DATA EMORY\JC_Analysis\JCVGAT15\vgat15_w10d7_z4400_VM_taskopto_optopost_G912_180721_vidY_250tr_24cel_11mW'
% MouseIDfull =  'vgat15w10d7'
% MouseID = 'Mouse#4'
% cd (DATA_folder)
% load('vid15raw.mat'); vid=vid15raw; clear vid15raw;

% VGAT17-w10d7
DATA_folder = 'D:\DATA EMORY\JC_Analysis\JCVGAT17\vgat17_w10d7_z4300_VM_taskopto_optopost_G912_180729_vidY_120tr_38cel_11mW'
MouseIDfull =  'vgat17w10d7'
MouseID = 'Mouse#5'
cd (DATA_folder)
load('vid17raw.mat'); vid=vid17raw; clear vid17raw;

% Load trial info
load('Ntrial_type.mat')

%% edit each trial for Number and Type
MANUAL = 'ON'

trtype = {'correct' ;'impulsive'; 'omission' } %;'opto stim'}

for tt=size(trtype):-1:1
    TRIAL_TYPE_ID = trtype{tt}
    mkdir(['.\VID_annotated\' TRIAL_TYPE_ID(1:3)])
    
    if TRIAL_TYPE_ID(1:3) == 'cor'        
        if MANUAL == 'ON';
            idxtrial = sort([4 65 117 147 196]);
        else
            idxtrial = sort([trial.idx_errorDelay_PL_CL trial.idx_errorDelay_PL_CR trial.idx_errorDelay_PR_CL trial.idx_errorDelay_PR_CR])
        end 
        
    elseif  TRIAL_TYPE_ID(1:3) == 'imp'       
        if MANUAL == 'ON';
            idxtrial = sort([15 145]);
        else
            idxtrial = sort([trial.idx_correct_L trial.idx_correct_R])
        end
        
    elseif  TRIAL_TYPE_ID(1:3) == 'omi'  
        if MANUAL == 'ON';
            idxtrial = sort([44 91]);
        else
            idxtrial = sort([trial.idx_NoLick])
        end
        
    elseif  TRIAL_TYPE_ID(1:3) == 'opt'
        idxtrial = sort([trial.idx_NoLick_opto trial.idx_correct_L_opto trial.idx_correct_R_opto])
    end
    
    ntrial = size(idxtrial,2)
    nframe = size(vid,2)
    
    % set to 10 trials example vid
    if  MANUAL ~= 'ON'
        idxtrial = idxtrial(1:floor(ntrial/30):end)
    end
    %     idxtrial = [48 86 106 141 165]
    ntrial = size(idxtrial,2)
    
    % loop
    close all
    for ii=nframe:-1:1;
        for jj=ntrial:-1:1
            if ii <= idxtrial(jj)*100 & ii > (idxtrial(jj)*100)-100
                if ii == idxtrial(jj)*100
                    TT=4040;
                end
                TT=TT-40;
                
                if TT<750
                    EVTID = 'Air Puff epoch';
                elseif TT>750 & TT<1500
                    EVTID = 'DELAY epoch';
                elseif TT>1500 & TT<1600
                    EVTID = 'GO CUE';
                elseif TT>1600
                    EVTID = 'Response epoch';
                end
                
                if  TRIAL_TYPE_ID(1:3) == 'opt'
                    if TT>1000 & TT<2000;
                        OPTO = 'Blue Laser ON'
                    else
                        OPTO = ''
                    end
                end
                
                figure(ii),
                imagesc(vid(ii).cdata);
%                 PositionTXT = [0.13 0.69 0.24 0.24]
                A= annotation(gcf,'textbox',...
                    [0.13+0.26 0.69 0.24 0.24],...% [X,Y,Widht, height ]
                    'String',['trial #' num2str(idxtrial(jj)) ' (' TRIAL_TYPE_ID ')'] ,...
                    'Color',[1 0.7 0.7],...
                    'FontSize',13,...
                    'FontName','Arial',...
                    'LineStyle','none',...
                    'FitBoxToText','on');
                
                B= annotation(gcf,'textbox',...
                    [0.13+0.26 0.62 0.24 0.24],...% [X,Y,Widht, height ]
                    'String',['time= ' num2str(TT)  'ms  ' ],...
                    'Color',[0.7 1 0.7],...
                    'FontSize',13,...
                    'FontName','Arial',...
                    'LineStyle','none',...
                    'FitBoxToText','on');
          
                C= annotation(gcf,'textbox',...
                    [0.13+0.26 0.55 0.24 0.24],...% [X,Y,Widht, height ]
                    'String',[ EVTID ' '] ,...
                    'Color',[0.7 0.7 1],...
                    'FontSize',13,...
                    'FontName','Arial',...
                    'EdgeColor',[1 1 1],...
                    'LineStyle','none',...
                    'FitBoxToText','on');
                
                
                D= annotation(gcf,'textbox',...
                    [0.13+0.26 0.025 0.14 0.14],...% [X,Y,Widht, height ]
                    'String',[ MouseID ' '] ,...
                    'Color',[0.9 0.9 0.9],...
                    'FontSize',11,...
                    'FontName','Arial',...
                    'EdgeColor',[1 1 1],...
                    'LineStyle','none',...
                    'FitBoxToText','on');
                
                if  TRIAL_TYPE_ID(1:3) == 'opt'
                    E= annotation(gcf,'textbox',...
                        [0.13+0.26 0.25 0.14 0.14],...% [X,Y,Widht, height ]
                        'String',[ OPTO ' '] ,...
                        'Color','c',...
                        'FontSize',14,...
                        'FontName','Arial',...
                        'EdgeColor',[1 1 1],...
                        'LineStyle','none',...
                        'FitBoxToText','on');
                end
                
                
                
                
                pause(0.01);
                
                if MANUAL == 'ON'
                    mkdir('.\VID_annotated\Article\')
                    saveas(gcf,['.\VID_annotated\Article\frame' num2str(ii) '_' TRIAL_TYPE_ID(1:3)],'tif');
                else
                    saveas(gcf,['.\VID_annotated\' TRIAL_TYPE_ID(1:3) '\frame' num2str(ii)],'tif');
                    %                 if TRIAL_TYPE_ID(1:3) == 'opt'
                    %                     saveas(gcf,['.\VID_annotated\' TRIAL_TYPE_ID(1:3) '\opt-5tr\frame' num2str(ii)],'tif');
                    %                 end
                end
                pause(0.01);
                close (ii);
                
            end
        end
    end
end