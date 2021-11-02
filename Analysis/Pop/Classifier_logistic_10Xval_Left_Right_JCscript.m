% Classifier_logistic_10Xval_Left_Right_JCscript
% Julien Catanese 10/31/2018

close all,

tr_type_comp1 = {'cCL', 'cCR'};

Comp_list=[{tr_type_comp1}];
trial_type_list =[];
trial_type_list = Comp_list{1}


for Ev=1:3
    if Ev==1
        psth_center_evtID = 'APuff'
    elseif Ev==2
        psth_center_evtID = 'Delay'
    elseif Ev==3
        psth_center_evtID = 'GoCue'
    end
    
    % Other_trial_type = []
    clearvars -except Trtyp Ev psth_center_evtID trial_type_list Comp_list FolderID MouseID
    rng('shuffle');
    
    %%
    typeClass='FoldXVal'  %     typeClass='%HoldOut'
    learner= {'logistic'}%, 'svm'}
    
    %% Define Center Evt
    disp(['psth Center = ' psth_center_evtID]);
    disp('NO Other trial type, plotting LvR only');
    
    psth_trig_evt =  psth_center_evtID; % 'Delay', % 'APuff' % 'GoCue' % 'Valve' %'Licks'
    NbTtrialType = max(size(trial_type_list));
    %% load data
    load('info.mat');
    load('evt.mat');
    load('time.mat');
    load('Ntrial_type.mat');
    sr=info.info_freq_parameters.board_dig_in_sample_rate;
    
    %% GET trigtimes : times of the psth_trig_evt is sec   (1 vector of times)
    % define trial idx
    trig_end = find(diff(evt_trial)<0);
    trig_st = find(diff(evt_trial)>0);
    
    idx_trial_start = trig_end(1:end-1); %-(1*sr);
    idx_trial_end = trig_st(2:end) ;     %-(1*sr);
    
    Ntrials=trial.Ntrial
    
    % define common events
    evt_valve = evt_rwd_L + evt_rwd_R;
    evt_lick = evt_lick_L + evt_lick_R;
    evt_puff =evt_puff_L + evt_puff_R;
    evt_prelick = evt_delay + evt_puff;
    
    
    trigtimes=[];
    for tr=1:Ntrials
        time_tr = time(idx_trial_start(tr):idx_trial_end(tr));
        if psth_trig_evt=='APuff'
            puff_tr=  evt_puff(idx_trial_start(tr):idx_trial_end(tr));
            idx_puff_st = find(diff(puff_tr)>0); % start of the puff
            time_puff_st = time_tr(idx_puff_st);
            trigtimes = [trigtimes time_puff_st]; % in sec
            
        elseif psth_trig_evt=='Delay'
            delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
            idx_delay_st = find(diff(delay_tr)>0) ; % start of the delay
            time_delay_st = time_tr(idx_delay_st);
            trigtimes = [trigtimes time_delay_st]; % in sec
            
        elseif psth_trig_evt=='GoCue'
            delay_tr= evt_delay(idx_trial_start(tr):idx_trial_end(tr));
            idx_GO_st = find(diff(delay_tr)<0) ; % end of the delay
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
            
            disp(' WRONG INPUT ARG pick one: psth_trig_evt = APuff, Delay, Licks, Gocue or Valve' )
            
        end
        
    end
    
    %%
    trigtimes_all=trigtimes';
    trigtimes_cCL = trigtimes(trial.idx_correct_L);
    trigtimes_cCR = trigtimes(trial.idx_correct_R);
    trigtimes_iCL = trigtimes([trial.idx_errorDelay_PL_CL trial.idx_errorDelay_PR_CL]);
    trigtimes_iCR = trigtimes([trial.idx_errorDelay_PL_CR trial.idx_errorDelay_PR_CR]);
    trigtimes_oCR = trigtimes([trial.idx_correct_R_opto]);
    trigtimes_oCL = trigtimes([trial.idx_correct_L_opto]);
    trigtimes_ocC = trigtimes([trial.idx_correct_L_opto trial.idx_correct_R_opto]);
    trigtimes_oNO = trigtimes([trial.idx_NoLick_opto]);
    trigtimes_eNO = trigtimes([trial.idx_NoLick]);
    trigtimes_ePL = trigtimes([trial.idx_errorResp_PL]);
    trigtimes_ePR = trigtimes([trial.idx_errorResp_PR]);
    
    trigtimes_cor = sort([trigtimes_cCR trigtimes_cCL]);
    trigtimes_imp = sort([trigtimes_iCL trigtimes_iCR]);
    trigtimes_opt = sort([trigtimes_oNO trigtimes_ocC]);
    
    %% GET spxtimes = SPIKE times in Sec (1 vector of times)
    clust_file = dir('times_*S*Ch*_sub.mat');
    all_nspx_A = []; all_nspx_B = [];
    ncell=0;
    for ncluf = 1:max(size(clust_file)) %  channel loop
        load(clust_file(ncluf).name); %  e.g.: load('times_S2Ch6_man.mat')
        Nclust = max(cluster_class(:,1));
        
        for CLUST= 1:Nclust %Cluster loop within each channel
            idx_spk = find(cluster_class(:,1)==CLUST);
            spxtimes = cluster_class(idx_spk,2)/10^3;
            ncell=ncell+1;
            for ii=1:NbTtrialType
                if trial_type_list{ii} == 'cCL'
                    trigtimes = trigtimes_cCL ;
                elseif trial_type_list{ii} == 'cCR'
                    trigtimes = trigtimes_cCR ; %disp('Correct Left trials');
                elseif trial_type_list{ii} == 'eNO'
                    trigtimes = trigtimes_eNO ; %disp('Correct Right trials');
                end
                
                % DELAY : Count Spikes within Delay Epoch for each trial ([delaystart+0] to [delaystart+750ms] )
                if ii==1
                    nspx_A = mnspx(spxtimes,trigtimes,0,750);
                    all_nspx_A = [all_nspx_A nspx_A];
                elseif ii==2
                    nspx_B = mnspx(spxtimes,trigtimes,0,750);
                    all_nspx_B = [all_nspx_B nspx_B];
                end
            end
        end
    end
    %%
    Nrepet=20;
    Lnall=[];
    
    for ll=1:max(size(learner))
        
        for Nr=1:Nrepet
            
            NA= size(all_nspx_A,1)
            NB= size(all_nspx_B,1)
            
            Mintr=min(NA,NB)
            %             Ntr = NL+NR;
            Ntr = Mintr*2;
            Y=ones(Ntr,1);
            Y(1:Mintr)=0; % trials Label 1 or 0
            X0 = [all_nspx_A(1:Mintr,:); all_nspx_B(1:Mintr,:)]; %Nspk for each trials
            
            ordered_tr =[X0 Y];
            
            Ln=[];
            for icell=1:ncell
                RandCell=randperm(ncell);
                X=X0(:,RandCell(randperm(icell)));
                
                % Create a ClassificationLinear object by using fitclinear.
                [Mdl,FitInfo] = fitclinear(X,Y,'Learner', learner{ll}, 'Lambda', 1/2);
                %                 [Mdl,FitInfo] = fitclinear(X,Y,'Learner', learner{ll}, 'Lambda', 'auto', 'Regularization', 'ridge');
                %                 [Mdl,FitInfo] = fitclinear(X,Y,'Learner', learner{ll}, 'Lambda', 'auto', 'Regularization', 'lasso');
                
                
                rng('shuffle'); %random function
                n = size(X,1);
                Nfold_cv = 10;
                if typeClass== '%HoldOut'
                    cvp = cvpartition(n,'Holdout',Nfold_cv/100);
                    % cvp is a CVPartition object that defines the random partition of n data into training and test sets.
                    % Hold-out cross validation partition
                    %    NumObservations: 31572
                    %        NumTestSets: 1
                    %          TrainSize: 29994
                    %           TestSize: 1578
                    
                    idxTrain = training(cvp); % Extract training set indices
                    X2 = X';
                    Mdl = fitclinear(X2(:,idxTrain),Y(idxTrain),'ObservationsIn','columns');
                    % Predict observations and classification error for the hold out sample.
                    
                    idxTest = test(cvp); % Extract test set indices
                    labels = predict(Mdl,X2(:,idxTest),'ObservationsIn','columns')
                    L = loss(Mdl,X2(:,idxTest),Y(idxTest),'ObservationsIn','columns');
                    Ln=[Ln L];
                elseif typeClass== 'FoldXVal'
                    cvp = cvpartition(n,'KFold',Nfold_cv);
                    %           K-fold cross validation partition
                    %    NumObservations: 74
                    %        NumTestSets: 10
                    %          TrainSize: 67  66  66  66  66  67  67  67  67  67
                    %           TestSize: 7  8  8  8  8  7  7  7  7  7
                    %
                    for incv=1:Nfold_cv
                        idxTrain = training(cvp,incv); % Extract training set indices
                        X2 = X';
                        Mdl = fitclinear(X2(:,idxTrain),Y(idxTrain),'ObservationsIn','columns');
                        % Predict observations and classification error for the hold out sample.
                        
                        idxTest = test(cvp,incv); % Extract test set indices
                        labels = predict(Mdl,X2(:,idxTest),'ObservationsIn','columns');
                        L = loss(Mdl,X2(:,idxTest),Y(idxTest),'ObservationsIn','columns');
                        
                        LnKF(incv)=L;
                    end
                    Ln=[Ln mean(LnKF)];
                end
            end
            %                 figure(1),
            %                 hold on,
            %                 NCellAxis=[1:1:size(Ln,2)];
            %                 plot(NCellAxis, Ln);
            Lnall=[Lnall;Ln];
        end
        %%
        Mouse=info.info_notes.MouseID;
        Day=info.info_notes.Day;
        M=mean(Lnall);
        S=std(Lnall);
        
        figure(2), hold on,
        subplot(1,3,Ev), hold on,
        
        
        
        plot( [1:1:ncell] , ones(1,ncell)/2, 'k--')
        col = [0.5 0.5 0.5]
        plot( [1:1:ncell] , M, 'Color',col, 'LineWidth', 2)
        M1=M; save(['ClassM1_' psth_center_evtID '_' trial_type_list{1} '_' trial_type_list{2}],'M1','trial_type_list')
        
    end
    %%
    
    ylim([0 1])
    xlabel('#Neurons', 'FontSize', 11)
    ylabel('Proba errors', 'FontSize', 11)
    title(psth_center_evtID, 'FontWeight','bold' ,'FontSize', 12)
end

%%
legend(['chance level  (Mouse: ' Mouse ' ' Day ')'], ' Left-Right Mean LogClass10X Nrepet20')
saveas(gcf, [Mouse '_' Day '_' trial_type_list{1} '_vs_' trial_type_list{2} '_Log_Lin_Classif_' num2str(Nfold_cv) 'FoldCV'], 'png' )
saveas(gcf, ['D:\JC_Figures\' Mouse '_' Day '_' trial_type_list{1} '_vs_' trial_type_list{2} '_Log_Lin_Classif_' num2str(Nfold_cv) 'FoldCV'], 'png' )
disp('saved')