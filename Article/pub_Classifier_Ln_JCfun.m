function [Ln] = pub_Classifier_Ln_JCfun(ncell, nspx_A, nspx_B, parfig)
% [Ln] = pub_Classifier_Ln_JCfun(ncell, nspx_A, nspx_B, parfig)
% nspx_A = nb spike in trial type A
% new classifier
% JC 4/23/2019

rng('shuffle');
NA= size(nspx_A,1);
NB= size(nspx_B,1);
Mintr=min(NA,NB);
Ntr = Mintr*2;

if parfig.ControlShuffle == 0
    Y=ones(Ntr,1);
    Y(1:Mintr)=0; % trials Label 1 or 0
elseif parfig.ControlShuffle == 1
    A = randperm(Ntr);
    Y = A' >(Ntr/2);
end


X0 = [nspx_A(1:Mintr,:); nspx_B(1:Mintr,:)]; %Nspk for each trials and each cells (X0 = ntrials x ncells)

ordered_tr =[X0 Y];

Ln=[];
for icell=1:ncell
    RandCell=randperm(ncell);
    X=X0(:,RandCell(randperm(icell)));
    n = size(X,1);
    
    % Create a ClassificationLinear object by using fitclinear.
    [Mdl,FitInfo] = fitclinear(X,Y,'Learner', parfig.learner, 'Lambda', 1/2); %'Regularization', 'ridge' or 'lasso');
    
    if parfig.typeClass== '%HoldOut'
        cvp = cvpartition(n,'Holdout', parfig.Nfold/100);
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
        labels = predict(Mdl,X2(:,idxTest),'ObservationsIn','columns');
        L = loss(Mdl,X2(:,idxTest),Y(idxTest),'ObservationsIn','columns');
        Ln=[Ln L];
        
    elseif parfig.typeClass== 'FoldXVal'
        cvp = cvpartition(n,'KFold', parfig.Nfold);
        %           K-fold cross validation partition
        %    NumObservations: 74
        %        NumTestSets: 10
        %          TrainSize: 67  66  66  66  66  67  67  67  67  67
        %           TestSize: 7  8  8  8  8  7  7  7  7  7
        %
        for incv=1:parfig.Nfold
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
