% pub_table2_coordinate_AP_ML_DV_JCscript
% create and save a table ('Tcoord') that contain all cells informations
% Fields: ncell, nMouse, nDay, nChan, nShank, AP, ML, DV, VM, VL.
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tcoord')
% written by Julien Catanese 12/05/2018
% last updated JC 2/7/2019.

clear all;
AP_measured=[]; ML_measured=[]; DV_measured=[];
AP_corr=[]; ML_corr=[]; DV_corr=[]; C=[]; B=[]; A=[];

cd('D:\JC_Analysis');
load('listcell.mat');  T1= listcell;

%% Add nSess (Session number from 1:nSess instead of the string Day)
listSess=unique(T1.SessID)
nSess = max(size(listSess))
for ii= 1:nSess
    A(:,ii)= strcmp(T1.SessID, listSess{ii});
end
B=A.*[1:1:nSess];
C=sum(B')'; % vector with Session Number from 1 to 11.
T1= addvars(T1, C, 'Before', 2, 'NewVariableNames', {'nSess'});
T1(1,:)

%% Add nMouse (Mice number from 1:nMouse instead of the string MouseID)
ncell=T1.ncell(end)
C=[]; B=[]; A=[];
for ii=1:ncell;
    MID{ii}= T1.MouseID(ii,:);
end
listMice=unique(MID')
nMice = max(size(listMice))
for ii= 1:nMice
    A(:,ii)= strcmp(MID, listMice{ii})
end
B=A.*[1:1:nMice];
C=sum(B')';% vector with Session Number from 1 to 11.
T1= addvars(T1, C, 'Before', 3, 'NewVariableNames', {'nMouse'});
T1(1:10,:)

%% Manual Attribution of Coordinate (based on Stereotaxic Measurement and Histo reconstruction and Ephys Markers)
% 'Sh1M_' = The coordinate in the table corespond to channel #8 (always tip) on Shank #1 (the most Medial M_)

SessionNAME = ['11_w10d4'; '12_w11d5'; '14_w14d2'; '14_w14d5'; '14_w14d7'; '14_w14d8'; '15_w08d5'; '15_w10d3'; '17_w10d3'; '17_w10d4'; '17_w10d7']
Ch_Origin   = ['S1Ch8M_' ; 'S2Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8MP' ; 'S1Ch8M_' ; 'S2Ch8_L' ; 'S1Ch8M_' ];
Ch_Orient   = ['ML'      ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ; 'DG'      ; 'ML'      ;  'ML'     ; 'ML'      ];
AP_measured = [1.30      ;  1.58     ;  1.10     ;  1.55     ;  1.50     ;  1.55     ;  1.90     ;   1.50    ; 1.20      ; 1.00      ; 1.45      ];
ML_measured = [0.30      ;  0.90     ;  0.75     ;  0.90     ;  0.80     ;  0.85     ;  0.70     ;   0.65    ; 0.50      ; 1.10      ; 0.65      ];
DV_measured = [4.30      ;  4.30     ;  4.65     ;  4.30     ;  4.35     ;  4.90     ;  4.70     ;   4.25    ; 4.30      ; 4.20      ; 4.30      ];

APAdj       = [+0.5      ;  +0.0     ;  -0.0     ;  +0.0     ;  +0.0     ;  -0.0     ;  -0.0     ;   +0.0    ; +0.0      ; +0.0      ; +0.0      ];
MLAdj       = [+0.1      ;  +0.0     ;  -0.0     ;  +0.0     ;  +0.0     ;  -0.0     ;  -0.0     ;   +0.0    ; +0.1      ; +0.0      ; +0.0      ];
DVAdj       = [+0.0      ;  +0.0     ;  -0.0     ;  +0.0     ;  +0.0     ;  -0.0     ;  -0.0     ;   +0.0    ; +0.0      ; +0.0      ; +0.0      ];

AP_corr = AP_measured + APAdj;
ML_corr = ML_measured + MLAdj;
DV_corr = DV_measured + DVAdj ;

Tst= table( Ch_Origin, Ch_Orient, AP_corr, ML_corr, DV_corr, listSess);

%%

T2=[];
T2=[T1(:,7) T1(:,1:3) ]; T2(1,:)
nShank = str2num(T1.ChanID(:,2));
nChan  = str2num(T1.ChanID(:,5));

T2=addvars(T2, nShank, nChan);
T2(1:5,:)
%%
% StereoT.Sh_Origin
% T2(1:5,:)

AP=[]; ML=[];DV=[];
for ii=1:ncell
    %   coln =
    Orig=Tst.Ch_Origin(T1.nSess(ii),:)
    ShOrig=str2num(Orig(2))
    ChOrig=str2num(Orig(5))
    Orient=Orig(6:7);
    
    if Orient=='M_'
        AP(ii) = Tst.AP_corr(T1.nSess(ii))
        ML(ii) = Tst.ML_corr(T1.nSess(ii))-(0.2*(ShOrig-T2.nShank(ii)))
        DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)))
        
    elseif Orient=='_L'
        AP(ii) = Tst.AP_corr(T1.nSess(ii))
        ML(ii) = Tst.ML_corr(T1.nSess(ii))+(0.2*(ShOrig-T2.nShank(ii)))
        DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)))
        
    elseif Orient=='_P'
        AP(ii) = Tst.AP_corr(T1.nSess(ii))+(0.2*(ShOrig-T2.nShank(ii)))
        ML(ii) = Tst.ML_corr(T1.nSess(ii))
        DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)))
        
    elseif Orient=='MP'
        AP(ii) = Tst.AP_corr(T1.nSess(ii))+(0.14*(ShOrig-T2.nShank(ii)))
        ML(ii) = Tst.ML_corr(T1.nSess(ii))-(0.14*(ShOrig-T2.nShank(ii)))
        DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)))
        
    end
    
end
%% Add Coordinates to Table 
AP=-AP'; ML=-ML'; DV=-DV';
T2=addvars(T2, AP, ML, DV); T2(T2.nMouse==1,:)

%% Select the Axis with less variability 
NML=max(size(unique(round(T2.ML,1))))
NAP=max(size(unique(round(T2.AP,1))))
if NAP<NML
    display('NAP<NML')
    AP_slice =unique(round(T2.AP,1))
else
    display('NML<NAP')
    ML_slice = unique(round(T2.ML,1))
end

%% Manual definition of the limite of VM and VL/VA in the Mouse Atlas (Paxinos)
% Tlim = table(ML_slice)
ML_Slice  = [-0.4 ; -0.6 ; -0.7 ; -0.8 ; -0.9 ; -1.0 ; -1.1 ; -1.2 ; -1.3 ]
APmax_VM    = [-1.7 ; -1.0 ; -1.0 ; -1.0 ; -1.0 ; -1.0 ; -1.1 ; -1.2 ; -1.3 ]
APmin_VM    = [-2.3 ; -2.2 ; -2.2 ; -2.0 ; -2.2 ; -2.2 ; -2.2 ; -2.0 ; -1.7 ]
Zmax_VM     = [-4.1 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.9 ]
Zmin_VM     = [-4.3 ; -4.5 ; -4.5 ; -4.5 ; -4.4 ; -4.4 ; -4.3 ; -4.2 ; -4.1 ]
Tlim=table(ML_Slice, APmax_VM, APmin_VM, Zmax_VM, Zmin_VM)

ML_Slic_VL  = [-0.4 ; -0.6 ; -0.7 ; -0.8 ; -0.9 ; -1.0 ; -1.1 ; -1.2 ; -1.3 ]
APmax_VL    = [-0   ; -1.0 ; -1.2 ; -1.0 ; -1.0 ; -1.0 ; -1.0 ; -1.0 ; -1.0 ]
APmin_VL    = [-0   ; -1.5 ; -2.0 ; -1.9 ; -1.9 ; -1.8 ; -1.8 ; -1.5 ; -1.4 ]
Zmax_VL     = [-0   ; -3.6 ; -3.5 ; -3.4 ; -3.4 ; -3.3 ; -3.1 ; -3.1 ; -3.0 ]
Zmin_VL     = [-0   ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.9 ]
Tlim=addvars(Tlim, APmax_VL, APmin_VL, Zmin_VL, Zmax_VL)

%%
for ii=1:ncell
    MLcol= find(ML_Slice == round(T2.ML(ii),1))
    VM(ii)= T2.AP(ii)>=APmin_VM(MLcol) & T2.AP(ii)<=APmax_VM(MLcol) & T2.DV(ii)>=Zmin_VM(MLcol) & T2.DV(ii)<Zmax_VM(MLcol);
    VL(ii)= T2.AP(ii)>=APmin_VL(MLcol) & T2.AP(ii)<=APmax_VL(MLcol) & T2.DV(ii)>=Zmin_VL(MLcol) & T2.DV(ii)<=Zmax_VL(MLcol);
    nVL=sum(VL)
    nVM=sum(VM)
    nVMVL=nVL+nVM
    ntot=max(size(VM))
end
%%
% VM=str2num(num2str(VM))';
VM=VM';
VL=VL';
T3=[];
T3=addvars(T2,VM,VL)
Tcoord= T3;
save('D:\JC_Analysis\listcell.mat','Tcoord', '-append')

%% Plot 3-D ALL Channels ALL sessions in color (to indicate which is VM/VL or Others)
figure, scatter3(T3.AP, T3.ML, T3.DV, 'k')
hold on, scatter3(T3.AP(T3.VM), T3.ML(T3.VM), T3.DV(T3.VM),'r')
hold on, scatter3(T3.AP(T3.VL), T3.ML(T3.VL), T3.DV(T3.VL),'b')
xlabel('AP'), ylabel('ML'), zlabel('DV')
legend('other', 'VM', 'VL')
title(['nVM=' num2str(nVM) '  nVL=' num2str(nVL) '  nVM/VL=' num2str(nVMVL) '  ntot=' num2str(ntot)])

%% Plot 2-D Channels for each single Sessions in color (to indicate which is VM/VL or Others)
close all
for ii = 1:11
    figure,
    scatter3(T3.AP(T3.nSess==ii), T3.ML(T3.nSess==ii), T3.DV(T3.nSess==ii), 'k')
    hold on, scatter3(T3.AP(T3.VM & T3.nSess==ii), T3.ML(T3.VM & T3.nSess==ii), T3.DV(T3.VM & T3.nSess==ii),'r')
    hold on, scatter3(T3.AP(T3.VL & T3.nSess==ii), T3.ML(T3.VL & T3.nSess==ii), T3.DV(T3.VL & T3.nSess==ii),'b')
    xlabel('AP'), ylabel('ML'), zlabel('DV')
    legend('other', 'VM', 'VL')
    title([ T1.MouseID(min(find(T1.nSess==ii)),:) ' ' T1.Day(min(find(T1.nSess==ii)),:)])
end
%%
% close all

for ii = 2% 1:11
    hold on,
    scatter(-T3.AP(T3.nSess==ii), T3.DV(T3.nSess==ii), 'k')
    hold on, scatter(-T3.AP(T3.VM & T3.nSess==ii),  T3.DV(T3.VM & T3.nSess==ii),'r')
    hold on, scatter(-T3.AP(T3.VL & T3.nSess==ii),  T3.DV(T3.VL & T3.nSess==ii),'b')
    xlabel('AP'), zlabel('DV')
    legend('other', 'VM', 'VL')
    title([ T1.MouseID(min(find(T1.nSess==ii)),:) ' ' T1.Day(min(find(T1.nSess==ii)),:)])
end



