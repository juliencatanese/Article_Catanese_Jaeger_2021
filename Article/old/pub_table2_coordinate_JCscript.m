% pub_table2_coordinate_JCscript
% create and save a table ('Tcoord') that contain all cells informations
% Fields: ncell, nMouse, nDay, nChan, nShank, AP, ML, DV, VM, VL.
% save/load('D:\JC_Analysis\listcell.mat','listcell', 'Tcoord')
% written by Julien Catanese 12/05/2018
% last updated JC 2/7/2019.

clear all;
AP_measured=[]; ML_measured=[]; DV_measured=[];
AP_corr=[]; ML_corr=[]; DV_corr=[]; C=[]; B=[]; A=[];

cd('D:\JC_Analysis');
listcell_File = 'listcell.mat'
load(listcell_File);  T1= listcell; T1=removevars(T1,[5,6,7,8])

%% Add nSess (Session number from 1:nSess instead of the string Day)
listSess=unique(listcell.SessID)
nSess = max(size(listSess))
for ii= 1:nSess
    A(:,ii)= strcmp(listcell.SessID, listSess{ii});
end
B=A.*[1:1:nSess];
C=sum(B')'; % vector with Session Number from 1 to 11.
T1= addvars(T1, C, 'Before', 5, 'NewVariableNames', {'nSess'});
T1(1,:);

%% Add nMouse (Mice number from 1:nMouse instead of the string MouseID)
ncell=T1.ncell(end);
C=[]; B=[]; A=[];
for ii=1:ncell;
    MID{ii}= T1.MouseID(ii,:);
end
listMice=unique(MID');
nMice = max(size(listMice))
for ii= 1:nMice
    A(:,ii)= strcmp(MID, listMice{ii});
end
B=A.*[1:1:nMice];
C=sum(B')'; % vector with Session Number from 1 to 11.
T1= addvars(T1, C, 'Before', 6, 'NewVariableNames', {'nMouse'});
T1(1,:)

%% Manual Attribution of Coordinate (based on Stereotaxic Measurement and Histo reconstruction and Ephys Markers)
% 'Sh1M_' = The coordinate in the table corespond to channel #8 (always tip) on Shank #1 (the most Medial M_)


SessionNAME = ['11_w10d4'; '12_w11d5'; '14_w14d2'; '14_w14d5'; '14_w14d7'; '14_w14d8'; '15_w08d5'; '15_w10d3'; '15_w10d7'; '15_w10d8'; '17_w10d3'; '17_w10d4'; '17_w10d5'; '17_w10d7'; '17_w10d8'  ];
Ch_Origin   = ['S1Ch8M_' ; 'S2Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S1Ch8_P' ; 'S4Ch8MA' ; 'S4Ch8MA' ;  'S1Ch8LP'; 'S1Ch8M_' ; 'S2Ch8_L' ; 'S1Ch4_Q' ; 'S1Ch8M_' ;  'S1Ch8M_'  ];
Ch_Orient   = ['ML'      ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ;  'AP'     ; 'DG'      ;  'DG'     ;  'DG'     ;  'ML'     ;  'ML'     ;  'Q4'     ;    'ML'   ;   'ML'      ];
ElectrodeID = ['CC4F_n4S';'CCAE_o3S' ;'G912_n4S' ;'G912_o4S' ;'CC4F_n4S' ;'G912_o4S' ;'G912_o4S' ; 'CAED_n4S';'G912_o4S' ;'G912_o4S' ;'CAED_n4S' ; 'CCAE_o3S'; 'EOBD_o1S'; 'G912_o4S'; 'G912_o4S'  ];
AP_measured = [1.30      ;  1.58     ;  1.10     ;  1.55     ;  1.50     ;  1.55     ;  1.90     ;  1.25     ;   1.25    ;   1.5     ;  1.00     ;  1.00     ;   1.3     ;     1.45  ;     1.2     ];
ML_measured = [0.30      ;  0.90     ;  0.75     ;  0.90     ;  0.80     ;  0.85     ;  0.70     ;  0.9      ;   0.85    ;   0.7     ;  0.50     ;  1.10     ;   0.65    ;     0.65  ;     0.65    ];
DV_measured = [4.30      ;  4.30     ;  4.65     ;  4.30     ;  4.35     ;  4.90     ;  4.70     ;  4.25     ;   4.4     ;   4.6     ;  4.30     ;  4.20     ;  4.30    ;     4.30  ;     4.5     ];
APAdj       = [+0.5      ;  +0.0     ;  +0.1     ;  +0.1     ;  +0.1     ;  +0.1     ;  +0.0     ;   +0.0    ;   +0.0    ;   +0.0    ;  +0.20    ;  +0.0     ;   +0.0    ;     +0.0  ;     +0.0    ];
MLAdj       = [+0.1      ;  +0.0     ;  +0.0     ;  +0.0     ;  +0.0     ;  +0.0     ;  +0.0     ;   +0.0    ;   +0.0    ;   +0.0    ;  +0.1     ;  +0.0     ;   +0.1    ;     +0.0  ;     +0.0    ];
DVAdj       = [+0.0      ;  +0.0     ;  -0.35     ;  +0.0     ;  +0.0    ;  -0.4     ;  +0.0     ;   +0.0    ;   +0.0    ;   +0.0    ;  +0.0     ;  +0.0     ;   +0.0    ;     +0.0  ;     -0.3    ];


AP_corr = AP_measured + APAdj;
ML_corr = ML_measured + MLAdj;
DV_corr = DV_measured + DVAdj ;

Tst= table( Ch_Origin, Ch_Orient, AP_corr, ML_corr, DV_corr, listSess)

%%
nShank = str2num(T1.ChanID(:,2));
nChan  = str2num(T1.ChanID(:,5));
T0=[]; T0=T1;
T2=addvars(T0,nShank,nChan); T2(1,:)

%%
AP=[]; ML=[];DV=[];
for ii=1:ncell;
    %   coln =;
    Orig=Tst.Ch_Origin(T1.nSess(ii),:);
    ShOrig=str2num(Orig(2));
    ChOrig=str2num(Orig(5));
    Orient=Orig(6:7);
    if ShOrig~=4
        if Orient=='M_'
            AP(ii) = Tst.AP_corr(T1.nSess(ii));
            ML(ii) = Tst.ML_corr(T1.nSess(ii))+(0.2*(-ShOrig+T2.nShank(ii)));
            DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
            
        elseif Orient=='_L'
            AP(ii) = Tst.AP_corr(T1.nSess(ii));
            ML(ii) = Tst.ML_corr(T1.nSess(ii))-(0.2*(-ShOrig+T2.nShank(ii)));
            DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
            
        elseif Orient=='_P'
            AP(ii) = Tst.AP_corr(T1.nSess(ii))-(0.2*(-ShOrig+T2.nShank(ii)));
            ML(ii) = Tst.ML_corr(T1.nSess(ii));
            DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
            
        elseif Orient=='_Q'
            AP(ii) = Tst.AP_corr(T1.nSess(ii));
            ML(ii) = Tst.ML_corr(T1.nSess(ii));
            DV(ii) = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
            
        elseif Orient=='LP'
            AP(ii) = Tst.AP_corr(T1.nSess(ii))-(((600*cos(pi/6))/3000)*(-ShOrig+T2.nShank(ii)));
            ML(ii) = Tst.ML_corr(T1.nSess(ii))-(((600*sin(pi/6))/3000)*(-ShOrig+T2.nShank(ii)));
            DV(ii) = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
        else
            lvkjsdlvjsdlvkj
        end
        
    elseif ShOrig==4
        if Orient=='MA'
            AP(ii) = Tst.AP_corr(T1.nSess(ii))+(((600*cos(pi/6))/3000)*(ShOrig-T2.nShank(ii)));
            ML(ii) = Tst.ML_corr(T1.nSess(ii))+(((600*sin(pi/6))/3000)*(ShOrig-T2.nShank(ii)));
            DV(ii)  = Tst.DV_corr(T1.nSess(ii))-(0.1*(ChOrig-T2.nChan(ii)));
            
        end
    end
end
AP=-AP'; ML=-ML'; DV=-DV';
%% Add Coordinates to Table
nShank = str2num(T1.ChanID(:,2));
nChan  = str2num(T1.ChanID(:,5));
T0=[]; T0=T1; T2=[];
T2=addvars(T0, nShank, nChan, AP, ML, DV); T2(1,:)

%% Select the Axis with less variability
NML=max(size(unique(round(T2.ML,1))))
NAP=max(size(unique(round(T2.AP,1))))
if NAP<NML
    display('NAP<NML')
    AP_slice =unique(round(T2.AP,1));
else
    display('NML<NAP')
    ML_slice = unique(round(T2.ML,1))
end

%% Manual definition of the limite of VM and VL/VA in the Mouse Atlas (Paxinos)
% Tlim = table(ML_slice);
ML_Slice    = [-0.4 ; -0.5 ; -0.6 ; -0.7 ; -0.8 ; -0.9 ; -1.0 ; -1.1 ; -1.2 ; -1.3 ];
APmax_VM    = [-1.7 ; -1.7 ; -1.0 ; -1.0 ; -1.0 ; -1.0 ; -1.0 ; -1.1 ; -1.2 ; -1.3 ];
APmin_VM    = [-2.3 ; -2.3 ; -2.2 ; -2.2 ; -2.0 ; -2.2 ; -2.2 ; -2.2 ; -2.0 ; -1.7 ];
Zmax_VM     = [-4.1 ; -4.0 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.9 ];
Zmin_VM     = [-4.5 ; -4.5 ; -4.45 ; -4.40 ; -4.35 ; -4.30 ; -4.3 ; -4.3 ; -4.2 ; -4.1 ];
Tlim=table(ML_Slice, APmax_VM, APmin_VM, Zmax_VM, Zmin_VM);

ML_Slic_VL  = [-0.4 ; -0.5 ; -0.6 ; -0.7 ; -0.8 ; -0.9 ; -1.0 ; -1.1 ; -1.2 ; -1.3 ];
APmax_VL    = [-0   ; -0   ; -1.0 ; -1.2 ; -1.0 ; -1.0 ; -1.0 ; -1.0 ; -1.0 ; -1.0 ];
APmin_VL    = [-0   ; -0   ; -1.5 ; -2.0 ; -1.9 ; -1.9 ; -1.8 ; -1.8 ; -1.5 ; -1.4 ];
Zmax_VL     = [-0   ; -0   ; -3.6 ; -3.5 ; -3.4 ; -3.4 ; -3.3 ; -3.1 ; -3.1 ; -3.0 ];
Zmin_VL     = [-0   ; -0   ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.8 ; -3.9 ];
Tlim=addvars(Tlim, APmax_VL, APmin_VL, Zmin_VL, Zmax_VL)

%%
VM=[];  VL=[];
for ii=1:ncell
    MLcol= find(ML_Slice == round(T2.ML(ii),1));
    VM(ii)= T2.AP(ii)>=APmin_VM(MLcol) & T2.AP(ii)<=APmax_VM(MLcol) & T2.DV(ii)>=Zmin_VM(MLcol) & T2.DV(ii)<Zmax_VM(MLcol);
    VL(ii)= T2.AP(ii)>=APmin_VL(MLcol) & T2.AP(ii)<=APmax_VL(MLcol) & T2.DV(ii)>=Zmin_VL(MLcol) & T2.DV(ii)<=Zmax_VL(MLcol);
    nVL=sum(VL)
    nVM=sum(VM)
    nVMVL=nVL+nVM
    ntot=max(size(VM))
end
%%
% VM=str2num(num2str(VM))';
VM=logical(VM');
VL=logical(VL');
T3=[];
T3=addvars(T2,VM,VL);
Tcoord= T3;
save(['D:\JC_Analysis\' listcell_File],'Tcoord', '-append')

%% Plot 2-D Channels for each single Sessions in color (to indicate which is VM/VL or Others)
close all
for ii = 1:size(Tst,1)
    figure,
    pause(0.2)
    scatter3(Tcoord.AP(Tcoord.nSess==ii), Tcoord.ML(Tcoord.nSess==ii), Tcoord.DV(Tcoord.nSess==ii), 'k')
    hold on, scatter3(Tcoord.AP(Tcoord.VM & Tcoord.nSess==ii), Tcoord.ML(Tcoord.VM & Tcoord.nSess==ii), Tcoord.DV(Tcoord.VM & Tcoord.nSess==ii),'r')
    hold on, scatter3(Tcoord.AP(Tcoord.VL & Tcoord.nSess==ii), Tcoord.ML(Tcoord.VL & Tcoord.nSess==ii), Tcoord.DV(Tcoord.VL & Tcoord.nSess==ii),'b')
    xlabel('AP'), ylabel('ML'), zlabel('DV')
    legend('other', 'VM', 'VL')
    title([ T1.MouseID(min(find(T1.nSess==ii)),:) ' ' T1.Day(min(find(T1.nSess==ii)),:)])
end

%% Plot 3-D ALL Channels ALL sessions in color (to indicate which is VM/VL or Others)
figure, scatter3(Tcoord.AP, Tcoord.ML, Tcoord.DV, 'k')
hold on, scatter3(Tcoord.AP(Tcoord.VM), Tcoord.ML(Tcoord.VM), Tcoord.DV(Tcoord.VM),'r')
hold on, scatter3(Tcoord.AP(Tcoord.VL), Tcoord.ML(Tcoord.VL), Tcoord.DV(Tcoord.VL),'b')
xlabel('AP'), ylabel('ML'), zlabel('DV')
legend('other', 'VM', 'VL')
title(['nVM=' num2str(nVM) '  nVL=' num2str(nVL) '  nVM/VL=' num2str(nVMVL) '  ntot=' num2str(ntot)])

