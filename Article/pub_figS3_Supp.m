% pub_figS3_Supp
% by JC 8/25/2019

%% SUPP FIG S3A: TASK Exct vs Inhib
% 3D PLOT for each cell types 1, 2 , 3
SaveFolder= 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\supplementary'
load('Tfig1_VMopto.mat')
load('Tfig2_cor.mat')

for nc= 1:Tfig1_VMopto.ncell(end) %545
    VAL(nc,1)=  logical(strcmp(Tfig1_VMopto.AreaID{nc},'VAL'));
end

Xvmvl= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL));
Yvmvl= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL));
Zvmvl= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL));

Xvl= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & VAL));
Yvl= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & VAL));
Zvl= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & VAL));

X1= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & Tfig2_cor.z_exct));
Y1= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & Tfig2_cor.z_exct));
Z1= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & Tfig2_cor.z_exct));

X2= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & Tfig2_cor.z_inib));
Y2= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & Tfig2_cor.z_inib));
Z2= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & Tfig2_cor.z_inib));

figure,

scatter3(X1(:),Y1(:),Z1(:),'o','filled'), view(-60,60)
hold on,
scatter3(X2(:),Y2(:),Z2(:),'o','filled'), view(-60,60)
hold on
scatter3(Xvmvl(:),Yvmvl(:),Zvmvl(:),'r'), view(-60,60)
hold on,
scatter3(Xvl(:),Yvl(:),Zvl(:),'b'), view(-60,60)

xlabel('AP')
ylabel('ML')
zlabel('DV')
legend ('task+','task-', 'VM', 'VAL')

saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3A.fig')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3A.png')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3A.emf')

saveas(gcf, [SaveFolder '\S3A.fig'])
saveas(gcf, [SaveFolder '\S3A.png'])
saveas(gcf, [SaveFolder '\S3A.emf'])

% Distrib and stats 
disp('AP')
[Px1 Hx1] = ranksum(X1,X2); if Hx1; Px1,  [mean(X1), std(X1); mean(X2), std(X2)], end  % Px1 = 0.0364
disp('ML')
[Py1 Hy1] = ranksum(Y1,Y2); if Hy1; Py1, [mean(Y1), std(Y1); mean(Y2), std(Y2)],end
disp('DV')
[Pz1 Hz1] = ranksum(Z1,Z2); if Hz1; Pz1, [mean(Z1), std(Z1); mean(Z2), std(Z2)],end 

figure, 
subplot(3,1,1)
hold on
BMIN = max([min(X1) min(X2)]); 
BMAX = min([max(X1) max(X2)]); 
Hy1 =histogram(X1,'Normalization','count', 'BinWidth',0.1,'BinLimits',[BMIN,BMAX]), 
Hy2 =histogram(X2,'Normalization','count', 'BinWidth',0.1,'BinLimits',[BMIN,BMAX]),
xlabel('AP')
legend('Task+ cells','Task- cells')
subplot(3,1,2)
hold on
BMIN = max([min(Y1) min(Y2) ]); 
BMAX = min([max(Y1) max(Y2) ]); 
Hy1 =histogram(Y1,'Normalization','count', 'BinWidth',0.1,'BinLimits',[BMIN,BMAX]), 
Hy2 =histogram(Y2,'Normalization','count', 'BinWidth',0.1,'BinLimits',[BMIN,BMAX]),
xlabel('ML')
legend('Task+ cells','Task- cells')
subplot(3,1,3)
hold on
histogram(Z1,10, 'Normalization','count'), 
histogram(Z2,10,'Normalization','count')
xlabel('DV')
legend('Task+ cells','Task- cells')

saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3B.fig')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3B.png')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3B.emf')
saveas(gcf, [SaveFolder '\S3B.fig'])
saveas(gcf, [SaveFolder '\S3B.png'])
saveas(gcf, [SaveFolder '\S3B.emf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUPP FIG S3B: TASK Type1 vs Type2 vs Type3
PUF_idx = Tfig2_cor.z_exct_UniMod_1puf;
DEL_idx = Tfig2_cor.z_exct_UniMod_2del;
RES_idx = Tfig2_cor.z_exct_UniMod_3res;

X1= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & PUF_idx));
Y1= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & PUF_idx));
Z1= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & PUF_idx));

X2= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & DEL_idx));
Y2= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & DEL_idx));
Z2= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & DEL_idx));

X3= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & RES_idx));
Y3= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & RES_idx));
Z3= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & RES_idx));

figure,
scatter3(X1(:),Y1(:),Z1(:),'o','filled'), view(-60,60)
hold on
scatter3(X2(:),Y2(:),Z2(:),'o','filled'), view(-60,60)
hold on
scatter3(X3(:),Y3(:),Z3(:),'y','filled'), view(-60,60)
hold on
scatter3(Xvmvl(:),Yvmvl(:),Zvmvl(:),'r'), view(-60,60)
hold on,
scatter3(Xvl(:),Yvl(:),Zvl(:),'b'), view(-60,60)

xlabel('AP')
ylabel('ML')
zlabel('DV')
legend ('sample','delay','response', 'VM', 'VAL')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3C.fig')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3C.png')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3C.emf')

saveas(gcf, [SaveFolder '\S3C.fig'])
saveas(gcf, [SaveFolder '\S3C.png'])
saveas(gcf, [SaveFolder '\S3C.emf'])

% Distrib and stats 
disp('AP')
[Px1 Hx1] = ranksum(X1,[X2; X3]); if Hx1; Px1, end % p=0.0037
[Px2 Hx2] = ranksum(X2,[X1; X3]); if Hx2; Px2, end 
[Px3 Hx3] = ranksum(X3,[X1; X2]); if Hx3; Px3, end 
disp('ML')
[Py1 Hy1] = ranksum(Y1,[Y2; Y3]); if Hy1; Py1, end % p=0.0257
[Py2 Hy2] = ranksum(Y2,[Y1; Y3]); if Hy2; Py2, end % p=0.0086
[Py3 Hy3] = ranksum(Y3,[Y2; Y1]); if Hy3; Py3, end 
disp('DV')
[Pz1 Hz1] = ranksum(Z1,[Z2; Z3]); if Hz1; Pz1, end 
[Pz2 Hz2] = ranksum(Z2,[Z1; Z3]); if Hz2; Pz2, end 
[Pz3 Hz3] = ranksum(Z3,[Z1; Z2]); if Hz3; Pz3, end 

figure, 
subplot(3,1,1)
hold on
histogram(X1,'Normalization','count'), histogram(X2,'Normalization','count'), histogram(X3,'Normalization','count')
xlabel('AP')
legend('sample','delay','response')
subplot(3,1,2)
hold on
histogram(Y1,'Normalization','count'), histogram(Y2,'Normalization','count'), histogram(Y3,'Normalization','count')
xlabel('ML')
legend('sample','delay','response')
subplot(3,1,3)
hold on
histogram(Z1,'Normalization','count'), histogram(Z2,'Normalization','count'), histogram(Z3,'Normalization','count')
xlabel('DV')
legend('sample','delay','response')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3D.fig')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3D.png')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3D.emf')

saveas(gcf, [SaveFolder '\S3D.fig'])
saveas(gcf, [SaveFolder '\S3D.png'])
saveas(gcf, [SaveFolder '\S3D.emf'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUPP FIG S3C: TASK UNI vs DUO vs TRI MODAL
TRI_idx = Tfig2_cor.z_exct_triMod_123;
DUO_idx = Tfig2_cor.z_exct_biMod_12 | Tfig2_cor.z_exct_biMod_13 | Tfig2_cor.z_exct_biMod_23;
UNI_idx = Tfig2_cor.z_exct_UniMod_1puf | Tfig2_cor.z_exct_UniMod_2del | Tfig2_cor.z_exct_UniMod_3res;


X1= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & UNI_idx));
Y1= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & UNI_idx));
Z1= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & UNI_idx));

X2= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & DUO_idx));
Y2= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & DUO_idx));
Z2= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & DUO_idx));

X3= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & TRI_idx));
Y3= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & TRI_idx));
Z3= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & TRI_idx));

figure,
scatter3(X1(:),Y1(:),Z1(:),'o','filled'), view(-60,60)
hold on
scatter3(X2(:),Y2(:),Z2(:),'o','filled'), view(-60,60)
hold on
scatter3(X3(:),Y3(:),Z3(:),'g','filled'), view(-60,60)
hold on
scatter3(Xvmvl(:),Yvmvl(:),Zvmvl(:),'r'), view(-60,60)
hold on,
scatter3(Xvl(:),Yvl(:),Zvl(:),'b'), view(-60,60)

xlabel('AP')
ylabel('ML')
zlabel('DV')
legend ('Uni','Bi','Tri', 'VM', 'VAL')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3E.fig')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3E.png')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3E.emf')
saveas(gcf, [SaveFolder '\S3E.fig'])
saveas(gcf, [SaveFolder '\S3E.png'])
saveas(gcf, [SaveFolder '\S3E.emf'])


% Distrib and stats 
disp('AP')
[Px1 Hx1] = ranksum(X1,[X2; X3]); if Hx1; Px1, end 
[Px2 Hx2] = ranksum(X2,[X1; X3]); if Hx2; Px2, end 
[Px3 Hx3] = ranksum(X3,[X1; X2]); if Hx3; Px3, end 
disp('ML')
[Py1 Hy1] = ranksum(Y1,[Y2; Y3]); if Hy1; Py1, end 
[Py2 Hy2] = ranksum(Y2,[Y1; Y3]); if Hy2; Py2, end 
[Py3 Hy3] = ranksum(Y3,[Y2; Y1]); if Hy3; Py3, end 
disp('DV')
[Pz1 Hz1] = ranksum(Z1,[Z2; Z3]); if Hz1; Pz1, end % P=0.0013
[Pz2 Hz2] = ranksum(Z2,[Z1; Z3]); if Hz2; Pz2, end 
[Pz3 Hz3] = ranksum(Z3,[Z1; Z2]); if Hz3; Pz3, end % P=0.0238

figure, 
subplot(3,1,1)
hold on
histogram(X1,'Normalization','count'), histogram(X2,'Normalization','count'), histogram(X3,'Normalization','count')
xlabel('AP')
legend('Uni','Bi','Tri')
subplot(3,1,2)
hold on
histogram(Y1,'Normalization','count'), histogram(Y2,'Normalization','count'), histogram(Y3,'Normalization','count')
xlabel('ML')
legend('Uni','Bi','Tri')
subplot(3,1,3)
hold on
histogram(Z1,'Normalization','count'), histogram(Z2,'Normalization','count'), histogram(Z3,'Normalization','count')
xlabel('DV')
legend('Uni','Bi','Tri')

saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3F.fig')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3F.png')
saveas(gcf, 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\Fig3_cor\SupplementaryS3\S3F.emf')

saveas(gcf, [SaveFolder '\S3F.fig'])
saveas(gcf, [SaveFolder '\S3F.png'])
saveas(gcf, [SaveFolder '\S3F.emf'])

