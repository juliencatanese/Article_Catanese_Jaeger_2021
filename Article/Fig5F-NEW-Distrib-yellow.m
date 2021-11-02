bool_sel(~bool_sel)=1

SaveFolder= 'C:\Users\catan\Documents\EMORY\ARTICLE_JC_DJ\supplementary'

load('Tfig1_VMopto.mat')

for nc= 1:Tfig1_VMopto.ncell(end) %545
    VAL(nc,1)=  logical(strcmp(Tfig1_VMopto.AreaID{nc},'VAL'));
end



Xvmvl= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & bool_sel));
Yvmvl= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & bool_sel));
Zvmvl= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & bool_sel));

Xvl= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & VAL & bool_sel));
Yvl= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & VAL & bool_sel));
Zvl= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & VAL & bool_sel));

X1= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib & bool_sel));
Y1= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib & bool_sel));
Z1= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_inib & bool_sel));

X2= Tfig1_VMopto.AP(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_exct & bool_sel));
Y2= Tfig1_VMopto.ML(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_exct & bool_sel));
Z2= Tfig1_VMopto.DV(logical(Tfig1_VMopto.VMVL & Tfig1_VMopto.Opto_exct & bool_sel));

X3= Tfig1_VMopto.AP(logical(Tfig1_VMopto.Opto_post_sess & Tfig1_VMopto.VMVL & ~Tfig1_VMopto.Opto_exct & ~Tfig1_VMopto.Opto_inib & bool_sel));
Y3= Tfig1_VMopto.ML(logical(Tfig1_VMopto.Opto_post_sess & Tfig1_VMopto.VMVL & ~Tfig1_VMopto.Opto_exct & ~Tfig1_VMopto.Opto_inib & bool_sel));
Z3= Tfig1_VMopto.DV(logical(Tfig1_VMopto.Opto_post_sess & Tfig1_VMopto.VMVL & ~Tfig1_VMopto.Opto_exct & ~Tfig1_VMopto.Opto_inib & bool_sel));

figure, hold on
scatter3(Xvmvl(:),Yvmvl(:),Zvmvl(:),'k','filled'), view(-60,60),hold on
scatter3(X2(:),Y2(:),Z2(:),'y','filled'), view(-60,60),hold on
scatter3(X1(:),Y1(:),Z1(:),'c','filled'), view(-60,60),hold on
% hold on
% scatter3(Xvmvl(:),Yvmvl(:),Zvmvl(:),'r'), view(-60,60)
% hold on,
% scatter3(Xvl(:),Yvl(:),Zvl(:),'b'), view(-60,60)
xlabel('AP')
ylabel('ML')
zlabel('DV')
% legend ('VM/VAL','inib','exct', 'VM', 'VAL')

% saveas(gcf, [SaveFolder '\S1A.fig'])
% saveas(gcf, [SaveFolder '\S1A.png'])
% saveas(gcf, [SaveFolder '\S1A.emf'])

% Distrib and stats 
disp('AP')
[Px1 Hx1] = ranksum(X1,X2); if Hx1; Px1, end  
disp('ML')
[Py1 Hy1] = ranksum(Y1,Y2); if Hy1; Py1, end % Px1 =  0.0016 
disp('DV')
[Pz1 Hz1] = ranksum(Z1,Z2); if Hz1; Pz1, end 
