% fig_bar_RvL_3epochDZ_opto_JCscript

contra_type1 = Tdz_LR.cont_puf & Tdz_LR.VMVL & Tcombo_z.z_exct
contra_type2 = Tdz_LR.cont_del & Tdz_LR.VMVL & Tcombo_z.z_exct
contra_type3 = Tdz_LR.cont_res & Tdz_LR.VMVL & Tcombo_z.z_exct

ipsi_type1 = Tdz_LR.ipsi_puf & Tdz_LR.VMVL & Tcombo_z.z_exct
ipsi_type2 = Tdz_LR.ipsi_del & Tdz_LR.VMVL & Tcombo_z.z_exct
ipsi_type3 = Tdz_LR.ipsi_res & Tdz_LR.VMVL & Tcombo_z.z_exct

% opto cells
opto_exct_contra_type1 = Tdz_LR.cont_puf & Tdz_LR.VMVL & Tcombo_z.z_exct  & Tcombo_z.Opto_exct
opto_exct_contra_type2 = Tdz_LR.cont_del & Tdz_LR.VMVL & Tcombo_z.z_exct & Tcombo_z.Opto_exct
opto_exct_contra_type3 = Tdz_LR.cont_res & Tdz_LR.VMVL & Tcombo_z.z_exct & Tcombo_z.Opto_exct
opto_exct_ipsi_type1 = Tdz_LR.ipsi_puf & Tdz_LR.VMVL & Tcombo_z.z_exct & Tcombo_z.Opto_exct
opto_exct_ipsi_type2 = Tdz_LR.ipsi_del & Tdz_LR.VMVL & Tcombo_z.z_exct & Tcombo_z.Opto_exct
opto_exct_ipsi_type3 = Tdz_LR.ipsi_res & Tdz_LR.VMVL & Tcombo_z.z_exct & Tcombo_z.Opto_exct

opto_inib_contra_type1 = Tdz_LR.cont_puf & Tdz_LR.VMVL & Tcombo_z.z_exct  & Tcombo_z.Opto_inib
opto_inib_contra_type2 = Tdz_LR.cont_del & Tdz_LR.VMVL & Tcombo_z.z_exct & Tcombo_z.Opto_inib
opto_inib_contra_type3 = Tdz_LR.cont_res & Tdz_LR.VMVL & Tcombo_z.z_exct & Tcombo_z.Opto_inib
opto_inib_ipsi_type1 = Tdz_LR.ipsi_puf & Tdz_LR.VMVL & Tcombo_z.z_exct & Tcombo_z.Opto_inib
opto_inib_ipsi_type2 = Tdz_LR.ipsi_del & Tdz_LR.VMVL & Tcombo_z.z_exct & Tcombo_z.Opto_inib
opto_inib_ipsi_type3 = Tdz_LR.ipsi_res & Tdz_LR.VMVL & Tcombo_z.z_exct & Tcombo_z.Opto_inib


% count 
Ncont1 = sum(contra_type1)
Ncont2 = sum(contra_type2)
Ncont3 = sum(contra_type3)
Nipsi1 = sum(ipsi_type1)
Nipsi2 = sum(ipsi_type2)
Nipsi3 = sum(ipsi_type3)
Ncont1_opto_inib = sum(opto_inib_contra_type1)
Ncont2_opto_inib = sum(opto_inib_contra_type2)
Ncont3_opto_inib = sum(opto_inib_contra_type3)
Nipsi1_opto_inib = sum(opto_inib_ipsi_type1)
Nipsi2_opto_inib = sum(opto_inib_ipsi_type2)
Nipsi3_opto_inib = sum(opto_inib_ipsi_type3)
Ncont1_opto_exct = sum(opto_exct_contra_type1)
Ncont2_opto_exct = sum(opto_exct_contra_type2)
Ncont3_opto_exct = sum(opto_exct_contra_type3)
Nipsi1_opto_exct = sum(opto_exct_ipsi_type1)
Nipsi2_opto_exct = sum(opto_exct_ipsi_type2)
Nipsi3_opto_exct = sum(opto_exct_ipsi_type3)


%
data = [sum(contra_type1),  sum(ipsi_type1), sum(opto_inib_contra_type1), sum(opto_inib_ipsi_type1), sum(opto_exct_contra_type1), sum(opto_exct_ipsi_type1);...
        sum(contra_type2) sum(ipsi_type2), sum(opto_inib_contra_type2), sum(opto_inib_ipsi_type2), sum(opto_exct_contra_type2), sum(opto_exct_ipsi_type2); ...
        sum(contra_type3) sum(ipsi_type3), sum(opto_inib_contra_type3), sum(opto_inib_ipsi_type3), sum(opto_exct_contra_type3), sum(opto_exct_ipsi_type3)]
figure, 
bar(data, 'grouped');

data = [(Ncont1/(Ncont1+Nipsi1))- (Nipsi1/(Ncont1+Nipsi1)), ... 
        (Ncont1_opto_inib/(Ncont1_opto_inib+Nipsi1_opto_inib))-(Nipsi1_opto_inib/(Ncont1_opto_inib+Nipsi1_opto_inib)),...
        (Ncont1_opto_exct/(Ncont1_opto_exct+Nipsi1_opto_exct))-(Nipsi1_opto_exct/(Ncont1_opto_exct+Nipsi1_opto_exct));...
        (Ncont2/(Ncont2+Nipsi2))- (Nipsi2/(Ncont2+Nipsi2)),...
        (Ncont2_opto_inib/(Ncont2_opto_inib+Nipsi2_opto_inib))-(Nipsi2_opto_inib/(Ncont2_opto_inib+Nipsi2_opto_inib)),...          
        (Ncont2_opto_exct/(Ncont2_opto_exct+Nipsi2_opto_exct))-(Nipsi2_opto_exct/(Ncont2_opto_exct+Nipsi2_opto_exct));...     
        (Ncont3/(Ncont3+Nipsi3))- (Nipsi3/(Ncont3+Nipsi3)),...
        (Ncont3_opto_inib/(Ncont3_opto_inib+Nipsi3_opto_inib))-(Nipsi3_opto_inib/(Ncont3_opto_inib+Nipsi3_opto_inib)),...        
        (Ncont3_opto_exct/(Ncont3_opto_exct+Nipsi3_opto_exct))-(Nipsi3_opto_exct/(Ncont3_opto_exct+Nipsi3_opto_exct))]; 

figure, 
colormap('white')
b=bar(data*100, 'grouped', 'FaceColor', 'flat'); title('Exct only');
for k = 1:size(data,2)
    b(k).CData = k;
end
b(2).EdgeColor = 'c'; b(2).LineWidth = 2;
b(3).EdgeColor = 'm'; b(3).LineWidth = 2;
% legend('all'; 'opto-'; 'opto+') ; ylabel('delta %'); ylim([-50 50]);   
% bar(data, 'grouped'); title('Exct only');  legend('all'; 'opto-'; ) 



