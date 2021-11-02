% fig_bar_RvL_3epochDZ_Nosig_JCscript 

contra_type1 = Tdz_LR.cont_puf & Tdz_LR.VMVL & Tcombo_z.z_exct
ipsi_type1 = Tdz_LR.ipsi_puf & Tdz_LR.VMVL & Tcombo_z.z_exct
nosig_type1 = ~Tdz_LR.cont_puf & ~Tdz_LR.ipsi_puf & Tdz_LR.VMVL & Tcombo_z.z_exct

contra_type2 = Tdz_LR.cont_del & Tdz_LR.VMVL & Tcombo_z.z_exct
ipsi_type2 = Tdz_LR.ipsi_del & Tdz_LR.VMVL & Tcombo_z.z_exct
nosig_type2 = ~Tdz_LR.cont_del & ~Tdz_LR.ipsi_del & Tdz_LR.VMVL & Tcombo_z.z_exct

contra_type3 = Tdz_LR.cont_res & Tdz_LR.VMVL & Tcombo_z.z_exct
ipsi_type3 = Tdz_LR.ipsi_res & Tdz_LR.VMVL & Tcombo_z.z_exct
nosig_type3 = ~Tdz_LR.cont_res & ~Tdz_LR.ipsi_res & Tdz_LR.VMVL & Tcombo_z.z_exct

% plot data into a bar chart
data = [sum(contra_type1),  sum(ipsi_type1), sum(nosig_type1);...
        sum(contra_type2),  sum(ipsi_type2), sum(nosig_type2);...
        sum(contra_type3),  sum(ipsi_type3), sum(nosig_type3);]
figure, 
bar(data, 'grouped');

%%%%%%
contra_type1 = Tdz_LR.cont_puf   & ~Tdz_LR.ipsi_puf & Tdz_LR.VMVL & Tcombo_z.z_exct
ipsi_type1   = Tdz_LR.ipsi_puf   & ~Tdz_LR.cont_puf & Tdz_LR.VMVL & Tcombo_z.z_exct
nosig_type1  = ~Tdz_LR.cont_puf  & ~Tdz_LR.ipsi_puf & Tdz_LR.VMVL & Tcombo_z.z_exct
cmplx_type1  = Tdz_LR.cont_puf   &  Tdz_LR.ipsi_puf & Tdz_LR.VMVL & Tcombo_z.z_exct

contra_type2 =  Tdz_LR.cont_del  & ~Tdz_LR.ipsi_del & Tdz_LR.VMVL & Tcombo_z.z_exct
ipsi_type2   =  Tdz_LR.ipsi_del  & ~Tdz_LR.cont_del & Tdz_LR.VMVL & Tcombo_z.z_exct
nosig_type2  = ~Tdz_LR.cont_del  & ~Tdz_LR.ipsi_del & Tdz_LR.VMVL & Tcombo_z.z_exct
cmplx_type2  =  Tdz_LR.cont_del  &  Tdz_LR.ipsi_del & Tdz_LR.VMVL & Tcombo_z.z_exct

contra_type3 = Tdz_LR.cont_res   & ~Tdz_LR.ipsi_res & Tdz_LR.VMVL & Tcombo_z.z_exct
ipsi_type3   = Tdz_LR.ipsi_res   & ~Tdz_LR.cont_res & Tdz_LR.VMVL & Tcombo_z.z_exct
nosig_type3  = ~Tdz_LR.cont_res  & ~Tdz_LR.ipsi_res & Tdz_LR.VMVL & Tcombo_z.z_exct
cmplx_type3  = Tdz_LR.cont_res   &  Tdz_LR.ipsi_res & Tdz_LR.VMVL & Tcombo_z.z_exct


data = [sum(contra_type1),  sum(ipsi_type1), sum(nosig_type1), sum(cmplx_type1) ;...
        sum(contra_type2),  sum(ipsi_type2), sum(nosig_type2), sum(cmplx_type2);...
        sum(contra_type3),  sum(ipsi_type3), sum(nosig_type3), sum(cmplx_type3);]
figure, 
bar(data, 'grouped');




% figure, 
% colormap('white')
% b=bar(data, 'grouped', 'FaceColor', 'flat'); title('Exct only');
% for k = 1:size(data,2)
%     b(k).CData = k;
% end
% b(2).EdgeColor = 'c'; b(2).LineWidth = 2;
% b(3).EdgeColor = 'm'; b(3).LineWidth = 2;
% 



