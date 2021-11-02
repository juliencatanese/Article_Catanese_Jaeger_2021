% fig_barPerc_ipsicontra_7types_JCscript

contra_type1 = contra_dz & ~ipsi_dz & Tcombo_z.z_exct_UniMod_1puf; 
contra_type2 = contra_dz & ~ipsi_dz & Tcombo_z.z_exct_UniMod_2del;
contra_type3 = contra_dz & ~ipsi_dz & Tcombo_z.z_exct_UniMod_3res;

ipsi_type1 = ~contra_dz & ipsi_dz & Tcombo_z.z_exct_UniMod_1puf; 
ipsi_type2 = ~contra_dz & ipsi_dz & Tcombo_z.z_exct_UniMod_2del;
ipsi_type3 = ~contra_dz & ipsi_dz & Tcombo_z.z_exct_UniMod_3res;

contra_type4 = contra_dz & ~ipsi_dz & Tcombo_z.z_exct_biMod_12pd; 
contra_type5 = contra_dz & ~ipsi_dz & Tcombo_z.z_exct_biMod_23dr;
contra_type6 = contra_dz & ~ipsi_dz & Tcombo_z.z_exct_biMod_13pr;

ipsi_type4 = ~contra_dz & ipsi_dz & Tcombo_z.z_exct_biMod_12pd; 
ipsi_type5 = ~contra_dz & ipsi_dz & Tcombo_z.z_exct_biMod_23dr;
ipsi_type6 = ~contra_dz & ipsi_dz & Tcombo_z.z_exct_biMod_13pr;

contra_type7 = contra_dz & ~ipsi_dz & Tcombo_z.z_exct_triMod_123; 
ipsi_type7 = ~contra_dz & ipsi_dz & Tcombo_z.z_exct_triMod_123;

close all,  figure, 
data = [sum(contra_type1),  sum(ipsi_type1) ;...
     sum(contra_type2) sum(ipsi_type2); ...
     sum(contra_type3) sum(ipsi_type3);...
     sum(contra_type4),  sum(ipsi_type4) ;...
     sum(contra_type5) sum(ipsi_type5); ...
     sum(contra_type6) sum(ipsi_type6);...
     sum(contra_type7) sum(ipsi_type7)]
bar(data./sum(data)*100, 'grouped')

