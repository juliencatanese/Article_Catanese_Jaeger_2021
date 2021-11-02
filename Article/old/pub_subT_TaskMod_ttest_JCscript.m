% pub_comp_EphysType_JCscript

puf_Exc = Tcombo.VMVL &...
    (Tephys.H_puf & (Tephys.Fr_puf_Hz > Tephys.Fr_BLE_Hz))    &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz))) ; 

del_Exc = Tcombo.VMVL &...
    (Tephys.H_del & (Tephys.Fr_del_Hz > Tephys.Fr_BLE_Hz))    &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz))) ; 

res_Exc  = Tcombo.VMVL &...
    (Tephys.H_res & (Tephys.Fr_res_Hz > Tephys.Fr_BLE_Hz))    &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz))) ; 

type7_Exc = Tcombo.VMVL & puf_Exc & del_Exc & res_Exc;
type6_Exc = Tcombo.VMVL & puf_Exc & res_Exc & ~type7_Exc;
type5_Exc = Tcombo.VMVL & del_Exc & res_Exc & ~type7_Exc;
type4_Exc = Tcombo.VMVL & puf_Exc & del_Exc & ~type7_Exc;
type3_Exc = Tcombo.VMVL & res_Exc & ~del_Exc & ~puf_Exc; 
type2_Exc = Tcombo.VMVL & del_Exc & ~res_Exc & ~puf_Exc; 
type1_Exc = Tcombo.VMVL & puf_Exc & ~del_Exc & ~res_Exc; 
type0     = Tcombo.VMVL & Tephys.H_puf==0 &  Tephys.H_del==0 & Tephys.H_res==0; 


Ntype0 = sum(type0)
Ntype1_Exc = sum(type1_Exc)
Ntype2_Exc = sum(type2_Exc)
Ntype3_Exc = sum(type3_Exc)
Ntype4_Exc = sum(type4_Exc)
Ntype5_Exc = sum(type5_Exc)
Ntype6_Exc = sum(type6_Exc)
Ntype7_Exc = sum(type7_Exc)
