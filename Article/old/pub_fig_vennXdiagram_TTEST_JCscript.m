% pub_fig_vennX_Diagram_JCscript 
% load ('listcell.mat');
% function error = vennX( data, resolution )
% vennX - draws an area proportional venn diagram
% INPUT data = 
% |A|
% |A and B|
% |B|
% |B and C|
% |C|
% |C and A|
% |A and B and C|
% resolution is between 0.1 and 0.001 (small value better resolution but longer time) 

close all
load ('listcell.mat');

%% 1 - VENNX DIAGRAM = ALL 

Npuf = sum(Tephys.H_puf & Tcombo.VMVL)
Ndel = sum(Tephys.H_del & Tcombo.VMVL)
Nres = sum(Tephys.H_res & Tcombo.VMVL)

type1 = sum(Tcombo.Type1Cell & Tcombo.VMVL);
type2 = sum(Tcombo.Type2Cell & Tcombo.VMVL);
type3 = sum(Tcombo.Type3Cell & Tcombo.VMVL);

i12 = sum(Tephys.H_puf & Tephys.H_del & Tcombo.VMVL);
i13 = sum(Tephys.H_puf & Tephys.H_res & Tcombo.VMVL);
i23 = sum(Tephys.H_del & Tephys.H_res & Tcombo.VMVL);
i123 = sum(Tephys.H_puf & Tephys.H_del & Tephys.H_res & Tcombo.VMVL);

data = [type1, i12-i123 , type2, i23-i123, type3, i13-i123, i123] 
resolution = 0.3; 

vennX( data, resolution, colmap )
title(['ALL cell Sig in VM/VL  (#cells=' num2str( sum(data)) ')'])


%% 2 - VENNX DIAGRAM = positive only  
Hpuf_pos = (Tephys.H_puf & (Tephys.Fr_puf_Hz > Tephys.Fr_BLE_Hz)) &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz))); 

Hdel_pos = (Tephys.H_del & (Tephys.Fr_del_Hz > Tephys.Fr_BLE_Hz)) &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz))); 

Hres_pos = (Tephys.H_res & (Tephys.Fr_res_Hz > Tephys.Fr_BLE_Hz)) &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz))); 

Npuf = sum(Tcombo.VMVL & Hpuf_pos)
Ndel = sum(Tcombo.VMVL & Hdel_pos) 
Nres = sum(Tcombo.VMVL & Hres_pos)  

i12 = sum(Tcombo.VMVL & Hpuf_pos & Hdel_pos);
i13 = sum(Tcombo.VMVL & Hpuf_pos & Tephys.H_res);
i23 = sum(Tcombo.VMVL & Hdel_pos & Hres_pos );
i123= sum(Tcombo.VMVL & Hpuf_pos &  Hdel_pos & Hres_pos);

data = [Npuf-(i12-i123)-i123-(i13-i123) , ... 
    i12-i123, ...
    Ndel-(i12-i123)-i123-(i23-i123), ...
    i23-i123, ...
    Nres-(i23-i123)-i123-(i13-i123), ...
    i13-i123, ...
    i123]; 

resolution = parfig.resolution; 

vennX( data, resolution, colmap )
title(['POSITIVE only (#cells=' num2str( sum(data)) ')'])
colormapeditor

%% 3 - VENNX DIAGRAM = Negative only  
Hpuf_pos = (Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz)) &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz > Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz > Tephys.Fr_BLE_Hz))); 

Hdel_pos = (Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz)) &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz > Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_res & (Tephys.Fr_res_Hz > Tephys.Fr_BLE_Hz))); 

Hres_pos = (Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz)) &...
    (~(Tephys.H_del & (Tephys.Fr_del_Hz > Tephys.Fr_BLE_Hz))) &...
    (~(Tephys.H_puf & (Tephys.Fr_puf_Hz > Tephys.Fr_BLE_Hz))); 

Npuf = sum(Tcombo.VMVL & Hpuf_pos)
Ndel = sum(Tcombo.VMVL & Hdel_pos) 
Nres = sum(Tcombo.VMVL & Hres_pos)  

i12 = sum(Tcombo.VMVL & Hpuf_pos & Hdel_pos);
i13 = sum(Tcombo.VMVL & Hpuf_pos & Tephys.H_res);
i23 = sum(Tcombo.VMVL & Hdel_pos & Hres_pos );
i123= sum(Tcombo.VMVL & Hpuf_pos &  Hdel_pos & Hres_pos);

data = [Npuf-(i12-i123)-i123-(i13-i123), i12-i123, Ndel-(i12-i123)-i123-(i23-i123), i23-i123, Nres-(i23-i123)-i123-(i13-i123), i13-i123, i123] 
resolution = 0.015; 

vennX( data, resolution, colmap )
title(['NEGATIVE only (#cells=' num2str( sum(data)) ')'])


%% 4 - VENNX DIAGRAM = Complex only  

Hpuf_pos = (Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz)) &...
    (((Tephys.H_del & (Tephys.Fr_del_Hz > Tephys.Fr_BLE_Hz))) |...
    ((Tephys.H_res & (Tephys.Fr_res_Hz > Tephys.Fr_BLE_Hz)))) |... 
    (Tephys.H_puf & (Tephys.Fr_puf_Hz > Tephys.Fr_BLE_Hz))   &...
    (((Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) |...
    ((Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz)))) ; 

Hdel_pos = (Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz)) &...
    (((Tephys.H_puf & (Tephys.Fr_puf_Hz > Tephys.Fr_BLE_Hz))) |...
    ((Tephys.H_res & (Tephys.Fr_res_Hz > Tephys.Fr_BLE_Hz)))) |...
    (Tephys.H_del & (Tephys.Fr_del_Hz > Tephys.Fr_BLE_Hz))   &...
    (((Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz))) |...
    ((Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz)))) ; 

Hres_pos = (Tephys.H_res & (Tephys.Fr_res_Hz < Tephys.Fr_BLE_Hz)) &...
    (((Tephys.H_del & (Tephys.Fr_del_Hz > Tephys.Fr_BLE_Hz))) |...
    ((Tephys.H_puf & (Tephys.Fr_puf_Hz > Tephys.Fr_BLE_Hz)))) |... 
    (Tephys.H_res & (Tephys.Fr_res_Hz > Tephys.Fr_BLE_Hz))   &...
    (((Tephys.H_del & (Tephys.Fr_del_Hz < Tephys.Fr_BLE_Hz))) |...
    ((Tephys.H_puf & (Tephys.Fr_puf_Hz < Tephys.Fr_BLE_Hz)))) ; 

Npuf = sum(Tcombo.VMVL & Hpuf_pos)
Ndel = sum(Tcombo.VMVL & Hdel_pos) 
Nres = sum(Tcombo.VMVL & Hres_pos)  

i12 = sum(Tcombo.VMVL & Hpuf_pos & Hdel_pos);
i13 = sum(Tcombo.VMVL & Hpuf_pos & Tephys.H_res);
i23 = sum(Tcombo.VMVL & Hdel_pos & Hres_pos );
i123= sum(Tcombo.VMVL & Hpuf_pos &  Hdel_pos & Hres_pos);

data = [Npuf-(i12-i123)-i123-(i13-i123) , i12-i123, Ndel-(i12-i123)-i123-(i23-i123), i23-i123, Nres-(i23-i123)-i123-(i13-i123), i13-i123, i123] 
resolution = 0.015; 

vennX( data, resolution, colmap )
title(['COMPLEXE only (#cells=' num2str( sum(data)) ')'])

