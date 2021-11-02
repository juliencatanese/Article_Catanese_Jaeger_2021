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

%% 1 - VENNX DIAGRAM POSITIVE 

Npuf = sum(Tephys_z.zpuf_exct & Tcombo.VMVL)
Ndel = sum(Tephys_z.zdel_exct & Tcombo.VMVL)
Nres = sum(Tephys_z.zres_exct & Tcombo.VMVL)

i12 = sum(Tcombo.VMVL & Tephys_z.zpuf_exct & Tephys_z.zdel_exct);
i13 = sum(Tcombo.VMVL & Tephys_z.zpuf_exct & Tephys_z.zres_exct);
i23 = sum(Tcombo.VMVL & Tephys_z.zdel_exct & Tephys_z.zres_exct );
i123= sum(Tcombo.VMVL & Tephys_z.zpuf_exct &  Tephys_z.zdel_exct & Tephys_z.zres_exct);

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

