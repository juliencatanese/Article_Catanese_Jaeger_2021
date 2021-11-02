
        
    % TCOMBO
    Tcombo=[];
    Tcombo=listcell(:,1:5); 
    AreaID=Tfig1_VMopto.AreaID; 
    VMVL=Tfig1_VMopto.VMVL;    
    ncell =  Tfig1_VMopto.ncell; 
    nSess =  Tfig1_VMopto.nSess; 
    nMouse = Tfig1_VMopto.nMouse;
    
    Tcombo = addvars(Tcombo, nMouse, nSess ,VMVL, AreaID);   
    
    Opto_inib = Tfig1_VMopto.Opto_inib; 
    Opto_exct = Tfig1_VMopto.Opto_exct; 
    Opto_post_sess= Tfig1_VMopto.Opto_post_sess;
    
    Tcombo= addvars(Tcombo, Opto_inib, Opto_exct, Opto_post_sess);    
    
    z_exct=Tfig2_cor.z_exct; z_inib=Tfig2_cor.z_inib; z_complex=Tfig2_cor.z_complex;
    z_nosig=Tfig2_cor.z_nosig; z_exct_UniMod_1puf=Tfig2_cor.z_exct_UniMod_1puf;
    z_exct_UniMod_2del=Tfig2_cor.z_exct_UniMod_2del; z_exct_UniMod_3res=Tfig2_cor.z_exct_UniMod_3res;
    z_exct_biMod_12=Tfig2_cor.z_exct_biMod_12; z_exct_biMod_13=Tfig2_cor.z_exct_biMod_13;
    z_exct_biMod_23=Tfig2_cor.z_exct_biMod_23; z_exct_triMod_123=Tfig2_cor.z_exct_triMod_123;   
    
    Tcombo= addvars(Tcombo, z_exct, z_inib, z_complex, z_nosig, ...
        z_exct_UniMod_1puf, z_exct_UniMod_2del, z_exct_UniMod_3res, ...
        z_exct_biMod_12, z_exct_biMod_13, z_exct_biMod_23, ...
        z_exct_triMod_123);
    
    dzRL_exct=Tfig3_RvL.dzRL_exct; 
    dzRL_inib=Tfig3_RvL.dzRL_inib; 
    dzRL_complex=Tfig3_RvL.dzRL_complex; 
    dzRL_nosig=Tfig3_RvL.dzRL_nosig;  
    dzRL_exct_UniMod_1puf=Tfig3_RvL.dzRL_exct_UniMod_1puf;  
    dzRL_exct_UniMod_2del=Tfig3_RvL.dzRL_exct_UniMod_2del; 
    dzRL_exct_UniMod_3res=Tfig3_RvL.dzRL_exct_UniMod_3res; 
    dzRL_exct_biMod_12=Tfig3_RvL.dzRL_exct_biMod_12; 
    dzRL_exct_biMod_13=Tfig3_RvL.dzRL_exct_biMod_13; 
    dzRL_exct_biMod_23=Tfig3_RvL.dzRL_exct_biMod_23; 
    dzRL_exct_triMod_123=Tfig3_RvL.dzRL_exct_triMod_123;   
    
    Tcombo= addvars(Tcombo, dzRL_exct, dzRL_inib, dzRL_complex, dzRL_nosig, ...
        dzRL_exct_UniMod_1puf, dzRL_exct_UniMod_2del, dzRL_exct_UniMod_3res, ...
        dzRL_exct_biMod_12, dzRL_exct_biMod_13, dzRL_exct_biMod_23, ...
        dzRL_exct_triMod_123);
    
    save('Tcombo.mat','Tcombo')
    disp('Tcombo SAVED'); Tcombo(1,:)
    


%% NEXT GENERATION 

%     % TCOMBO
%     Tcombo=[];
%     Tcombo=listcell(:,1:5); AreaID=Tfig1_VMopto.AreaID; VMVL=Tfig1_VMopto.VMVL;    
%     ncell =  Tfig1_VMopto.ncell; nSess =  Tfig1_VMopto.nSess; nMouse = Tfig1_VMopto.nMouse;
%     Tcombo = addvars(Tcombo, nMouse, nSess ,VMVL, AreaID);   
%     Opto_inib = Tfig1_VMopto.Opto_inib; Opto_exct = Tfig1_VMopto.Opto_exct; Opto_post_sess= Tfig1_VMopto.Opto_post_sess;
%     Tcombo= addvars(Tcombo, Opto_inib, Opto_exct, Opto_post_sess);    
%     z_exct=Tfig2_cor.z_exct; z_inib=Tfig2_cor.z_inib; z_complex=Tfig2_cor.z_complex;
%     z_nosig=Tfig2_cor.z_nosig; z_exct_UniMod_1puf=Tfig2_cor.z_exct_UniMod_1puf;
%     z_exct_UniMod_2del=Tfig2_cor.z_exct_UniMod_2del; z_exct_UniMod_3res=Tfig2_cor.z_exct_UniMod_3res;
%     z_exct_biMod_12=Tfig2_cor.z_exct_biMod_12; z_exct_biMod_13=Tfig2_cor.z_exct_biMod_13;
%     z_exct_biMod_23=Tfig2_cor.z_exct_biMod_23; z_exct_triMod_123=Tfig2_cor.z_exct_triMod_123;   
%     Tcombo= addvars(Tcombo, z_exct, z_inib, z_complex, z_nosig, ...
%         z_exct_UniMod_1puf, z_exct_UniMod_2del, z_exct_UniMod_3res, ...
%         z_exct_biMod_12, z_exct_biMod_13, z_exct_biMod_23, ...
%         z_exct_triMod_123);
%     dzRL_exct=Tfig3_RvL.dzRL_exct; 
%     dzRL_inib=Tfig3_RvL.dzRL_inib; 
%     dzRL_complex=Tfig3_RvL.dzRL_complex; 
%     dzRL_nosig=Tfig3_RvL.dzRL_nosig;  
%     dzRL_exct_UniMod_1puf=Tfig3_RvL.dzRL_exct_UniMod_1puf;  
%     dzRL_exct_UniMod_2del=Tfig3_RvL.dzRL_exct_UniMod_2del; 
%     dzRL_exct_UniMod_3res=Tfig3_RvL.dzRL_exct_UniMod_3res; 
%     dzRL_exct_biMod_12=Tfig3_RvL.dzRL_exct_biMod_12; 
%     dzRL_exct_biMod_13=Tfig3_RvL.dzRL_exct_biMod_13; 
%     dzRL_exct_biMod_23=Tfig3_RvL.dzRL_exct_biMod_23; 
%     dzRL_exct_triMod_123=Tfig3_RvL.dzRL_exct_triMod_123;   
%     Tcombo= addvars(Tcombo, dzRL_exct, dzRL_inib, dzRL_complex, dzRL_nosig, ...
%         dzRL_exct_UniMod_1puf, dzRL_exct_UniMod_2del, dzRL_exct_UniMod_3res, ...
%         dzRL_exct_biMod_12, dzRL_exct_biMod_13, dzRL_exct_biMod_23, ...
%         dzRL_exct_triMod_123);
%     save('Tcombo.mat','Tcombo')
%     disp('Tcombo SAVED'); Tcombo(1,:)
%     