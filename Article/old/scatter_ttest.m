%% Scatters pop
close all;
parfig.XlimMAX = 20;

parfig.title = 'PUFF CELL'
Tephys_type= Tephys(logical(Tcombo.PuffCell),:)
pub_fig2_scatter_JCfun(Tephys_type, parfig)

parfig.title = 'DELAY CELL'
Tephys_type= Tephys(logical(Tcombo.Type2Cell),:)
pub_fig2_scatter_JCfun(Tephys_type, parfig)

parfig.title = 'RESPONSE CELL'
Tephys_type= Tephys(logical(Tcombo.RespCell),:)
pub_fig2_scatter_JCfun(Tephys_type, parfig)

parfig.title = 'BOTH CELL'
Tephys_type= Tephys(logical(Tcombo.BothCell),:)
pub_fig2_scatter_JCfun(Tephys_type, parfig)

