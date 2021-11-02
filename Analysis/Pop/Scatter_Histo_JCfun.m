function Scatter_Histo_JCfun(T2)

Fr_Delta_Up = diff([T2.Fr_BLE_Hz T2.SMA_MAX_Hz],1,2);
Fr_Delta_Down = diff([T2.Fr_BLE_Hz T2.SMA_MIN_Hz],1,2);
Fr_Delta_Diff =  diff([ abs(Fr_Delta_Down)  Fr_Delta_Up],1,2);
Fr_Delta_Sum =  sum([ abs(Fr_Delta_Down)  Fr_Delta_Up],2);

%% Fig2E
BMIN=1
BMAX=25
XMIN=0;
XMAX=20;
figure,
subplot(3,3,1), histogram(T2.Fr_BLE_Hz,50,'BinLimits',[BMIN,BMAX], 'EdgeColor','k','FaceColor','k'), title('VM BLE Fr'), xlabel('Hz'), ylabel('#cell'); xlim([0 20])
subplot(3,3,2), histogram(T2.SMA_MAX_Hz,50,'BinLimits',[BMIN,BMAX], 'EdgeColor','k','FaceColor','k'), title('VM MAX Fr'), xlabel('Hz'), ylabel('#cell'); xlim([0 20])
subplot(3,3,3), histogram(T2.SMA_MIN_Hz,50,'BinLimits',[BMIN,BMAX], 'EdgeColor','k','FaceColor','k'), title('VM MIN Fr'), xlabel('Hz'), ylabel('#cell'); xlim([0 20])

subplot(3,3,4), scatter(T2.zSMA_MAX_del, T2.zSMA_MAX_res, 5, 'k'), xlabel('delay (Zscore)'), ylabel('response (Zscore)'),
xlim([XMIN XMAX]),ylim([XMIN XMAX]), hold on, plot([0:1:100],[0:1:100]),

subplot(3,3,5), scatter(T2.zSMA_MAX_puf, T2.zSMA_MAX_del, 5, 'k'), xlabel('response (Zscore)'), ylabel('puff (Zscore)'),
xlim([XMIN XMAX]),ylim([XMIN XMAX]), hold on, plot([0:1:100],[0:1:100]),

subplot(3,3,6), scatter(T2.zSMA_MAX_res, T2.zSMA_MAX_puf, 5, 'k'), xlabel('delay (Zscore)'), ylabel('puff (Zscore)'),
xlim([XMIN XMAX]),ylim([XMIN XMAX]), hold on, plot([0:1:100],[0:1:100]),

subplot(3,3,7), scatter(Fr_Delta_Up, abs(Fr_Delta_Down), 5, 'k'),xlim([0 15]),ylim([0 15]), xlabel('Delta Up (Hz)'),ylabel('Delta Down (Hz)'),
hold on, plot([0:1:40],[0:1:40])

subplot(3,3,8), scatter(Fr_Delta_Diff, Fr_Delta_Sum, 5, 'k'),xlim([-15 15]),ylim([0 15]), xlabel('Delta Diff(Hz)'),ylabel('Delta Sum (Hz)'),
hold on, plot([0:1:40],[0:1:40]), hold on, plot(-[0:1:40],[0:1:40])

subplot(3,3,9), scatter(Fr_Delta_Diff, T2.Fr_BLE_Hz, 5, 'k'),xlim([-25 15]),ylim([0 60]), xlabel('Delta Diff (Hz)'),ylabel('Baseline (Hz)'),
hold on, plot([0:1:40],[0:1:40]), hold on, plot(-[0:1:40],[0:1:40])









