DKest = [(2*pi*10)^2 0; 0 25];
for j1 = 1:300
    DKest_new(1,1) = DmeasW * (DKest(1,1) + 2*DKest(2,1)*Tf + DKest(2,2)*Tf^2) ./ (DmeasW + DKest(1,1) + 2*DKest(2,1)*Tf + DKest(2,2)*Tf^2);
    DKest_new(2,1) = DmeasW * (DKest(2,1) + DKest(2,2)*Tf) ./  (DmeasW + DKest(1,1) + 2*DKest(2,1)*Tf + DKest(2,2)*Tf^2);
    DKest_new(1,2) = DKest(2,1);
    DKest_new(2,2) = (DmeasW*(DKest(2,2) + Dksi) + (DKest(1,1) + 2*DKest(2,1)*Tf + DKest(2,2)*Tf^2)*(DKest(2,2) + Dksi) - (DKest(2,1) + DKest(2,2)*Tf)^2) ./  (DmeasW + DKest(1,1) + 2*DKest(2,1)*Tf + DKest(2,2)*Tf^2);
    DKest = DKest_new;
end
KalmanFLL.setKoeff([DKest(1,1); DKest(2,1)] / DmeasW);
KalmanFLL.calcBand();
KResFLL.DteorW = DKest(1,1);
KResFLL.ud = zeros(1, K);
KResFLL.X = cell(1, 2);
KResFLL.X{1} = nan(1, K);
KResFLL.X{2} = nan(1, K);
KResFLL.ErrX1 = zeros(1, K);

