DKest = [(pi)^2 0 0; 0 (2*pi*1)^2 0; 0 0 1];
for j1 = 1:300
    DKextr = Fp*DKest*Fp' + [0 0 0; 0 0 0; 0 0 Dksi];
    DKest = inv(inv(DKextr) + [1/DmeasPLL3rdOrder 0 0;0 0 0;0 0 0]);
end
KalmanPLL3rdOrder.setKoeff([DKest(1,1); DKest(2,1); DKest(3,1) ] / DmeasPLL3rdOrder);
KalmanPLL3rdOrder.calcBand();
KResPLL3rdOrder.DteorP = DKest(1,1);
KResPLL3rdOrder.DteorW = DKest(2,2);
KResPLL3rdOrder.udPLL3rdOrder = zeros(1, K);
KResPLL3rdOrder.X = cell(1, 2);
KResPLL3rdOrder.X{1} = nan(1, K);
KResPLL3rdOrder.X{2} = nan(1, K);
KResPLL3rdOrder.ErrX1 = zeros(1, K);
KResPLL3rdOrder.ErrX2 = zeros(1, K);
