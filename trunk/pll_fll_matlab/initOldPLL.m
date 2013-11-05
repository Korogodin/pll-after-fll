DKest = [(pi)^2 0 0; 0 (2*pi*1)^2 0; 0 0 1];
for j1 = 1:300
    DKextr = Fp*DKest*Fp' + [0 0 0; 0 0 0; 0 0 Dksi];
    DKest = inv(inv(DKextr) + [1/DmeasOldPLL 0 0;0 0 0;0 0 0]);
end
KalmanOldPLL.setKoeff([DKest(1,1); DKest(2,1); DKest(3,1) ] / DmeasOldPLL);
KalmanOldPLL.calcBand();
KResOldPLL.DteorP = DKest(1,1);
KResOldPLL.DteorW = DKest(2,2);
KResOldPLL.udPhi = zeros(1, K);
KResOldPLL.X = cell(1, 2);
KResOldPLL.X{1} = nan(1, K);
KResOldPLL.X{2} = nan(1, K);
KResOldPLL.ErrX1 = zeros(1, K);
KResOldPLL.ErrX2 = zeros(1, K);
