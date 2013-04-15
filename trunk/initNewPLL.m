DKest = [(pi)^2 0 0; 0 (2*pi*10)^2 0; 0 0 25];
for j1 = 1:30
    DKextr = Fp*DKest*Fp' + [0 0 0; 0 0 0; 0 0 Dksi];
    DKest = inv(inv(DKextr) + [1/DmeasNewP 0 0;0 0 0;0 0 0]);
end
KalmanNewPLL.setKoeff([DKest(1,1); DKest(2,1); DKest(3,1) ] / DmeasNewP);
KalmanNewPLL.calcBand();
KResNewPLL.DteorP = DKest(1,1);
KResNewPLL.DteorW = DKest(2,2);
KResNewPLL.ud = zeros(1, K);
KResNewPLL.X = cell(1, 2);
KResNewPLL.X{1} = nan(1, K);
KResNewPLL.X{2} = nan(1, K);
KResNewPLL.ErrX1 = zeros(1, K);