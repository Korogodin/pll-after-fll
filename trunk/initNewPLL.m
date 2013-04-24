DKest = [(pi)^2 0 0; 0 (2*pi*1)^2 0; 0 0 1];
H = [1 0 0; 0 1 0];
D = [DmeasNewP 0; 0 DmeasW];
for j1 = 1:300
    DKextr = Fp*DKest*Fp' + [0 0 0; 0 0 0; 0 0 Dksi];
    DKest = inv(inv(DKextr) + H'*inv(D)*H);
end
KalmanNewPLL.setKoeff(DKest*H'*inv(D));
KalmanNewPLL.calcBand();
KResNewPLL.DteorP = DKest(1,1);
KResNewPLL.DteorW = DKest(2,2);
KResNewPLL.udPhi = zeros(1, K);
KResNewPLL.udW = zeros(1, K);
KResNewPLL.X = cell(1, 2);
KResNewPLL.X{1} = nan(1, K);
KResNewPLL.X{2} = nan(1, K);
KResNewPLL.ErrX1 = zeros(1, K);
KResNewPLL.ErrX2 = zeros(1, K);
