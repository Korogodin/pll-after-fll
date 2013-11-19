DKest = [(pi)^2 0 0; 0 (2*pi*1)^2 0; 0 0 1];
H = [1 0 0];
G = [0; 0; 1];
for j1 = 1:300
    DKextr = Fp*DKest*Fp' + G*G'*Dksi;
    DKest = inv(inv(DKextr) + H'*H*(1/DmeasPhi));
end
KalmanPLL.setKoeff(DKest*H'*(1/DmeasPhi));
% KalmanPLL.setKoeffFromSpectr(S_ksi, S_meas_phi)
KalmanPLL.calcBand();
KResPLL.DteorPhi = DKest(1,1);
KResPLL.DteorW = DKest(2,2);
% KResPLL.DteorPhi = 2*(S_ksi*S_meas_phi^5)^(1/6);
% KResPLL.DteorW = 3*(S_ksi*S_meas_phi)^(1/2);
KResPLL.udPhi = zeros(1, K);
KResPLL.X = cell(1, 2);
KResPLL.X{1} = nan(1, K);
KResPLL.X{2} = nan(1, K);
KResPLL.ErrX1 = zeros(1, K);
KResPLL.ErrX2 = zeros(1, K);
