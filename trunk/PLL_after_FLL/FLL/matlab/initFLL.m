% DKest = [(2*pi*10)^2 0; 0 1^2];
% H = [1 0]; % это вектор-строка
% G = [0; 1]; % а это вектор-столбец
% for j1 = 1:6000
%     DKextr = Ff*DKest*Ff' + G*G'*Dksi;
%     DKest = inv(inv(DKextr) + H'*H*(1/DekvW));
% end
% KalmanFLL.setKoeff(DKest*H'*(1/DekvW));
% KResFLL.DteorW = DKest(1,1);
% KResFLL.DteorV = DKest(2,2);

KalmanFLL.setKoeffFromSpectr(Sksi, Sw)
KalmanFLL.calcBand();
KResFLL.DteorW = (4*Sksi*Sw^3)^(1/4);
KResFLL.udW = nan(1, K);
KResFLL.X = cell(1, 2);
KResFLL.X{1} = nan(1, K);
KResFLL.X{2} = nan(1, K);
KResFLL.ErrX1 = nan(1, K);
KResFLL.ErrX2 = nan(1, K);

