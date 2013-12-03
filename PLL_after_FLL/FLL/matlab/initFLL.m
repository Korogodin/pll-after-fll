DKest = [(2*pi*1)^2 0; 0 1];
H = [1 0]; % это вектор-строка
G = [0; 1]; % а это вектор-столбец
for j1 = 1:300
    DKextr = Ff*DKest*Ff' + G*G'*Dksi;
    DKest = inv(inv(DKextr) + H'*H*(1/DekvW));
end
KalmanFLL.setKoeff(DKest*H'*(1/DekvW));
KalmanFLL.calcBand();
KResFLL.DteorW = DKest(1,1);
% KResFLL.DteorV = DKest(2,2);
KResFLL.udW = nan(1, K);
KResFLL.X = cell(1, 2);
KResFLL.X{1} = nan(1, K);
KResFLL.X{2} = nan(1, K);
KResFLL.ErrX1 = nan(1, K);
KResFLL.ErrX2 = nan(1, K);

