clear
close all
clc

Tmod = 30; % Время моделирования

L = 813; % число точек на интервале интегрирования в корреляторе
Tc = 0.02; % Период интегрирования в корреляторе
Td = Tc/L; % интервал дискретизации
f0 = 1/(3*Td); % промежуточная частота

Tf = 0.02; % Период работы фильтров


K = fix(Tmod/Tc); % Число интервалов накопления в корреляторе за время моделирования

stdn = 8; % СКО шума квадратурных сумм
qcno_ist = 45; % с/ш в дБ
qcno = 10^(qcno_ist / 10); % с/ш в разах

A = sqrt(2*qcno*Td)*stdn;

Xist = [0; 0; 0];
Xoporn = [0; 0; 0];
XestFLL = [0; 0];
XextrFLL = [0; 0];
XestNewPLL = [0; 0; 0];
XextrNewPLL = [0; 0; 0];

Ff = [1 Tf
      0 1];

Fc = [1 Tc 0
      0 1  Tc
      0 0  1]; % Переходная матрица для модели процесса
  
  
Fp = [1 Tf 0
      0 1  Tf
      0 0  1]; % Переходная матрица для модели процесса

% Расчет параметров формирующего шума. GLONASS, page 162
alpha = 0.1; % Ширина спектра ускорения, с^-1
std_a = 1; %СКЗ ускорения
S_ksi = 2*(33*std_a)^2 * alpha; %Спектральная плотность формирующего шума
stdIst = sqrt(S_ksi * Tf); %СКО формирующего шума
Dksi = stdIst^2;
nIst = randn(1,K);

DmeasW = (6/(qcno*Tc^3))*(1 + 1/(qcno*Tc));
DmeasNewP = (stdn^2*L/2)/(A*L/2)^2;
SdFLL = (A*L/2)^2 * Tc^2 / 12;
SdNewPLL = (A*L/2);

KalmanFLL = CKalmanEqMesConstK(XestFLL, Ff, Tf);
KalmanNewPLL = CKalmanEqMesConstK(XestNewPLL, Fp, Tf);


initFLL;
initNewPLL;

I = nan(1,K);
Q = nan(1,K);
omega_ist = nan(1,K);
phase_ist = nan(1,K);

IF_phaseEnd = 0;
for k = 1:K
    KalmanFLL.Extrapolate();
    KalmanNewPLL.Extrapolate();
    Xoporn(2) = KalmanFLL.Xextr(1);
    
    IF_phase = IF_phaseEnd + 2*pi*f0*(1:L)*Td;
    IF_phaseEnd = IF_phase(end);
    y = A * cos(IF_phase + Xist(1) + Xist(2)*((1:L)-1)*Td) + 0*stdn*randn(1,L);
    
    I_oporn = cos(IF_phase + Xoporn(1) + Xoporn(2)*((1:L)-1)*Td);
    Q_oporn = sin(IF_phase + Xoporn(1) + Xoporn(2)*((1:L)-1)*Td);
    
    I(k) = y * I_oporn';
    Q(k) = y * Q_oporn';
    
    Ih(k) = - (y .* ((1:L)-1)*Td) * Q_oporn';
    Qh(k) = (y .* ((1:L)-1)*Td) * I_oporn';
    
    KResFLL.ud(k) = I(k)*Ih(k) + Q(k)*Qh(k);
    dPhi = (KalmanNewPLL.Xextr(2) - Xoporn(2))*Tf/2*1 + (KalmanNewPLL.Xextr(1) - Xoporn(1));
    KResNewPLL.ud(k) = -(I(k)*sin(dPhi) + Q(k)*cos(dPhi));
    
    KalmanFLL.Estimate(KResFLL.ud(k)/SdFLL);
    KalmanNewPLL.Estimate(KResNewPLL.ud(k)/SdNewPLL);
    Xoporn = Fc*Xoporn;
    
    omega_ist(k) = Xist(2);
    phase_ist(k) = Xist(1);
    
%     Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % Модель изменения истинного вектора.
    if k > 200
    Xist(1) = pi/6;
    else
    Xist(1) = 0;
    end
    KResFLL.ErrX1(k) = Xist(2) - KalmanFLL.Xextr(1);
    KResNewPLL.ErrX2(k) = Xist(2) - KalmanNewPLL.Xextr(2);
    KResNewPLL.ErrX1(k) = Xist(1) - KalmanNewPLL.Xextr(1);
    
    KResFLL.X{1}(k) = KalmanFLL.Xest(1);
    KResNewPLL.X{1}(k) = KalmanNewPLL.Xest(1);
    KResNewPLL.X{2}(k) = KalmanNewPLL.Xest(2);

    if ~mod(k,fix(K/10))
        fprintf('Progress: %.0f%%\n', 100*k/K);
    end

end

t = (1:K)*Tc;
plot(t,KResFLL.X{1},t,omega_ist);
