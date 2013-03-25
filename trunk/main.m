clear
close all
clc

Tmod = 30; % Время моделирования
Tf = 0.001; % Период работы фильтров 
Tc = 0.001; % Период интегрирования в корреляторе
K = fix(Tmod/Tc); 

qcno_ist = 45*ones(1,K); % // SNR

H = 5.5; % Hz, полоса

Xextr = [0; 0; 0]; % Вектор экстраполяций
Xest = Xextr;

XextrP = [0; 0; 0]; % Вектор экстраполяций
XestP = XextrP;

Ff = [1 0  0
      0 1  Tf
      0 0  1]; % Переходная матрица для модели частоты (в темпе фильтра)
  
Fc = [1 Tc 0
      0 1  Tc
      0 0  1]; % Переходная матрица для модели частоты (в темпе коррелятора)
  
Fincorr =  [1 Tc 0
            0 1  0
            0 0  1]; % Переходная матрица набега фазы в корреляторе

% ЧАП
Ko = nan(3,1); % Вектор-столбец коэффициентов фильтра
Ko(3) = 2*16/9*H^2; % Коэффициенты непрерывной системы в установившемся режиме
Ko(2) = sqrt(2*Ko(3));
Ko(1) = 0;    % Это не баг, это фича! Из-за этого нуля система, на самом деле, - второго порядка
Ko = Ko*Tf; % Переход к коэффициентам дискретной системы

% ФАП
HP = 45; % Hz, полоса
KoP = nan(3,1); % Вектор-столбец коэффициентов фильтра
KoP(3) = (1.2*HP)^3; % Коэффициенты непрерывной системы в установившемся режиме
KoP(2) = 2*(KoP(3))^(2/3);
KoP(1) = 2*(KoP(3))^(1/3);
KoP = KoP*Tf; % Переход к коэффициентам дискретной системы

Xist = [0; 0; 0]; % Истинный вектор состояния
% S_ksi = 3.48e5; % See GLONASS, page 162
alpha = 0.1;
std_a = 40;
S_ksi = 2*(33*std_a)^2 * alpha;
stdIst = Tc*sqrt(S_ksi / Tc);
nIst = randn(1,K); 

stdn_IQ = ones(1,K)*8; % СКО шума квадратурных сумм

A_IQ = nan(1,K); % // Memory allocation
A_IQ_eff = nan(1,K);

I = nan(1,K); % // Memory allocation
Q = nan(1,K);

EpsPhi = nan(1, K); % // Memory allocation
EpsPhiP = nan(1,K);
EpsW = nan(1, K);

ErrW = nan(1,K);
ErrWP = nan(1,K);

omega = nan(1,K);
omega_ist = nan(1,K);
omegaP = nan(1,K);

phase = nan(1,K);
phase_ist = nan(1,K);
phaseP = nan(1,K);

Ud = nan(1,K);
UdP = nan(1,K);

nI = 1*stdn_IQ.*randn(1,K); % // I-comp noise
nQ = 1*stdn_IQ.*randn(1,K); % // Q-comp noise

w = 0; Isum = 0; Qsum = 0; Iold = 1; Qold = 0;
for k = 1:K

    % // Расчет стат.эквивалентов корреляционных сумм
    EpsPhi(k) = Xist(1) - Xextr(1);
    EpsW(k) = Xist(2) - Xextr(2);

    qcno = 10.^(qcno_ist(k)/10);
    A_IQ(k) = stdn_IQ(k) .* sqrt(2 * qcno * Tc);
    A_IQ_eff(k) = A_IQ(k)*sinc(EpsW(k)*Tc/2 /pi);

    mI = A_IQ_eff(k) * cos(EpsW(k)*Tc/2 + EpsPhi(k));
    mQ = - A_IQ_eff(k) * sin(EpsW(k)*Tc/2 + EpsPhi(k));
    I(k) = mI + nI(k);
    Q(k) = mQ + nQ(k);
    Isum = Isum + I(k);
    Qsum = Qsum + Q(k);

    Xextr = Fincorr * Xextr; % Набег фазы в корреляторе к концу накопления    
    
    w = w + 1;
    if w == fix(Tf/Tc)                   
       
        EpsPhiP(k) = Xist(1) - XextrP(1);
        UdP(k) = -(I(k)*sin((XextrP(2) - Xextr(2))*Tf/2*1 + (XextrP(1) - Xextr(1))) + Q(k)*cos((XextrP(2) - Xextr(2))*Tf/2*1 + (XextrP(1) - Xextr(1))));
        SdP = A_IQ(k);
%         KoP = [0.5; 0; 0];
        XestP = XextrP + KoP*UdP(k)/SdP;                 % Вектор оценок на очередной интервал фильтра
        XextrP = Fc*XestP;                                      % Экстраполяция на следующий интервал
%         XextrP = XestP;

        %     Ud = f(I(k), Q(k), I(k-1), Q(k-1), ...);      % Дискриминатор
        Ud(k) = (I(k)*Qold - Q(k)*Iold);
        %     Sd = f(A_IQ);             % Крутизна дискриминационной характеристики
        Sd = Tc*(A_IQ(k)*Tf/Tc)^2 * 1.3;
        Xest = Xextr + Ko*Ud(k)/Sd;                 % Вектор оценок на очередной интервал фильтра
        Xextr = Ff*Xest;                         % Экстраполяция на следующий интервал
                
        w = 0;
        Iold = Isum; Isum = 0;
        Qold = Qsum; Qsum = 0;
    end

    ErrW(k) = Xist(2) - Xest(2);
    ErrWP(k) = Xist(2) - XestP(2);
    
    omega(k) = Xest(2);
    omegaP(k) = XestP(2);
    omega_ist(k) = Xist(2);
    
    phase_ist(k) = Xist(1);
    phaseP(k) = XestP(1);
    
%     if k == 12000
%         Xist(2) = 100;
%     end
%     if k == 3000
%         Xist(1) = Xist(1) + 2*pi/3;
%     end
        Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % Здесь может быть любая другая модель изменения истинного вектора состояния

    if ~mod(k,fix(K/10))
        fprintf('Progress: %.0f%%\n', 100*k/K);
    end
end

t = (1:K)*Tc;

hF = 0;
hF = figure(hF + 1);
subplot(2,1,1)
plot(t, omega / 2 / pi, t, omegaP / 2 / pi, t, omega_ist / 2 / pi);
title('Freq')
ylabel('\omega, Hz');
xlabel('t, s');
subplot(2,1,2)
plot(t, ErrW / 2 / pi, t, ErrWP / 2 / pi);
ylabel('Err_{\omega}, Hz');
xlabel('t, s');

hF = figure(hF + 1);
subplot(2,1,1)
plot(EpsW, Ud, '.')
title('Diskriminators')
ylabel('Ud')
xlabel('\epsilon_{\omega}, rad/s');
subplot(2,1,2)
plot(mod(EpsPhiP, 2*pi), UdP, '.')
ylabel('UdP')
xlabel('\epsilon_{\phi}, rad');

hF = figure(hF + 1);
subplot(2,1,1)
plot(t, phase / 2 / pi, t, phaseP / 2 / pi, t, phase_ist / 2 / pi);
title('Phase')
ylabel('\phi, cycles');
xlabel('t, s');
subplot(2,1,2)
plot(t, phase_ist / 2 / pi - phaseP / 2 / pi);
ylabel('Err\phi, cycles');
xlabel('t, s');
