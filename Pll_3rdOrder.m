clear
close all
clc

Tmod = 30; % Время моделирования
Tf = 0.001; % Период работы фильтров
Tc = 0.001; % Период интегрирования в корреляторе
L = fix(Tmod/Tc); % Число интервалов накопления в корреляторе за время моделирования

Hpll = 245; % Hz, полоса ФАП
qcno_ist = 45*ones(1,L); % SNR на каждом интервале, дБГц

Xextr = [0; 0; 0]; % Вектор экстраполяций, начальное состояние
Xest = Xextr;

Fc = [1 Tc 0
    0 1  Tc
    0 0  1]; % Переходная матрица для модели частоты (в темпе коррелятора)

Fincorr = [1 Tc 0
    0 1  0
    0 0  1]; % Переходная матрица набега фазы в корреляторе

Ko = nan(3,1); % Вектор-столбец коэффициентов фильтра
Ko(3) = (1.2*Hpll)^3; % Коэффициенты непрерывной системы в установившемся режиме
Ko(2) = 2*(Ko(3))^(2/3);
Ko(1) = 2*(Ko(3))^(1/3);
Ko = Ko*Tf; % Переход к коэффициентам дискретной системы

Xist = [0; 0; 0]; % Истинный вектор состояния
% Расчет параметров формирующего шума. GLONASS, page 162
alpha = 0.1; % Ширина спектра ускорения, с^-1
std_a = 40; %СКЗ ускорения
S_ksi = 2*(33*std_a)^2 * alpha; %Спектральная плотность формирующего шума
stdIst = Tc*sqrt(S_ksi / Tc); %СКО формирующего шума
nIst = randn(1,L);

A_IQ = nan(1,L); % Амплитуда квадратурных компонент
A_IQ_eff = nan(1,L); % Онаже, включая функцию sinc
I = nan(1,L);
Q = nan(1,L);
EpsPhi = nan(1, L); % Рассогласование по фазе (ист - экстр)
EpsW = nan(1, L);   % Рассогласование по частоте (ист - экстр)
ErrW = nan(1,L);    % Ошибка оценки частоты в системе ЧАП
ErrPhi = nan(1,L);
omega = nan(1,L);   % Вектор оценки частоты в системе ЧАП
omega_ist = nan(1,L); % Реальная частота
phase = nan(1,L);   % ???
phase_ist = nan(1,L); % Реальная фаза
Ud = nan(1,L);      % Отсчеты дискриминатора в ЧАП

stdn_IQ = ones(1,L)*8; % СКО шума квадратурных сумм

nI = 1*stdn_IQ.*randn(1,L); % I-comp noise
nQ = 1*stdn_IQ.*randn(1,L); % Q-comp noise

w = 0; Isum = 0; Qsum = 0;
for k = 1:L
    
    % Расчет стат.эквивалентов корреляционных сумм
    EpsPhi(k) = Xist(1) - Xextr(1);
    EpsW(k) = Xist(2) - Xextr(2);
    
    qcno = 10.^(qcno_ist(k)/10);
    A_IQ(k) = stdn_IQ(k) .* sqrt(2 * qcno * Tc);
    A_IQ_eff(k) = A_IQ(k)*sinc(EpsW(k)*Tc/2 /pi);
    
    mI = A_IQ_eff(k) * cos(EpsW(k)*Tc/2 + EpsPhi(k));
    mQ = - A_IQ_eff(k) * sin(EpsW(k)*Tc/2 + EpsPhi(k));
    I(k) = mI + 0*nI(k);
    Q(k) = mQ + 0*nQ(k);
    
    Isum = Isum + I(k);
    Qsum = Qsum + Q(k);
    
    Xextr = Fincorr * Xextr; % Набег фазы в корреляторе к концу накопления
    w = w + 1;
    if w == fix(Tf/Tc)
        
        EpsPhi(k) = Xist(1) - Xextr(1);
        % Фазовый дискриминатор
        Ud(k) = -atan(Q(k)/I(k));
        Sd = 1; % Крутизна ФД
        Xest = Xextr + Ko*Ud(k)/Sd;  % Вектор оценок на очередной интервал
        Xextr = Fc*Xest;                % Экстраполяция на следующий интервал
        
        w = 0;
        Isum = 0;
        Qsum = 0;
    end
    
    ErrW(k) = Xist(2) - Xest(2);
    ErrPhi(k) = Xist(1) - Xest(1);
    omega(k) = Xest(2);
    omega_ist(k) = Xist(2);
    phase_ist(k) = Xist(1);
    phase(k) = Xest(1);
    
    if ~mod(k,fix(L/10))
        fprintf('Progress: %.0f%%\n', 100*k/L);
    end
    Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % Модель изменения истинного вектора.
end
t = (1:L)*Tc;
hF = 0;
hF = figure(hF + 1);
subplot(2,1,1)
plot(t, omega / 2 / pi, t, omega_ist / 2 / pi);
title('Freq')
ylabel('\omega, Hz');
xlabel('t, s');
subplot(2,1,2)
plot(t, ErrW / 2 / pi);
ylabel('Err_{\omega}, Hz');
xlabel('t, s');

hF = figure(hF + 1);
subplot(2,1,1)
plot(EpsW, Ud, '.')
title('Diskriminators')
ylabel('Ud')
xlabel('\epsilon_{\omega}, rad/s');
subplot(2,1,2)
plot(mod(EpsPhi, 2*pi), Ud, '.')
ylabel('Ud')
xlabel('\epsilon_{\phi}, rad');

hF = figure(hF + 1);
subplot(2,1,1)
plot(t, phase / 2 / pi, t, phase_ist / 2 / pi);
title('Phase')
ylabel('\phi, cycles');
xlabel('t, s');
subplot(2,1,2)
plot(t, phase_ist / 2 / pi - phase / 2 / pi);
ylabel('Err\phi, cycles');
xlabel('t, s');

