clear
close all
clc

Tmod = 30; % Время моделирования
Tf = 0.001; % Период работы фильтров
Tc = 0.001; % Период интегрирования в корреляторе
N = fix(Tmod/Tc); % Число интервалов накопления в корреляторе за время моделирования

qcno_ist = 37; % SNR на каждом интервале, дБГц
stdn_IQ = 8; % СКО шума квадратурных сумм
qcno = 10.^(qcno_ist/10); A_IQ = stdn_IQ .* sqrt(2 * qcno * Tc); % Амплитуда квадратурных компонент

Xist = [0; 0; 0]; % Истинный вектор состояния
Xextr = [0; 0; 0]; % Вектор экстраполяций, начальное состояние
Xest = Xextr;

Fc = [1 Tc 0
        0 1  Tc
        0 0  1]; % Переходная матрица для модели процесса фазы (в темпе коррелятора)

Fincorr = [1 Tc 0
              0 1  0
              0 0  1]; % Переходная матрица набега фазы в корреляторе

Hpll = 25; % Hz, полоса ФАП
Ko = nan(3,1); % Вектор-столбец коэффициентов фильтра
Ko(3) = (1.2*Hpll)^3; % Коэффициенты непрерывной системы в установившемся режиме
Ko(2) = 2*(Ko(3))^(2/3);
Ko(1) = 2*(Ko(3))^(1/3);
Ko = Ko*Tf; % Переход к коэффициентам дискретной системы

% Расчет параметров формирующего шума. GLONASS, page 162
alpha = 0.1; % Ширина спектра ускорения, с^-1
std_a = 20; %СКЗ ускорения
S_ksi = 2*(33*std_a)^2 * alpha; %Спектральная плотность формирующего шума
stdIst = Tc*sqrt(S_ksi / Tc); %СКО формирующего шума
nIst = randn(1,N);

A_IQ_eff = nan(1,N); % % Амплитуда квадратурных компонент с учетом функции sinc
I = nan(1,N);
Q = nan(1,N);
EpsPhi = nan(1, N); % Рассогласование по фазе (ист - экстр)
EpsW = nan(1, N);   % Рассогласование по частоте (ист - экстр)
phi_est = nan(1, N);
omega_est = nan(1, N);
phi_ist = nan(1, N);
omega_ist = nan(1, N);
Ud = nan(1,N);      % Отсчеты дискриминатора в ЧАП

nI = 1*stdn_IQ.*randn(1,N); % I-comp noise
nQ = 1*stdn_IQ.*randn(1,N); % Q-comp noise

w = 0; Isum = 0; Qsum = 0;
for n = 1:N
    
    % Расчет стат.эквивалентов корреляционных сумм
    EpsPhi(n) = Xist(1) - Xextr(1);
    EpsW(n) = Xist(2) - Xextr(2);
    A_IQ_eff(n) = A_IQ*sinc(EpsW(n)*Tc/2 /pi);
    
    mI = A_IQ_eff(n) * cos(EpsW(n)*Tc/2 + EpsPhi(n));
    mQ = - A_IQ_eff(n) * sin(EpsW(n)*Tc/2 + EpsPhi(n));
    I(n) = mI + nI(n);
    Q(n) = mQ + nQ(n);
    
    Isum = Isum + I(n);
    Qsum = Qsum + Q(n);
    
    w = w + 1;
    if w == fix(Tf/Tc)
        
        % Фазовый дискриминатор
        Ud(n) = -atan(Q(n)/I(n));
        Sd = 1; % Крутизна ФД
        Xest = Xextr + Ko*Ud(n)/Sd;  % Вектор оценок на очередной интервал
        Xextr = Fc*Xest;                % Экстраполяция на следующий интервал
        
        w = 0;
        Isum = 0;
        Qsum = 0;
    end
    
    omega_est(n) = Xest(2);
    omega_ist(n) = Xist(2);
    phi_ist(n) = Xist(1);
    phi_est(n) = Xest(1);
    
    if ~mod(n,fix(N/10))
        fprintf('Progress: %.0f%%\n', 100*n/N);
    end
    Xist = Fc*Xist + [0; 0; 1]*nIst(n)*stdIst; % Модель изменения истинного вектора.
end

t = (1:N)*Tc;
ErrW = omega_ist - omega_est;
ErrPhi = phi_ist - phi_est;
RMSE_phi =sqrt(mean(ErrPhi.^2)) / 2/pi;
RMSE_omega = sqrt(mean(ErrW.^2)) / 2/pi;
fprintf('RMSE_phi = %f deg, RMSE_omega = %f Hz \n', 360*RMSE_phi, RMSE_omega);

hF = 0;
hF = figure(hF + 1);
subplot(2,1,1)
plot(t, omega_ist/2/pi, t, omega_est/2/pi);
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
plot(t, phi_ist/ 2 / pi, t, phi_est / 2 / pi);
title('Phase')
ylabel('\phi, cycles');
xlabel('t, s');
subplot(2,1,2)
plot(t, ErrPhi/2/pi);
ylabel('Err\phi, cycles');
xlabel('t, s');


