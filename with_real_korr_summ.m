clear
close all
clc

Tmod = 0.01; % Время моделирования

L = 113; % число точек на интервале интегрирования в корреляторе
Tc = 0.001; % Период интегрирования в корреляторе
Td = Tc/L; % интервал дискретизации
f0 = 1/(3*Td); % промежуточная частота

Tf = 0.001; % Период работы фильтров
Tc = 0.001; % Период интегрирования в корреляторе

K = fix(Tmod/Tc); % Число интервалов накопления в корреляторе за время моделирования

stdn = 8; % СКО шума квадратурных сумм
qcno_ist = 45; % с/ш в дБ
qcno = 10^(qcno_ist / 10); % с/ш в разах

A = sqrt(2*qcno*Td)*stdn;

Xist = [0; 0; 0];
Xoporn = [0; 0; 0];

Fc = [1 Tc 0
    0 1  Tc
    0 0  1]; % Переходная матрица для модели процесса

% Расчет параметров формирующего шума. GLONASS, page 162
alpha = 0.1; % Ширина спектра ускорения, с^-1
std_a = 40; %СКЗ ускорения
S_ksi = 2*(33*std_a)^2 * alpha; %Спектральная плотность формирующего шума
stdIst = sqrt(S_ksi * Tf); %СКО формирующего шума
nIst = randn(1,K);


I = nan(1,K);
Q = nan(1,K);

IF_phaseEnd = 0;
for k = 1:K
    IF_phase = IF_phaseEnd + 2*pi*f0*(1:L)*Td;
    IF_phaseEnd = IF_phase(end);
    y = A * cos(IF_phase + Xist(1) + Xist(2)*((1:L)-1)*Td) + stdn*randn(1,L);
    
    I_oporn = cos(IF_phase + Xoporn(1) + Xoporn(2)*((1:L)-1)*Td);
    Q_oporn = sin(IF_phase + Xoporn(1) + Xoporn(2)*((1:L)-1)*Td);
    
    I(k) = y * I_oporn';
    Q(k) = y * Q_oporn';
    
    Ih(k) = - (y * ((1:L)-1)*Td) * Q_oporn';
    Qh(k) = (y * ((1:L)-1)*Td) * I_oporn';
    
    Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % Модель изменения истинного вектора.
    
end

