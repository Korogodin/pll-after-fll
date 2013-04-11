clear
close all
clc

global K

Tmod = 30; % Время моделирования
Tf = 0.001; % Период работы фильтров
Tc = 0.0005; % Период интегрирования в корреляторе
K = fix(Tmod/Tc); % Число интервалов накопления в корреляторе за время моделирования
stdn_IQ = 8; % СКО шума квадратурных сумм

Ff = [1  Tf
      0  1]; % Переходная матрица для модели частоты (в темпе фильтра)

Fp = [1 Tf 0
      0 1  Tf
      0 0  1]; % Переходная матрица для модели фазы (в темпе фильтра)

Fc = [1 Tc 0
      0 1  Tc
      0 0  1]; % Переходная матрица для модели фазы

Fincorr = [1 Tc
           0 1]; % Переходная матрица набега фазы в корреляторе

std_a = [10:10:100];
CKO_W_PLL = nan(1,length(std_a));
CKO_W_FLL = nan(1,length(std_a));
CKO_Phi_PLL = nan(1,length(std_a));
for j = 1:length(std_a)
    fprintf('SKZ = %.0f\n',std_a(j))
    MemoryAlloc; % Резервирование памяти
    
    qcno_ist = 45; % SNR на каждом интервале, дБГц
    qcno = 10.^(qcno_ist/10);
    A_IQ = stdn_IQ .* sqrt(2 * qcno * Tc); % Амплитуда квадратурных компонент
    
    Xist = [0; 0; 0]; % Истинный вектор состояния
    XextrFLL = [0; 0]; % Вектор экстраполяций, начальное состояние
    XestFLL = XextrFLL;
    XextrPLL = [0; 0; 0]; % Вектор экстраполяций в ФАП
    XestPLL = XextrPLL;
    Xcorr = [0; 0];
    
    % Расчет параметров формирующего шума. GLONASS, page 162
    alpha = 0.1; % Ширина спектра ускорения, с^-1
%     std_a = 40; %СКЗ ускорения
    S_ksi = 2*(33*std_a(j))^2 * alpha; %Спектральная плотность формирующего шума
    stdIst = sqrt(S_ksi * Tf); %СКО формирующего шума
    nIst = randn(1,K);
    
    % Расчет СКО эквивалентных шумов
    stdEqNoise = stdn_IQ/A_IQ;
    stdEqNoiseW = sqrt((1/(qcno*Tc^3))*(1+1/(2*qcno*Tc)));
    
    CKO_W_teor(j) = (4*S_ksi*(Tf*stdEqNoiseW^2)^3)^(1/8) / (2*pi);
    
    % ЧАП
%     HFLL = 5.5; % Hz, полоса ЧАП
    KFLL = nan(2,1); % Вектор-столбец коэффициентов фильтра
    KFLL(2) = (stdIst/(stdEqNoiseW*Tf));
%     KFLL(2) = 2*16/9*HFLL^2; % Коэффициенты непрерывной системы в установившемся режиме
    KFLL(1) = sqrt(2*KFLL(2));
    KFLL = KFLL*Tf; % Переход к коэффициентам дискретной системы
    
    % ФАП
    KPLL = nan(3,1); % Вектор-столбец коэффициентов фильтра
    %     KPLL(3) = (1.2*HPLL(j))^3; % Коэффициенты непрерывной системы в установившемся режиме
    KPLL(3) = (stdIst/(stdEqNoise*Tf));
    KPLL(2) = 2*(KPLL(3))^(2/3);
    KPLL(1) = 2*(KPLL(3))^(1/3);
    KPLL = KPLL*Tf; % Переход к коэффициентам дискретной системы
    
    nI = 1*stdn_IQ.*randn(1,K); % I-comp noise
    nQ = 1*stdn_IQ.*randn(1,K); % Q-comp noise
    
    w = 0; Isum = 0; Qsum = 0; Iold = 1; Qold = 0;
    Ik1 = 1; Qk1 = 0;
    for k = 1:K
        
        % Расчет стат.эквивалентов корреляционных сумм
        Xcorr = Fincorr * Xcorr; % Набег фазы в корреляторе к концу накопления
        EpsPhi(k) = Xist(1) - Xcorr(1);
        EpsW(k) = Xist(2) - Xcorr(2);
        
        A_IQ_eff(k) = A_IQ*sinc(EpsW(k)*Tc/2 /pi);
        
        mI = A_IQ_eff(k) * cos(EpsW(k)*Tc/2 + EpsPhi(k));
        mQ = - A_IQ_eff(k) * sin(EpsW(k)*Tc/2 + EpsPhi(k));
        I(k) = mI + nI(k);
        Q(k) = mQ + nQ(k);
        Isum = Isum + I(k);
        Qsum = Qsum + Q(k);
        
        w = w + 1;
        if w == fix(Tf/Tc)
            Ik2 = I(k); Qk2 = Q(k);
            % Фазовый дискриминатор
            UdPLL(k) = -(Isum*sin((XextrPLL(2) - Xcorr(2))*Tf/2*1 + (XextrPLL(1) - Xcorr(1))) + Qsum*cos((XextrPLL(2) - Xcorr(2))*Tf/2*1 + (XextrPLL(1) - Xcorr(1))));
            UdPLL_mean(k) = A_IQ*sinc((Xist(2) - Xcorr(2))*Tf/2 /pi)*sin((Xist(2) - XextrPLL(2))*Tf/2*1 + (Xist(1) - XextrPLL(1)));
            SdPLL = A_IQ; % Крутизна ФД
            XestPLL = XextrPLL + KPLL*UdPLL(k)/SdPLL;  % Вектор оценок на очередной интервал
            XextrPLL = Fp*XestPLL;                % Экстраполяция на следующий интервал
            
            UdFLL(k) = (Ik2*Qk1 - Qk2*Ik1); % Частотный дискриминатор
            SdFLL = Tc * A_IQ^2; % Крутизна ЧД
            XestFLL = XextrFLL + KFLL*UdFLL(k)/SdFLL;  % Вектор оценок на очередной интервал фильтра
            XextrFLL = Ff*XestFLL;             % Экстраполяция на следующий интервал
            
            Xcorr(2) = XextrFLL(1); % Управляем частотой в корреляторе из ЧАП
            w = 0;
            Iold = Isum; Isum = 0;
            Qold = Qsum; Qsum = 0;
            Ik1 = Ik2; Qk1 = Qk2;
        elseif w == 1
            Ik1 = I(k); Qk1 = Q(k);
        end
        
        Err_W_FLL(k) = Xist(2) - XestFLL(1);
        Err_W_PLL(k) = Xist(2) - XestPLL(2);
        Err_Phi_PLL(k) = Xist(1) - XestPLL(1);
        
        Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % Модель изменения истинного вектора.
        
        if ~mod(k,fix(K/10))
            fprintf('Progress: %.0f%%\n', 100*k/K);
        end
    end
    CKO_W_PLL(j) = sqrt(mean((Err_W_PLL./(2*pi)).^2));
    CKO_W_FLL(j) = sqrt(mean((Err_W_FLL./(2*pi)).^2));
    CKO_Phi_PLL(j) = sqrt(mean((Err_Phi_PLL./(2*pi)*360).^2));
end

D_etta_phi = mean((UdPLL - UdPLL_mean).^2); % флуктуационная характеристика дискриминатора ФАП

t = (1:K)*Tc;

plot(t, Err_W_FLL /2 /pi);
title('Freq error ')
ylabel('\deltaf, Hz');
xlabel('t, s');

% figure(2)
% subplot(2,1,1)
% plot(HPLL,CKO_Phi_PLL);
% title('\Phi RMSE, \Omega RMSE')
% ylabel('\Phi RMSE');
% xlabel('HPLL, Hz');
% subplot(2,1,2)
% plot(HPLL,CKO_W_PLL);
% ylabel('\Omega RMSE');
% xlabel('HPLL, Hz');


