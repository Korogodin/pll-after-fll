clear 
close all
clc

Tmod = 60; % Время моделирования
Np = 100; %  число прогонов схемы
f0 = 10.3e6; % промежуточная частота
Td = 1/(3.3712*f0); % интервал дискретизации
Tc = 0.005; % Период интегрирования в корреляторе
Tf = 0.005; % Период работы фильтров
L = fix(Tc/Td); % Число отсчетов на интервале накопления в корреляторе

K = fix(Tmod/Tc); % Число интервалов накопления в корреляторе за время моделирования

stdn = 8; % СКО шума
stdnIQ = sqrt((stdn^2)*L/2);

% Расчет параметров формирующего шума. GLONASS, page 162
alpha = 0.1; % Ширина спектра ускорения, с^-1
std_a = 1; %СКЗ ускорения
S_ksi = 2*(33*std_a)^2 * alpha; %Спектральная плотность формирующего шума
stdIst = sqrt(S_ksi * Tf); %СКО формирующего шума
Dksi = stdIst^2; % Дисперсия формирующего шума

qcno_dB = [10:0.5:50];
Nq = length(qcno_dB);

% KResPLL.CKO_W_TEOR = zeros(1, Nq);
% KResPLL.CKO_W = zeros(1, Nq);
% KResPLL.CKO_Phi_TEOR = zeros(1, Nq);
% KResPLL.CKO_Phi = zeros(1, Nq);
% KResPLL.CKO_PhiSredn = zeros(1, Nq);
% KResPLL.Band = zeros(1, Nq);

Fc = [1 Tc 0
    0 1  Tc
    0 0  1]; % Переходная матрица для модели процесса

Fp = [1 Tf 0
    0 1  Tf
    0 0  1]; % Переходная матрица для ССФ

for j = 1:Nq
    
    qcno = 10^(qcno_dB(j) / 10); % с/ш в разах
    
    A = sqrt(4*qcno*stdn^2*Td);
    Aiq = A*L/2;
    
    % Крутизна и флуктуационная характеристика дискриминатора
    SdPLL = Aiq*erf(sqrt(qcno*Tc)); % крутизна ДХ ФД
    DmeasPhi = stdnIQ^2 / SdPLL^2; % дисперсия экв. шумов ФД (на входе)
    S_meas_phi = DmeasPhi*Tf; % СПМ шумов дискриминатора
    
    DestPhi = 0;
    DestW = 0;
    DestPhiSredn = 0;
    
    for m = 1:Np
        
        Xist = [0; 0; 0];
        Xoporn = [0; 0; 0];
        XestPLL = [0; 0; 0];
        
        nIst = randn(1,K);
        
        KalmanPLL = CKalmanEqMesConstK(XestPLL, Fp, Tf); % вызываем класс фильтра Калмана
        KalmanPLL.Xextr = [0; 0; 0];
        
        % ининциализация переменных / массивов + там же рассчет коэффициентов для
        % фильтра Калмана
        initPLL;
        
%         tau = 20/KalmanPLL.Band; % примерное время установления переходного процесса
%         k_ust = find((1:K)*Tc >=tau, 1, 'first'); % индекс k, когда наступает установившийся режим
%         KResPLL.Band(j) = KalmanPLL.Band;

        I = nan(1,K);
        Q = nan(1,K);
        
        nI = stdnIQ*randn(1,K);
        nQ = stdnIQ*randn(1,K);
        
        omega_ist = nan(1,K);
        phase_ist = nan(1,K);
        
        for k = 1:K
            phaseExtrOld = KalmanPLL.Xextr(1);
            KalmanPLL.Extrapolate();
            Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % Модель изменения истинного вектора.
            %Xoporn(2) = (KalmanPLL.Xextr(1) - phaseExtrOld)/Tf;
            Xoporn(1) = KalmanPLL.Xextr(1);
            Xoporn(2) = KalmanPLL.Xextr(2);
            deltaPhi = Xist(1) - Xoporn(1);
            deltaW = Xist(2) - Xoporn(2);
            
            mI = Aiq*cos(deltaPhi+deltaW*Tc/2)*sinc(deltaW*Tc/2/pi);
            mQ = -Aiq*sin(deltaPhi+deltaW*Tc/2)*sinc(deltaW*Tc/2/pi);
            
            h = sign(randn(1));
            
            I(k) = h*mI + 1*nI(k);
            Q(k) = h*mQ + 1*nQ(k);
            
            KResPLL.udPhi(k) = -sign(I(k)) * Q(k);
%           KResPLL.udPhi(k) = SdPLL*(Xist(1)-KalmanPLL.Xextr(1)) + stdnIQ*randn(1,1);
            
            KalmanPLL.Estimate(KResPLL.udPhi(k)/SdPLL);
            
%             KResPLL.X{1}(k) = KalmanPLL.Xest(1);
%             KResPLL.X{2}(k) = KalmanPLL.Xest(2);
%             
%             phase_ist(k) = Xist(1);
%             omega_ist(k) = Xist(2);
            
%             subplot(3,1,1)
%             plot((1:K)*Tc, [phase_ist*180/pi; KResPLL.X{1}*180/pi]);
%             title('\phi_{istinnaya} \phi_{ocenka}')
%             subplot(3,1,2)
%             plot((1:K)*Tc, KResPLL.ErrX1*180/pi);
%             title('error \phi')
%             subplot(3,1,3)
%             plot((1:K)*Tc, KResPLL.udPhi/SdPLL);
%             title('ud\phi')
%             drawnow
            
%             if k*Tc >= tau
                KResPLL.ErrX1(k) = Xist(1) - KalmanPLL.Xest(1);
                KResPLL.ErrX2(k) = Xist(2) - KalmanPLL.Xest(2);
                KResPLL.ErrSrednPhase(k) = (2*(Xist(1)-KalmanPLL.Xest(1))+(Xist(2)-KalmanPLL.Xest(2))*Tf/2)/2;
%             end
            
            Xoporn(1) = Xoporn(1) + Xoporn(2)*Tc;
                                    
        end
        
        DestPhi = mean(KResPLL.ErrX1.^2) + 0*DestPhi;
        DestW = mean(KResPLL.ErrX2.^2) + 0*DestW;
        DestPhiSredn = mean(KResPLL.ErrSrednPhase.^2) + 0*DestPhiSredn;
        
        save_statistic;
        
        if ~mod(m,fix(Np/10))
            fprintf('Progress po Np: %.2f%%\n', m*100/Np);
        end
        clear('KalmanPLL');
    end
    
%     KResPLL.CKO_Phi_TEOR(j) = sqrt(KResPLL.DteorPhi)*180/pi;
%     KResPLL.CKO_W_TEOR(j) = sqrt(KResPLL.DteorW)/2 /pi;
%     KResPLL.CKO_Phi(j) = sqrt(DestPhi/Np)*180/pi;
%     KResPLL.CKO_PhiSredn(j) = sqrt(DestPhiSredn/Np)*180/pi;
%     KResPLL.CKO_W(j) = sqrt(DestW/Np)/2/pi;
   
    fprintf('Progress po q_cno: %.2f%%\n', j*100/Nq)
end

t = (1:K)*Tc;

% plot(qcno_dB, [KResPLL.CKO_Phi; KResPLL.CKO_PhiSredn; KResPLL.CKO_Phi_TEOR]);
% title('RMSE \phi');
% figure
% plot(qcno_dB, [KResPLL.CKO_W; KResPLL.CKO_W_TEOR]);
% title('RMSE \omega');
% figure
% plot(qcno_dB,(KResPLL.CKO_Phi./KResPLL.CKO_Phi_TEOR -1)*100);
% title('Procent');
% figure
% plot(qcno_dB, KResPLL.Band);
% title('Bandwidth');


