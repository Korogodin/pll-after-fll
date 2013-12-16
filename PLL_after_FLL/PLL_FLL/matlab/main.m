clear
close all
clc

Tmod = 60; % Время моделирования
Np = 1000; %  число прогонов схемы
f0 = 10.3e6; % промежуточная частота
Td = 1/(3.3712*f0); % интервал дискретизации
Tms = 0.001;
Tc = 0.005; % Период интегрирования в корреляторе
Tf = 0.005; % Период работы фильтров
L = fix(Tc/Td); % Число отсчетов на интервале накопления в корреляторе
Lms = fix(Tms/Td); % Число отсчетов на 1 мс

K = fix(Tmod/Tc); % Число интервалов накопления в корреляторе за время моделирования

stdn = 8; % СКО шума
stdnIQ = sqrt((stdn^2)*L/2);
stdnIQms = sqrt((stdn^2)*Lms/2);

% Расчет параметров формирующего шума. GLONASS, page 162
alpha = 0.1; % Ширина спектра ускорения, с^-1
std_a = 40; %СКЗ ускорения
Sksi = 2*(33*std_a)^2 * alpha; %Спектральная плотность формирующего шума
stdIst = sqrt(Sksi * Tf); %СКО формирующего шума
Dksi = stdIst^2; % Дисперсия формирующего шума

qcno_dB = [28:0.5:32];
Nq = length(qcno_dB);

Ff = [1 Tf
    0 1]; % Переходная матрица для ССЧ

Fp = [1 Tf 0
      0 1  Tf
      0 0  1]; % Переходная матрица для ССФ

Fc = [1 Tc 0
      0 1  Tc
      0 0  1]; % Переходная матрица для модели процесса

for j = 1:Nq
    
    qcno = 10^(qcno_dB(j) / 10); % с/ш в разах
    
    A = sqrt(4*qcno*stdn^2*Td);
    Aiq = A*L/2;
    Aiq_ms = A*Lms/2;
    
    % Дискриминаторы
    % Частотный дискриминатор
    SdFLL = 1/12 * (Aiq)^2  * Tc^2; % крутизна ДХ ЧД
    DekvW = (6/(qcno*Tc^3))*(1 + 1/(qcno*Tc)); % дисперсия экв. шумов ЧД 
    Sw = DekvW*Tf; % СПМ шумов дискриминатора
    
    % Фазовый дискриминатор
    SdPLL = Aiq*erf(sqrt(qcno*Tc)); % крутизна ДХ ФД
    DekvPhi = stdnIQ^2 / SdPLL^2; % дисперсия экв. шумов ФД (на входе) (в нуле дискр. х-ки)
    Sphi = DekvPhi*Tf; % СПМ шумов дискриминатора
    
    for m = 1:Np
        
        Xist = [0; 0; 0];
        XestPLL = [0; 0; 0];
        Xoporn = [0; 0];
        XestFLL = [0; 0];
        
        nIst = randn(1,K);
        
        % Создаем фильтры для систем
        KalmanPLL = CKalmanEqMesConstK(XestPLL, Fp, Tf);
        KalmanPLL.Xextr = [0; 0; 0];
        
        KalmanFLL = CKalmanEqMesConstK(XestFLL, Ff, Tf);
        KalmanFLL.Xextr = [0; 0];
        
        % ининциализация переменных / массивов + там же рассчет коэффициентов для
        % фильтра Калмана
        initFLL;
        initPLL;
        
        I = zeros(1,K);
        Q = zeros(1,K);
        Ih = zeros(1,K);
        Qh = zeros(1,K);
        
        nIms = stdnIQms*randn(5,K);
        nQms = stdnIQms*randn(5,K);
        
        phase_ist = nan(1,K);
        omega_ist = nan(1,K);
                
        for k = 1:K
            %Xoporn(2) = (KalmanPLL.Xextr(1) - phaseExtrOld)/Tf;
            
            Xoporn(2) = KalmanFLL.Xextr(1);
            
            dPhiOporn = Xist(1) - Xoporn(1);
            dWOporn = Xist(2) - Xoporn(2);
            dPhi = Xist(1) - KalmanPLL.Xextr(1);
            dW = Xist(2) - KalmanPLL.Xextr(2);
            dPhiextr = dPhiOporn - dPhi;
            dWextr = dWOporn - dW;
            
            dPhiDovorot = (dWextr)*Tc/2 + (dPhiextr);
            
            h = sign(randn(1));
            
            for n = 1:5
                mIms = Aiq_ms*h*cos(dPhiOporn+(n-1)*dWOporn*Tms+dWOporn*Tms/2)*sinc(dWOporn*Tms/2/pi);
                mQms = -Aiq_ms*h*sin(dPhiOporn+(n-1)*dWOporn*Tms+dWOporn*Tms/2)*sinc(dWOporn*Tms/2/pi);
                Ims = mIms + nIms(n,k);
                Qms = mQms + nQms(n,k);
                I(k) = Ims + I(k);
                Q(k) = Qms + Q(k);
                Ih(k) = -n*Tms*Qms + Ih(k);
                Qh(k) = n*Tms*Ims + Qh(k);
            end
                      
            KResFLL.udW(k) = I(k)*Ih(k) + Q(k)*Qh(k);
            KResPLL.udPhi(k) = - sign((I(k)*cos(dPhiDovorot) - Q(k)*sin(dPhiDovorot))) * (I(k)*sin(dPhiDovorot) + Q(k)*cos(dPhiDovorot));
            
            KalmanFLL.Estimate(KResFLL.udW(k)/SdFLL);
            KalmanPLL.Estimate(KResPLL.udPhi(k)/SdPLL);
            
            KResFLL.ErrX1(k) = Xist(2) - KalmanFLL.Xest(1);
            
            KResPLL.ErrX1(k) = Xist(1) - KalmanPLL.Xest(1);
            KResPLL.ErrX2(k) = Xist(2) - KalmanPLL.Xest(2);
            
%             plotInTime;
                        
            Xoporn(1) = Xoporn(1) + Xoporn(2)*Tc;
            
            KalmanFLL.Extrapolate();
            KalmanPLL.Extrapolate();
            
            Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % Модель изменения истинного вектора.
        end
       
        DestPhiPLL = mean(KResPLL.ErrX1.^2);
        if (DestPhiPLL >= pi^2)
            DestPhiPLL = pi^2;
        end
        DestWPLL = mean(KResPLL.ErrX2.^2);
        
        DestWFLL = mean(KResFLL.ErrX1.^2);
        
%         fprintf('CkoPhiPLL = %.3f CkoPhiTeorPLL = %.3f\n', sqrt(DestPhiPLL)*180/pi, sqrt(KResPLL.DteorPhi)*180/pi)
%         fprintf('CkoWPLL = %.3f CkoWTeorPLL = %.3f CkoWFLL = %.3f CkoWTeorFLL = %.3f\n!!!!!!!!\n',...
%             sqrt(DestWPLL)/2/pi, sqrt(KResPLL.DteorW)/2/pi, sqrt(DestWFLL)/2/pi, sqrt(KResFLL.DteorW)/2/pi)
        save_statistic;
        
        if ~mod(m,fix(Np/10))
            fprintf('Progress po Np: %.2f%%\n', m*100/Np);
        end
        clear('KalmanFLL');
        clear('KalmanPLL');
    end
    
    fprintf('Progress po q_cno: %.2f%%\n', j*100/Nq)
end



