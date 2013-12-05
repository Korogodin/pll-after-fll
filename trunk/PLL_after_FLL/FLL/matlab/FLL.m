clear
close all
clc

Tmod = 60; % Время моделирования
Np = 100; %  число прогонов схемы
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
std_a = 1; %СКЗ ускорения
Sksi = 2*(33*std_a)^2 * alpha; %Спектральная плотность формирующего шума
stdIst = sqrt(Sksi * Tf); %СКО формирующего шума
Dksi = stdIst^2; % Дисперсия формирующего шума

qcno_dB = [10:2:50];
Nq = length(qcno_dB);

Ff = [1 Tf
      0 1]; % Переходная матрица для ССЧ

Fc = [1 Tc 0
      0 1  Tc
      0 0  1]; % Переходная матрица для модели процесса

for j = 1:Nq
    
    qcno = 10^(qcno_dB(j) / 10); % с/ш в разах
    
    A = sqrt(4*qcno*stdn^2*Td);
    Aiq = A*L/2;
    Aiq_ms = A*Lms/2;
    
    % Крутизна и флуктуационная характеристика дискриминатора
    SdFLL = 1/12 * (Aiq)^2  * Tc^2; % крутизна ДХ ЧД
    DekvW = (6/(qcno*Tc^3))*(1 + 1/(qcno*Tc)); % дисперсия экв. шумов ЧД (на входе?)
    Sw = DekvW*Tf;
   
    for m = 1:Np
        
        Xist = [0; 0; 0];
        Xoporn = [0; 0];
        XestFLL = [0; 0];
        
        nIst = randn(1,K);
        
        % Создаем фильтры для систем
        KalmanFLL = CKalmanEqMesConstK(XestFLL, Ff, Tf);
        KalmanFLL.Xextr = [0; 0];
      
        % ининциализация переменных / массивов + там же рассчет коэффициентов для
        % фильтра Калмана
        initFLL;
        
        I = zeros(1,K);
        Q = zeros(1,K);
        Ih = zeros(1,K);
        Qh = zeros(1,K);
        
        nIms = stdnIQms*randn(5,K);
        nQms = stdnIQms*randn(5,K);
        
        omega_ist = nan(1,K);
        
        for k = 1:K
            %Xoporn(2) = (KalmanPLL.Xextr(1) - phaseExtrOld)/Tf;
            
            Xoporn(2) = KalmanFLL.Xextr(1);
            
            dPhiOporn = Xist(1) - Xoporn(1);
            dWOporn = Xist(2) - Xoporn(2);
            
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
                      
            KalmanFLL.Estimate(KResFLL.udW(k)/SdFLL);

            KResFLL.ErrX1(k) = Xist(2) - KalmanFLL.Xest(1);
            
%             KResFLL.X{1}(k) = KalmanFLL.Xest(1);
%             KResFLL.X{2}(k) = KalmanFLL.Xest(2);
%             
%             omega_ist(k) = Xist(2);
%             
%             subplot(3,1,1)
%             plot((1:K)*Tc, [omega_ist/2/pi; KResFLL.X{1}/2/pi]);
%             title('\omega_{istinnaya} \omega_{ocenka}')
%             subplot(3,1,2)
%             plot((1:K)*Tc, KResFLL.ErrX1/2/pi);
%             title('error \omega')
%             subplot(3,1,3)
%             plot((1:K)*Tc, KResFLL.udW/SdFLL);
%             title('ud\omega')
%             drawnow
            
            Xoporn(1) = Xoporn(1) + Xoporn(2)*Tc;
            
            KalmanFLL.Extrapolate();
            Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % Модель изменения истинного вектора.
        end
        
        DestW = mean(KResFLL.ErrX1.^2);
%         fprintf('CkoW = %.3f CkoWTeor = %.3f\n', sqrt(DestW)/2/pi, sqrt(KResFLL.DteorW)/2/pi)
        save_statistic;

        if ~mod(m,fix(Np/10))
            fprintf('Progress po Np: %.2f%%\n', m*100/Np);
        end
        clear('KalmanFLL');
    end

    fprintf('Progress po q_cno: %.2f%%\n', j*100/Nq)
end



