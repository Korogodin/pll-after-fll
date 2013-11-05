clear all
close all
clc


Tmod = 12; % Время моделирования
M = 1000;
L = 413; % число точек на интервале интегрирования в корреляторе
Tc = 0.02; % Период интегрирования в корреляторе
Tf = 0.02; % Период работы фильтров
Td = Tc/L; % интервал дискретизации
f0 = 1/(3.3712*Td); % промежуточная частота


K = fix(Tmod/Tc); % Число интервалов накопления в корреляторе за время моделирования

stdn = 8e-2; % СКО шума 

qcno_ist = [10:1:42];
N = length(qcno_ist);


KResFLL.CKO_W_TEOR = zeros(1, N);
KResFLL.CKO_W = zeros(1, N);
KResNewPLL.CKO_W_TEOR = zeros(1, N);
KResNewPLL.CKO_W = zeros(1, N);
KResOldPLL.CKO_W_TEOR = zeros(1, N);
KResOldPLL.CKO_W = zeros(1, N);

KResFLL.Band = zeros(1, N);

for j = 1:N
    
    fprintf('qcno = %.2f dBHz\n', qcno_ist(j))
    
    for m = 1:M
        
        % qcno_ist = 45; % с/ш в дБ
        qcno = 10^(qcno_ist(j) / 10); % с/ш в разах
        
        A = sqrt(4*qcno*Td)*stdn;
        
        Xist = [0; 0; 0];
        Xoporn = [0; 0; 0];
        XopornOldPLL = [0; 0; 0];
        XestFLL = [0; 0];
        XextrFLL = [0; 0];
        XestNewPLL = [0; 0; 0];
        XextrNewPLL = [0; 0; 0];
        XestOldPLL = [0; 0; 0];
        XextrOldPLL = [0; 0; 0];
        
        Ff = [1 Tf
              0 1]; % Переходная матрица для ЧАП
        
        Fc = [1 Tc 0
              0 1  Tc
              0 0  1]; % Переходная матрица для модели процесса
        
        
        Fp = [1 Tf 0
              0 1  Tf
              0 0  1]; % Переходная матрица для ФАП
        
        % Расчет параметров формирующего шума. GLONASS, page 162
        alpha = 0.1; % Ширина спектра ускорения, с^-1
        std_a = 1; %СКЗ ускорения
        S_ksi = 2*(33*std_a)^2 * alpha; %Спектральная плотность формирующего шума
        stdIst = sqrt(S_ksi * Tf); %СКО формирующего шума
        Dksi = stdIst^2;
        nIst = randn(1,K);
        
        DmeasOldPLL = (stdn^2*L/2)/(A*L/2)^2;
        DmeasW = (6/(qcno*Tc^3))*(1 + 1/(qcno*Tc)); % дисперсия экв. шумов ЧД
        DmeasNewP = (stdn^2*L/2)/(A*L/2)^2; % дисперсия экв. шумов ФД
       
        SdFLL = 1/12 * (A*L/2)^2  * Tc^2; % крутизна ДХ ЧД
        SdNewPLL = (A*L/2); % крутизна ДХ ФД(нового)
        SdOldPLL = (A*L/2);
        
        KalmanFLL = CKalmanEqMesConstK(XestFLL, Ff, Tf); % вызываем класс фильтра Калмана для ЧАП
        KalmanNewPLL = CKalmanEqMesConstK(XestNewPLL, Fp, Tf); % вызываем класс фильтра Калмана для ФАП(новой)
        KalmanOldPLL = CKalmanEqMesConstK(XestOldPLL, Fp, Tf);
        
        % ининциализация переменных / массивов + там же рассчет коэффициентов для
        % фильтра Калмана
        initFLL;
        initNewPLL;
        initOldPLL;
        
        I = nan(1,K);
        Q = nan(1,K);
        I_OldPLL = nan(1,K);
        Q_OldPLL = nan(1,K);
        Ih = nan(1,K);
        Qh = nan(1,K);
        
        omega_ist = nan(1,K);
        phase_ist = nan(1,K);
        
        
        IF_phaseEnd = 0;
        phaseExtrOld = 0;
        for k = 1:K
            KalmanFLL.Extrapolate();
            KalmanNewPLL.Extrapolate();
            KalmanOldPLL.Extrapolate();
            
            XopornOldPLL(2) = (KalmanOldPLL.Xextr(1) - phaseExtrOld)/Tf ;
            Xoporn(2) = KalmanFLL.Xextr(1);
            KalmanNewPLL.Xextr(2) = KalmanFLL.Xextr(1);
            
            IF_phase = IF_phaseEnd + 2*pi*f0*(1:L)*Td;
            IF_phaseEnd = IF_phase(end);
            y = sign(1/3 * randn(1)) * A * cos(IF_phase + Xist(1) + Xist(2)*((1:L)-1)*Td) + 1*stdn*randn(1,L);
            
            I_oporn = cos(IF_phase + Xoporn(1) + Xoporn(2)*((1:L)-1)*Td);
            Q_oporn = sin(IF_phase + Xoporn(1) + Xoporn(2)*((1:L)-1)*Td);
            
            I_opornOldPLL = cos(IF_phase + XopornOldPLL(1) + XopornOldPLL(2)*((1:L)-1)*Td);
            Q_opornOldPLL = sin(IF_phase + XopornOldPLL(1) + XopornOldPLL(2)*((1:L)-1)*Td);
                        
            I(k) = y * I_oporn';
            Q(k) = y * Q_oporn';
            
            I_OldPLL(k) = y * I_opornOldPLL';
            Q_OldPLL(k) = y * Q_opornOldPLL';
            
            Ih(k) = - (y .* ((1:L)-1)*Td) * Q_oporn';
            Qh(k) = (y .* ((1:L)-1)*Td) * I_oporn';
            
            KResFLL.ud(k) = I(k)*Ih(k) + Q(k)*Qh(k); % расчет отсчета ЧД (I*I' + Q*Q')
            
            dPhi = (KalmanNewPLL.Xextr(2) - Xoporn(2))*Tf/2 + (KalmanNewPLL.Xextr(1) - Xoporn(1)); % расчет отсчета ФД(нового)
            KResNewPLL.udPhi(k) = -sign((I(k)*cos(dPhi) - Q(k)*sin(dPhi))) * (I(k)*sin(dPhi) + Q(k)*cos(dPhi));
            KResNewPLL.udW(k) = KResFLL.ud(k) + SdFLL*(KalmanNewPLL.Xextr(2) - Xoporn(2));
            KResOldPLL.udPhi(k) = -Q_OldPLL(k)*sign(I_OldPLL(k));           


            KalmanFLL.Estimate(KResFLL.ud(k)/SdFLL);
            KalmanNewPLL.Estimate([KResNewPLL.udPhi(k)/SdNewPLL; KResNewPLL.udW(k)/SdFLL]);
            KalmanOldPLL.Estimate(KResOldPLL.udPhi(k)/SdOldPLL);
            
            phaseExtrOld = KalmanOldPLL.Xextr(1);
            Xoporn = Fc*Xoporn;
            XopornOldPLL = Fc*XopornOldPLL;
            
            omega_ist(k) = Xist(2);
            phase_ist(k) = Xist(1);
            
            
            KResFLL.ErrX1(k) = Xist(2) - KalmanFLL.Xest(1);
            KResNewPLL.ErrX1(k) = Xist(1) - KalmanNewPLL.Xest(1);
            KResNewPLL.ErrX2(k) = Xist(2) - KalmanNewPLL.Xest(2);
            KResOldPLL.ErrX1(k) = Xist(1) - KalmanOldPLL.Xest(1);
            KResOldPLL.ErrX2(k) = Xist(2) - KalmanOldPLL.Xest(2);
            
            
            Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % Модель изменения истинного вектора.

            KResFLL.X{1}(k) = KalmanFLL.Xest(1);
            KResNewPLL.X{1}(k) = KalmanNewPLL.Xest(1);
            KResNewPLL.X{2}(k) = KalmanNewPLL.Xest(2);
            KResOldPLL.X{1}(k) = KalmanOldPLL.Xest(1);
            KResOldPLL.X{2}(k) = KalmanOldPLL.Xest(2);
            
        end
        
        
        KResFLL.CKO_W_TEOR(j) = sqrt(KResFLL.DteorW)/2 /pi + KResFLL.CKO_W_TEOR(j);
        KResFLL.CKO_W(j) = sqrt(mean(KResFLL.ErrX1.^2))/2 /pi + KResFLL.CKO_W(j);
        KResNewPLL.CKO_W_TEOR(j) = sqrt(KResNewPLL.DteorW)/2 /pi + KResNewPLL.CKO_W_TEOR(j);
        KResNewPLL.CKO_W(j) = sqrt(mean(KResNewPLL.ErrX2.^2))/2 /pi + KResNewPLL.CKO_W(j);
        KResOldPLL.CKO_W_TEOR(j) = sqrt(KResOldPLL.DteorW)/2 /pi + KResOldPLL.CKO_W_TEOR(j);
        KResOldPLL.CKO_W(j) = sqrt(mean(KResOldPLL.ErrX2.^2))/2 /pi + KResOldPLL.CKO_W(j);
        
        KResFLL.Band(j) = KalmanFLL.Band;
        
        % CKO_Phi_NEW_PLL(j) = sqrt(mean((KResNewPLL.ErrX1./(2*pi)*360).^2));
    end
end
KResFLL.CKO_W_TEOR = KResFLL.CKO_W_TEOR / M;
KResFLL.CKO_W = KResFLL.CKO_W / M;
KResNewPLL.CKO_W_TEOR = KResNewPLL.CKO_W_TEOR / M;
KResNewPLL.CKO_W = KResNewPLL.CKO_W / M;
KResOldPLL.CKO_W_TEOR = KResOldPLL.CKO_W_TEOR / M;
KResOldPLL.CKO_W = KResOldPLL.CKO_W / M;
t = (1:K)*Tc;

% savefile = ['stats\3systemi_uprPLLizFLL_skz' num2str(std_a) '_Tmod' num2str(Tmod) '_Tochek' num2str(M) '.mat'];
% save(savefile);

