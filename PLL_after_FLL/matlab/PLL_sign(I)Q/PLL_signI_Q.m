clear all
close all



Tmod = 5; % ����� �������������
Np = 100; %  ����� �������� ����� (���������� ??)
f0 = 10.3e6; % ������������� �������
Td = 1/(3.3712*f0); % �������� �������������
Tc = 0.005; % ������ �������������� � �����������
Tf = 0.005; % ������ ������ ��������
L = fix(Tc/Td); % ����� �������� �� ��������� ���������� � �����������

K = fix(Tmod/Tc); % ����� ���������� ���������� � ����������� �� ����� �������������

stdn = 8; % ��� ����
stdnIQ = sqrt((stdn^2)*L/2);

qcno_dB = [45];
N = length(qcno_dB);

KResPLL.CKO_W_TEOR = zeros(1, N);
KResPLL.CKO_W = zeros(1, N);
KResPLL.CKO_Phi_TEOR = zeros(1, N);
KResPLL.CKO_Phi = zeros(1, N);

KResPLL.Band = zeros(1, N);

for j = 1:N
    
    fprintf('qcno = %.2f dBHz\n', qcno_dB(j))
    
    for m = 1:Np
        
        fprintf('Progress: %.2f%%\n', m*100/Np)
        
        qcno = 10^(qcno_dB(j) / 10); % �/� � �����
        
        A = sqrt(4*qcno*stdn^2*Td);
        Aiq = A*L/2;
        
        Xist = [0; 0; 0];
        Xoporn = [0; 0; 0];

        XestPLL = [0; 0; 0];
        XextrPLL = [0; 0; 0];

        Fc = [1 Tc 0
              0 1  Tc
              0 0  1]; % ���������� ������� ��� ������ ��������
                
        Fp = [1 Tf 0
              0 1  Tf
              0 0  1]; % ���������� ������� ��� ���
        
        % ������ ���������� ������������ ����. GLONASS, page 162
        alpha = 0.1; % ������ ������� ���������, �^-1
        std_a = 1; %��� ���������
        S_ksi = 2*(33*std_a)^2 * alpha; %������������ ��������� ������������ ����
        stdIst = sqrt(S_ksi * Tf); %��� ������������ ����
        Dksi = stdIst^2;
        nIst = randn(1,K);
        
        SdPLL = Aiq*erf(sqrt(qcno*Tc)); % �������� �� ��
        DmeasPhi = stdnIQ^2 / SdPLL^2; % ��������� ���. ����� �� (�� �����)
        S_meas_phi = DmeasPhi*Tf;
        
        KalmanPLL = CKalmanEqMesConstK(XestPLL, Fp, Tf); % �������� ����� ������� �������

        % �������������� ���������� / �������� + ��� �� ������� ������������� ���
        % ������� �������
        initPLL;
        
        I = nan(1,K);
        Q = nan(1,K);
        
        nI = stdnIQ*randn(1,K);
        nQ = stdnIQ*randn(1,K);
        
        omega_ist = nan(1,K);
        phase_ist = nan(1,K);
        
        
        IF_phaseEnd = 0;
        phaseExtrOld = 0;
        for k = 1:K
            KalmanPLL.Extrapolate();
            
            Xoporn(2) = (KalmanPLL.Xextr(1) - phaseExtrOld)/Tf ;
            
            IF_phase = IF_phaseEnd + 2*pi*f0*(1:L)*Td;
            IF_phaseEnd = IF_phase(end);
            
            
            mI = Aiq*cos(deltaphi+deltaW*Tc/2)*sinc(deltaW*Tc/2/pi);
            mQ = -Aiq*sin(deltaPhi+deltaW*Tc/2)*sinc(deltaW*Tc/2/pi);
            
            h = sign(1/3*randn(1));
            
            I(k) = h*mI + nI(k);
            Q(k) = h*mQ + nQ(k);
            
            y = sign(1/3 * randn(1)) * A * cos(IF_phase + Xist(1) + Xist(2)*((1:L)-1)*Td) + 1*stdn*randn(1,L);
            
            I_oporn = cos(IF_phase + Xoporn(1) + Xoporn(2)*((1:L)-1)*Td);
            Q_oporn = sin(IF_phase + Xoporn(1) + Xoporn(2)*((1:L)-1)*Td);
                                 
            I(k) = y * I_oporn';
            Q(k) = y * Q_oporn';

            KResPLL.udPhi(k) = -sign(I(k)) * Q(k);
           
            KalmanPLL.Estimate(KResPLL.udPhi(k)/SdPLL);
            
            KResPLL.X{1}(k) = KalmanPLL.Xest(1);
            KResPLL.X{2}(k) = KalmanPLL.Xest(2);
            
            phaseExtrOld = KalmanPLL.Xextr(1);
            Xoporn = Fc*Xoporn;
            
            phase_ist(k) = Xist(1);
            omega_ist(k) = Xist(2);
                        
            subplot(2,1,1)
            plot((1:K)*Tc, [phase_ist*180/pi; KResPLL.X{1}*180/pi]);
            title('\phi_{istinnaya} \phi_{ocenka}')
            subplot(2,1,2)
            plot((1:K)*Tc, (phase_ist - KResPLL.X{1})*180/pi);
            title('\phi_{error}')
            drawnow
                        
            KResPLL.ErrX1(k) = Xist(1) - KalmanPLL.Xest(1);
            KResPLL.ErrX2(k) = Xist(2) - KalmanPLL.Xest(2);
            
            Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % ������ ��������� ��������� �������.

        end
        
        
        KResPLL.CKO_W_TEOR(j) = sqrt(KResPLL.DteorW)/2 /pi + KResPLL.CKO_W_TEOR(j);
        KResPLL.CKO_W(j) = sqrt(mean(KResPLL.ErrX2.^2))/2 /pi + KResPLL.CKO_W(j);
        KResPLL.CKO_Phi_TEOR = sqrt(KResPLL.DteorPhi)*180/pi + KResPLL.CKO_Phi_TEOR(j);
        KResPLL.CKO_Phi(j) = sqrt(mean(KResPLL.ErrX1.^2))*180/pi + KResPLL.CKO_Phi(j);
                
        KResPLL.Band(j) = KalmanPLL.Band;
      
    end
end
KResPLL.CKO_W_TEOR = KResPLL.CKO_W_TEOR / Np;
KResPLL.CKO_W = KResPLL.CKO_W / Np;
KResPLL.CKO_Phi_TEOR = KResPLL.CKO_Phi_TEOR / Np;
KResPLL.CKO_Phi = KResPLL.CKO_Phi / Np;

t = (1:K)*Tc;

