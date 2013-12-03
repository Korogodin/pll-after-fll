clear
close all
clc

Tmod = 60; % ����� �������������
Np = 100; %  ����� �������� �����
f0 = 10.3e6; % ������������� �������
Td = 1/(3.3712*f0); % �������� �������������
Tc = 0.005; % ������ �������������� � �����������
Tf = 0.005; % ������ ������ ��������
L = fix(Tc/Td); % ����� �������� �� ��������� ���������� � �����������

K = fix(Tmod/Tc); % ����� ���������� ���������� � ����������� �� ����� �������������

stdn = 8; % ��� ����
stdnIQ = sqrt((stdn^2)*L/2);

% ������ ���������� ������������ ����. GLONASS, page 162
alpha = 0.1; % ������ ������� ���������, �^-1
std_a = 1; %��� ���������
S_ksi = 2*(33*std_a)^2 * alpha; %������������ ��������� ������������ ����
stdIst = sqrt(S_ksi * Tf); %��� ������������ ����
Dksi = stdIst^2; % ��������� ������������ ����

qcno_dB = [45];
Nq = length(qcno_dB);

Ff = [1 Tf
      0 1]; % ���������� ������� ��� ���

Fc = [1 Tc 0
      0 1  Tc
      0 0  1]; % ���������� ������� ��� ������ ��������

Fp = [1 Tf 0
      0 1  Tf
      0 0  1]; % ���������� ������� ��� ���

for j = 1:Nq
    
    qcno = 10^(qcno_dB(j) / 10); % �/� � �����
    
    A = sqrt(4*qcno*stdn^2*Td);
    Aiq = A*L/2;
    
    % �������� � �������������� �������������� ��������������
    SdPLL = Aiq*erf(sqrt(qcno*Tc)); % �������� �� ��
    DekvPhi = stdnIQ^2 / SdPLL^2; % ��������� ���. ����� �� (�� �����) (� ���� �����. �-��)
    S_meas_phi = DekvPhi*Tf; % ��� ����� ��������������
    
    SdFLL = 1/12 * (Aiq)^2  * Tc^2; % �������� �� ��
    DekvW = (6/(qcno*Tc^3))*(1 + 1/(qcno*Tc)); % ��������� ���. ����� �� (�� �����?)
    
    %     DestPhi = 0;
    %     DestW = 0;
    %     DestPhiSredn = 0;
    
    for m = 1:Np
        
        Xist = [0; 0; 0];
        Xoporn = [0; 0; 0];
        XestPLL = [0; 0; 0];
        XestFLL = [0; 0];
        
        nIst = randn(1,K);
        
        % ������� ������� ��� ������
        KalmanFLL = CKalmanEqMesConstK(XestFLL, Ff, Tf);
        KalmanFLL.Xextr = [0; 0];
        
        KalmanPLL = CKalmanEqMesConstK(XestPLL, Fp, Tf); 
        KalmanPLL.Xextr = [0; 0; 0];
        
        % �������������� ���������� / �������� + ��� �� ������� ������������� ���
        % ������� �������
        initPLL;
        initFLL;
        
        I = nan(1,K);
        Q = nan(1,K);
        Ih = nan(1,K);
        Qh = nan(1,K);
        
        nI = stdnIQ*randn(1,K);
        nQ = stdnIQ*randn(1,K);
        
        omega_ist = nan(1,K);
        phase_ist = nan(1,K);
        phaseExtrOld = 0;
        
        for k = 1:K
            %Xoporn(2) = (KalmanPLL.Xextr(1) - phaseExtrOld)/Tf;
            
            Xoporn(2) = KalmanFLL.Xextr(2);
            
            deltaPhiOporn = Xist(1) - Xoporn(1);
            deltaWOporn = Xist(2) - Xoporn(2);
            
            mI = Aiq*cos(deltaPhiOporn+deltaWOporn*Tc/2)*sinc(deltaWOporn*Tc/2/pi);
            mQ = -Aiq*sin(deltaPhiOporn+deltaWOporn*Tc/2)*sinc(deltaWOporn*Tc/2/pi);
            
            h = sign(randn(1));
            
            I(k) = h*mI + 1*nI(k);
            Q(k) = h*mQ + 1*nQ(k);
            
            
            KResFLL.udW(k) = I(k)*Ih(k) + Q(k)*Qh(k);
            
            KResPLL.udPhi(k) = -sign(I(k)) * Q(k);
            % KResPLL.udPhi(k) = SdPLL*(Xist(1)-KalmanPLL.Xextr(1)) + stdnIQ*randn(1,1);
            
            KalmanPLL.Estimate(KResPLL.udPhi(k)/SdPLL);
            
            % if k*Tc >= tau
            KResPLL.ErrX1(k) = Xist(1) - KalmanPLL.Xest(1);
            KResPLL.ErrX2(k) = Xist(2) - KalmanPLL.Xest(2);
            KResPLL.ErrSrednPhase(k) = (2*(Xist(1)-KalmanPLL.Xest(1))+(Xist(2)-KalmanPLL.Xest(2))*Tf/2)/2;
            % end
            
            %             KResPLL.X{1}(k) = KalmanPLL.Xest(1);
            %             KResPLL.X{2}(k) = KalmanPLL.Xest(2);
            %
            %             phase_ist(k) = Xist(1);
            %             omega_ist(k) = Xist(2);
            %
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
            
            Xoporn(1) = Xoporn(1) + Xoporn(2)*Tc;
            
            phaseExtrOld = KalmanPLL.Xextr(1);
            KalmanPLL.Extrapolate();
            Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % ������ ��������� ��������� �������.
        end
        
        DestPhi = mean(KResPLL.ErrX1.^2);
        if (DestPhi >= pi^2)
            DestPhi = pi^2;
        end
        DestW = mean(KResPLL.ErrX2.^2);
        DestPhiSredn = mean(KResPLL.ErrSrednPhase.^2);
        fprintf('CkoPhi = %.3f CkoPhiTeor = %.3f\n', sqrt(DestPhi)*180/pi, sqrt(KResPLL.DteorPhi)*180/pi)
%         save_statistic;
        
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


