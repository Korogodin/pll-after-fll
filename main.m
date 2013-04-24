clear
close all
clc

% global K

f0 = 10.3e6;
Td = 1/(5*f0);

Tmod = 300; % ����� �������������
Tf = 0.001; % ������ ������ ��������
Tc = 0.001; % ������ �������������� � �����������
L = fix(Tc/Td)
K = fix(Tmod/Tc); % ����� ���������� ���������� � ����������� �� ����� �������������
stdn_IQ = 8; % ��� ���� ������������ ����

Ff = [1  Tf
      0  1]; % ���������� ������� ��� ������ ������� (� ����� �������)

Fp = [1 Tf 0
      0 1  Tf
      0 0  1]; % ���������� ������� ��� ������ ���� (� ����� �������)

Fc = [1 Tc 0
      0 1  Tc
      0 0  1]; % ���������� ������� ��� ������ ����

Fincorr = [1 Tc
           0 1]; % ���������� ������� ������ ���� � �����������

qcno_ist = [10:50];

CKO_W_PLL = nan(1,length(qcno_ist));
CKO_W_FLL = nan(1,length(qcno_ist));
CKO_Phi_PLL = nan(1,length(qcno_ist));
CKO_W_teor = nan(1,length(qcno_ist));

for j = 1:length(qcno_ist)
    fprintf('qcno = %.3f dB\n',qcno_ist(j))
MemoryAlloc; % �������������� ������

% qcno_ist = 45; % SNR �� ������ ���������, ����
qcno = 10.^(qcno_ist(j)/10);
A_IQ = stdn_IQ .* sqrt(2 * qcno * Tc); % ��������� ������������ ���������

Xist = [0; 0; 0]; % �������� ������ ���������
XextrFLL = [0; 0]; % ������ �������������, ��������� ���������
XestFLL = XextrFLL;
XextrPLL = [0; 0; 0]; % ������ ������������� � ���
XestPLL = XextrPLL;
XestPLL_old = [0; 0; 0];
Xcorr = [0; 0];

% ������ ���������� ������������ ����. GLONASS, page 162
alpha = 0.1; % ������ ������� ���������, �^-1
std_a = 40; %��� ���������
S_ksi = 2*(33*std_a)^2 * alpha; %������������ ��������� ������������ ����
stdIst = sqrt(S_ksi * Tf); %��� ������������ ����
nIst = randn(1,K);

% ������ ��� ������������� �����
stdEqNoiseP = stdn_IQ/A_IQ; % ����
% stdEqNoiseW = sqrt((1/(qcno*Tc^3))*(1+1/(2*qcno*Tc))); % �������
stdEqNoiseW = sqrt((6/(qcno*Tc^3))*(1+1/(2*qcno*Tc))); % �������


CKO_W_teor(j) = (4*S_ksi*(Tf*stdEqNoiseW^2)^3)^(1/8) / (2*pi);


% ���
KFLL = nan(2,1); % ������-������� ������������� �������
KFLL(2) = (stdIst/(stdEqNoiseW*Tf));
KFLL(1) = sqrt(2*KFLL(2));
KFLL = KFLL*Tf; % ������� � ������������� ���������� �������

% ���
KPLL = nan(3,1); % ������-������� ������������� �������
KPLL(3) = (stdIst/(stdEqNoiseP*Tf));
KPLL(2) = 2*(KPLL(3))^(2/3);
KPLL(1) = 2*(KPLL(3))^(1/3);
KPLL = KPLL*Tf; % ������� � ������������� ���������� �������

nI = 1*stdn_IQ.*randn(1,K); % I-comp noise
nQ = 1*stdn_IQ.*randn(1,K); % Q-comp noise

w = 0; Isum = 0; Qsum = 0; Iold = 1; Qold = 0;

for k = 1:K
    XestPLL_old = XestPLL;
    % ������ ����.������������ �������������� ����
    Xcorr = Fincorr * Xcorr; % ����� ���� � ����������� � ����� ����������
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
        % ������� �������������
        UdPLL(k) = -(Isum*sin((XextrPLL(2) - Xcorr(2))*Tf/2*1 + (XextrPLL(1) - Xcorr(1))) + Qsum*cos((XextrPLL(2) - Xcorr(2))*Tf/2*1 + (XextrPLL(1) - Xcorr(1))));
        UdPLL_mean(k) = A_IQ*sinc((Xist(2) - Xcorr(2))*Tf/2 /pi)*sin((Xist(2) - XextrPLL(2))*Tf/2*1 + (Xist(1) - XextrPLL(1)));
        SdPLL = A_IQ; % �������� ��
        XestPLL = XextrPLL + KPLL*UdPLL(k)/SdPLL;  % ������ ������ �� ��������� ��������
        XextrPLL = Fp*XestPLL;                % ������������� �� ��������� ��������
        
        UdFLL(k) = (Xist(2) - XextrFLL(1)) + stdEqNoiseW*randn(1,1); % ��������� �������������
        SdFLL = 1; % �������� ��
        XestFLL = XextrFLL + KFLL*UdFLL(k)/SdFLL;  % ������ ������ �� ��������� �������� �������
        XextrFLL = Ff*XestFLL;             % ������������� �� ��������� ��������
        
        Xcorr(2) = XextrFLL(1); % ��������� �������� � ����������� �� ���
        XextrPLL(2) = XextrFLL(1);
        w = 0;
        Iold = Isum; Isum = 0;
        Qold = Qsum; Qsum = 0;
    end
    
    Err_W_FLL(k) = Xist(2) - XestFLL(1);
    Err_W_PLL(k) = Xist(2) - (XestPLL(1) - XestPLL_old(1))/Tf;
    Err_Phi_PLL(k) = Xist(1) - XestPLL(1);
    
    Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % ������ ��������� ��������� �������.
    
%     omega_est_PLL(k) = XestPLL(2);
%     omega_est_FLL(k) = XestFLL(1);
%     omega_ist(k) = Xist(2);
    
    if ~mod(k,fix(K/10))
        fprintf('Progress: %.0f%%\n', 100*k/K);
    end
end
CKO_W_PLL(j) = sqrt(mean((Err_W_PLL./(2*pi)).^2));
CKO_W_FLL(j) = sqrt(mean((Err_W_FLL./(2*pi)).^2));
CKO_Phi_PLL(j) = sqrt(mean((Err_Phi_PLL./(2*pi)*360).^2));
end

t = (1:K)*Tc;

% plot(t, Err_W_FLL /2 /pi);
% title('Freq error ')
% ylabel('\deltaf, Hz');
% xlabel('t, s');

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


