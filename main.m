clear
close all
clc

global K

Tmod = 30; % ����� �������������
Tf = 0.001; % ������ ������ ��������
Tc = 0.001; % ������ �������������� � �����������
K = fix(Tmod/Tc); % ����� ���������� ���������� � ����������� �� ����� �������������

MemoryAlloc; % �������������� ������

qcno_ist = 45*ones(1,K); % SNR �� ������ ���������, ����

Hfll = 5.5; % Hz, ������ ���

Xextr = [0; 0; 0]; % ������ �������������, ��������� ���������
Xest = Xextr;

XextrP = [0; 0; 0]; % ������ ������������� � ���
XestP = XextrP;

Ff = [1 0  0
    0 1  Tf
    0 0  1]; % ���������� ������� ��� ������ ������� (� ����� �������)

Fc = [1 Tc 0
    0 1  Tc
    0 0  1]; % ���������� ������� ��� ������ ������� (� ����� �����������)

Fincorr = [1 Tc 0
    0 1  0
    0 0  1]; % ���������� ������� ������ ���� � �����������

% ���
Ko = nan(3,1); % ������-������� ������������� �������
Ko(3) = 2*16/9*Hfll^2; % ������������ ����������� ������� � �������������� ������
Ko(2) = sqrt(2*Ko(3));
Ko(1) = 0;    % ��� �� ���, ��� ����! ��-�� ����� ���� �������, �� ����� ����, - ������� �������
Ko = Ko*Tf; % ������� � ������������� ���������� �������

% ���
HP = 45; % Hz, ������ ���
KoP = nan(3,1); % ������-������� ������������� �������
KoP(3) = (1.2*HP)^3; % ������������ ����������� ������� � �������������� ������
KoP(2) = 2*(KoP(3))^(2/3);
KoP(1) = 2*(KoP(3))^(1/3);
KoP = KoP*Tf; % ������� � ������������� ���������� �������

Xist = [0; 0; 0]; % �������� ������ ���������
% ������ ���������� ������������ ����. GLONASS, page 162
alpha = 0.1; % ������ ������� ���������, �^-1
std_a = 40; %��� ���������
S_ksi = 2*(33*std_a)^2 * alpha; %������������ ��������� ������������ ����
stdIst = Tc*sqrt(S_ksi / Tc); %��� ������������ ����
nIst = randn(1,K); 

stdn_IQ = ones(1,K)*8; % ��� ���� ������������ ����

nI = 1*stdn_IQ.*randn(1,K); % I-comp noise
nQ = 1*stdn_IQ.*randn(1,K); % Q-comp noise

w = 0; Isum = 0; Qsum = 0; Iold = 1; Qold = 0;
for k = 1:K
    
    % ������ ����.������������ �������������� ����
    EpsPhi(k) = Xist(1) - Xextr(1);
    EpsW(k) = Xist(2) - Xextr(2);
    
    qcno = 10.^(qcno_ist(k)/10);
    A_IQ(k) = stdn_IQ(k) .* sqrt(2 * qcno * Tc);
    A_IQ_eff(k) = A_IQ(k)*sinc(EpsW(k)*Tc/2 /pi);
    
    mI = A_IQ_eff(k) * cos(EpsW(k)*Tc/2 + EpsPhi(k));
    mQ = - A_IQ_eff(k) * sin(EpsW(k)*Tc/2 + EpsPhi(k));
    I(k) = mI + nI(k);
    Q(k) = mQ + nQ(k);
    Isum = Isum + I(k);
    Qsum = Qsum + Q(k);
    
    Xextr = Fincorr * Xextr; % ����� ���� � ����������� � ����� ����������
    
    w = w + 1;
    if w == fix(Tf/Tc)
        
        EpsPhiP(k) = Xist(1) - XextrP(1);
        % ������� �������������
        UdP(k) = -(I(k)*sin((XextrP(2) - Xextr(2))*Tf/2*1 + (XextrP(1) - Xextr(1))) + Q(k)*cos((XextrP(2) - Xextr(2))*Tf/2*1 + (XextrP(1) - Xextr(1))));
        SdP = A_IQ(k); % �������� ��
        XestP = XextrP + KoP*UdP(k)/SdP;  % ������ ������ �� ��������� �������� 
        XextrP = Fc*XestP;                % ������������� �� ��������� ��������
        
        Ud(k) = (I(k)*Qold - Q(k)*Iold); % ��������� �������������
        Sd = Tc*(A_IQ(k)*Tf/Tc)^2 * 1.3; % �������� ��
        Xest = Xextr + Ko*Ud(k)/Sd;  % ������ ������ �� ��������� �������� �������
        Xextr = Ff*Xest;             % ������������� �� ��������� ��������
        
        w = 0;
        Iold = Isum; Isum = 0;
        Qold = Qsum; Qsum = 0;
    end
    
    ErrW(k) = Xist(2) - Xest(2);
    ErrWP(k) = Xist(2) - XestP(2);
    
    omega(k) = Xest(2);
    omegaP(k) = XestP(2);
    omega_ist(k) = Xist(2);
    
    phase_ist(k) = Xist(1);
    phaseP(k) = XestP(1);
    
    %     if k == 12000
    %         Xist(2) = 100;
    %     end
    %     if k == 3000
    %         Xist(1) = Xist(1) + 2*pi/3;
    %     end
    Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % ������ ��������� ��������� �������.
    
    if ~mod(k,fix(K/10))
        fprintf('Progress: %.0f%%\n', 100*k/K);
    end
end

t = (1:K)*Tc;

hF = 0;
hF = figure(hF + 1);
subplot(2,1,1)
plot(t, omega / 2 / pi, t, omegaP / 2 / pi, t, omega_ist / 2 / pi);
title('Freq')
ylabel('\omega, Hz');
xlabel('t, s');
subplot(2,1,2)
plot(t, ErrW / 2 / pi, t, ErrWP / 2 / pi);
ylabel('Err_{\omega}, Hz');
xlabel('t, s');

hF = figure(hF + 1);
subplot(2,1,1)
plot(EpsW, Ud, '.')
title('Diskriminators')
ylabel('Ud')
xlabel('\epsilon_{\omega}, rad/s');
subplot(2,1,2)
plot(mod(EpsPhiP, 2*pi), UdP, '.')
ylabel('UdP')
xlabel('\epsilon_{\phi}, rad');

hF = figure(hF + 1);
subplot(2,1,1)
plot(t, phase / 2 / pi, t, phaseP / 2 / pi, t, phase_ist / 2 / pi);
title('Phase')
ylabel('\phi, cycles');
xlabel('t, s');
subplot(2,1,2)
plot(t, phase_ist / 2 / pi - phaseP / 2 / pi);
ylabel('Err\phi, cycles');
xlabel('t, s');
