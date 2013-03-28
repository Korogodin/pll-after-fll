clear
close all
clc

Tmod = 30; % ����� �������������
Tf = 0.001; % ������ ������ ��������
Tc = 0.001; % ������ �������������� � �����������
L = fix(Tmod/Tc); % ����� ���������� ���������� � ����������� �� ����� �������������

Hpll = 245; % Hz, ������ ���
qcno_ist = 45*ones(1,L); % SNR �� ������ ���������, ����

Xextr = [0; 0; 0]; % ������ �������������, ��������� ���������
Xest = Xextr;

Fc = [1 Tc 0
    0 1  Tc
    0 0  1]; % ���������� ������� ��� ������ ������� (� ����� �����������)

Fincorr = [1 Tc 0
    0 1  0
    0 0  1]; % ���������� ������� ������ ���� � �����������

Ko = nan(3,1); % ������-������� ������������� �������
Ko(3) = (1.2*Hpll)^3; % ������������ ����������� ������� � �������������� ������
Ko(2) = 2*(Ko(3))^(2/3);
Ko(1) = 2*(Ko(3))^(1/3);
Ko = Ko*Tf; % ������� � ������������� ���������� �������

Xist = [0; 0; 0]; % �������� ������ ���������
% ������ ���������� ������������ ����. GLONASS, page 162
alpha = 0.1; % ������ ������� ���������, �^-1
std_a = 40; %��� ���������
S_ksi = 2*(33*std_a)^2 * alpha; %������������ ��������� ������������ ����
stdIst = Tc*sqrt(S_ksi / Tc); %��� ������������ ����
nIst = randn(1,L);

A_IQ = nan(1,L); % ��������� ������������ ���������
A_IQ_eff = nan(1,L); % �����, ������� ������� sinc
I = nan(1,L);
Q = nan(1,L);
EpsPhi = nan(1, L); % ��������������� �� ���� (��� - �����)
EpsW = nan(1, L);   % ��������������� �� ������� (��� - �����)
ErrW = nan(1,L);    % ������ ������ ������� � ������� ���
ErrPhi = nan(1,L);
omega = nan(1,L);   % ������ ������ ������� � ������� ���
omega_ist = nan(1,L); % �������� �������
phase = nan(1,L);   % ???
phase_ist = nan(1,L); % �������� ����
Ud = nan(1,L);      % ������� �������������� � ���

stdn_IQ = ones(1,L)*8; % ��� ���� ������������ ����

nI = 1*stdn_IQ.*randn(1,L); % I-comp noise
nQ = 1*stdn_IQ.*randn(1,L); % Q-comp noise

w = 0; Isum = 0; Qsum = 0;
for k = 1:L
    
    % ������ ����.������������ �������������� ����
    EpsPhi(k) = Xist(1) - Xextr(1);
    EpsW(k) = Xist(2) - Xextr(2);
    
    qcno = 10.^(qcno_ist(k)/10);
    A_IQ(k) = stdn_IQ(k) .* sqrt(2 * qcno * Tc);
    A_IQ_eff(k) = A_IQ(k)*sinc(EpsW(k)*Tc/2 /pi);
    
    mI = A_IQ_eff(k) * cos(EpsW(k)*Tc/2 + EpsPhi(k));
    mQ = - A_IQ_eff(k) * sin(EpsW(k)*Tc/2 + EpsPhi(k));
    I(k) = mI + 0*nI(k);
    Q(k) = mQ + 0*nQ(k);
    
    Isum = Isum + I(k);
    Qsum = Qsum + Q(k);
    
    Xextr = Fincorr * Xextr; % ����� ���� � ����������� � ����� ����������
    w = w + 1;
    if w == fix(Tf/Tc)
        
        EpsPhi(k) = Xist(1) - Xextr(1);
        % ������� �������������
        Ud(k) = -atan(Q(k)/I(k));
        Sd = 1; % �������� ��
        Xest = Xextr + Ko*Ud(k)/Sd;  % ������ ������ �� ��������� ��������
        Xextr = Fc*Xest;                % ������������� �� ��������� ��������
        
        w = 0;
        Isum = 0;
        Qsum = 0;
    end
    
    ErrW(k) = Xist(2) - Xest(2);
    ErrPhi(k) = Xist(1) - Xest(1);
    omega(k) = Xest(2);
    omega_ist(k) = Xist(2);
    phase_ist(k) = Xist(1);
    phase(k) = Xest(1);
    
    if ~mod(k,fix(L/10))
        fprintf('Progress: %.0f%%\n', 100*k/L);
    end
    Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % ������ ��������� ��������� �������.
end
t = (1:L)*Tc;
hF = 0;
hF = figure(hF + 1);
subplot(2,1,1)
plot(t, omega / 2 / pi, t, omega_ist / 2 / pi);
title('Freq')
ylabel('\omega, Hz');
xlabel('t, s');
subplot(2,1,2)
plot(t, ErrW / 2 / pi);
ylabel('Err_{\omega}, Hz');
xlabel('t, s');

hF = figure(hF + 1);
subplot(2,1,1)
plot(EpsW, Ud, '.')
title('Diskriminators')
ylabel('Ud')
xlabel('\epsilon_{\omega}, rad/s');
subplot(2,1,2)
plot(mod(EpsPhi, 2*pi), Ud, '.')
ylabel('Ud')
xlabel('\epsilon_{\phi}, rad');

hF = figure(hF + 1);
subplot(2,1,1)
plot(t, phase / 2 / pi, t, phase_ist / 2 / pi);
title('Phase')
ylabel('\phi, cycles');
xlabel('t, s');
subplot(2,1,2)
plot(t, phase_ist / 2 / pi - phase / 2 / pi);
ylabel('Err\phi, cycles');
xlabel('t, s');

