clear
close all
clc

Tmod = 0.01; % ����� �������������

L = 113; % ����� ����� �� ��������� �������������� � �����������
Tc = 0.001; % ������ �������������� � �����������
Td = Tc/L; % �������� �������������
f0 = 1/(3*Td); % ������������� �������

Tf = 0.001; % ������ ������ ��������
Tc = 0.001; % ������ �������������� � �����������

K = fix(Tmod/Tc); % ����� ���������� ���������� � ����������� �� ����� �������������

stdn = 8; % ��� ���� ������������ ����
qcno_ist = 45; % �/� � ��
qcno = 10^(qcno_ist / 10); % �/� � �����

A = sqrt(2*qcno*Td)*stdn;

Xist = [0; 0; 0];
Xoporn = [0; 0; 0];

Fc = [1 Tc 0
    0 1  Tc
    0 0  1]; % ���������� ������� ��� ������ ��������

% ������ ���������� ������������ ����. GLONASS, page 162
alpha = 0.1; % ������ ������� ���������, �^-1
std_a = 40; %��� ���������
S_ksi = 2*(33*std_a)^2 * alpha; %������������ ��������� ������������ ����
stdIst = sqrt(S_ksi * Tf); %��� ������������ ����
nIst = randn(1,K);


I = nan(1,K);
Q = nan(1,K);

IF_phaseEnd = 0;
for k = 1:K
    IF_phase = IF_phaseEnd + 2*pi*f0*(1:L)*Td;
    IF_phaseEnd = IF_phase(end);
    y = A * cos(IF_phase + Xist(1) + Xist(2)*((1:L)-1)*Td) + stdn*randn(1,L);
    
    I_oporn = cos(IF_phase + Xoporn(1) + Xoporn(2)*((1:L)-1)*Td);
    Q_oporn = sin(IF_phase + Xoporn(1) + Xoporn(2)*((1:L)-1)*Td);
    
    I(k) = y * I_oporn';
    Q(k) = y * Q_oporn';
    
    Ih(k) = - (y * ((1:L)-1)*Td) * Q_oporn';
    Qh(k) = (y * ((1:L)-1)*Td) * I_oporn';
    
    Xist = Fc*Xist + [0; 0; 1]*nIst(k)*stdIst; % ������ ��������� ��������� �������.
    
end

