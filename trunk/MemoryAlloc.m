%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%      Memory      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%    Allocation    %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_IQ_eff = nan(1,K); % ��������� ������������ ���������, ������� ������� sinc
I = nan(1,K);
Q = nan(1,K);
EpsPhi = nan(1, K); % ��������������� �� ���� (��� - �����)
EpsW = nan(1, K);   % ��������������� �� ������� (��� - �����)
Err_W_FLL = nan(1,K);    % ������ ������ ������� � ������� ���
Err_W_PLL = nan(1,K); % ������ ������ ������� � ������� ���
Err_Phi_PLL = nan(1,K);
UdFLL = nan(1,K);      % ������� �������������� � ���
UdPLL = nan(1,K);     % ������� �������������� � ���