clear all
clc
close all
rng('shuffle'); % Reinit for randn

%settings

doScurve = 1; %������ ����������������� ��������������
doEQnoise = 0; %������ �������������� ��������������

Tc = 20e-3; %����� ���������� � �����������
L = 51000; %����� �������� �� ��������� ����������
Td = Tc/L; %��� �������������
f0 = 1/(3.3712*Td); %������� �������
stdn = 8; %��� ���� ���������
stdnIQ = sqrt((stdn^2)*L/2);

%diskriminator settings
Xist = [pi/3 1e3]; %�������� �������� ������� \lambda

deltaPhiextr = 0; %������� ����� �������������� ��� � ������ ��� (�������� deltaPhi � ������)
deltaWextr = 0; %������� ����� �������������� ��� � ������ ��� (�������� deltaW � ������)
deltaWoporn = 0;
deltaPhioporn = 0;


%������ ����������������� ��������������
if doScurve
    
    Np = 15000;
    
    qcno_dB = [10:2:57];
    Nq = length(qcno_dB);
    
    diskretPhiExtr = 0.1;
    phi_extr = [-pi:diskretPhiExtr:pi];
    Nphi = length(phi_extr);
    
    Scurve = zeros(1,Nphi);
    ScurveTeor = nan(1,Nphi);
    SdPhi = nan(1,Nq);
    SdPhiTeor = nan(1,Nq);
    
    hF = 0;
    
    for j = 1:Nq
        
        fprintf('*** udPhi processing *** Progress: %.2f%%\n', j*100/Nq)
        
        qcno = 10^(qcno_dB(j)/10);
        A = sqrt(4*qcno*Td*stdn^2);
        Aiq = A*L/2;

        
        SdPhiTeor(j) = Aiq * sinc(deltaWoporn*Tc/2/pi);
              
        for m = 1:Np
            if (~mod(m,100))
                fprintf('*** Nakoplenie: %.2f%%\n', m*100/Np)
            end
            
            udPhi = nan(1,Nphi);
%             n = stdn*randn(1,L);
                        
            for k = 1:Nphi
                Xextr = [phi_extr(k) Xist(2)];
                Xoporn = Xextr - [deltaPhiextr deltaWextr];
                
                dPhi = (deltaWextr)*Tc/2 + (deltaPhiextr);
                
                I = Aiq*cos((Xist(1)-Xoporn(1))+(Xist(2)-Xoporn(2))*Tc/2)*sinc((Xist(2)-Xoporn(2))*Tc/2/pi) + stdnIQ*randn(1,1);
                Q = -Aiq*sin((Xist(1)-Xoporn(1))+(Xist(2)-Xoporn(2))*Tc/2)*sinc((Xist(2)-Xoporn(2))*Tc/2/pi) + stdnIQ*randn(1,1);
%                 y = A * cos(2*pi*f0*(1:L)*Td + 2*pi*Xist(2)*(1:L)*Td + Xist(1)) + 1*n;

%                 Iop = cos(2*pi*f0*(1:L)*Td + 2*pi*Xoporn(2)*(1:L)*Td + Xoporn(1));
%                 Qop = sin(2*pi*f0*(1:L)*Td + 2*pi*Xoporn(2)*(1:L)*Td + Xoporn(1));
%                 
%                 I = y*Iop';
%                 Q = y*Qop';
                
                %  udPhi(k) = -sign((I*cos(dPhi) - Q*sin(dPhi))) * (I*sin(dPhi) + Q*cos(dPhi));
                udPhi(k) = -(I*sin(dPhi) + Q*cos(dPhi));
            end
            Scurve = Scurve + udPhi;
        end
        Scurve = Scurve / Np;
        ScurveTeor = (A*L/2) * sinc(deltaWoporn*Tc/2/pi) * sin(0*Tc/2 + (Xist(1) - phi_extr));
        
        index = find(phi_extr >= Xist(1), 1, 'first');
        SdPhi(j) = (Scurve(index-1)-Scurve(index+1))/(2*diskretPhiExtr);
%         plot((Xist(1) - phi_extr)*180/pi, [Scurve; ScurveTeor; SdPhiTeor(j)*(Xist(1) - phi_extr); SdPhi(j)*(Xist(1) - phi_extr)]);
%         ylim([min(Scurve) max(Scurve)]);

    end
        plot(qcno_dB, [SdPhi; SdPhiTeor]);
        filename = ['Tc=' num2str(Tc) '_Sd(q).mat'];
        save(filename, 'qcno_dB','SdPhi', 'SdPhiTeor');
end

%������ �������������� ��������������
if doEQnoise
    
    qcno_dB = [10:5:60];
    Nq = length(qcno_dB);
    DmeasPhi = nan(1,Nq);
    DmeasPhi_Teor = nan(1,Nq);
    
    for j = 1:length(qcno_dB)
        
        Np = 10000;
        udPhi = nan(1,Np);
        
        qcno = 10^(qcno_dB(j)/10);
        A = sqrt(4*qcno*Td*stdn^2);
        
        DmeasPhi_Teor(j) = (stdn^2*L/2)/(A*L/2)^2;
        
        Xextr = Xist;
        Xoporn = Xextr;
        
        S = A * cos(2*pi*f0*(1:L)*Td + 2*pi*Xist(2)*(1:L)*Td + Xist(1));
        
        dPhi = (Xextr(1) - Xoporn(1))*Tc/2 + (Xextr(1) - Xoporn(1));
        
        Iop = cos(2*pi*f0*(1:L)*Td + 2*pi*Xoporn(2)*(1:L)*Td + Xoporn(1));
        Qop = sin(2*pi*f0*(1:L)*Td + 2*pi*Xoporn(2)*(1:L)*Td + Xoporn(1));
        
        for k = 1:Np
            y = sign(1/3 * randn(1)) * S  + 1*stdn*randn(1,L);
            I = y*Iop';
            Q = y*Qop';
            udPhi(k) = -sign((I*cos(dPhi) - Q*sin(dPhi))) * (I*sin(dPhi) + Q*cos(dPhi));
        end
        DmeasPhi(j) = mean(udPhi.^2)/(A*L/2)^2;
    end
    figure(2)
    plot(qcno_dB, (180/pi)*sqrt([DmeasPhi_Teor; DmeasPhi]));
end



