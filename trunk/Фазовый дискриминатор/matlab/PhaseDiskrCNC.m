clear all
clc
close all
rng('shuffle'); % Reinit for randn

%settings

doScurve = 1; %������ ����������������� ��������������
doEQnoise = 0; %������ �������������� ��������������

f0 = 10.3e6; %������� �������
Td = 1/(3.3712*f0); %��� �������������
Tc = 20e-3; %����� ���������� � �����������
L = round(Tc/Td); %����� �������� �� ��������� ����������

% Td = Tc/L; %��� �������������
% f0 = 1/(3.3712*Td); %������� �������

stdn = 8; %��� ���� ���������
stdnIQ = sqrt((stdn^2)*L/2);

%diskriminator settings
Xist = [pi/3 2*pi*1e3]; %�������� �������� ������� \lambda

deltaPhioporn = 0; % !!!
deltaWoporn = 0;

qcno_dB = [45];
Nq = length(qcno_dB);
Np = 300000;

%������ ����������������� ��������������
if doScurve
    
    diskretPhiExtr = pi/180*2;
    phi_extr = [pi/3-1.8*pi/3:diskretPhiExtr:pi/3+1.8*pi/3];
    Nphi = length(phi_extr);
    
    deltaPhi = Xist(1) - phi_extr; %!!!
    
    ScurveTeor = nan(1,Nphi);
    SdPhi = nan(1,Nq);
    SdPhiTeor = nan(1,Nq);
    
    udPhi = zeros(1,Nphi);
    
    
    for j = 1:Nq
        
        fprintf('*** udPhi processing *** Progress: %.2f%%\n', j*100/Nq)
        
        qcno = 10^(qcno_dB(j)/10);
        A = sqrt(4*qcno*Td*stdn^2);
        Aiq = A*L/2;
        
        SdPhiTeor(j) = Aiq * erf(sqrt(qcno*Tc));
        
        for k = 1:Nphi
            
            fprintf('*** deltaPhi Progress: %.2f%%\n', k*100/Nphi)
            
            Xextr = [phi_extr(k) Xist(2)];
            
            deltaPhiextr = deltaPhioporn - deltaPhi(k);
            deltaWextr = deltaWoporn - (Xist(2)-Xextr(2));
            
            dPhi = (deltaWextr)*Tc/2 + (deltaPhiextr);
            
            mI = Aiq*cos(deltaPhioporn+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi);
            mQ = -Aiq*sin(deltaPhioporn+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi);
            
            nI = stdnIQ*randn(1,Np);
            nQ = stdnIQ*randn(1,Np);
            
            for m = 1:Np
                
                h = sign(1/3 * randn(1));
                
                I = h * mI + nI(m);
                Q = h * mQ + nQ(m);
                
                udPhi(k) = udPhi(k) - sign((I*cos(dPhi) - Q*sin(dPhi))) * (I*sin(dPhi) + Q*cos(dPhi));
            end
        end
        udPhi = udPhi / Np;
        ScurveTeor = Aiq * sin(deltaWoporn*Tc/2 + deltaPhi).*erf(sqrt(qcno*Tc)*cos(deltaPhi));
        
        index = find(phi_extr >= Xist(1), 1, 'first');
        SdPhi(j) = (udPhi(index-1)-udPhi(index+1))/(2*diskretPhiExtr);
        plot_ud(deltaPhi*180/pi, [udPhi; ScurveTeor; SdPhi(j)*deltaPhi; SdPhiTeor(j)*deltaPhi]);
        ylim([min(ScurveTeor) max(ScurveTeor)]);
        
    end
%     filename = ['deltaPhioporn=' num2str(deltaPhioporn*180/pi) '_ud(deltaPhioporn).mat'];
%     filename = ['ud_Tc=' num2str(Tc) ' .mat'];
%     save(filename, 'deltaPhi', 'udPhi', 'ScurveTeor', 'SdPhiTeor', 'SdPhi');
%     plot(qcno_dB, [SdPhi; SdPhiTeor]);
%     filename = ['Tc=' num2str(Tc) '_Sd(q).mat'];
%     save(filename, 'qcno_dB','SdPhi', 'SdPhiTeor');
end

%������ �������������� ��������������
if doEQnoise
    
    DmeasPhiTeor = nan(1,Nq);
    DmeasPhi = nan(1,Nq);
    udPhi = nan(1,Np);
    
    for j = 1:Nq
        
        fprintf('*** EQnoise processing *** Progress: %.2f%%\n', j*100/Nq)
        
        qcno = 10^(qcno_dB(j)/10);
        A = sqrt(4*qcno*Td*stdn^2);
        Aiq = A*L/2;
        SdPhi = Aiq;
        
        DmeasPhiTeor(j) = stdnIQ^2/SdPhi^2;
        
        for m = 1:Np
            if (~mod(m,5000))
                fprintf('*** Nakoplenie: %.2f%%\n', m*100/Np)
            end
            
            Xextr = Xist;
            deltaPhi = 0;
            
            deltaPhiextr = deltaPhioporn - deltaPhi;
            deltaWextr = deltaWoporn - (Xist(2)-Xextr(2));
            
            dPhi = (deltaWextr)*Tc/2 + (deltaPhiextr);
            
            I = Aiq*cos(deltaPhioporn+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi) + stdnIQ*randn(1,1);
            Q = -Aiq*sin(deltaPhioporn+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi) + stdnIQ*randn(1,1);
            %                 y = A * cos(2*pi*f0*(1:L)*Td + 2*pi*Xist(2)*(1:L)*Td + Xist(1)) + 1*n;
            
            %                 Iop = cos(2*pi*f0*(1:L)*Td + 2*pi*Xoporn(2)*(1:L)*Td + Xoporn(1));
            %                 Qop = sin(2*pi*f0*(1:L)*Td + 2*pi*Xoporn(2)*(1:L)*Td + Xoporn(1));
            %
            %                 I = y*Iop';
            %                 Q = y*Qop';
            
            %  udPhi(k) = -sign((I*cos(dPhi) - Q*sin(dPhi))) * (I*sin(dPhi) + Q*cos(dPhi));
            udPhi(1,m) = -(I*sin(dPhi) + Q*cos(dPhi));
        end
        DmeasPhi(j) = mean(udPhi.^2)/SdPhi^2;
    end
    %     filename = ['DmeasPhi(q)_Tc=' num2str(Tc) '.mat'];
    %     save(filename, 'qcno_dB', 'DmeasPhi', 'DmeasPhiTeor');
    plot(qcno_dB, (180/pi)*sqrt([DmeasPhi; DmeasPhiTeor]));
    
end




