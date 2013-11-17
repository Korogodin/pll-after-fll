clear all
clc
close all
rng('shuffle'); % Reinit for randn

%settings

doScurve = 1; %расчет дискриминационной характеристики
doEQnoise = 0; %расчет флуктуационной характеристики

f0 = 10.3e6; %промежуточна€ частота
Td = 1/(3.3712*f0); %шаг дискретизации
Tc = 5e-3; %врем€ накоплени€ в коррел€торе
L = round(Tc/Td); %число отсчетов на интервале накоплени€

stdn = 8; %— ќ шума приемника
stdnIQ = sqrt((stdn^2)*L/2);

%diskriminator settings
Xist = [pi/3 2*pi*1e3]; %истинные значени€ вектора \lambda

deltaPhioporn = 0; % !!!
deltaWoporn = 0;

qcno_dB = [45];
Nq = length(qcno_dB);
Np = 300000;

%расчет дискриминационной характеристики
if doScurve
    
    diskretPhiExtr = pi/180*2;
    phi_extr = [pi/3-3.5*pi/3:diskretPhiExtr:pi/3+3.5*pi/3];
    Nphi = length(phi_extr);
    
    deltaPhi = Xist(1) - phi_extr; %!!!
    
    ScurveTeor = zeros(1,Nphi);
    SdPhi = zeros(1,Nq);
    SdPhiTeor = zeros(1,Nq);
    
    udPhi = zeros(1,Nphi);
    
    
    for j = 1:Nq
        
        qcno = 10^(qcno_dB(j)/10);
        A = sqrt(4*qcno*Td*stdn^2);
        Aiq = A*L/2;
        
        SdPhiTeor(j) = Aiq * sinc(deltaWoporn*Tc/2/pi);
        
        for k = 1:Nphi
            
            fprintf('*** Current progress: %.2f%% ... General progress: %.2f%%\n', k*100/Nphi, j*100/Nq)
            
            Xextr = [phi_extr(k) Xist(2)];
            
            deltaPhiextr = deltaPhioporn - deltaPhi(k);
            deltaWextr = deltaWoporn - (Xist(2)-Xextr(2));
            
            dPhi = (deltaWextr)*Tc/2 + (deltaPhiextr);
            
            mI = Aiq*cos(deltaPhioporn+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi);
            mQ = -Aiq*sin(deltaPhioporn+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi);
            
            nI = stdnIQ*randn(1,Np);
            nQ = stdnIQ*randn(1,Np);
            
            for m = 1:Np
              
                I = mI + nI(m);
                Q = mQ + nQ(m);
                
                udPhi(k) = udPhi(k) - (I*sin(dPhi) + Q*cos(dPhi));
            end
        end
        udPhi = udPhi / Np;
        ScurveTeor = Aiq * sin(0*Tc/2 + deltaPhi);
        
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


%расчет флуктуационной характеристики
if doEQnoise
    
    phi_extr = [Xist(1)+0];
    deltaPhi = Xist(1) - phi_extr; %!!!
    
    DPhiTeor1 = nan(1,Nq);
    DPhiTeor2 = nan(1,Nq);
    DPhi = nan(1,Nq);
    SdPhiTeor = nan(1,Nq);
    
    udPhi = zeros(1,Np);
    
    for j = 1:Nq
        
        qcno = 10^(qcno_dB(j)/10);
        A = sqrt(4*qcno*Td*stdn^2);
        Aiq = A*L/2;
        
        SdPhiTeor = Aiq;
        
        fprintf('*** Progress: %.2f%%\n', j*100/Nq)
        
        Xextr = [phi_extr Xist(2)];
        
        deltaPhiextr = deltaPhioporn - deltaPhi;
        deltaWextr = deltaWoporn - (Xist(2)-Xextr(2));
        
        dPhi = (deltaWextr)*Tc/2 + (deltaPhiextr);
        
        mI = Aiq*cos(deltaPhioporn+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi);
        mQ = -Aiq*sin(deltaPhioporn+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi);
        
        
        DPhiTeor1(j) = stdnIQ^2/SdPhiTeor^2;
        DPhiTeor2(j) = stdnIQ^2*sinc((deltaWextr)*Tc/2)^2/SdPhiTeor^2;
        
        nI = stdnIQ*randn(1,Np);
        nQ = stdnIQ*randn(1,Np);
        
        for m = 1:Np
                       
            I = mI + nI(m);
            Q = mQ + nQ(m);
            
            udPhi(1,m) = -(I*sin(dPhi) + Q*cos(dPhi));
        end
        DPhi(j) = mean((udPhi - mean(udPhi)).^2)/SdPhiTeor^2;
    end
    
%     filename = ['Dphi(q)_Tc=' num2str(Tc) '_deltaphi=' num2str(180/pi*deltaPhi) '.mat'];
%     save(filename, 'qcno_dB', 'DPhi', 'DPhiTeor1', 'DPhiTeor2');
    plot(qcno_dB, 180/pi * [sqrt(DPhi); sqrt(DPhiTeor1); sqrt(DPhiTeor2)]);
    
end




