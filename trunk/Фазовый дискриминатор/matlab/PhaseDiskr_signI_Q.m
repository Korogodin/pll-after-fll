clear all
clc
close all
rng('shuffle'); % Reinit for randn

%settings

doScurve = 0; %расчет дискриминационной характеристики
doEQnoise = 1; %расчет флуктуационной характеристики

f0 = 10.3e6; %промежуточна€ частота
Td = 1/(3.3712*f0); %шаг дискретизации
Tc = 1e-3; %врем€ накоплени€ в коррел€торе
L = round(Tc/Td); %число отсчетов на интервале накоплени€

stdn = 8; %— ќ шума приемника
stdnIQ = sqrt((stdn^2)*L/2);

%diskriminator settings
Xist = [pi/3 2*pi*1e3]; %истинные значени€ вектора \lambda

deltaWoporn = 0;

qcno_dB = [10:2:57];
Nq = length(qcno_dB);
Np = 300000;

%расчет дискриминационной характеристики
if doScurve
    
    diskretPhiOporn = pi/180*2;
    phi_oporn = [pi/3-2*pi/3:diskretPhiOporn:pi/3+2*pi/3];
    Nphi = length(phi_oporn);
    
    deltaPhiOporn = Xist(1) - phi_oporn; %!!!
    
    ScurveTeor = zeros(1,Nphi);
    SdPhi = zeros(1,Nq);
    SdPhiTeor = zeros(1,Nq);
    
    udPhi = zeros(1,Nphi);
    
    
    for j = 1:Nq
        
        qcno = 10^(qcno_dB(j)/10);
        A = sqrt(4*qcno*Td*stdn^2);
        Aiq = A*L/2;
        
        SdPhiTeor(j) = Aiq * erf(sqrt(qcno*Tc));
        
        for k = 1:Nphi
            
            fprintf('*** Current progress: %.2f%% ... General progress: %.2f%%\n', k*100/Nphi, j*100/Nq)
          
            mI = Aiq*cos(deltaPhiOporn(k)+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi);
            mQ = -Aiq*sin(deltaPhiOporn(k)+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi);
            
            nI = stdnIQ*randn(1,Np);
            nQ = stdnIQ*randn(1,Np);
            
            for m = 1:Np
                
                h = sign(1/3 * randn(1));
                
                I = h * mI + nI(m);
                Q = h * mQ + nQ(m);
                
                udPhi(k) = udPhi(k) - sign(I)*Q;
            end
        end
        udPhi = udPhi / Np;
        ScurveTeor = Aiq * sin(0*Tc/2 + deltaPhiOporn).*erf(sqrt(qcno*Tc)*cos(deltaPhiOporn));
        
        index = find(phi_oporn >= Xist(1), 1, 'first');
        SdPhi(j) = (udPhi(index-1)-udPhi(index+1))/(2*diskretPhiOporn);
        plot_ud(deltaPhiOporn*180/pi, [udPhi; ScurveTeor; SdPhi(j)*deltaPhiOporn; SdPhiTeor(j)*deltaPhiOporn]);
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
    
    deltaPhiOporn = 0; %!!!
    
    DPhiTeor1 = nan(1,Nq);
    DPhiTeor2 = nan(1,Nq);
    DPhi = nan(1,Nq);
    SdPhiTeor = nan(1,Nq);
    
    udPhi = zeros(1,Np);
    
    for j = 1:Nq
        
        qcno = 10^(qcno_dB(j)/10);
        A = sqrt(4*qcno*Td*stdn^2);
        Aiq = A*L/2;
        
        SdPhiTeor = Aiq * erf(sqrt(qcno*Tc));
        
        fprintf('*** Progress: %.2f%%\n', j*100/Nq)
                
        mI = Aiq*cos(deltaPhiOporn+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi);
        mQ = -Aiq*sin(deltaPhiOporn+deltaWoporn*Tc/2)*sinc(deltaWoporn*Tc/2/pi);
        
        ScurveTeor = Aiq * sin(0*Tc/2 + deltaPhiOporn).*erf(sqrt(qcno*Tc)*cos(deltaPhiOporn));
        
        DPhiTeor1(j) = (stdnIQ^2 + mQ^2 - ScurveTeor^2) / SdPhiTeor^2 ;
        DPhiTeor2(j) = stdnIQ^2 / SdPhiTeor^2;
        
        nI = stdnIQ*randn(1,Np);
        nQ = stdnIQ*randn(1,Np);
        
        for m = 1:Np
            
            h = sign(1/3 * randn(1));
            
            I = h * mI + nI(m);
            Q = h * mQ + nQ(m);
            
            udPhi(1,m) = - sign(I)*Q;
        end
        DPhi(j) = mean((udPhi - mean(udPhi)).^2)/SdPhiTeor^2;
    end
    
%     filename = ['Dphi(q)_Tc=' num2str(Tc) '_deltaphi=' num2str(180/pi*deltaPhiOporn) '.mat'];
%     save(filename, 'qcno_dB', 'DPhi', 'DPhiTeor1', 'DPhiTeor2', 'SdPhiTeor');
    plot(qcno_dB, 180/pi * [sqrt(DPhi); sqrt(DPhiTeor1); sqrt(DPhiTeor2)])
    
end
