clear all
clc
close all

Tc = 1e-3;
L = 413;
Td = Tc/L;
f0 = 1/(3.3712*Td);


Xist = [pi/3 1e3];

% charact = char('Scurve', 'EQ noise');

for charact = {'Scurve', 'EQ noise'}
    if strcmp(charact, 'Scurve')
        
        qcno = 10^(45/10);
        stdn = 8;
        
        A = sqrt(4*qcno*Td*stdn^2);
        
        phi_extr = [-pi:2*pi/500:pi];
        
        for k = 1:length(phi_extr)
            Xextr = [phi_extr(k) 1e3];
            Xoporn = Xextr;
            
            y = sign(1/3 * randn(1)) * A * cos(2*pi*f0*(1:L)*Td + 2*pi*Xist(2)*(1:L)*Td + Xist(1)) + 0*stdn*randn(1,L);
                        
            dPhi = (Xextr(1) - Xoporn(1))*Tc/2 + (Xextr(1) - Xoporn(1));
            Iop = cos(2*pi*f0*(1:L)*Td + 2*pi*Xoporn(2)*(1:L)*Td + Xoporn(1));
            Qop = sin(2*pi*f0*(1:L)*Td + 2*pi*Xoporn(2)*(1:L)*Td + Xoporn(1));
            
            I = y*Iop';
            Q = y*Qop';
            
            udPhi(k) = -sign((I*cos(dPhi) - Q*sin(dPhi))) * (I*sin(dPhi) + Q*cos(dPhi));
        end
        plot((Xist(1) - phi_extr)*180/pi, [udPhi; A*L/2*(Xist(1) - phi_extr)]);
        ylim([min(udPhi) max(udPhi)]);    
    
    elseif strcmp(charact, 'EQ noise')
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
    end


