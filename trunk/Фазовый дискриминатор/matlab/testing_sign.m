close all

stdnIQ = 8;
qcno_dB = 45;
qcno = 10^(qcno_dB/10);
T = 5e-3;
Aiq = stdnIQ*sqrt(2*T*qcno);

Np = 30000;



deltaPhiOporn = [-pi:pi/180*1:pi];
Nphi = length(deltaPhiOporn);

udPlus = zeros(1,Nphi);
udMinus = zeros(1,Nphi);

Pimore = zeros(1,Nphi);
Piless = zeros(1,Nphi);

for i = 1:Nphi
         
    mI = Aiq*cos(deltaPhiOporn(i));
    mQ = -Aiq*sin(deltaPhiOporn(i));
    
    nI = stdnIQ*randn(1,Np);
    nQ = stdnIQ*randn(1,Np);
    
    Pimore(i) = 0.5 + 0.5*erf(mI/(stdnIQ*sqrt(2)));
    Piless(i) = 0.5 - 0.5*erf(mI/(stdnIQ*sqrt(2)));
    
    for k = 1:Np
        
        Iplus = mI + nI(k);
        Qplus = mQ + nQ(k);
        
        Iminus = -mI + nI(k);
        Qminus = -mQ + nQ(k);
        
        udPlus(i) = udPlus(i) - sign(Iplus)*Qplus;
        udMinus(i) = udMinus(i) - sign(Iminus)*Qminus;
        
    end
end
udPlus = udPlus/Np;
udMinus = udMinus/Np;


plot(deltaPhiOporn*180/pi, [udPlus; Aiq*(Pimore - Piless).*sin(deltaPhiOporn)]);

