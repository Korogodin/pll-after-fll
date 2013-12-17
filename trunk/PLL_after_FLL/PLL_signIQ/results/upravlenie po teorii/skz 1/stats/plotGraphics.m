load('Stat.mat');

u1 = 0;
for u = 1:length(StatFile.qcno_dB)
    if (StatFile.Np(u) > 0)
        u1 = u1 + 1;
        Graph.Np(u1) = StatFile.Np(u);
        Graph.qcno_dB(u1) = StatFile.qcno_dB(u);
        
        Graph.DestPhi(u1) = StatFile.DestPhi(u)/StatFile.Np(u);
        Graph.DestW(u1) = StatFile.DestW(u)/StatFile.Np(u);
        Graph.DestPhiTeor(u1) = StatFile.DestPhiTeor(u);
        Graph.DestWTeor(u1) = StatFile.DestWTeor(u);
        Graph.PLLBand(u1) = StatFile.PLLBand(u);
     
    end
    
end

plot(Graph.qcno_dB, [sqrt(Graph.DestPhi)*180/pi; 
    sqrt(Graph.DestPhiTeor)*180/pi]);
xlabel('q_{c/n0}, dBHz');
ylabel('RMSE_{\phi}, deg');
% ylim([0 100]);

figure
plot(Graph.qcno_dB, [sqrt(Graph.DestW)/2/pi; 
    sqrt(Graph.DestWTeor)/2/pi])
xlabel('q_{c/n0}, dBHz');
ylabel('RMSE_{\omega}, Hz');

figure
plot(Graph.qcno_dB, Graph.PLLBand)
title('\DeltaF, Hz')
xlabel('q_{c/n0}, dBHz');
