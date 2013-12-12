load('Stat.mat');

u1 = 0;
for u = 1:length(StatFile.qcno_dB)
    if (StatFile.Np(u) > 0)
        u1 = u1 + 1;
        Graph.Np(u1) = StatFile.Np(u);
        Graph.qcno_dB(u1) = StatFile.qcno_dB(u);
        
        Graph.DestPhiPLL(u1) = StatFile.DestPhiPLL(u)/StatFile.Np(u);
        Graph.DestWPLL(u1) = StatFile.DestWPLL(u)/StatFile.Np(u);
        Graph.DestPhiTeorPLL(u1) = StatFile.DestPhiTeorPLL(u);
        Graph.DestWTeorPLL(u1) = StatFile.DestWTeorPLL(u);
        Graph.PLLBand(u1) = StatFile.PLLBand(u);
        
        Graph.DestWFLL(u1) = StatFile.DestWFLL(u)/StatFile.Np(u);
        Graph.DestWTeorFLL(u1) = StatFile.DestWTeorFLL(u);
        Graph.FLLBand(u1) = StatFile.FLLBand(u);
       
    end
    
end

plot(Graph.qcno_dB, [sqrt(Graph.DestPhiPLL)*180/pi; 
    sqrt(Graph.DestPhiTeorPLL)*180/pi]);
xlabel('q_{c/n0}, dBHz');
ylabel('RMSE_{\phi}, deg');
% ylim([0 100]);
figure
plot(Graph.qcno_dB, [sqrt(Graph.DestWFLL)/2/pi; 
    sqrt(Graph.DestWTeorFLL)/2/pi])
xlabel('q_{c/n0}, dBHz');
ylabel('RMSE_{\omega} FLL, Hz');

figure
plot(Graph.qcno_dB, [sqrt(Graph.DestWPLL)/2/pi; 
    sqrt(Graph.DestWTeorPLL)/2/pi])
xlabel('q_{c/n0}, dBHz');
ylabel('RMSE_{\omega} PLL, Hz');

figure
plot(Graph.qcno_dB, [Graph.PLLBand; Graph.FLLBand])
title('\DeltaF, Hz')


