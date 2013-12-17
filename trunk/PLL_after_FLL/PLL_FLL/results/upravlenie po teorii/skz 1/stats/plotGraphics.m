load('StatSignIQ.mat');
StatFileSign = StatFile;
clear('StatFile');

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

u2 = 0;
for h = 1:length(StatFileSign.qcno_dB)
    if (StatFileSign.Np(h) > 0)
        u2 = u2 + 1;
        GraphSign.Np(u2) = StatFileSign.Np(h);
        GraphSign.qcno_dB(u2) = StatFileSign.qcno_dB(h);
        GraphSign.DestPhi(u2) = StatFileSign.DestPhi(h)/StatFileSign.Np(h);
        GraphSign.DestW(u2) = StatFileSign.DestW(h)/StatFileSign.Np(h);
        GraphSign.DestPhiTeor(u2) = StatFileSign.DestPhiTeor(h);
        GraphSign.DestWTeor(u2) = StatFileSign.DestWTeor(h);
        GraphSign.PLLBand(u2) = StatFileSign.PLLBand(h);
        
    end
end

plot(Graph.qcno_dB, [sqrt(Graph.DestPhiPLL)*180/pi; sqrt(Graph.DestPhiTeorPLL)*180/pi],...
    GraphSign.qcno_dB, [sqrt(GraphSign.DestPhi)*180/pi; sqrt(GraphSign.DestPhiTeor)*180/pi]);
xlabel('q_{c/n0}, dBHz');
title('RMSE_{\phi}, deg');
% ylim([0 100]);
figure
plot(Graph.qcno_dB, [sqrt(Graph.DestWFLL)/2/pi; sqrt(Graph.DestWTeorFLL)/2/pi])
xlabel('q_{c/n0}, dBHz');
title('RMSE_{\omega} FLL, Hz');

figure
plot(Graph.qcno_dB, [sqrt(Graph.DestWPLL)/2/pi; sqrt(Graph.DestWTeorPLL)/2/pi])
xlabel('q_{c/n0}, dBHz');
title('RMSE_{\omega} PLL, Hz');

figure
plot(Graph.qcno_dB, [Graph.PLLBand; Graph.FLLBand])
title('\DeltaF, Hz')
title('q_{c/n0}, dBHz');

figure
plot(Graph.qcno_dB, [sqrt(Graph.DestWFLL)/2/pi; sqrt(Graph.DestWTeorFLL)/2/pi; sqrt(Graph.DestWPLL)/2/pi; sqrt(Graph.DestWTeorPLL)/2/pi],...
    GraphSign.qcno_dB, [sqrt(GraphSign.DestW)/2/pi; sqrt(GraphSign.DestWTeor)/2/pi]);
xlabel('q_{c/n0}, dBHz');
title('RMSE_{\omega} FLL, Hz; RMSE_{\omega} PLL, Hz');
