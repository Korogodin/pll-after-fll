load('Stat.mat', 'StatFile');
save('StatBackup.mat', 'StatFile');

StatFile.Std_a = std_a;
StatFile.Tc = Tc;
StatFile.type = 'PLL+FLL';

for p = 1:StatFile.len_qcno_dB
    if qcno_dB(j) == StatFile.qcno_dB(p)
        
        StatFile.Np(p) = StatFile.Np(p) + 1;
        
        StatFile.DestPhiPLL(p) = StatFile.DestPhiPLL(p) + DestPhiPLL; 
        StatFile.DestWPLL(p) = StatFile.DestWPLL(p) + DestWPLL;
        StatFile.DestPhiTeorPLL(p) = KResPLL.DteorPhi;
        StatFile.DestWTeorPLL(p) = KResPLL.DteorW;
        StatFile.PLLBand(p) = KalmanPLL.Band;
        
        StatFile.DestWFLL(p) = StatFile.DestWFLL(p) + DestWFLL;
        StatFile.DestWTeorFLL(p) = KResFLL.DteorW;
        StatFile.FLLBand(p) = KalmanFLL.Band;
        
    end
end

save('Stat.mat', 'StatFile')
