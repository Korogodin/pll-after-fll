load('Stat.mat', 'StatFile');
save('StatBackup.mat', 'StatFile');

StatFile.Std_a = std_a;
StatFile.Tc = Tc;

for p = 1:StatFile.len_qcno_dB
    if qcno_dB(j) == StatFile.qcno_dB(p)
        
        StatFile.Np(p) = StatFile.Np(p) + 1;
        StatFile.DestW(p) = StatFile.DestW(p) + DestW;
        StatFile.DestW_Theor(p) = KResFLL.DteorW;
        StatFile.FLL_Band(p) = KalmanFLL.Band;
        
    end
end

save('Stat.mat', 'StatFile')
