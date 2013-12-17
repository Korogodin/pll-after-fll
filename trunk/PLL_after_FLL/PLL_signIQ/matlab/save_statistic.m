load('Stat.mat', 'StatFile');
save('StatBackup.mat', 'StatFile');

StatFile.Std_a = std_a;
StatFile.Tc = Tc;
StatFile.type = 'PLL sign(I)Q';

for p = 1:StatFile.len_qcno_dB
    if qcno_dB(j) == StatFile.qcno_dB(p)
        
        StatFile.Np(p) = StatFile.Np(p) + 1;
        
        StatFile.DestPhi(p) = StatFile.DestPhi(p) + DestPhi; 
        StatFile.DestW(p) = StatFile.DestW(p) + DestW;
        StatFile.DestPhiTeor(p) = KResPLL.DteorPhi;
        StatFile.DestWTeor(p) = KResPLL.DteorW;
        StatFile.PLLBand(p) = KalmanPLL.Band;
                        
    end
end

save('Stat.mat', 'StatFile')
