
load('Stat.mat', 'StatFile');
save('StatBackup.mat', 'StatFile');

for h = 1:StatFile.len_qcno_dB
    if qcno_dB(j) == StatFile.qcno_dB(h)
        
        StatFile.Np(h) = StatFile.Np(h) + 1;
        StatFile.DestPhi(h) = StatFile.DestPhi(h) + DestPhi;
        StatFile.DestW(h) = StatFile.DestW(h) + DestW;
        StatFile.DestPhiSredn(h) = StatFile.DestPhiSredn(h) + DestPhiSredn;
        StatFile.DestPhi_Theor(h) = KResPLL.DteorPhi;
        StatFile.DestW_Theor(h) = KResPLL.DteorW;
        
        if ~isnan(KalmanPLL.Band)
            StatFile.PLL_Band(h) = KalmanPLL.Band;
        end
        
    end
end

save('Stat.mat', 'StatFile')
