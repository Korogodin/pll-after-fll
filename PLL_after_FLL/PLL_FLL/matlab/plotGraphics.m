load('Stat.mat');

plot(StatFile.qcno_dB, [sqrt(StatFile.DestPhi./StatFile.Np)*180/pi; 
    sqrt(StatFile.DestPhiSredn./StatFile.Np)*180/pi; sqrt(StatFile.DestPhi_Theor)*180/pi]);
title('RMSE \phi')
ylim([0 100]);
figure
plot(StatFile.qcno_dB, [sqrt(StatFile.DestW./StatFile.Np)/2/pi; 
    sqrt(StatFile.DestW_Theor)/2/pi])
title('RMSE \omega')
ylim([0 5]);
figure
plot(StatFile.qcno_dB, StatFile.PLL_Band)
title('\DeltaF, Hz')


