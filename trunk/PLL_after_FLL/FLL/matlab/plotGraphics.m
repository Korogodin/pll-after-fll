load('Stat.mat');

plot(StatFile.qcno_dB, [sqrt(StatFile.DestW./StatFile.Np)/2/pi; 
    sqrt(StatFile.DestW_Theor)/2/pi]);
title('RMSE \omega')
% ylim([0 max(sqrt(StatFile.DestW_Theor)/2/pi)]);
figure
plot(StatFile.qcno_dB, StatFile.FLL_Band)
title('\DeltaF, Hz')


