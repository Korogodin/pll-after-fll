KResFLL.X{1}(k) = KalmanFLL.Xest(1);
KResFLL.X{2}(k) = KalmanFLL.Xest(2);

KResPLL.X{1}(k) = KalmanPLL.Xest(1);
KResPLL.X{2}(k) = KalmanPLL.Xest(2);

phase_ist(k) = Xist(1);
omega_ist(k) = Xist(2);

subplot(3,2,1)
plot((1:K)*Tc, [omega_ist/2/pi; KResFLL.X{1}/2/pi; KResPLL.X{2}/2/pi]);
title('\omega_{istinnaya} \omega_{ocenka}FLL \omega_{ocenka}PLL')
subplot(3,2,3)
plot((1:K)*Tc, [KResFLL.ErrX1/2/pi; KResPLL.ErrX2/2/pi]);
title('error_{\omega FLL} error_{\omega PLL}')
subplot(3,2,5)
plot((1:K)*Tc, KResFLL.udW/SdFLL);
title('ud_{\omega}FLL ud_{\omega}PLL')

subplot(3,2,2)
plot((1:K)*Tc, [phase_ist*180/pi; KResPLL.X{1}*180/pi]);
title('phase_{istinnaya} phase_{ocenka}')
subplot(3,2,4)
plot((1:K)*Tc, KResPLL.ErrX1*180/pi);
title('error phase')
subplot(3,2,6)
plot((1:K)*Tc, KResPLL.udPhi/SdPLL);
title('ud_{phase}')
drawnow