StatFile.qcno_dB = 10:0.5:50;
StatFile.len_qcno_dB = length(StatFile.qcno_dB);
StatFile.Std_a = zeros(1,1);
StatFile.Tc = zeros(1,1);

StatFile.Np = zeros(1, StatFile.len_qcno_dB);
StatFile.DestW = zeros(1, StatFile.len_qcno_dB);
StatFile.DestW_Theor = zeros(1, StatFile.len_qcno_dB);
StatFile.FLL_Band = zeros(1, StatFile.len_qcno_dB);


save('Stat.mat', 'StatFile');