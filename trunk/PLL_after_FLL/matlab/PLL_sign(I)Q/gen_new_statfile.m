StatFile.qcno_dB = 10:0.5:50;
StatFile.len_qcno_dB = length(StatFile.qcno_dB);

StatFile.Np = zeros(1, StatFile.len_qcno_dB);
StatFile.DestPhi = zeros(1, StatFile.len_qcno_dB);
StatFile.DestW = zeros(1, StatFile.len_qcno_dB);
StatFile.DestPhiSredn = zeros(1, StatFile.len_qcno_dB);
StatFile.DestPhi_Theor = zeros(1, StatFile.len_qcno_dB);
StatFile.DestW_Theor = zeros(1, StatFile.len_qcno_dB);

StatFile.PLL_Band = zeros(1, StatFile.len_qcno_dB);


save('Stat.mat', 'StatFile');