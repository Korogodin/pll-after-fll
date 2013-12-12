StatFile.qcno_dB = 10:0.5:50;
StatFile.len_qcno_dB = length(StatFile.qcno_dB);

StatFile.Np = zeros(1, StatFile.len_qcno_dB);

StatFile.DestPhiPLL = zeros(1, StatFile.len_qcno_dB);
StatFile.DestWPLL = zeros(1, StatFile.len_qcno_dB);
StatFile.DestPhiTeorPLL = zeros(1, StatFile.len_qcno_dB);
StatFile.DestWTeorPLL = zeros(1, StatFile.len_qcno_dB);
StatFile.PLLBand = zeros(1, StatFile.len_qcno_dB);

StatFile.DestWFLL = zeros(1, StatFile.len_qcno_dB);
StatFile.DestWTeorFLL = zeros(1, StatFile.len_qcno_dB);
StatFile.FLLBand = zeros(1, StatFile.len_qcno_dB);

save('Stat.mat', 'StatFile');