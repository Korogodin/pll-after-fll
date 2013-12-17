StatFile.qcno_dB = 10:0.5:50;
StatFile.len_qcno_dB = length(StatFile.qcno_dB);

StatFile.Np = zeros(1, StatFile.len_qcno_dB);
StatFile.DestPhi = zeros(1, StatFile.len_qcno_dB);
StatFile.DestW = zeros(1, StatFile.len_qcno_dB);
StatFile.DestPhiSredn = zeros(1, StatFile.len_qcno_dB);
StatFile.DestPhiTeor = zeros(1, StatFile.len_qcno_dB);
StatFile.DestWTeor = zeros(1, StatFile.len_qcno_dB);

StatFile.PLLBand = zeros(1, StatFile.len_qcno_dB);

save('Stat.mat', 'StatFile');