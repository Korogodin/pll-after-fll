load('Stat.mat');
u1 = 0;
for u = 1:length(StatFile.qcno_dB)
    if (StatFile.Np(u) > 0)
        u1 = u1 + 1;
        Graph.Np(u1) = StatFile.Np(u);
        Graph.qcno_dB(u1) = StatFile.qcno_dB(u);
        Graph.DestW(u1) = StatFile.DestW(u)/StatFile.Np(u);
        Graph.DestW_Theor(u1) = StatFile.DestW_Theor(u);
        Graph.FLL_Band(u1) = StatFile.FLL_Band(u);
       
    end
    
end
plot(Graph.qcno_dB, [sqrt(Graph.DestW)/2/pi; sqrt(Graph.DestW_Theor)/2/pi]);
ylabel('— ќш оценки частоты, √ц')
xlabel('q_{c/n0}, дЅ√ц')
xlim([15 50]);
set(gca,'box','off');
grid on
ylim([0 60]);
figure
plot(Graph.qcno_dB, Graph.FLL_Band)
ylabel('Ўумова€ полоса, √ц')
xlabel('q_{c/n0}, дЅ√ц')
xlim([15 50]);
set(gca,'box','off');
grid on


