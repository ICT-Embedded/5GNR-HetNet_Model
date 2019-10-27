%Set line styles
Macro_all_DL = 'r--';
Micro_all_DL = 'b--';
Macro_all_UL = 'r-*';
Micro_all_UL = 'b-*';

%Plot DL spectral efficiency for macrocell vs small cell UEs
figure(100);
hold on; box on; grid on;
plot(SE_DL_active,...
    linspace(0,1,Ktotal),Macro_all_DL,'LineWidth',1);
plot(SE_DL_SC_active,...
    linspace(0,1,Ktotal_SC),Micro_all_DL,'LineWidth',1);
plot(SE_UL_active,...
    linspace(0,1,Ktotal),Macro_all_UL,'LineWidth',1);
plot(SE_UL_SC_active,...
    linspace(0,1,Ktotal_SC),Micro_all_UL,'LineWidth',1);
legend('Macro UE, DL','Micro UE, DL','Macro UE, UL','Micro UE, UL',...
    'Location','SouthEast','AutoUpdate','off');
xlabel('UE spectral efficiencies [bits/s/Hz]');
ylabel('CDF');
xlim([0 20]);