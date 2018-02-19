Br10ns = load('Figures\Data for plotting\Br10ns.mat');
Br100ps = load('Figures\Data for plotting\Br100ps.mat');
figure(1);clf
subplot(2,1,1); hold on;
plot(Br100ps.final_x_values,.03*Br100ps.TFY_total_ONandOFF,'Linewidth',2);
plot(Br100ps.final_x_values,Br100ps.TFY_pump_probe_avg,'Linewidth',2);
plot(Br10ns.final_x_values,Br10ns.TFY_pump_probe_avg,'Linewidth',2);
legend('.03 * Br','Br 100ps','Br 10ns');
xlabel('Energy, keV')
title('Br K-edge spectra, CsPbBr3')
ylabel('Pump-Probe Intensities')
subplot(2,1,2); hold on;
plot(Br100ps.final_x_values,Br100ps.TFY_pump_probe_avg,'Linewidth',2);
plot(Br10ns.final_x_values,3 * Br10ns.TFY_pump_probe_avg,'Linewidth',2);
legend('Br 100ps','3 * Br 10ns');
xlabel('Energy, keV')
title('Br K-edge TFY spectra, CsPbBr3')
ylabel('Pump-Probe Intensities')


Pb10ns = load('Figures\Data for plotting\Pb10ns.mat');
Pb100ps = load('Figures\Data for plotting\Pb100ps.mat');
figure(2);clf
subplot(2,1,1); hold on;
plot(Pb100ps.final_x_values,.01*Pb100ps.TFY_total_ONandOFF,'Linewidth',2);
plot(Pb100ps.final_x_values,Pb100ps.TFY_pump_probe_avg,'Linewidth',2);
plot(Pb10ns.final_x_values,Pb10ns.TFY_pump_probe_avg,'Linewidth',2);
legend('.01 * Pb','Pb 100ps','Pb 10ns');
xlabel('Energy, keV')
title('Pb L3-edge spectra, CsPbBr3')
ylabel('Pump-Probe Intensities')
subplot(2,1,2); hold on;
plot(Pb100ps.final_x_values,Pb100ps.TFY_pump_probe_avg,'Linewidth',2);
plot(Pb10ns.final_x_values,5 * Pb10ns.TFY_pump_probe_avg,'Linewidth',2);
legend('Pb 100ps','5 * Pb 10ns');
xlabel('Energy, keV')
title('Pb L3-edge TFY spectra, CsPbBr3')
ylabel('Pump-Probe Intensities')

figure(3);clf;hold on;
HERFD_sig_ON = Pb100ps.average_data_table(:,Pb100ps.detectorindices.HERFD_LaserON);
HERFD_sig_OFF = Pb100ps.average_data_table(:,Pb100ps.detectorindices.HERFD_LaserOFF);
plot(Pb100ps.final_x_values(Pb100ps.HERFD_pump_probe_avg~=0),nonzeros(Pb100ps.HERFD_pump_probe_avg),'Linewidth',2);
plot(Pb100ps.final_x_values(HERFD_sig_ON~=0),nonzeros(HERFD_sig_ON),'Linewidth',2);
plot(Pb100ps.final_x_values(HERFD_sig_OFF~=0),nonzeros(HERFD_sig_OFF),'Linewidth',2);
title([num2str(Pb100ps.scandetails),' Average of ',num2str(Pb100ps.counter),' scans']); xlabel('Energy (keV)'); ylabel('Partial Fluorescence Yield (A.U.)');
legend('Time resolved HERFD','HERFD, Laser ON','HERFD, Laser OFF');
figure(4);clf;
hold off;
tfy_scaler=10;
plot(Pb100ps.final_x_values(Pb100ps.HERFD_pump_probe_avg~=0),smooth(nonzeros(Pb100ps.HERFD_pump_probe_avg),1),Pb100ps.final_x_values,tfy_scaler*Pb100ps.TFY_pump_probe_avg,'Linewidth',2);
legend('pump-probe HERFD',[num2str(tfy_scaler),' * pump-probe TFY']);title([Pb100ps.scandetails,' ',sprintf('Average of %1.0f scans',Pb100ps.counter)]); xlabel('Energy (keV)');