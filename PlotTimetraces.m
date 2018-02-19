clear;figure(3);clf;
Br134675 = load('Figures\Data for plotting\BrTimeScan13.4675.mat');
Br13472 = load('Figures\Data for plotting\BrTimeScan13.472.mat');
Br134765 = load('Figures\Data for plotting\BrTimeScan13.4765.mat');
Pb = load('Figures\Data for plotting\PbTimeScan.mat');

% Brline1 = -(5/3)*(Br134675.average_data_table(:,57)-Br134675.average_data_table(:,56)./Br134675.average_data_table(:,5));
% Brline2 = Br13472.average_data_table(:,57)-Br13472.average_data_table(:,56)./Br13472.average_data_table(:,5);
% Brline3 = -(5/3)*(Br134765.average_data_table(:,57)-Br134765.average_data_table(:,56)./Br134765.average_data_table(:,5));
% Pbline1 = .00015+12*(Pb.average_data_table(:,57)-Pb.average_data_table(:,56))./Pb.average_data_table(:,5);

figure(3); hold on;
plot(Br134675.plottable_x,Br134675.plottable_y,'Linewidth',1);
plot(Br13472.plottable_x,Br13472.plottable_y,'Linewidth',1);
plot(Br134765.plottable_x,Br134765.plottable_y,'Linewidth',1);
plot(Pb.plottable_x,Pb.plottable_y,'Linewidth',1);
legend('Br Pre-edge','Br White line','Br Post-edge','Pb');
xlabel('Time, nanoseconds')
title('Time Traces for CsPbBr3')
ylabel('Pump-Probe Intensities, individually scaled')
% xlim([-1.71e-7 -1.65e-7])