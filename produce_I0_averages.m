function [ final_data_table ] = produce_I0_averages( final_data_table, detectorindices, scantype, final_x_values, file_list )
%PRODUCE_I0_AVERAGES Takes complete data table and computes and applies i0
%   To improve signal-to-noise, similar scans (with matching x-values) are
%   corrected with a common i_0.

%% Set up initial variables
I_zero = squeeze(final_data_table(:,detectorindices.I0,:));
run_set={}; %n_remaining_runs_to_sort = size(I_zero,3); groupnumber=1;
if strcmp(scantype,'Timescan')
    all_x_values = squeeze(final_data_table(:,detectorindices.TimeIndex,:));
elseif strcmp(scantype,'Spectrum')
    all_x_values = squeeze(final_data_table(:,detectorindices.EnergyIndex,:));
end

%% For each scan, make a list of matching scan that share identical x-values.
for i=1:size(final_data_table,3)
    run_set{i}=find(ismember(all_x_values',all_x_values(:,i)','rows')); % Identify runs with identical sets of x values
end

%% Produce averaged I_zero for each scan:
for i=1:size(final_data_table,3) % For every run...
    % Throw out data where I_0 has dropped with respect to similar scans
%     if numel(run_set{i})>2 % Skip this if there aren't enough scans to make a distribution.
        for current_x=1:size(final_data_table,1) % For every x value...
            current_i0 = squeeze(final_data_table(current_x,detectorindices.I0,run_set{i}));
            if sum(current_i0)>0
%                 i0_fit = fitdist(current_i0(current_i0~=0),'Normal'); % Make Normal distribution of the I_zero values at current x value
%                 lower_i0_limit{i}(j) = i0_fit.mu-2*i0_fit.sigma;
%                 upper_i0_limit{i}(j) = i0_fit.mu+2*i0_fit.sigma;
                lower_i0_limit{i}(current_x) = mean(current_i0,'omitnan')*.97;
                upper_i0_limit{i}(current_x) = mean(current_i0,'omitnan')*1.03;
                extremei0_points{current_x} = find(current_i0<lower_i0_limit{i}(current_x) | current_i0>upper_i0_limit{i}(current_x)); % Find runs with extremely high or low I_0 values (within the current list)
                if not(isempty(extremei0_points{current_x}))
                    extremei0_runs{current_x} = run_set{i}(extremei0_points{current_x});% Convert that bad run/point list to global scan number
                else
                    extremei0_runs{current_x} = [];
                end
            end
        end
        if length(lower_i0_limit{i})<length(final_x_values)
            lower_i0_limit{i}(length(final_x_values))= NaN(1); % Add a little trailing data so that this variable doesn't terminate prematurely.
            upper_i0_limit{i}(length(final_x_values))= NaN(1); % Add a little trailing data so that this variable doesn't terminate prematurely.
            extremei0_runs{size(final_data_table,1)}=[]; % Add a little trailing data so that this variable doesn't terminate prematurely.
        end
        lower_i0_limit{i}(lower_i0_limit{i}==0) = NaN(1);
        upper_i0_limit{i}(upper_i0_limit{i}==0) = NaN(1);
%     end
end

%% Plot I_zeros before any bad points are removed
figure(8);clf;
for i=1:size(final_data_table,3)
    hold on;subplot(1,3,1);
    plot(final_x_values(final_data_table(:,detectorindices.I0,i)~=0),nonzeros(final_data_table(:,detectorindices.I0,i))); title('Individual I_{0} before removal of bad points');
    plot(final_x_values(~isnan(lower_i0_limit{i})),lower_i0_limit{i}(~isnan(lower_i0_limit{i})),final_x_values(~isnan(lower_i0_limit{i})),upper_i0_limit{i}(~isnan(lower_i0_limit{i})),'Linewidth',2);
    legend(file_list);
end

figure(9);clf;
for i=1:size(final_data_table,3)
    hold on;
    plot(final_x_values(final_data_table(:,detectorindices.I0,i)~=0),nonzeros(final_data_table(:,detectorindices.I0,i))); title('Individual I_{0}');
end

%% Set data for outlying points to zero
for current_x=1:size(final_data_table,1) % For every x value...
    final_data_table(current_x,:,extremei0_runs{current_x})=NaN(1); % ...clear out entire row corresponding to bad I_zero values
end

final_data_table(final_data_table==0)=NaN(1);

%% Produce average I_zero for each run
for i=1:size(final_data_table,3) % For every run...
%     averaged_I_zero{i} = sum(final_data_table( :, detectorindices.I0, run_set{i} ),3,'omitnan')./sum(final_data_table( :, detectorindices.I0, run_set{i} )~=0,3);
    averaged_I_zero{i} = mean(final_data_table( :, detectorindices.I0, run_set{i} ),3,'omitnan');
    averaged_I_zero{i}(averaged_I_zero{i}==0) = NaN(1);
end

%% Plot I_0
figure(8);
for i=1:size(final_data_table,3)
    subplot(1,3,2);hold on;
    current_plot = nonzeros(final_data_table(:,detectorindices.I0,i));
    plot(final_x_values(~isnan(lower_i0_limit{i})),current_plot(~isnan(lower_i0_limit{i}))); title('Individual I_{0}');
    plot(final_x_values(~isnan(lower_i0_limit{i})),lower_i0_limit{i}(~isnan(lower_i0_limit{i})),final_x_values(~isnan(lower_i0_limit{i})),upper_i0_limit{i}(~isnan(lower_i0_limit{i})),'Linewidth',2);
    subplot(1,3,3);hold on;
    plot(final_x_values(averaged_I_zero{i}~=0),nonzeros(averaged_I_zero{i})); title('Average I_{0} per set');
end
legend(file_list);

%% Apply those averaged I_zeros
for i=1:size(final_data_table,3)
    final_data_table(:,detectorindices.TFY_LaserON,i)=final_data_table(:,detectorindices.TFY_LaserON,i)./averaged_I_zero{i};
    final_data_table(:,detectorindices.TFY_LaserOFF,i)=final_data_table(:,detectorindices.TFY_LaserOFF,i)./averaged_I_zero{i};
    final_data_table(:,detectorindices.HERFD_LaserON,i)=final_data_table(:,detectorindices.HERFD_LaserON,i)./averaged_I_zero{i};
    final_data_table(:,detectorindices.HERFD_LaserOFF,i)=final_data_table(:,detectorindices.HERFD_LaserOFF,i)./averaged_I_zero{i};
    final_data_table(:,detectorindices.I0,i)=I_zero(:,i);
end

end

