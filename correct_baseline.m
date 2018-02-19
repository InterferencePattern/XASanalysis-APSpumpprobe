function [ final_data_table,discardedruns ] = correct_baseline( baselinecorrection_on, scantype, final_data_table,detectorindices,final_x_values, baseline_indices )
%normalize_data Normalizes baselines 3D matrix of spectral/timescan data
%   Detailed explanation goes here
discardedruns=0;
corrected_runs=[];

if strcmp(scantype,'Timescan')
    for i=1:size(final_data_table,3)
        % For each scan, determine which x_values exist:
        extant_rows = find(final_data_table(:,detectorindices.TFY_LaserON,i)~=0 & isfinite(final_data_table(:,detectorindices.TFY_LaserON,i)));
        beforetimezero = find(final_x_values>80e-12);
        x_baseline_TFYon = mean(final_data_table(beforetimezero,detectorindices.TFY_LaserON,i),1,'omitnan');
        x_baseline_TFYoff = mean(final_data_table(beforetimezero,detectorindices.TFY_LaserOFF,i),1,'omitnan');
        x_baseline_HERFDon = mean(final_data_table(beforetimezero,detectorindices.HERFD_LaserON,i),1,'omitnan');
        x_baseline_HERFDoff = mean(final_data_table(beforetimezero,detectorindices.HERFD_LaserOFF,i),1,'omitnan');
        
        if baselinecorrection_on
            for j=1:length(extant_rows)
                final_data_table(extant_rows(j),detectorindices.TFY_LaserON,i) = final_data_table(extant_rows(j),detectorindices.TFY_LaserON,i)-x_baseline_TFYon;
                final_data_table(extant_rows(j),detectorindices.TFY_LaserOFF,i) = final_data_table(extant_rows(j),detectorindices.TFY_LaserOFF,i)-x_baseline_TFYoff;
                final_data_table(extant_rows(j),detectorindices.HERFD_LaserON,i) = final_data_table(extant_rows(j),detectorindices.HERFD_LaserON,i)-x_baseline_HERFDon;
                final_data_table(extant_rows(j),detectorindices.HERFD_LaserOFF,i) = final_data_table(extant_rows(j),detectorindices.HERFD_LaserOFF,i)-x_baseline_HERFDoff;
            end
        end
    end
elseif strcmp(scantype,'Spectrum')
    for i=1:size(final_data_table,3)
        extant_rows = find(final_data_table(:,detectorindices.TFY_LaserON,i)~=0 & isfinite(final_data_table(:,detectorindices.TFY_LaserON,i)));
        if not(isempty(extant_rows)) && extant_rows(5)<baseline_indices
            % Take the average of the first 5 data points as the baseline:
            x_baseline_TFYon = mean(final_data_table(extant_rows(1:5),detectorindices.TFY_LaserON,i),'omitnan');
            x_baseline_TFYoff = mean(final_data_table(extant_rows(1:5),detectorindices.TFY_LaserOFF,i),'omitnan');
            x_baseline_HERFDon = mean(final_data_table(extant_rows(1:5),detectorindices.HERFD_LaserON,i),'omitnan');
            x_baseline_HERFDoff = mean(final_data_table(extant_rows(1:5),detectorindices.HERFD_LaserOFF,i),'omitnan');
            
            % Actually subtract the baseline from this scan:
            if baselinecorrection_on
                for j=1:length(extant_rows)
                    final_data_table(extant_rows(j),detectorindices.TFY_LaserON,i) = final_data_table(extant_rows(j),detectorindices.TFY_LaserON,i) - mean([x_baseline_TFYon x_baseline_TFYoff]);
                    final_data_table(extant_rows(j),detectorindices.TFY_LaserOFF,i) = final_data_table(extant_rows(j),detectorindices.TFY_LaserOFF,i) - mean([x_baseline_TFYon x_baseline_TFYoff]);
                    final_data_table(extant_rows(j),detectorindices.HERFD_LaserON,i) = final_data_table(extant_rows(j),detectorindices.HERFD_LaserON,i) - mean([x_baseline_HERFDon x_baseline_HERFDoff]);
                    final_data_table(extant_rows(j),detectorindices.HERFD_LaserOFF,i) = final_data_table(extant_rows(j),detectorindices.HERFD_LaserOFF,i) - mean([x_baseline_HERFDon x_baseline_HERFDoff]);
                end
            end
            corrected_runs = cat(1,corrected_runs,i);
        else
            discardedruns = discardedruns+1;
        end
    end
    %% Do alternative baseline correction if data does not reach into the pre-edge
    for i=1:size(final_data_table,3)
        extant_rows = find(final_data_table(:,detectorindices.TFY_LaserON,i)~=0 & isfinite(final_data_table(:,detectorindices.TFY_LaserON,i)));
        if not(isempty(extant_rows)) && not(extant_rows(5)<baseline_indices)
            % Determine the index of the point to be taken:
            first_x_point=extant_rows(1);
            
            % Set the averages determined previously as the baseline:
            x_baseline_TFYon = final_data_table(first_x_point,detectorindices.TFY_LaserON,i) - mean(mean(final_data_table([first_x_point-1 first_x_point+1],detectorindices.TFY_LaserON,corrected_runs),'omitnan'));
            x_baseline_TFYoff = final_data_table(first_x_point,detectorindices.TFY_LaserOFF,i) - mean(mean(final_data_table([first_x_point-1 first_x_point+1],detectorindices.TFY_LaserOFF,corrected_runs),'omitnan'));
            x_baseline_HERFDon = final_data_table(first_x_point,detectorindices.HERFD_LaserON,i) - mean(mean(final_data_table([first_x_point-1 first_x_point+1],detectorindices.HERFD_LaserON,corrected_runs),'omitnan'));
            x_baseline_HERFDoff = final_data_table(first_x_point,detectorindices.HERFD_LaserOFF,i) - mean(mean(final_data_table([first_x_point-1 first_x_point+1],detectorindices.HERFD_LaserOFF,corrected_runs),'omitnan'));
            
            % Actually subtract the baseline from this scan:
            if baselinecorrection_on
                for j=1:length(extant_rows)
                    final_data_table(extant_rows(j),detectorindices.TFY_LaserON,i) = final_data_table(extant_rows(j),detectorindices.TFY_LaserON,i)-x_baseline_TFYon;
                    final_data_table(extant_rows(j),detectorindices.TFY_LaserOFF,i) = final_data_table(extant_rows(j),detectorindices.TFY_LaserOFF,i)-x_baseline_TFYoff;
                    final_data_table(extant_rows(j),detectorindices.HERFD_LaserON,i) = final_data_table(extant_rows(j),detectorindices.HERFD_LaserON,i)-x_baseline_HERFDon;
                    final_data_table(extant_rows(j),detectorindices.HERFD_LaserOFF,i) = final_data_table(extant_rows(j),detectorindices.HERFD_LaserOFF,i)-x_baseline_HERFDoff;
                end
            end
        end
    end
else
    error('Scan type does not match Spectrum or Timescan.')
end

end

