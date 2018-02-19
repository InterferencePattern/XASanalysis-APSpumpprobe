function [ final_data_table ] = normalize_data( final_data_table, detectorindices, normalization_start_index, normalization_stop_index )
%normalize_data Normalizes 3D matrix of spectral/timescan data
%   Detailed explanation goes here

for i=1:size(final_data_table,3)
    %% Prepare normalization factors for TFY and HERFD data:
    TFY_ROI = nonzeros((final_data_table(normalization_start_index:normalization_stop_index,detectorindices.TFY_LaserOFF,i)+final_data_table(normalization_start_index:normalization_stop_index,detectorindices.TFY_LaserON,i)));
    HERFD_ROI = nonzeros((final_data_table(normalization_start_index:normalization_stop_index,detectorindices.HERFD_LaserOFF,i)+final_data_table(normalization_start_index:normalization_stop_index,detectorindices.HERFD_LaserON,i)));
    mean_TFY_ROI = mean(TFY_ROI,'omitnan');
    mean_HERFD_ROI = mean(HERFD_ROI,'omitnan');
    
    if isfinite(mean_TFY_ROI)
        normalization_factor_TFY(i) = abs(mean_TFY_ROI);
    else
        normalization_factor_TFY(i) = 0;
    end
    
    if isfinite(mean_HERFD_ROI)
        normalization_factor_HERFD(i) = abs(mean_HERFD_ROI);
    else
        normalization_factor_HERFD(i) = 0;
    end
    
    %% Apply normalization factors to TFY and HERFD data:
    if normalization_factor_TFY(i)~=0
        final_data_table(:,detectorindices.TFY_LaserON,i) = final_data_table(:,detectorindices.TFY_LaserON,i)./normalization_factor_TFY(i);
        final_data_table(:,detectorindices.TFY_LaserOFF,i) = final_data_table(:,detectorindices.TFY_LaserOFF,i)./normalization_factor_TFY(i);
    end
    if normalization_factor_HERFD(i)~=0
        final_data_table(:,detectorindices.HERFD_LaserON,i) = final_data_table(:,detectorindices.HERFD_LaserON,i)./normalization_factor_HERFD(i);
        final_data_table(:,detectorindices.HERFD_LaserOFF,i) = final_data_table(:,detectorindices.HERFD_LaserOFF,i)./normalization_factor_HERFD(i);
    end
    
    % If you get NaNs for the normalization factors, it means there is a
    % run somewhere with no data in the region of interest.
end

end