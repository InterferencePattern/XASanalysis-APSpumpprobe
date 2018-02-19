function [ data_table ] = move_values( data_import, x_values, detectorlist_minimal, final_detectornames , x_rounding_factor, final_x_values, num_detectors)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

data_table = zeros([length(final_x_values) num_detectors]);
% data_are_strings = ischar(data_import{1}(1));

%% Convert cell matrix into 2D array:
data_import_nums = [];
for current_x_index=1:length(data_import)
    data_import_nums = cat(1,data_import_nums,str2double(data_import{current_x_index}));
    final_x_index = find(final_x_values==x_values(current_x_index),1); % Find the index of this x value in the final dataset
    for current_detector=1:numel(data_import{1}) % Loop over all available detectors in current dataset
        final_detector_index = find(strcmp(final_detectornames,detectorlist_minimal(current_detector)),1);
        % Find the appropriate "matching x" row in the existing
        % matrix to add the current data to.
        data_table(final_x_index, final_detector_index) = data_import_nums(current_x_index, current_detector); % Insert current single value into its place in the table.
    end
end
data_table(:,2)=round(x_rounding_factor*data_table(:,2))/x_rounding_factor;

end

