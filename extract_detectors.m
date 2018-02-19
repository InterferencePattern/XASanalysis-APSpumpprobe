function [detectorlist_minimal_array,NHEADERLINES] = extract_detectors(data, filename)
% This function makes a list of the detectors in use for a given run.
%   This function was written by Thomas Rossi.
%Find the number of headerlines
        for i=1:numel(data(:))
            if strcmp(data{i},'')
            elseif strcmp(data{i},'# Column Descriptions:')
                firstdetector_index = i+1;
            elseif data{i}(1)~='#'
                NHEADERLINES=i-1;
                disp(filename);
%                 disp(NHEADERLINES);
                break;
            else
            end
        end
        
        % Use header lines to fetch detector names:
        detectorlist_fullnames = data(firstdetector_index:NHEADERLINES-1);
        detectorlist_split = cell([numel(detectorlist_fullnames) 1]);
        for i=1:numel(detectorlist_fullnames)
            detectorlist_split{i} = strsplit(detectorlist_fullnames{i},']');
            detectorlist_minimal(i) = detectorlist_split{i,1}(2);
        end
        
        detectorlist_minimal_array = string(detectorlist_minimal);
        
%        detectorlist_minimal_array = [];
%        for i = 2:length(detectorlist_minimal)
%            detectorlist_minimal_array(i)=detectorlist_minimal{i};
%        end

end

