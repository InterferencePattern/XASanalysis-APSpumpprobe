function [ data_table ] = loaddatascans(foldername, file_list,k, scantype, edgetype )
%LOADDATASCANS For a given file, data is parsed and loaded into a 2D matrix
%   Detailed explanation goes here

    %% Set limits on spectrum:
    Br_lower_spectral_limit = 13.452;
    Pb_lower_spectral_limit = 13.005;
    x_rounding_factor = 4000; %4000 corresponds to .25eV
    % Just to figure out how many data files have been converted ADU to
    % photons:

    if file_list{k} > 279
        disp('Detectors 62-63 have been converted to number of photons.')
        if file_list{k} > 321
            disp('Undulator might not be scanning.')
        end
    else
        disp('Detectors 62-63 are in ADU.')
    end
    
    filename=strcat(foldername,file_list{k});
    
    %Import the ascii data
    fid=fopen(filename);
    if fid>0
%         counter=counter+1; %counter for the number of scans
%         fprintf('Processing scan number %1.0f.\n',k);
        data=importdata(filename);
        
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
        detectorlist_minimal=string(detectorlist_minimal);
        
        %Import the data into a cell array, ignoring header lines:
        clear data_import x_values_to_add;
        for i=(NHEADERLINES+1):numel(data)
            data_import{i-NHEADERLINES}=strsplit(data{i});
        end
        
        % Remove last 4 data points for motor glitch:
        if strcmp(scantype,'Spectrum')
            data_import = data_import(4:end);
        end
                
        %% Put the values (for a single scan) into a table:
        clear x_values;
        for n=1:length(data_import)
            x_values(n,1) = str2num(cell2mat(data_import{1,n}(2)));
        end
        if strcmp(scantype,'Spectrum')
            x_values = round(x_rounding_factor*x_values)/x_rounding_factor;
            if strcmp(edgetype,'Files/Br') || strcmp(edgetype,'Files\Br')
                x_values(x_values<Br_lower_spectral_limit)=[];
                data_import=data_import(x_values>=Br_lower_spectral_limit);
            elseif strcmp(edgetype,'Files/Pb') || strcmp(edgetype,'Files\Pb')
                x_values(x_values<Pb_lower_spectral_limit)=[];
                data_import=data_import(x_values>=Pb_lower_spectral_limit);
            end
        elseif strcmp(scantype,'Timescan')
            %             x_values = round(10000*x_values)/10000;
        else
            disp('Unknown Scan Type');
        end
        %x_values=unique(x_values);
        
        %% Just make the complete list of x-values to insert data later
        if not(isequal(x_values,final_x_values))
            % Find new values
            x_values_to_add = x_values(not(ismember(x_values,final_x_values)));
            n_new_x_values = sum(not(ismember(x_values,final_x_values)));
            % Tack new x_values onto the end of the previously-loaded scans, then sort new values into final_x_values
            [final_x_values,x_sort_index] = sort(cat(1,final_x_values,x_values_to_add));
            % Tack zeros onto the end of the x-values, then sort them into the list to expand and accomodate .
            final_data_table = cat(1,final_data_table,zeros([n_new_x_values size(final_data_table,2) size(final_data_table,3)]));
            if not(isempty(final_data_table))
                final_data_table = final_data_table(x_sort_index,:,:); % Expanding existing final_data_table to accomodate current run.
            end
        end
        
        %% Insert new columns for added detectors:
        if not(isequal(detectorlist_minimal,final_detectornames))
            if not(isempty(final_detectornames))
                detectornames_to_add = detectorlist_minimal(not(ismember(detectorlist_minimal,final_detectornames)));
                n_new_detectors = sum(not(ismember(detectorlist_minimal,final_detectornames)));
                % Insert zeros into final_data_table for all previous scans
                final_data_table = cat(2,final_data_table,zeros([size(final_data_table,1) n_new_detectors size(final_data_table,3)]));
            else
                detectornames_to_add = detectorlist_minimal;
                n_new_detectors = numel(detectorlist_minimal);
            end
            % Insert new detector names into list
            [final_detectornames,detector_sort_index] = sort(cat(2,final_detectornames,detectornames_to_add));
            
            if not(isempty(final_data_table))
                final_data_table = final_data_table(:,detector_sort_index,:);
            end
        end
        %% Using that complete x*y space, move the data to the appropriate location:
        data_table = zeros([length(final_x_values) size(final_data_table,2)]);
        for current_x_index=1:length(data_import) % Loop over each x value in the current scan
            final_x_index = find(final_x_values==x_values(current_x_index),1); % Find the index of this x value in the final dataset
            for current_detector=1:numel(data_import{1}) % Loop over all available detectors
                final_detector_index = find(strcmp(final_detectornames,detectorlist_minimal(current_detector)),1);
                % Find the appropriate "matching x" row in the existing
                % matrix to add the current data to.
                if ischar(cell2mat(data_import{current_x_index}(current_detector)))
                    data_table(final_x_index,final_detector_index)=str2num(cell2mat(data_import{current_x_index}(current_detector)));
                else
                    data_table(final_x_index,final_detector_index)=cell2mat(data_import{current_x_index}(current_detector));
                end
            end
        end
        data_table(:,2)=round(x_rounding_factor*data_table(:,2))/x_rounding_factor;
        if isempty(final_data_table)
            final_data_table = data_table;
        else
            final_data_table(:,:,counter)=data_table(:,:);
        end

        fclose(fid);
    else
        %execute that part in case of multiple scans
        filename=strcat(sprintf('%s%04.0f_001.asc',basename,k));
        fid=fopen(filename);
        if fid>0
            counter_multiple=0;
            while 1>0
                counter_multiple=counter_multiple+1;
                filename=strcat(sprintf('%s%04.0f_%03.0f.asc',basename,k,counter_multiple));
                fid=fopen(filename);
                if fid>0
                    fclose(fid);
                else
                    break;
                end
            end
            counter_multiple=counter_multiple-1;
            
            for l=1:counter_multiple
                counter=counter+1; %counter for the number of scans
                fprintf('Processing scan number %1.0f_%03.0f.\n',k,l);
                filename=strcat(sprintf('%s%04.0f_%03.0f.asc',basename,k,l));
                data=importdata(filename);
                
                %Find the number of headerlines
                for i=1:numel(data(:))
                    if strcmp(data{i},'')
                    elseif strcmp(data{i},'# Column Descriptions:')
                        firstdetector_index = i+1;
                    elseif data{i}(1)~='#'
                        NHEADERLINES=i-1;
                        disp(filename);
                        disp(NHEADERLINES);
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
                detectorlist_minimal=string(detectorlist_minimal);
                
                % Check if the current run has 6-second integration in
                % addition to 3-second integration:
                find(strcmp(detectorlist_minimal,'  7idb:userCalc13.VAL, boxcarOFF (6s), '))
                
                %Import the data into a cell array, ignoring header lines:
                clear data_import x_values_to_add;
                for i=(NHEADERLINES+1):numel(data)
                    data_import{i-NHEADERLINES}=strsplit(data{i});
                end
                
                % Remove last 4 data points for motor glitch:
                if strcmp(scantype,'Spectrum')
                    data_import = data_import(4:end);
                end
                
                %% Put the values (for a single scan) into a table:
                clear x_values;
                for n=1:length(data_import)
                    x_values(n,1) = str2num(cell2mat(data_import{1,n}(2)));
                end
                if strcmp(scantype,'Spectrum')
                    x_values = round(x_rounding_factor*x_values)/x_rounding_factor;
                    if strcmp(edgetype,'Files/Br') || strcmp(edgetype,'Files\Br')
                        x_values(x_values<Br_lower_spectral_limit)=[];
                        data_import=data_import(x_values>=Br_lower_spectral_limit);
                    elseif strcmp(edgetype,'Files/Pb') || strcmp(edgetype,'Files\Pb')
                        x_values(x_values<Pb_lower_spectral_limit)=[];
                        data_import=data_import(x_values>=Pb_lower_spectral_limit);
                    end
                elseif strcmp(scantype,'Timescan')
                    %             x_values = round(10000*x_values)/10000;
                else
                    disp('Unknown Scan Type');
                end
                %x_values=unique(x_values);
                
                %% Insert new x values and shift previous final_data_table
                if not(isequal(x_values,final_x_values))
                    % Find new values
                    x_values_to_add = x_values(not(ismember(x_values,final_x_values)));
                    n_new_x_values = sum(not(ismember(x_values,final_x_values)));
                    % Insert new values into final_x_values
                    [final_x_values,x_sort_index] = sort(cat(1,final_x_values,x_values_to_add));
                    % Insert zeros into final_data_table for all previous scans.
                    final_data_table = cat(1,final_data_table,zeros([n_new_x_values size(final_data_table,2) size(final_data_table,3)]));
                    if not(isempty(final_data_table))
                        final_data_table = final_data_table(x_sort_index,:,:); % This sorts final_data_table by x-value
                    end
                end
                
                %% Insert new columns for added detectors:
                if not(isequal(detectorlist_minimal,final_detectornames))
                    if not(isempty(final_detectornames))
                        detectornames_to_add = detectorlist_minimal(not(ismember(detectorlist_minimal,final_detectornames)));
                        n_new_detectors = sum(not(ismember(detectorlist_minimal,final_detectornames)));
                        % Insert zeros into final_data_table for all previous scans
                        final_data_table = cat(2,final_data_table,zeros([size(final_data_table,1) n_new_detectors size(final_data_table,3)]));
                    else
                        detectornames_to_add = detectorlist_minimal;
                        n_new_detectors = numel(detectorlist_minimal);
                    end
                    % Insert new detector names into list
                    [final_detectornames,detector_sort_index] = sort(cat(2,final_detectornames,detectornames_to_add));
                    
                    if not(isempty(final_data_table))
                        final_data_table = final_data_table(:,detector_sort_index,:);
                    end
                end
                
                data_table = zeros([length(final_x_values) size(final_data_table,2)]);
                for current_x_index=1:length(data_import)
                    final_x_index = find(final_x_values==x_values(current_x_index));
                    for current_detector=1:numel(data_import{1})
                        final_detector_index = find(strcmp(final_detectornames,detectorlist_minimal(current_detector)));
                        % Find the appropriate "matching x" row in the existing
                        % matrix to add the current data to.
                        if ischar(cell2mat(data_import{current_x_index}(current_detector)))
                            data_table(final_x_index,final_detector_index)=str2num(cell2mat(data_import{current_x_index}(current_detector)));
                        else
                            data_table(final_x_index,final_detector_index)=cell2mat(data_import{current_x_index}(current_detector));
                        end
                    end
                end
                data_table(:,2)=round(x_rounding_factor*data_table(:,2))/x_rounding_factor;
            end
        else
        end
        if fid>0
            fclose(fid);
        else
        end
    end
end

