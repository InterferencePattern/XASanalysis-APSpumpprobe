% ProcessAPSdata.m processes the raw data from APS X-ray Absorption
% experiments performed in June 2017 by the LSU at EPFL.

clear
%% Choose the folder of all runs you want to average:
foldername=input('Give the name of the folder you want to load, no apostrophes, and slash at the end.','s');
cd(foldername); folderparse = strsplit(foldername);edgetype=lower(folderparse{2});foldernameparts = strsplit(foldername,'/');scandetails=foldernameparts{2};
if not(isempty(findstr(foldernameparts{2},'Spectra')));scantype='Spectrum';elseif not(isempty(findstr(foldernameparts{2},'Timescan')));scantype='Timescan';end;
foldercontents = dir('*.asc'); % Fetch all the scan numbers from available files
file_list = {foldercontents.name}';cd ../;cd ../;
xray_duration_FWHM = 80e-12; %seconds
monochromator_res_eV = .7; 
% xray_duration = xray_duration_FWHM / 2.355; %seconds

%% Decide whether to use integrated signal (Br) or photon counting (Pb)
if strcmp(scantype,'Spectrum')
    % Set limits on spectrum:
    Br_lower_spectral_limit = 13.452;
    Pb_lower_spectral_limit = 13.005;
    x_rounding_factor = 4000; %4000 corresponds to .25eV, units of 1/keV
elseif strcmp(scantype,'Timescan')
    x_rounding_factor = 1/(1e12); % Units are 1/seconds so this rounds to 1ps.
else
    error('Unknown scan type, does not conform to Spectrum or Timescan.')
end

baselinecorrection_on=true; % Revert this to true for final analysis
i0_correction_on=true; % Revert this to true for final analysis
normalization_on=true; % Revert this to true for final analysis
dropbadframes=true; % Revert this to true for final analysis
normalize_pp=false; % Revert this to true for final analysis

if strcmp('files/pb',edgetype)
    photoncounting=true;
elseif strcmp('files/br',edgetype)
    photoncounting=false;
else
    error('Unknown edge type. Is it Pb or Br?')
end

plotindividually=true;

%% Import ASCII data and make single matrix of all scans:
counter=0;
final_data_table = []; final_x_values = [];final_detectornames = {};
for k = 1:length(file_list)
    %Import the ascii data
    filename=strcat(foldername,file_list{k});
    fid=fopen(filename);
    if fid>0
        counter=counter+1; %counter for the number of scans
        fprintf('Processing scan number %1.0f.\n',k);
        data=importdata(filename);
        [detectorlist_minimal, NHEADERLINES] = extract_detectors(data, filename);
        
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
            if strcmp(edgetype,'files/br') || strcmp(edgetype,'files/br')
                x_values(x_values<Br_lower_spectral_limit)=[];
                data_import=data_import(x_values>=Br_lower_spectral_limit);
            elseif strcmp(edgetype,'files/pb') || strcmp(edgetype,'files/pb')
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
            if isempty(final_detectornames)
                detectornames_to_add = detectorlist_minimal;
                n_new_detectors = numel(detectorlist_minimal);
            else
                detectornames_to_add = detectorlist_minimal(not(ismember(detectorlist_minimal,final_detectornames)));
                n_new_detectors = sum(not(ismember(detectorlist_minimal,final_detectornames)));
                % Insert zeros into final_data_table for all previous scans
                final_data_table = cat(2,final_data_table,zeros([size(final_data_table,1) n_new_detectors size(final_data_table,3)]));
            end
            % Insert new detector names into list
            [final_detectornames,detector_sort_index] = sort(cat(2,final_detectornames,detectornames_to_add));
            
            if not(isempty(final_data_table))
                final_data_table = final_data_table(:,detector_sort_index,:);
            end
        end
        
        % Take this scan's data matrix and rearrange values to match prior
        % scans' format:
        num_detectors = numel(data_import{1});
        [ data_table ] = move_values( data_import, x_values, detectorlist_minimal, final_detectornames , x_rounding_factor, final_x_values, num_detectors);
        
        % Tack this scan onto prior scans: 
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
                [detectorlist_minimal, NHEADERLINES] = extract_detectors(data,filename);
                
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
                    if strcmp(edgetype,'files/br') || strcmp(edgetype,'files/br')
                        x_values(x_values<Br_lower_spectral_limit)=[];
                        data_import=data_import(x_values>=Br_lower_spectral_limit);
                    elseif strcmp(edgetype,'files/pb') || strcmp(edgetype,'files/pb')
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
                
                % Take this scan's data matrix and rearrange values to match prior
                % scans' format:
                num_detectors = numel(data_import{1});
                [ data_table ] = move_values( data_import, x_values, detectorlist_minimal, final_detectornames , x_rounding_factor, final_x_values, num_detectors);
                
                if isempty(final_data_table)
                    final_data_table = data_table;
                else
                    final_data_table(:,:,counter)=data_table(:,:);
                end
            end
        else
        end
        if fid>0
            fclose(fid);
        else
        end
    end
end

%% Find I0 and detector indexes:
[ detectorindices ] = find_detector_indices( final_detectornames, photoncounting );

%% Split 6-second data into 2 separate scans of 3-second data:
if not(photoncounting)
    % Find runs that have 6-second data:
    scanswith6 = find(squeeze(sum(final_data_table(:,detectorindices.TFY_LaserOFF6sec,:)~=0 & isfinite(final_data_table(:,detectorindices.TFY_LaserOFF6sec,:)))));
    % Concatenate the data set with another copy of the 6-second data scan:
    final_data_table=cat(3,final_data_table,final_data_table(:,:,scanswith6));
    % Overwrite original 3-second data with 2*6-second minus 3-second data:
    final_data_table(:,detectorindices.TFY_LaserOFF,scanswith6) = 2*final_data_table(:,detectorindices.TFY_LaserOFF6sec,scanswith6)-final_data_table(:,detectorindices.TFY_LaserOFF3sec,scanswith6);
    final_data_table(:,detectorindices.TFY_LaserON,scanswith6) = 2*final_data_table(:,detectorindices.TFY_LaserON6sec,scanswith6)-final_data_table(:,detectorindices.TFY_LaserON3sec,scanswith6);
    final_data_table(:,detectorindices.HERFD_LaserOFF,scanswith6) = 2*final_data_table(:,detectorindices.HERFD_LaserOFF6sec,scanswith6)-final_data_table(:,detectorindices.HERFD_LaserOFF3sec,scanswith6);
    final_data_table(:,detectorindices.HERFD_LaserON,scanswith6) = 2*final_data_table(:,detectorindices.HERFD_LaserON6sec,scanswith6)-final_data_table(:,detectorindices.HERFD_LaserON3sec,scanswith6);
end

%% Reverse sign of the detector for integrating detectors:
if not(photoncounting)
    final_data_table(:,detectorindices.TFY_LaserOFF,:) = -(final_data_table(:,detectorindices.TFY_LaserOFF,:));
    final_data_table(:,detectorindices.TFY_LaserON,:)  = -(final_data_table(:,detectorindices.TFY_LaserON,:));
elseif photoncounting
    final_data_table(:,detectorindices.TFY_LaserOFF_bak,:)=-final_data_table(:,detectorindices.TFY_LaserOFF_bak,:); % Only flips sign of the backup data used when photon counting data is missing.
    final_data_table(:,detectorindices.TFY_LaserON_bak,:)=-final_data_table(:,detectorindices.TFY_LaserON_bak,:);
    
    %% If photon counting is missing for some runs, scale and copy integrated data instead:
    for i=1:size(final_data_table,3) % Loop over all scans...
        if sum(final_data_table(:,detectorindices.TFY_LaserOFF,i))==0
            final_data_table(:,detectorindices.TFY_LaserOFF,i) = final_data_table(:,detectorindices.TFY_LaserOFF_bak,i); % Replace photon-counting data with integrated data and scale.
            final_data_table(:,detectorindices.TFY_LaserON,i) = final_data_table(:,detectorindices.TFY_LaserON_bak,i);
        end
    end
end

%% Complicated I_zero averaging: Average together I_zeros for runs with the exact same x-values and apply averaged I_zero to those runs
if not(i0_correction_on)
    disp('Warning: I0 correction is off')
else
    [ final_data_table ] = produce_I0_averages( final_data_table, detectorindices, scantype, final_x_values, file_list );
end

%% Normalize data
if strcmp(scantype,'Timescan')
    % Correction for glitch in the phase shifter!
    PhaseShift = squeeze(final_data_table(:,detectorindices.PhaseShifter,1));
    final_x_values(PhaseShift<2560) = final_x_values(PhaseShift<2560) + 40e-12;
    
    %Correct for time-zero
    timezero = 1.6978e-7; % First reading: 169.76ns, second reading: 169.76ns
    final_x_values = final_x_values - timezero;

    % Do baseline correction:
    [ final_data_table, discardedruns ] = correct_baseline( baselinecorrection_on, scantype, final_data_table, detectorindices,final_x_values, 0 );
    
    % Then scale the data
    if normalization_on
        normalization_start_index=1; normalization_stop_index=length(final_x_values);
        [ final_data_table ] = normalize_data( final_data_table, detectorindices, normalization_start_index, normalization_stop_index);
    end
    
elseif strcmp(scantype,'Spectrum')
    if strcmp(edgetype,'files/br') || strcmp(edgetype,'files/br')
        normalization_start_index = find(final_x_values>13.52,1);
        normalization_stop_index = find(final_x_values>13.53,1);
        baseline_indices = find(final_x_values>13.46,1);
    elseif strcmp(edgetype,'files/pb') || strcmp(edgetype,'files/pb')
        normalization_start_index = find(final_x_values>13.045,1);
        normalization_stop_index = find(final_x_values>13.055,1);
        baseline_indices = find(final_x_values>13.02,1);
    end
    
    % Do baseline correction first:
    if not(isempty(baseline_indices))
        [ final_data_table,~ ] = correct_baseline( baselinecorrection_on, scantype, final_data_table,detectorindices,final_x_values, baseline_indices);
    end
    
    % Then scale the data:
    if normalization_on
        [ final_data_table ] = normalize_data( final_data_table, detectorindices, normalization_start_index, normalization_stop_index );
    end
end

%% Find and remove data points where the X-ray intensity dropped:
% Make a histogram at each x-value and delete remarkably low points:
if dropbadframes
    for current_x = 1:length(final_x_values)
        current_vals = squeeze(final_data_table(current_x,detectorindices.TFY_LaserOFF,:));
        % First drop the lowest value:
        %     [~,lowest_current] = min(current_vals); current_vals(lowest_current) = mean(current_vals,'omitnan');
        %     figure(12);subplot(2,1,1);hist(current_vals,100);
        x_ray_dropped_cutoff = mean(current_vals,'omitnan') - 3*std(current_vals,'omitnan');
        dropped_points{current_x} = find(current_vals<x_ray_dropped_cutoff);
        final_data_table(current_x,:,dropped_points{current_x})=NaN(1);
        %     subplot(2,1,2);hist(squeeze(final_data_table(current_x,detectorindices.TFY_LaserOFF,:)),100);
        %     pause(1)
    end
end

%% Adjust pump-probe magnitude to counteract laser drift (or at least prepare to use it later:
if normalize_pp
    pp_norm_term = squeeze(mean(abs(final_data_table(:,detectorindices.TFY_LaserON,:)-final_data_table(:,detectorindices.TFY_LaserOFF,:)),'omitnan'));
    pp_norm_term = pp_norm_term/mean(pp_norm_term,'omitnan');
    for i = 1:size(final_data_table,3)
        final_data_table(:,detectorindices.TFY_LaserOFF,i) = final_data_table(:,detectorindices.TFY_LaserOFF,i)./pp_norm_term(i);
        final_data_table(:,detectorindices.TFY_LaserON,i) = final_data_table(:,detectorindices.TFY_LaserON,i)./pp_norm_term(i);
    end
end

%% Examine the noise level of each scan:
% noise_laserON = sum(abs(diff(squeeze(final_data_table(:,detectorindices.TFY_LaserOFF,:)))),'omitnan');
% noise_laserOFF = sum(abs(diff(squeeze(final_data_table(:,detectorindices.TFY_LaserON,:)))),'omitnan');
% noise_normalizer = max(abs(squeeze(final_data_table(:,detectorindices.TFY_LaserON,:)-final_data_table(:,detectorindices.TFY_LaserOFF,:))));
% noisiness = (noise_laserON+noise_laserOFF)./noise_normalizer;
% % final_data_table(:,:,noisiness>300) = [];
% figure(5);hist(noisiness,20);title('Histogram of run noisiness');

%% Averaging of all available data:
% Remove x-values where there's no TFY data due to a detector or naming glitch:
TFY_present = find(sum(final_data_table(:,detectorindices.TFY_LaserOFF,:),3));
final_data_table = final_data_table(TFY_present,:,:);
final_x_values = final_x_values(TFY_present);

% average_data_table=mean(final_data_table,3);
n_3s_scans_per_point = squeeze(sum(final_data_table(:,detectorindices.TFY_LaserOFF3sec,:)~=0 & isfinite(final_data_table(:,detectorindices.TFY_LaserOFF3sec,:)),3)); % Only count nonzero values.
n_6s_scans_per_point = squeeze(sum(final_data_table(:,detectorindices.TFY_LaserOFF6sec,:)~=0 & isfinite(final_data_table(:,detectorindices.TFY_LaserOFF6sec,:)),3)); % Only count nonzero values.
disp([num2str(max(n_3s_scans_per_point)),' scans include 3s integration, ',num2str(max(n_6s_scans_per_point/2)),' of which are split 6s integration.'])

n_points_per_scan = squeeze(sum(final_data_table(:,detectorindices.TFY_LaserOFF,:)~=0 & isfinite(final_data_table(:,detectorindices.TFY_LaserOFF,:)),1)); % Only count nonzero values.

if ismember(0,n_points_per_scan)
    disp('Some runs have no data being used! Please check data for NaNs')
end

% For runs that have 6-second integration, double-count it and use that instead:
n_scans_per_point = sum(final_data_table(:,detectorindices.TFY_LaserOFF,:)~=0 & isfinite(final_data_table(:,detectorindices.TFY_LaserOFF,:)),3); % Only count nonzero values.

% Average all scans equally:
average_data_table=sum(final_data_table,3,'omitnan')./repmat(n_scans_per_point,[1 size(final_data_table,2)]);

TFY_total_ONandOFF = .5*nonzeros(average_data_table(:,detectorindices.TFY_LaserON)+average_data_table(:,detectorindices.TFY_LaserOFF));
TFY_pump_probe_avg = (average_data_table(:,detectorindices.TFY_LaserON)-average_data_table(:,detectorindices.TFY_LaserOFF));
HERFD_pump_probe_avg = average_data_table(:,detectorindices.HERFD_LaserON)-average_data_table(:,detectorindices.HERFD_LaserOFF);

%% For timescans, ensure signs are flipped when observing a bleach feature.
if strcmp(scantype,'Timescan')
    TFY_pump_probe_avg = abs(TFY_pump_probe_avg);
end
    
%% Filter out frequencies in the data higher than monochromator resolution:
try 
    if strcmp(scantype,'Spectrum')
        [TFY_pump_probe_avg] = monochromator_lowpass(TFY_pump_probe_avg,final_x_values,monochromator_res_eV);
    end
catch
    disp('Sampling rate was not high enough to remove high-frequency noise.');
end

%% Also get error for error bars:
finaltable_forerrors = final_data_table; finaltable_forerrors(final_data_table==0) = NaN(1);
if photoncounting
    TFY_pump_probe_error = std( finaltable_forerrors(:,detectorindices.TFY_LaserON,:)-finaltable_forerrors(:,detectorindices.TFY_LaserOFF,:),0,3,'omitnan' );
    HERFD_pump_probe_error = std( finaltable_forerrors(:,detectorindices.HERFD_LaserON,:)-finaltable_forerrors(:,detectorindices.HERFD_LaserOFF,:),0,3,'omitnan' );
else
    TFY_pump_probe_error = std( finaltable_forerrors(:,detectorindices.TFY_LaserON,:)-finaltable_forerrors(:,detectorindices.TFY_LaserOFF,:),0,3,'omitnan' );
    HERFD_pump_probe_error = std( finaltable_forerrors(:,detectorindices.HERFD_LaserON,:)-finaltable_forerrors(:,detectorindices.HERFD_LaserOFF,:),0,3,'omitnan' );
end

%% Do a final test of all data to make sure TFY exists:
discardedruns=0;
for i=1:size(final_data_table,3)
    totalTFY(i) = sum(final_data_table(:,detectorindices.TFY_LaserOFF,i),'omitnan');
    if totalTFY(i)==0;discardedruns = discardedruns+1;end
end

%% Plot all the individual scans, TFY and HERFD:
legend_plot = {};
if strcmp(scantype,'Spectrum')
    % Individual TFY plots
    if plotindividually
        figure(1);clf;hold off;
        for i=1:size(final_data_table,3)
            if totalTFY(i)>0 % Don't bother to plot runs without data.
                x_for_plotting = final_x_values(final_data_table(:,detectorindices.TFY_LaserON,i)~=0);
                x_for_plotting_pp = final_x_values((final_data_table(:,detectorindices.TFY_LaserON,i)-final_data_table(:,detectorindices.TFY_LaserOFF,i))~=0);
                avgONandOFF = .5*nonzeros(final_data_table(:,detectorindices.TFY_LaserON,i)+final_data_table(:,detectorindices.TFY_LaserOFF,i));
                transient_spectrum_y = nonzeros(final_data_table(:,detectorindices.TFY_LaserON,i)-final_data_table(:,detectorindices.TFY_LaserOFF,i));
                subplot(2,2,1); hold on; plot(x_for_plotting(~isnan(avgONandOFF)),avgONandOFF(~isnan(avgONandOFF)),'Linewidth',2);
                subplot(2,2,3); hold on; plot(x_for_plotting_pp(~isnan(transient_spectrum_y)),transient_spectrum_y(~isnan(transient_spectrum_y)),'Linewidth',2);
                legend_plot{length(legend_plot)+1}=sprintf('scan%1.0f',i);
            end
        end
        subplot(2,2,1);title(['Single scans, TFY, ',scandetails]); xlabel('Energy (keV)');ylabel('Fluorescence (A.U.)'); %legend(legend_plot);
        subplot(2,2,3);title(['Single scans, TFY, ',scandetails,' pump-probe transient']); xlabel('Energy (keV)'); ylabel('\Delta Fluorescence (A.U.)'); %legend(legend_plot); 
    end
elseif strcmp(scantype,'Timescan')
    if plotindividually
        figure(1);clf;hold off;
        for i=1:size(final_data_table,3)
            %             plot(-final_x_values(final_data_table(:,detectorindices.TFY_LaserON,i)~=0),final_data_table(:,detectorindices.TFY_LaserON,i)-final_data_table(:,detectorindices.TFY_LaserOFF,i),'Linewidth',2);
            plot(-final_x_values,final_data_table(:,detectorindices.TFY_LaserON,i)-final_data_table(:,detectorindices.TFY_LaserOFF,i),'Linewidth',2);
            legend_plot{i}=sprintf('scan%1.0f',i); hold on;
        end
        title(['Single scans, TFY, ',num2str(scandetails)]); xlabel('Time (sec)'); ylabel('\Delta Fluorescence (A.U.)'); %legend(legend_plot);
    end
else
    error('Unidentified scan type. Does not match "Timescan" or "Spectrum".');
end

%% Plot the averaged data
if strcmp(scantype,'Spectrum')
    figure(1);
    subplot(2,2,2);
    plot(nonzeros(final_x_values),TFY_total_ONandOFF,'Linewidth',2);
    title(['Average of ',num2str(counter),' scans, ',num2str(scandetails)]); xlabel('Energy (keV)');ylabel('Fluorescence (A.U.)');
    subplot(2,2,4);
    errorbar(final_x_values,TFY_pump_probe_avg,TFY_pump_probe_error,'Linewidth',2);
%     plot(final_x_values,-TFY_pump_probe_avg,'Linewidth',2);
    xlabel('Energy (keV)'); ylabel('\Delta Fluorescence (A.U.)');
elseif strcmp(scantype,'Timescan')
    figure(2);subplot(1,1,1); % Fitting has to be done in picoseconds, since MATLAB doesn't like dealing with tiny numbers.
    plottable_x = -(final_x_values(~isnan(TFY_pump_probe_avg)))*1e12;% Note that data is now in units of picoseconds because that makes optimization better (stopping conditions)
    [~,y_max] = max(abs(TFY_pump_probe_avg));
    plottable_y = TFY_pump_probe_avg(~isnan(TFY_pump_probe_avg))./TFY_pump_probe_avg(y_max);
    %% Fit timescans with biexponential decay (Results from MicroXAS were published as [64percent 542ps arbitrary_scaling_factor 104ns timezero_adjustment]):
    CsPbBr3_starting_fit_coeffss = [.64 542 1 104e3 20];
%     bestfit = fit(plottable_x,plottable_y,'c*(a*exp(-(x)/b) + (1-a)*exp(-(x)/d))*heaviside(x)','StartPoint',CsPbBr3_starting_fit_coeffss)
    bestfit = fit(plottable_x,plottable_y,'c*(a*exp(-(x+e)/b) + (1-a)*exp(-(x+e)/d))*(.5)*(1+erf((x+e)/(sqrt(2)*(80 / 2.355))))','StartPoint',CsPbBr3_starting_fit_coeffss);
    bestfit_pluserror = fit(plottable_x,plottable_y+TFY_pump_probe_error./TFY_pump_probe_avg(y_max),'c*(a*exp(-(x+e)/b) + (1-a)*exp(-(x+e)/d))*(1/2)*(1+erf((x+e)/(sqrt(2)*(80 / 2.355))))','StartPoint',CsPbBr3_starting_fit_coeffss);
    bestfit_minuserror = fit(plottable_x,plottable_y-TFY_pump_probe_error./TFY_pump_probe_avg(y_max),'c*(a*exp(-(x+e)/b) + (1-a)*exp(-(x+e)/d))*(1/2)*(1+erf((x+e)/(sqrt(2)*(80 / 2.355))))','StartPoint',CsPbBr3_starting_fit_coeffss);
    disp(['Time constants are ',num2str(bestfit.b),'ps and ',num2str(bestfit.d*1e-3),'ns.']);
    plot(bestfit,plottable_x,plottable_y,'.'); legend('TFY Laser ON - Laser OFF','Fitted Biexponential Decay');
    hold on;
    errorbar(plottable_x,plottable_y,TFY_pump_probe_error./TFY_pump_probe_avg(y_max),'.');%./TFY_pump_probe_avg(y_max)
%     plot(plottable_x,plottable_y,'.');
    hold off;
    title([num2str(scandetails),' Average of ',num2str(counter),' scans']); xlabel('Time (picoseconds)'); ylabel('\Delta Fluorescence (A.U.)');    
end

%% Plot the average time-resolved HERFD
if strcmp(scantype,'Spectrum') && (strcmp(edgetype,'files\pb') || strcmp(edgetype,'files/pb'))
    figure(3);clf;hold on;
    HERFD_sig_ON = average_data_table(:,detectorindices.HERFD_LaserON);
    HERFD_sig_OFF = average_data_table(:,detectorindices.HERFD_LaserOFF);
    %         errorbar(nonzeros(final_x_values),-nonzeros(average_data_table(:,detectorindices.HERFD_LaserON)-average_data_table(:,detectorindices.HERFD_LaserOFF)),HERFD_pump_probe_error,'Linewidth',2);
    plot(final_x_values(HERFD_pump_probe_avg~=0),nonzeros(HERFD_pump_probe_avg),'Linewidth',2);
    plot(final_x_values(HERFD_sig_ON~=0),nonzeros(HERFD_sig_ON),'Linewidth',2);
    plot(final_x_values(HERFD_sig_OFF~=0),nonzeros(HERFD_sig_OFF),'Linewidth',2);
    title([num2str(scandetails),' Average of ',num2str(counter),' scans']); xlabel('Energy (keV)'); ylabel('Partial Fluorescence Yield (A.U.)');
    legend('Time resolved HERFD','HERFD, Laser ON','HERFD, Laser OFF');
    figure(4);clf;
    hold off;
    tfy_scaler=10;
    plot(final_x_values(HERFD_pump_probe_avg~=0),smooth(nonzeros(HERFD_pump_probe_avg),1),final_x_values,tfy_scaler*TFY_pump_probe_avg,'Linewidth',2);
%     errorbar(final_x_values(HERFD_pump_probe_avg~=0),nonzeros(HERFD_pump_probe_avg),HERFD_pump_probe_error(HERFD_pump_probe_avg~=0),'Linewidth',2);
    legend('pump-probe HERFD',[num2str(tfy_scaler),' * pump-probe TFY']);title([scandetails,' ',sprintf('Average of %1.0f scans',counter)]); xlabel('Energy (keV)');
end

if discardedruns~=0 && photoncounting
    disp(['Warning! ',num2str(discardedruns),' out of ',num2str(size(final_data_table,3)),' runs have been discarded in photon counting mode!'])
end