%% Script to pull Ocean Color data
% Ocean Color Data Download with Progress Bar and Time Estimation
% Downloads NOAA VIIRS/OLCI DINEOF gap-filled ocean color data
% Variables: Chlorophyll-a, Kd490, and Suspended Particulate Matter
% Daily data via NOAA CoastWatch ERDDAP
%   VARIABLE-BY-VARIABLE approach (server-friendly)
%
% https://coastwatch.noaa.gov/erddap/griddap/noaacwNPPN20S3AkdSCIDINEOF2kmDaily.graph?kd_490%5B(2024-01-19T12:00:00Z)%5D%5B(0.0)%5D%5B(18.880205):(17.786455)%5D%5B(-66.005215):(-64.911465)%5D&.draw=surface&.vars=longitude%7Clatitude%7Ckd_490&.colorBar=%7C%7CLog%7C0.001%7C100%7C&.bgColor=0xffccccff

clear;clc

%% ========== CONFIGURATION SECTION ==========
% Define your date range here
START_YEAR = 2018;
START_MONTH = 1;        % 1-12
END_YEAR = 2018;
END_MONTH = 11;         % 1-12

% File handling option
KEEP_FILES = true;     % Set to true to keep downloaded files, false to delete them

% Puerto Rico/USVI region coordinates (matching other datasets)
lat_min = 17.0;   % Southern boundary
lat_max = 19.5;   % Northern boundary  
lon_min = -68.0;  % Western boundary
lon_max = -64.0;  % Eastern boundary

% Depth (surface level)
depth_str = '(0.0)';

% Define the three datasets - reordered to download Kd490 first as requested
datasets = {
    'noaacwNPPN20S3AkdSCIDINEOF2kmDaily', 'kd_490', 'Diffuse Attenuation Coefficient at 490nm';
    'noaacwNPPN20S3ASCIDINEOF2kmDaily', 'chlor_a', 'Chlorophyll-a Concentration';
    'noaacwNPPN20S3AspmSCIDINEOF2kmDaily', 'spm', 'Suspended Particulate Matter'
};
%% =============================================

% Get the project root directory
projectPath = matlab.project.rootProject().RootFolder;

% Define paths relative to the project root
dataPath = fullfile(projectPath, 'data');
srcPath = fullfile(projectPath, 'src');
outputPath = fullfile(projectPath, 'output/ocean_color');
tempPath = fullfile(projectPath, 'temp');

%% Calculate total months and generate year/month list
% Calculate total number of months in the date range
total_months = (END_YEAR - START_YEAR) * 12 + (END_MONTH - START_MONTH + 1);

% Generate list of all year-month combinations
year_month_list = [];
for year = START_YEAR:END_YEAR
    start_month_for_year = (year == START_YEAR) * START_MONTH + (year ~= START_YEAR) * 1;
    end_month_for_year = (year == END_YEAR) * END_MONTH + (year ~= END_YEAR) * 12;
    
    for month = start_month_for_year:end_month_for_year
        year_month_list = [year_month_list; year, month];
    end
end

% Initialize data arrays for each variable
kd490_data = [];
chlor_a_data = [];
spm_data = [];
times_data = [];

% Initialize timing variables
start_time = tic;
operation_times = [];

% Calculate total operations for progress tracking
total_operations = size(datasets, 1) * total_months;
current_operation = 0;

% Create progress bar
h = waitbar(0, 'Starting Ocean Color data download...', 'Name', 'Ocean Color Download Progress');

fprintf('Downloading Ocean Color data for Puerto Rico/USVI region:\n');
fprintf('Lat: %.1f to %.1f, Lon: %.1f to %.1f\n', lat_min, lat_max, lon_min, lon_max);
fprintf('Time period: %d-%02d to %d-%02d\n', START_YEAR, START_MONTH, END_YEAR, END_MONTH);
fprintf('Total months to download: %d\n', total_months);
fprintf('Variables: Kd490, Chlorophyll-a, Suspended Particulate Matter\n');
fprintf('Strategy: Complete each variable before moving to next (server-friendly)\n');
if KEEP_FILES
    fprintf('File handling: Keeping files\n');
else
    fprintf('File handling: Deleting files after processing\n');
end

%% Main download loop - VARIABLE BY VARIABLE
try
    fprintf('\n=== VARIABLE-BY-VARIABLE DOWNLOAD APPROACH ===\n');
    
    for var_idx = 1:size(datasets, 1)
        dataset_id = datasets{var_idx, 1};
        variable = datasets{var_idx, 2};
        var_name = datasets{var_idx, 3};
        
        fprintf('\n*** STARTING COMPLETE DOWNLOAD OF %s (%s) ***\n', upper(var_name), variable);
        fprintf('Will download %d months of %s data before moving to next variable\n', total_months, variable);
        
        % Initialize storage for this variable
        current_var_data = [];
        current_var_times = [];
        
        %% Download ALL months for this specific variable
        for i = 1:size(year_month_list, 1)
            year = year_month_list(i, 1);
            month = year_month_list(i, 2);
            current_operation = current_operation + 1;
            
            operation_start_time = tic;
            
            % Update progress bar
            progress = current_operation / total_operations;
            status_msg = sprintf('Variable %d/3: %s\nMonth %d/%d: %d-%02d\nOverall: %d/%d operations', ...
                var_idx, variable, i, total_months, year, month, current_operation, total_operations);
            
            % Add time estimation
            if current_operation > 1
                avg_time_per_operation = mean(operation_times);
                remaining_operations = total_operations - current_operation;
                est_remaining_time = avg_time_per_operation * remaining_operations;
                
                status_msg = sprintf('%s\nEst. remaining: %.1f minutes', ...
                    status_msg, est_remaining_time / 60);
            end
            
            waitbar(progress, h, status_msg);
            
            % Construct date range for the full month (like SST script)
            start_date = sprintf('%d-%02d-01T12:00:00Z', year, month);

            % Calculate end date for the month
            days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
            % Handle leap year for February
            if month == 2 && mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)
                days_in_month(2) = 29;
            end

            last_day = days_in_month(month);
            end_date = sprintf('%d-%02d-%02dT12:00:00Z', year, month, last_day);

            fprintf('  Downloading %s for %d-%02d (operation %d/%d)...\n', ...
                variable, year, month, current_operation, total_operations);
            
            % Construct URL for ocean color data
            base_url = 'https://coastwatch.noaa.gov/erddap/griddap/';
            query_str = sprintf('%s[(%s):(%s)][%s][(%.1f):(%.1f)][(%.1f):(%.1f)]', ...
                variable, start_date, end_date, depth_str, lat_min, lat_max, lon_min, lon_max);
            url = sprintf('%s%s.nc?%s', base_url, dataset_id, query_str);
            
            % Download with enhanced retry logic
            filename = fullfile(outputPath, sprintf('%s_%d_%02d.nc', variable, year, month));
            download_success = false;
            retry_count = 0;
            max_retries = 15;
            
            options = weboptions('Timeout', 60, 'RequestMethod', 'get');
            
            while ~download_success && retry_count < max_retries
                retry_count = retry_count + 1;
                
                try
                    fprintf('    Attempt %d/%d...\n', retry_count, max_retries);
                    websave(filename, url, options);
                    download_success = true;
                    fprintf('    Downloaded: %s\n', filename);
                    
                catch ME
                    if retry_count < max_retries
                        fprintf('    Download failed (attempt %d/%d): %s\n', retry_count, max_retries, ME.message);
                        
                        if contains(ME.message, '503') || contains(ME.message, 'Service Unavailable')
                            fprintf('    Server temporarily unavailable (503 error)\n');
                        elseif contains(ME.message, 'timeout') || contains(ME.message, 'Timeout')
                            fprintf('    Request timed out\n');
                        elseif contains(ME.message, '404') || contains(ME.message, 'Not Found')
                            fprintf('    Data not found for this date (404 error)\n');
                        else
                            fprintf('    Network error: %s\n', ME.message);
                        end
                        
                        fprintf('    Waiting 15 seconds before retry...\n');
                        pause(15);
                    else
                        warning('Failed to download %s for %d-%02d after %d attempts: %s', ...
                            variable, year, month, max_retries, ME.message);
                    end
                end
            end
            
            if ~download_success
                fprintf('    SKIPPING %s for %d-%02d (download failed)\n', variable, year, month);
                operation_time = toc(operation_start_time);
                operation_times(end+1) = operation_time;
                continue; % Skip to next month for this variable
            end
            
            % Read and process the data
            try
                temp_data = ncread(filename, variable);
                temp_time_data = ncread(filename, 'time');
                
                % Check for fill values
                try
                    fill_value = ncreadatt(filename, variable, '_FillValue');
                    fprintf('    Fill value: %.2f\n', fill_value);
                catch
                    fill_value = NaN;
                end
                
                % Apply fill value handling
                temp_data(temp_data == fill_value) = NaN;
                
                % Data quality check
                sample_val = temp_data(~isnan(temp_data));
                if ~isempty(sample_val)
                    sample_mean = mean(sample_val(1:min(100, length(sample_val))));
                    fprintf('    Sample mean %s: %.4f\n', variable, sample_mean);
                else
                    fprintf('    Warning: No valid %s data found for this month\n', variable);
                end
 
                % Check actual dimensions received
                fprintf('    Data dimensions: %s\n', mat2str(size(temp_data)));
                
                % Squeeze out the altitude dimension (dimension 3)
                temp_data = squeeze(temp_data);
                fprintf('    After squeeze: %s\n', mat2str(size(temp_data)));
                
                % Concatenate data for this variable
                if i == 1
                    % First month for this variable - initialize
                    current_var_data = temp_data;
                    current_var_times = temp_time_data(:);
                else
                    % Subsequent months - concatenate along time dimension (3rd dimension after squeeze)
                    current_var_data = cat(3, current_var_data, temp_data);
                    current_var_times = [current_var_times; temp_time_data(:)];
                end
                
                fprintf('    Month %d-%02d: Added %d time steps (Total so far: %d)\n', ...
                    year, month, size(temp_data, 3), size(current_var_data, 3));
                
                % Read coordinates on very first successful download
                if var_idx == 1 && i == 1
                    latitude = ncread(filename, 'latitude');
                    longitude = ncread(filename, 'longitude');
                    altitude = ncread(filename, 'altitude');
                    fprintf('    Grid dimensions: %d x %d\n', length(longitude), length(latitude));
                    fprintf('    Altitude: %.1f m\n', altitude);
                end
                
            catch ME
                warning('Failed to read data from %s: %s', filename, ME.message);
            end
            
            % File handling
            if KEEP_FILES
                fprintf('    Keeping file: %s\n', filename);
            else
                if exist(filename, 'file')
                    delete(filename);
                end
            end
            
            % Record timing
            operation_time = toc(operation_start_time);
            operation_times(end+1) = operation_time;
            
            fprintf('    Completed in %.1f seconds\n', operation_time);
        end
        
        %% Store the completed variable data
        if ~isempty(current_var_data)
            if strcmp(variable, 'kd_490')
                kd490_data = current_var_data;
                if isempty(times_data)
                    times_data = current_var_times;
                end
                fprintf('*** COMPLETED KD490: %d time steps ***\n', size(kd490_data, 3));
            elseif strcmp(variable, 'chlor_a')
                chlor_a_data = current_var_data;
                if isempty(times_data)
                    times_data = current_var_times;
                end
                fprintf('*** COMPLETED CHLOROPHYLL-A: %d time steps ***\n', size(chlor_a_data, 3));
            elseif strcmp(variable, 'spm')
                spm_data = current_var_data;
                if isempty(times_data)
                    times_data = current_var_times;
                end
                fprintf('*** COMPLETED SPM: %d time steps ***\n', size(spm_data, 3));
            end
        else
            fprintf('*** WARNING: No data collected for %s ***\n', variable);
        end
        
        % Brief pause between variables to be server-friendly
        if var_idx < size(datasets, 1)
            fprintf('\nPausing 30 seconds before next variable (server-friendly)...\n');
            for countdown = 30:-1:1
                fprintf('  %d...', countdown);
                if mod(countdown, 10) == 1 || countdown <= 5
                    fprintf('\n');
                end
                pause(1);
            end
            fprintf('Resuming downloads...\n');
        end
    end
    
    % Final progress update
    total_elapsed = toc(start_time);
    waitbar(1, h, sprintf('Download complete! Total time: %.1f minutes', total_elapsed/60));
    pause(2);
    close(h);
    
catch ME
    if ishandle(h)
        close(h);
    end
    rethrow(ME);
end



%% Clean up zero-byte files
fprintf('\n=== Checking for zero-byte files ===\n');
if exist(outputPath, 'dir')
    % Get all .nc files in the output directory
    nc_files = dir(fullfile(outputPath, '*.nc'));
    
    for i = 1:length(nc_files)
        file_path = fullfile(outputPath, nc_files(i).name);
        file_size = nc_files(i).bytes;
        
        if file_size == 0
            fprintf('Deleting zero-byte file: %s\n', nc_files(i).name);
            delete(file_path);
        end
    end
    fprintf('Zero-byte file check complete\n');
else
    fprintf('Output directory not found, skipping zero-byte check\n');
end



%% Data Processing and Statistics
fprintf('\n=== Processing Results ===\n');
fprintf('Total download time: %.2f minutes\n', total_elapsed/60);

% Check what data we successfully downloaded
variables_downloaded = {};
if ~isempty(kd490_data)
    fprintf('Kd490 data shape: %s\n', mat2str(size(kd490_data)));
    variables_downloaded{end+1} = 'kd490';
end
if ~isempty(chlor_a_data)
    fprintf('Chlorophyll-a data shape: %s\n', mat2str(size(chlor_a_data)));
    variables_downloaded{end+1} = 'chlor_a';
end
if ~isempty(spm_data)
    fprintf('SPM data shape: %s\n', mat2str(size(spm_data)));
    variables_downloaded{end+1} = 'spm';
end

fprintf('Time range: %d data points\n', length(times_data));
fprintf('Successfully downloaded variables: %s\n', strjoin(variables_downloaded, ', '));

% Create date range string for titles and filenames
date_range_str = sprintf('%d-%02d to %d-%02d', START_YEAR, START_MONTH, END_YEAR, END_MONTH);

%% Calculate statistics and create maps for each variable
fprintf('\n=== Creating Maps and Statistics ===\n');

% Process Kd490 if available
if ~isempty(kd490_data)
    fprintf('Processing Kd490 statistics...\n');
    
    % Calculate statistics (no squeeze needed - altitude already removed)
    mean_kd490 = mean(kd490_data, 3, 'omitnan');
    max_kd490 = max(kd490_data, [], 3, 'omitnan');
    min_kd490 = min(kd490_data, [], 3, 'omitnan');
    std_kd490 = std(kd490_data, 0, 3, 'omitnan');

    % Overall statistics
    overall_mean_kd490 = mean(kd490_data(:), 'omitnan');
    overall_max_kd490 = max(kd490_data(:), [], 'omitnan');
    overall_min_kd490 = min(kd490_data(:), [], 'omitnan');

    fprintf('Kd490: Mean: %.4f m-1, Max: %.4f m-1, Min: %.4f m-1\n', ...
        overall_mean_kd490, overall_max_kd490, overall_min_kd490);
    
    % Create maps
    figure('Position', [100, 100, 600, 500]);
    imagesc(longitude, latitude, mean_kd490');
    set(gca, 'YDir', 'normal');
    colorbar;
    colormap('jet');
    title(sprintf('Mean Diffuse Attenuation Coefficient at 490nm (%s)', date_range_str), 'FontSize', 14);
    xlabel('Longitude (°W)', 'FontSize', 12);
    ylabel('Latitude (°N)', 'FontSize', 12);
    clim([min(mean_kd490(:)) max(mean_kd490(:))]);
    grid on;
    saveas(gcf, fullfile(outputPath, 'mean_kd490_map.png'));
    
    fprintf('Kd490 map saved\n');
end

% Process Chlorophyll-a if available
if ~isempty(chlor_a_data)
    fprintf('Processing Chlorophyll-a statistics...\n');
    
    % Calculate statistics (no squeeze needed - altitude already removed)
    mean_chlor_a = mean(chlor_a_data, 3, 'omitnan');
    max_chlor_a = max(chlor_a_data, [], 3, 'omitnan');
    min_chlor_a = min(chlor_a_data, [], 3, 'omitnan');
    std_chlor_a = std(chlor_a_data, 0, 3, 'omitnan');    
    
    % Overall statistics
    overall_mean_chlor_a = mean(chlor_a_data(:), 'omitnan');
    overall_max_chlor_a = max(chlor_a_data(:), [], 'omitnan');
    overall_min_chlor_a = min(chlor_a_data(:), [], 'omitnan');

    fprintf('Chlorophyll-a: Mean: %.3f mg m-3, Max: %.3f mg m-3, Min: %.3f mg m-3\n', ...
        overall_mean_chlor_a, overall_max_chlor_a, overall_min_chlor_a);
    
    % Create maps
    figure('Position', [200, 200, 600, 500]);
    imagesc(longitude, latitude, mean_chlor_a');
    set(gca, 'YDir', 'normal');
    colorbar;
    colormap('jet');
    title(sprintf('Mean Chlorophyll-a Concentration (%s)', date_range_str), 'FontSize', 14);
    xlabel('Longitude (°W)', 'FontSize', 12);
    ylabel('Latitude (°N)', 'FontSize', 12);
    clim([min(mean_chlor_a(:)) max(mean_chlor_a(:))]);
    grid on;
    saveas(gcf, fullfile(outputPath, 'mean_chlor_a_map.png'));
    
    fprintf('Chlorophyll-a map saved\n');
end

% Process SPM if available
if ~isempty(spm_data)
    fprintf('Processing SPM statistics...\n');
    
    % Calculate statistics (no squeeze needed - altitude already removed)
    mean_spm = mean(spm_data, 3, 'omitnan');
    max_spm = max(spm_data, [], 3, 'omitnan');
    min_spm = min(spm_data, [], 3, 'omitnan');
    std_spm = std(spm_data, 0, 3, 'omitnan');

    % Overall statistics
    overall_mean_spm = mean(spm_data(:), 'omitnan');
    overall_max_spm = max(spm_data(:), [], 'omitnan');
    overall_min_spm = min(spm_data(:), [], 'omitnan');
    
    fprintf('SPM: Mean: %.3f g m-3, Max: %.3f g m-3, Min: %.3f g m-3\n', ...
        overall_mean_spm, overall_max_spm, overall_min_spm);
    
    % Create maps
    figure('Position', [300, 300, 600, 500]);
    imagesc(longitude, latitude, mean_spm');
    set(gca, 'YDir', 'normal');
    colorbar;
    colormap('jet');
    title(sprintf('Mean Suspended Particulate Matter (%s)', date_range_str), 'FontSize', 14);
    xlabel('Longitude (°W)', 'FontSize', 12);
    ylabel('Latitude (°N)', 'FontSize', 12);
    clim([min(mean_spm(:)) max(mean_spm(:))]);
    grid on;
    saveas(gcf, fullfile(outputPath, 'mean_spm_map.png'));
    
    fprintf('SPM map saved\n');
end

%% Export Data for R Analysis
fprintf('\n=== Exporting Data for R Analysis ===\n');

% Only export if we have at least one successful variable and coordinates
if ~isempty(variables_downloaded) && exist('longitude', 'var') && exist('latitude', 'var')
    
    % Create coordinate grids
    [lon_grid, lat_grid] = meshgrid(longitude, latitude);
    
    % Flatten grids for long format
    lon_vec = lon_grid(:);
    lat_vec = lat_grid(:);
    
    % Initialize summary table with coordinates
    summary_data = table(lon_vec, lat_vec, 'VariableNames', {'longitude', 'latitude'});
    
    % Add each available variable to the summary table
    if ~isempty(kd490_data)
        mean_kd490_t = mean_kd490';
        mean_kd490_vec = mean_kd490_t(:);
        max_kd490_t = max_kd490';
        max_kd490_vec = max_kd490_t(:);
        min_kd490_t = min_kd490';
        min_kd490_vec = min_kd490_t(:);
        std_kd490_t = std_kd490';
        std_kd490_vec = std_kd490_t(:);
        
        summary_data.mean_kd490 = mean_kd490_vec;
        summary_data.max_kd490 = max_kd490_vec;
        summary_data.min_kd490 = min_kd490_vec;
        summary_data.std_kd490 = std_kd490_vec;
    end
    
    if ~isempty(chlor_a_data)
        mean_chlor_a_t = mean_chlor_a';
        mean_chlor_a_vec = mean_chlor_a_t(:);
        max_chlor_a_t = max_chlor_a';
        max_chlor_a_vec = max_chlor_a_t(:);
        min_chlor_a_t = min_chlor_a';
        min_chlor_a_vec = min_chlor_a_t(:);
        std_chlor_a_t = std_chlor_a';
        std_chlor_a_vec = std_chlor_a_t(:);
        
        summary_data.mean_chlor_a = mean_chlor_a_vec;
        summary_data.max_chlor_a = max_chlor_a_vec;
        summary_data.min_chlor_a = min_chlor_a_vec;
        summary_data.std_chlor_a = std_chlor_a_vec;
    end
    
    if ~isempty(spm_data)
        mean_spm_t = mean_spm';
        mean_spm_vec = mean_spm_t(:);
        max_spm_t = max_spm';
        max_spm_vec = max_spm_t(:);
        min_spm_t = min_spm';
        min_spm_vec = min_spm_t(:);
        std_spm_t = std_spm';
        std_spm_vec = std_spm_t(:);
        
        summary_data.mean_spm = mean_spm_vec;
        summary_data.max_spm = max_spm_vec;
        summary_data.min_spm = min_spm_vec;
        summary_data.std_spm = std_spm_vec;
    end
    
    % Export summary CSV
    writetable(summary_data, fullfile(outputPath, 'ocean_color_summary_for_R.csv'));
    fprintf('Ocean color summary saved: ocean_color_summary_for_R.csv\n');
    
    % Export coordinate information
    lon_table = table(longitude(:), 'VariableNames', {'longitude'});
    lat_table = table(latitude(:), 'VariableNames', {'latitude'});
    writetable(lon_table, fullfile(outputPath, 'ocean_color_longitude.csv'));
    writetable(lat_table, fullfile(outputPath, 'ocean_color_latitude.csv'));
    
    % Grid info
    grid_info = table({'longitude'; 'latitude'; 'nlon'; 'nlat'}, ...
        {min(longitude); min(latitude); length(longitude); length(latitude)}, ...
        {max(longitude); max(latitude); length(longitude); length(latitude)}, ...
        'VariableNames', {'dimension', 'min_value', 'max_value'});
    writetable(grid_info, fullfile(outputPath, 'ocean_color_grid_info.csv'));
    
    fprintf('Coordinate files saved\n');
    
    %% 3. Export Individual Summary Matrices (for direct raster creation in R)
    fprintf('Exporting individual summary matrices...\n');
    
    % Create a function to save matrices in R-readable format
    save_matrix_for_R = @(matrix, filename) ...
        writematrix(matrix', fullfile(outputPath, [filename '.csv']));
    
    % Export Kd490 matrices if available
    if ~isempty(kd490_data)
        save_matrix_for_R(mean_kd490, 'mean_kd490_matrix');
        save_matrix_for_R(max_kd490, 'max_kd490_matrix');
        save_matrix_for_R(min_kd490, 'min_kd490_matrix');
        save_matrix_for_R(std_kd490, 'std_kd490_matrix');
        fprintf('Kd490 matrices saved (4 files)\n');
    end
    
    % Export Chlorophyll-a matrices if available
    if ~isempty(chlor_a_data)
        save_matrix_for_R(mean_chlor_a, 'mean_chlor_a_matrix');
        save_matrix_for_R(max_chlor_a, 'max_chlor_a_matrix');
        save_matrix_for_R(min_chlor_a, 'min_chlor_a_matrix');
        save_matrix_for_R(std_chlor_a, 'std_chlor_a_matrix');
        fprintf('Chlorophyll-a matrices saved (4 files)\n');
    end
    
    % Export SPM matrices if available
    if ~isempty(spm_data)
        save_matrix_for_R(mean_spm, 'mean_spm_matrix');
        save_matrix_for_R(max_spm, 'max_spm_matrix');
        save_matrix_for_R(min_spm, 'min_spm_matrix');
        save_matrix_for_R(std_spm, 'std_spm_matrix');
        fprintf('SPM matrices saved (4 files)\n');
    end
    
    %% 4. Export Time Series Data (Sample - first grid point)
    export_summary = {
        sprintf('Ocean Color Data Export Summary - %s', datestr(now))
        '=================================================='
        ''
        sprintf('Variables successfully downloaded: %s', strjoin(variables_downloaded, ', '))
        'Files Created:'
        '1. ocean_color_summary_for_R.csv - All summary statistics in long format'
        '2. ocean_color_longitude.csv, ocean_color_latitude.csv - Coordinate vectors'
        '3. ocean_color_grid_info.csv - Grid dimension and extent information'
        ''
        sprintf('Data Characteristics:')
        sprintf('- Grid size: %d x %d', length(longitude), length(latitude))
        sprintf('- Longitude range: %.3f to %.3f°W', min(longitude), max(longitude))
        sprintf('- Latitude range: %.3f to %.3f°N', min(latitude), max(latitude))
        sprintf('- Time period: %s', date_range_str)
        sprintf('- Total time steps: %d', length(times_data))
        sprintf('- Variables: Kd490 (m-1), Chlorophyll-a (mg m-3), SPM (g m-3)')
        ''
        'Ready for R analysis using exported CSV files'
    };
    
    % Write summary
    fileID = fopen(fullfile(outputPath, 'ocean_color_export_summary.txt'), 'w');
    for j = 1:length(export_summary)
        fprintf(fileID, '%s\n', export_summary{j});
    end
    fclose(fileID);
    
    fprintf('Export summary saved\n');
    
else
    fprintf('No successful downloads - skipping R export\n');
end

fprintf('\n=== Script Complete ===\n');
fprintf('Total processing time: %.2f minutes\n', total_elapsed/60);
if ~isempty(variables_downloaded)
    fprintf('Successfully processed: %s\n', strjoin(variables_downloaded, ', '));
    fprintf('All files saved to: %s\n', outputPath);
else
    fprintf('No variables were successfully downloaded\n');
end