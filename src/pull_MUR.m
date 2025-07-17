%% Script to pull SST data
% MUR SST Data Download with Progress Bar and Time Estimation
% Downloads 2017-2018 NASA JPL MUR Sea Surface Temperature data
% 0.01° resolution (~1km) global data via NOAA CoastWatch ERDDAP
%   16 July 2025

clear;clc

% Get the project root directory
projectPath = matlab.project.rootProject().RootFolder;

% Define paths relative to the project root
dataPath = fullfile(projectPath, 'data');
srcPath = fullfile(projectPath, 'src');
outputPath = fullfile(projectPath, 'output');
tempPath = fullfile(projectPath, 'temp');

%% Setup
years = [2017, 2018];
% Initialize data arrays for SST variables
analysed_sst_data = [];    % Analysed sea surface temperature
times_data = [];

% Calculate total number of months
total_months = length(years) * 12;
current_month = 0;

% Initialize timing variables
start_time = tic;
month_times = [];

% Create progress bar
h = waitbar(0, 'Starting MUR SST data download...', 'Name', 'MUR SST Download Progress');

% Puerto Rico/USVI region coordinates (adjust as needed)
lat_min = 17.0;   % Southern boundary
lat_max = 19.5;   % Northern boundary  
lon_min = -68.0;  % Western boundary
lon_max = -64.0;  % Eastern boundary

fprintf('Downloading MUR SST data for Puerto Rico/USVI region:\n');
fprintf('Lat: %.1f to %.1f, Lon: %.1f to %.1f\n', lat_min, lat_max, lon_min, lon_max);
fprintf('Time period: %d-%d\n', years(1), years(end));

%% Main download loop
try
    for year = years
        for month = 1:12
            current_month = current_month + 1;
            month_start_time = tic;
            
            % Update progress bar with current status
            progress = current_month / total_months;
            status_msg = sprintf('Downloading %d-%02d (%d/%d months)', year, month, current_month, total_months);
            
            % Add time estimation if we have enough data
            if current_month > 1
                avg_time_per_month = mean(month_times);
                remaining_months = total_months - current_month;
                est_remaining_time = avg_time_per_month * remaining_months;
                
                status_msg = sprintf('%s\nEstimated time remaining: %.1f minutes', ...
                    status_msg, est_remaining_time / 60);
            end
            
            waitbar(progress, h, status_msg);
            
            % Construct monthly URL for MUR SST data with no overlap
            start_date = sprintf('%d-%02d-01T09:00:00Z', year, month);
            if month == 12
                % December: go to Dec 31st, not Jan 1st of next year
                end_date = sprintf('%d-12-31T09:00:00Z', year);
            else
                % Other months: go to last day of current month
                days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
                % Handle leap year for February
                if month == 2 && mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)
                    days_in_month(2) = 29;
                end
                last_day = days_in_month(month);
                end_date = sprintf('%d-%02d-%02dT09:00:00Z', year, month, last_day);
            end
            
            % ERDDAP URL for MUR SST (daily data)
            base_url = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.nc';
            url = sprintf('%s?analysed_sst[(%s):(%s)][(%.1f):(%.1f)][(%.1f):(%.1f)]', ...
                         base_url, start_date, end_date, lat_min, lat_max, lon_min, lon_max);
            
            % Download with increased timeout and better error handling
            filename = sprintf('mur_sst_%d_%02d.nc', year, month);
            download_success = false;
            retry_count = 0;
            max_retries = 10;
            
            % Set longer timeout for large data requests
            options = weboptions('Timeout', 30, 'RequestMethod', 'get');
            
            while ~download_success && retry_count < max_retries
                try
                    websave(filename, url, options);
                    fprintf('Downloaded: %s\n', filename);
                    download_success = true;
                catch ME
                    retry_count = retry_count + 1;
                    if retry_count < max_retries
                        fprintf('Download failed, retrying (%d/%d): %s\n', retry_count, max_retries, filename);
                        fprintf('Error: %s\n', ME.message);
                        pause(5); % Longer wait before retry
                    else
                        warning('Failed to download %s after %d attempts: %s', filename, max_retries, ME.message);
                    end
                end
            end
            
            if ~download_success
                continue; % Skip to next month
            end
            
            % Read monthly data with proper data handling
            try
                monthly_sst = ncread(filename, 'analysed_sst');
                monthly_time_data = ncread(filename, 'time');
                
                % Check for fill values and scale factors
                try
                    fill_value = ncreadatt(filename, 'analysed_sst', '_FillValue');
                    fprintf('Fill value: %.2f\n', fill_value);
                catch
                    fill_value = NaN;
                end
                
                try
                    scale_factor = ncreadatt(filename, 'analysed_sst', 'scale_factor');
                    fprintf('Scale factor: %.6f\n', scale_factor);
                catch
                    scale_factor = 1;
                end
                
                try
                    add_offset = ncreadatt(filename, 'analysed_sst', 'add_offset');
                    fprintf('Add offset: %.6f\n', add_offset);
                catch
                    add_offset = 0;
                end
                
                % Apply proper scaling and fill value handling
                monthly_sst(monthly_sst == fill_value) = NaN;
                monthly_sst = monthly_sst * scale_factor + add_offset;
                
                % Check if data is in Kelvin (reasonable SST should be 270-310K or 15-35°C)
                sample_val = monthly_sst(~isnan(monthly_sst));
                if ~isempty(sample_val)
                    sample_mean = mean(sample_val(1:min(100, length(sample_val))));
                    fprintf('Sample mean before conversion: %.2f\n', sample_mean);
                    
                    % If values are in reasonable Kelvin range (270-310), convert to Celsius
                    if sample_mean > 250 && sample_mean < 350
                        monthly_sst = monthly_sst - 273.15;
                        fprintf('Converted from Kelvin to Celsius\n');
                    elseif sample_mean > -50 && sample_mean < 50
                        fprintf('Data appears to already be in Celsius\n');
                    else
                        warning('Unexpected temperature range - check data quality');
                    end
                end
                
                % Concatenate data along time dimension with proper handling
                if current_month == 1
                    % First month - initialize the arrays
                    analysed_sst_data = monthly_sst;
                    times_data = monthly_time_data(:);  % Ensure column vector
                else
                    % Subsequent months - concatenate along 3rd dimension (time)
                    analysed_sst_data = cat(3, analysed_sst_data, monthly_sst);
                    times_data = [times_data; monthly_time_data(:)];  % Ensure column vector
                end
                
                fprintf('Month %d-%02d: Added %d time steps (Total so far: %d)\n', ...
                    year, month, size(monthly_sst, 3), size(analysed_sst_data, 3));
                
                % Read coordinates on first iteration with metadata check
                if current_month == 1
                    latitude = ncread(filename, 'latitude');
                    longitude = ncread(filename, 'longitude');
                    fprintf('Grid dimensions: %d x %d\n', length(longitude), length(latitude));
                    
                    % Display metadata for debugging
                    try
                        sst_units = ncreadatt(filename, 'analysed_sst', 'units');
                        fprintf('SST units: %s\n', sst_units);
                    catch
                        fprintf('SST units: not specified\n');
                    end
                    
                    try
                        sst_long_name = ncreadatt(filename, 'analysed_sst', 'long_name');
                        fprintf('SST description: %s\n', sst_long_name);
                    catch
                    end
                    
                    % Show actual data range for verification
                    valid_sst = monthly_sst(~isnan(monthly_sst));
                    if ~isempty(valid_sst)
                        fprintf('Actual SST range: %.2f to %.2f\n', min(valid_sst), max(valid_sst));
                    end
                end
                
            catch ME
                warning('Failed to read data from %s: %s', filename, ME.message);
            end
            
            % Clean up file
            if exist(filename, 'file')
                delete(filename);
            end
            
            % Record timing for this month
            month_time = toc(month_start_time);
            month_times(end+1) = month_time;
            
            % Update console with detailed progress
            elapsed_time = toc(start_time);
            if current_month > 1
                avg_time_per_month = mean(month_times);
                est_total_time = avg_time_per_month * total_months;
                est_remaining_time = est_total_time - elapsed_time;
                
                fprintf('Progress: %d/%d (%.1f%%) | Elapsed: %.1f min | Remaining: %.1f min | Current file: %.1f sec\n', ...
                    current_month, total_months, progress*100, elapsed_time/60, est_remaining_time/60, month_time);
            else
                fprintf('Progress: %d/%d (%.1f%%) | Elapsed: %.1f min | Current file: %.1f sec\n', ...
                    current_month, total_months, progress*100, elapsed_time/60, month_time);
            end
        end
    end
    
    % Final progress update
    total_elapsed = toc(start_time);
    waitbar(1, h, sprintf('Download complete! Total time: %.1f minutes', total_elapsed/60));
    pause(2); % Show completion message briefly
    close(h);
    
catch ME
    % Close progress bar on error
    if ishandle(h)
        close(h);
    end
    rethrow(ME);
end

%% Data Processing and Statistics
fprintf('\n=== Processing Results ===\n');
fprintf('Total download time: %.2f minutes\n', total_elapsed/60);
fprintf('SST data shape: %s\n', mat2str(size(analysed_sst_data)));
fprintf('Time range: %d data points\n', length(times_data));

% Calculate SST statistics
fprintf('Calculating SST statistics...\n');

% Basic statistics
mean_sst = mean(analysed_sst_data, 3, 'omitnan');
max_sst = max(analysed_sst_data, [], 3, 'omitnan');
min_sst = min(analysed_sst_data, [], 3, 'omitnan');
std_sst = std(analysed_sst_data, 0, 3, 'omitnan');

% Overall statistics
overall_mean_sst = mean(analysed_sst_data(:), 'omitnan');
overall_max_sst = max(analysed_sst_data(:), [], 'omitnan');
overall_min_sst = min(analysed_sst_data(:), [], 'omitnan');
overall_std_sst = std(analysed_sst_data(:), 'omitnan');

fprintf('\n=== Sea Surface Temperature Statistics (2017-2018) ===\n');
fprintf('Sea Surface Temperature (SST):\n');
fprintf('  Mean: %.2f °C, Max: %.2f °C, Min: %.2f °C\n', overall_mean_sst, overall_max_sst, overall_min_sst);

%% Create Maps of SST Statistics
fprintf('\n=== Creating Maps ===\n');

%% Sea Surface Temperature Maps
% Mean SST
figure('Position', [100, 100, 600, 500]);
imagesc(longitude, latitude, mean_sst');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Mean Sea Surface Temperature (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([min(mean_sst(:)) max(mean_sst(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'mean_sst_map.png'));

% Max SST
figure('Position', [200, 200, 600, 500]);
imagesc(longitude, latitude, max_sst');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Maximum Sea Surface Temperature (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([min(max_sst(:)) max(max_sst(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'max_sst_map.png'));

% Min SST
figure('Position', [300, 300, 600, 500]);
imagesc(longitude, latitude, min_sst');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Minimum Sea Surface Temperature (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([min(min_sst(:)) max(min_sst(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'min_sst_map.png'));

% Std SST
figure('Position', [400, 400, 600, 500]);
imagesc(longitude, latitude, std_sst');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Standard Deviation Sea Surface Temperature (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(std_sst(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'std_sst_map.png'));

fprintf('All SST summary maps saved to output directory\n');

%% Export MUR SST Data for R Analysis
%   Exports processed SST data in R-friendly formats
%   Following the same format as SWAN wave data export

fprintf('\n=== Exporting Data for R Analysis ===\n');

%% 1. Export Summary Statistics as CSV (Long Format for Easy R Plotting)
fprintf('Creating summary statistics CSV for R...\n');

% Create coordinate grids
[lon_grid, lat_grid] = meshgrid(longitude, latitude);

% Flatten grids for long format
lon_vec = lon_grid(:);
lat_vec = lat_grid(:);

% Flatten summary statistics (transpose first to match coordinate orientation)
mean_sst_t = mean_sst';
mean_sst_vec = mean_sst_t(:);
max_sst_t = max_sst';
max_sst_vec = max_sst_t(:);
min_sst_t = min_sst';
min_sst_vec = min_sst_t(:);
std_sst_t = std_sst';
std_sst_vec = std_sst_t(:);

% Create summary table
summary_table = table(lon_vec, lat_vec, ...
    mean_sst_vec, max_sst_vec, min_sst_vec, std_sst_vec, ...
    'VariableNames', {'longitude', 'latitude', ...
    'mean_sst', 'max_sst', 'min_sst', 'std_sst'});

% Export summary CSV
writetable(summary_table, fullfile(outputPath, 'mur_sst_summary_for_R.csv'));
fprintf('Summary statistics saved: mur_sst_summary_for_R.csv\n');

%% 2. Export Coordinate Information
fprintf('Creating coordinate reference files...\n');

% Export coordinate vectors separately
lon_table = table(longitude(:), 'VariableNames', {'longitude'});
lat_table = table(latitude(:), 'VariableNames', {'latitude'});
writetable(lon_table, fullfile(outputPath, 'mur_sst_longitude.csv'));
writetable(lat_table, fullfile(outputPath, 'mur_sst_latitude.csv'));

% Also save grid dimensions for R raster creation
grid_info = table({'longitude'; 'latitude'; 'nlon'; 'nlat'}, ...
    {min(longitude); min(latitude); length(longitude); length(latitude)}, ...
    {max(longitude); max(latitude); length(longitude); length(latitude)}, ...
    'VariableNames', {'dimension', 'min_value', 'max_value'});
writetable(grid_info, fullfile(outputPath, 'mur_sst_grid_info.csv'));

fprintf('Coordinate files saved: mur_sst_longitude.csv, mur_sst_latitude.csv, mur_sst_grid_info.csv\n');

%% 3. Export Individual Summary Matrices (for direct raster creation in R)
fprintf('Exporting individual summary matrices...\n');

% Create a function to save matrices in R-readable format
save_matrix_for_R = @(matrix, filename) ...
    writematrix(matrix', fullfile(outputPath, [filename '.csv']));

% Export all summary statistics as individual matrices
save_matrix_for_R(mean_sst, 'mean_sst_matrix');
save_matrix_for_R(max_sst, 'max_sst_matrix');
save_matrix_for_R(min_sst, 'min_sst_matrix');
save_matrix_for_R(std_sst, 'std_sst_matrix');

fprintf('Individual matrices saved (4 files)\n');

%% 4. Export Time Series Data (Sample - first grid point)
fprintf('Creating time series sample data...\n');

% Convert MATLAB time to datetime
time_datetime = datetime(times_data, 'ConvertFrom', 'datenum');

% Extract time series for first valid grid point (avoid NaN locations)
[valid_i, valid_j] = find(~isnan(mean_sst), 1, 'first');

% Extract time series for this point
sst_ts = squeeze(analysed_sst_data(valid_i, valid_j, :));

% Create time series table
ts_table = table(time_datetime, sst_ts, ...
    'VariableNames', {'datetime', 'sst_celsius'});

% Add location info as metadata
sample_lon = longitude(valid_i);
sample_lat = latitude(valid_j);

writetable(ts_table, fullfile(outputPath, 'mur_sst_timeseries_sample.csv'));
fprintf('Time series sample saved: mur_sst_timeseries_sample.csv\n');
fprintf('Sample location: %.3f°W, %.3f°N\n', sample_lon, sample_lat);

%% 5. Export Summary Report
fprintf('Creating export summary...\n');

export_summary = {
    sprintf('MUR SST Data Export Summary - %s', datestr(now))
    '=================================================='
    ''
    'Files Created:'
    '1. mur_sst_summary_for_R.csv - All summary statistics in long format'
    '2. mur_sst_longitude.csv, mur_sst_latitude.csv - Coordinate vectors'
    '3. mur_sst_grid_info.csv - Grid dimension and extent information'
    '4. Individual matrix files (4 total):'
    '   - mean_sst_matrix.csv, max_sst_matrix.csv, etc.'
    '5. mur_sst_timeseries_sample.csv - Time series data for one location'
    ''
    sprintf('Data Characteristics:')
    sprintf('- Grid size: %d x %d', length(longitude), length(latitude))
    sprintf('- Longitude range: %.3f to %.3f°W', min(longitude), max(longitude))
    sprintf('- Latitude range: %.3f to %.3f°N', min(latitude), max(latitude))
    sprintf('- Time period: 2017-2018')
    sprintf('- Total time steps: %d', length(times_data))
    sprintf('- Variable: Sea Surface Temperature (°C)')
    ''
    'Ready for R analysis using exported CSV files'
};

% Write summary
fileID = fopen(fullfile(outputPath, 'sst_export_summary.txt'), 'w');
for i = 1:length(export_summary)
    fprintf(fileID, '%s\n', export_summary{i});
end
fclose(fileID);

fprintf('\n=== Export Complete ===\n');
fprintf('Total files created: %d\n', 9); % Updated count
fprintf('All files saved to: %s\n', outputPath);