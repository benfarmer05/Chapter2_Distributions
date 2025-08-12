%% Script to pull CARICOOS SWAN wave data
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
% Initialize data arrays for all variables
hsig_data = [];      % Significant wave height
hswell_data = [];    % Swell wave height  
dir_data = [];       % Wave direction
per_data = [];       % Wave period
times_data = [];

% Calculate total number of months
total_months = length(years) * 12;
current_month = 0;

% Initialize timing variables
start_time = tic;
month_times = [];

% Create progress bar
h = waitbar(0, 'Starting SWAN data download...', 'Name', 'SWAN Data Download Progress');

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
            
            % Construct monthly URL with no overlap - FIXED VERSION
            start_date = sprintf('%d-%02d-01T00:00:00Z', year, month);
            if month == 12
                % December: go to Dec 31st, not Jan 1st of next year
                end_date = sprintf('%d-12-31T23:59:59Z', year);
            else
                % Other months: go to last day of current month
                days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
                % Handle leap year for February
                if month == 2 && mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)
                    days_in_month(2) = 29;
                end
                last_day = days_in_month(month);
                end_date = sprintf('%d-%02d-%02dT23:59:59Z', year, month, last_day);
            end
            
            % Request all four wave variables in the URL
            % NOTE - was 'http://dm3.caricoos.org:8003/erddap' before and
            % ran much faster, but that port doesn't work now. this version
            % without a specified port does run, but VERY slowly. takes ~30
            % minutes per year of data
            url = sprintf('http://dm3.caricoos.org:8003/erddap/griddap/caricoos_dm3_d848_78a5_baa6.nc?Hsig[(%s):(%s)][(17.0):(19.5)][(-68.0):(-64.0)],Hswell[(%s):(%s)][(17.0):(19.5)][(-68.0):(-64.0)],Dir[(%s):(%s)][(17.0):(19.5)][(-68.0):(-64.0)],Per[(%s):(%s)][(17.0):(19.5)][(-68.0):(-64.0)]', ...
                          start_date, end_date, start_date, end_date, start_date, end_date, start_date, end_date);
            
            % Download with error handling and retry logic
            filename = sprintf('swan_%d_%02d.nc', year, month);
            download_success = false;
            retry_count = 0;
            max_retries = 10; %was 3; occasionally not enough
            
            while ~download_success && retry_count < max_retries
                try
                    websave(filename, url);
                    fprintf('Downloaded: %s\n', filename);
                    download_success = true;
                catch ME
                    retry_count = retry_count + 1;
                    if retry_count < max_retries
                        fprintf('Download failed, retrying (%d/%d): %s\n', retry_count, max_retries, filename);
                        pause(5); % Wait before retry
                    else
                        warning('Failed to download %s after %d attempts: %s', filename, max_retries, ME.message);
                    end
                end
            end
            
            if ~download_success
                continue; % Skip to next month
            end
            
            % Read monthly data for all variables
            try
                monthly_hsig = ncread(filename, 'Hsig');
                monthly_hswell = ncread(filename, 'Hswell');
                monthly_dir = ncread(filename, 'Dir');
                monthly_per = ncread(filename, 'Per');
                monthly_time_data = ncread(filename, 'time');
                
                % Concatenate data along time dimension with proper handling
                if current_month == 1
                    % First month - initialize the arrays
                    hsig_data = monthly_hsig;
                    hswell_data = monthly_hswell;
                    dir_data = monthly_dir;
                    per_data = monthly_per;
                    times_data = monthly_time_data(:);  % Ensure column vector
                else
                    % Subsequent months - concatenate along 3rd dimension (time)
                    hsig_data = cat(3, hsig_data, monthly_hsig);
                    hswell_data = cat(3, hswell_data, monthly_hswell);
                    dir_data = cat(3, dir_data, monthly_dir);
                    per_data = cat(3, per_data, monthly_per);
                    times_data = [times_data; monthly_time_data(:)];  % Ensure column vector
                end
                
                fprintf('Month %d-%02d: Added %d time steps (Total so far: %d)\n', ...
                    year, month, size(monthly_hsig, 3), size(hsig_data, 3));
                
                % Read coordinates on first iteration
                if current_month == 1
                    latitude = ncread(filename, 'latitude');
                    longitude = ncread(filename, 'longitude');
                    fprintf('Grid dimensions: %d x %d\n', length(longitude), length(latitude));
                    fprintf('Variables: Hsig, Hswell, Dir, Per\n');
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
fprintf('Data shape: %s\n', mat2str(size(hsig_data)));
fprintf('Time range: %d data points\n', length(times_data));

% Calculate wave statistics for all variables
fprintf('Calculating wave statistics...\n');

% Significant wave height (Hsig)
mean_hsig = mean(hsig_data, 3, 'omitnan');
max_hsig = max(hsig_data, [], 3, 'omitnan');
min_hsig = min(hsig_data, [], 3, 'omitnan');
std_hsig = std(hsig_data, 0, 3, 'omitnan');

% Swell wave height (Hswell)  
mean_hswell = mean(hswell_data, 3, 'omitnan');
max_hswell = max(hswell_data, [], 3, 'omitnan');
min_hswell = min(hswell_data, [], 3, 'omitnan');
std_hswell = std(hswell_data, 0, 3, 'omitnan');

% Wave direction (Dir) - use circular statistics for angles
mean_dir = atan2d(mean(sind(dir_data), 3, 'omitnan'), mean(cosd(dir_data), 3, 'omitnan'));
mean_dir(mean_dir < 0) = mean_dir(mean_dir < 0) + 360; % Convert to 0-360 range

% Wave period (Per)
mean_per = mean(per_data, 3, 'omitnan');
max_per = max(per_data, [], 3, 'omitnan');
min_per = min(per_data, [], 3, 'omitnan');
std_per = std(per_data, 0, 3, 'omitnan');

% Overall statistics
overall_mean_hsig = mean(hsig_data(:), 'omitnan');
overall_max_hsig = max(hsig_data(:), [], 'omitnan');
overall_min_hsig = min(hsig_data(:), [], 'omitnan');
overall_mean_hswell = mean(hswell_data(:), 'omitnan');
overall_max_hswell = max(hswell_data(:), [], 'omitnan');
overall_min_hswell = min(hswell_data(:), [], 'omitnan');
overall_mean_per = mean(per_data(:), 'omitnan');
overall_max_per = max(per_data(:), [], 'omitnan');
overall_min_per = min(per_data(:), [], 'omitnan');

fprintf('\n=== Wave Statistics Summary (2017-2018) ===\n');
fprintf('Significant Wave Height (Hsig):\n');
fprintf('  Mean: %.2f m, Max: %.2f m, Min: %.2f m\n', overall_mean_hsig, overall_max_hsig, overall_min_hsig);
fprintf('Swell Wave Height (Hswell):\n');
fprintf('  Mean: %.2f m, Max: %.2f m, Min: %.2f m\n', overall_mean_hswell, overall_max_hswell, overall_min_hswell);
fprintf('Wave Period (Per):\n');
fprintf('  Mean: %.2f s, Max: %.2f s, Min: %.2f s\n', overall_mean_per, overall_max_per, overall_min_per);



%% Add this section after the "Data Processing and Statistics" section
%% and before the final fprintf statements

%% Create Maps of Wave Statistics
fprintf('\n=== Creating Maps ===\n');

%% Significant Wave Height (Hsig) Maps
% Mean Hsig
figure('Position', [100, 100, 600, 500]);
imagesc(longitude, latitude, mean_hsig');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Mean Significant Wave Height (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(mean_hsig(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'mean_hsig_map.png'));

% Max Hsig
figure('Position', [200, 200, 600, 500]);
imagesc(longitude, latitude, max_hsig');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Maximum Significant Wave Height (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(max_hsig(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'max_hsig_map.png'));

% Min Hsig
figure('Position', [300, 300, 600, 500]);
imagesc(longitude, latitude, min_hsig');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Minimum Significant Wave Height (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(min_hsig(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'min_hsig_map.png'));

% Std Hsig
figure('Position', [400, 400, 600, 500]);
imagesc(longitude, latitude, std_hsig');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Standard Deviation Significant Wave Height (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(std_hsig(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'std_hsig_map.png'));

%% Swell Wave Height (Hswell) Maps
% Mean Hswell
figure('Position', [500, 500, 600, 500]);
imagesc(longitude, latitude, mean_hswell');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Mean Swell Wave Height (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(mean_hswell(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'mean_hswell_map.png'));

% Max Hswell
figure('Position', [600, 600, 600, 500]);
imagesc(longitude, latitude, max_hswell');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Maximum Swell Wave Height (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(max_hswell(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'max_hswell_map.png'));

% Min Hswell
figure('Position', [700, 700, 600, 500]);
imagesc(longitude, latitude, min_hswell');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Minimum Swell Wave Height (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(min_hswell(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'min_hswell_map.png'));

% Std Hswell
figure('Position', [800, 800, 600, 500]);
imagesc(longitude, latitude, std_hswell');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Standard Deviation Swell Wave Height (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(std_hswell(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'std_hswell_map.png'));

%% Wave Direction (Dir) Maps
% Mean Dir (use circular colormap for direction)
figure('Position', [100, 200, 600, 500]);
imagesc(longitude, latitude, mean_dir');
set(gca, 'YDir', 'normal');
colorbar;
colormap('hsv'); % Better for circular/directional data
title('Mean Wave Direction (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 360]);
grid on;
saveas(gcf, fullfile(outputPath, 'mean_dir_map.png'));

%% Wave Period (Per) Maps
% Mean Per
figure('Position', [200, 300, 600, 500]);
imagesc(longitude, latitude, mean_per');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Mean Wave Period (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(mean_per(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'mean_per_map.png'));

% Max Per
figure('Position', [300, 400, 600, 500]);
imagesc(longitude, latitude, max_per');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Maximum Wave Period (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(max_per(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'max_per_map.png'));

% Min Per
figure('Position', [400, 500, 600, 500]);
imagesc(longitude, latitude, min_per');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Minimum Wave Period (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(min_per(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'min_per_map.png'));

% Std Per
figure('Position', [500, 600, 600, 500]);
imagesc(longitude, latitude, std_per');
set(gca, 'YDir', 'normal');
colorbar;
colormap('jet');
title('Standard Deviation Wave Period (2017-2018)', 'FontSize', 14);
xlabel('Longitude (°W)', 'FontSize', 12);
ylabel('Latitude (°N)', 'FontSize', 12);
clim([0 max(std_per(:))]);
grid on;
saveas(gcf, fullfile(outputPath, 'std_per_map.png'));

fprintf('All wave summary maps saved to output directory\n');




%% Export SWAN Wave Data for R Analysis
%   Exports processed wave data in R-friendly formats
%   Run this after the main SWAN data processing script

fprintf('\n=== Exporting Data for R Analysis ===\n');

%% 1. Export Summary Statistics as CSV (Long Format for Easy R Plotting)
fprintf('Creating summary statistics CSV for R...\n');

% Create coordinate grids
[lon_grid, lat_grid] = meshgrid(longitude, latitude);

% Flatten grids for long format
lon_vec = lon_grid(:);
lat_vec = lat_grid(:);

% Flatten summary statistics (transpose first to match coordinate orientation)
mean_hsig_t = mean_hsig';
mean_hsig_vec = mean_hsig_t(:);
max_hsig_t = max_hsig';
max_hsig_vec = max_hsig_t(:);
min_hsig_t = min_hsig';
min_hsig_vec = min_hsig_t(:);
std_hsig_t = std_hsig';
std_hsig_vec = std_hsig_t(:);

mean_hswell_t = mean_hswell';
mean_hswell_vec = mean_hswell_t(:);
max_hswell_t = max_hswell';
max_hswell_vec = max_hswell_t(:);
min_hswell_t = min_hswell';
min_hswell_vec = min_hswell_t(:);
std_hswell_t = std_hswell';
std_hswell_vec = std_hswell_t(:);

mean_dir_t = mean_dir';
mean_dir_vec = mean_dir_t(:);
mean_per_t = mean_per';
mean_per_vec = mean_per_t(:);
max_per_t = max_per';
max_per_vec = max_per_t(:);
min_per_t = min_per';
min_per_vec = min_per_t(:);
std_per_t = std_per';
std_per_vec = std_per_t(:);

% Create summary table
summary_table = table(lon_vec, lat_vec, ...
    mean_hsig_vec, max_hsig_vec, min_hsig_vec, std_hsig_vec, ...
    mean_hswell_vec, max_hswell_vec, min_hswell_vec, std_hswell_vec, ...
    mean_dir_vec, mean_per_vec, max_per_vec, min_per_vec, std_per_vec, ...
    'VariableNames', {'longitude', 'latitude', ...
    'mean_hsig', 'max_hsig', 'min_hsig', 'std_hsig', ...
    'mean_hswell', 'max_hswell', 'min_hswell', 'std_hswell', ...
    'mean_dir', 'mean_per', 'max_per', 'min_per', 'std_per'});

% Export summary CSV
writetable(summary_table, fullfile(outputPath, 'swan_wave_summary_for_R.csv'));
fprintf('Summary statistics saved: swan_wave_summary_for_R.csv\n');

%% 2. Export Coordinate Information
fprintf('Creating coordinate reference files...\n');

% Export coordinate vectors separately
lon_table = table(longitude(:), 'VariableNames', {'longitude'});
lat_table = table(latitude(:), 'VariableNames', {'latitude'});
writetable(lon_table, fullfile(outputPath, 'swan_longitude.csv'));
writetable(lat_table, fullfile(outputPath, 'swan_latitude.csv'));

% Also save grid dimensions for R raster creation
grid_info = table({'longitude'; 'latitude'; 'nlon'; 'nlat'}, ...
    {min(longitude); min(latitude); length(longitude); length(latitude)}, ...
    {max(longitude); max(latitude); length(longitude); length(latitude)}, ...
    'VariableNames', {'dimension', 'min_value', 'max_value'});
writetable(grid_info, fullfile(outputPath, 'swan_grid_info.csv'));

fprintf('Coordinate files saved: swan_longitude.csv, swan_latitude.csv, swan_grid_info.csv\n');

%% 3. Export Individual Summary Matrices (for direct raster creation in R)
fprintf('Exporting individual summary matrices...\n');

% Create a function to save matrices in R-readable format
save_matrix_for_R = @(matrix, filename) ...
    writematrix(matrix', fullfile(outputPath, [filename '.csv']));

% Export all summary statistics as individual matrices
save_matrix_for_R(mean_hsig, 'mean_hsig_matrix');
save_matrix_for_R(max_hsig, 'max_hsig_matrix');
save_matrix_for_R(min_hsig, 'min_hsig_matrix');
save_matrix_for_R(std_hsig, 'std_hsig_matrix');

save_matrix_for_R(mean_hswell, 'mean_hswell_matrix');
save_matrix_for_R(max_hswell, 'max_hswell_matrix');
save_matrix_for_R(min_hswell, 'min_hswell_matrix');
save_matrix_for_R(std_hswell, 'std_hswell_matrix');

save_matrix_for_R(mean_dir, 'mean_dir_matrix');
save_matrix_for_R(mean_per, 'mean_per_matrix');
save_matrix_for_R(max_per, 'max_per_matrix');
save_matrix_for_R(min_per, 'min_per_matrix');
save_matrix_for_R(std_per, 'std_per_matrix');

fprintf('Individual matrices saved (13 files)\n');

%% 4. Export Time Series Data (Sample - first grid point)
fprintf('Creating time series sample data...\n');

% Convert MATLAB time to datetime
time_datetime = datetime(times_data, 'ConvertFrom', 'datenum');

% Extract time series for first valid grid point (avoid NaN locations)
[valid_i, valid_j] = find(~isnan(mean_hsig), 1, 'first');

% Extract time series for this point
hsig_ts = squeeze(hsig_data(valid_i, valid_j, :));
hswell_ts = squeeze(hswell_data(valid_i, valid_j, :));
dir_ts = squeeze(dir_data(valid_i, valid_j, :));
per_ts = squeeze(per_data(valid_i, valid_j, :));

% Create time series table
ts_table = table(time_datetime, hsig_ts, hswell_ts, dir_ts, per_ts, ...
    'VariableNames', {'datetime', 'hsig', 'hswell', 'direction', 'period'});

% Add location info as metadata in first row comment
sample_lon = longitude(valid_i);
sample_lat = latitude(valid_j);

writetable(ts_table, fullfile(outputPath, 'swan_timeseries_sample.csv'));
fprintf('Time series sample saved: swan_timeseries_sample.csv\n');
fprintf('Sample location: %.3f°W, %.3f°N\n', sample_lon, sample_lat);

%% 5. Export Summary Report
fprintf('Creating export summary...\n');

export_summary = {
    sprintf('SWAN Wave Data Export Summary - %s', datestr(now))
    '=================================================='
    ''
    'Files Created:'
    '1. swan_wave_summary_for_R.csv - All summary statistics in long format'
    '2. swan_longitude.csv, swan_latitude.csv - Coordinate vectors'
    '3. swan_grid_info.csv - Grid dimension and extent information'
    '4. Individual matrix files (13 total):'
    '   - mean_hsig_matrix.csv, max_hsig_matrix.csv, etc.'
    '5. swan_timeseries_sample.csv - Time series data for one location'
    ''
    sprintf('Data Characteristics:')
    sprintf('- Grid size: %d x %d', length(longitude), length(latitude))
    sprintf('- Longitude range: %.3f to %.3f°W', min(longitude), max(longitude))
    sprintf('- Latitude range: %.3f to %.3f°N', min(latitude), max(latitude))
    sprintf('- Time period: 2017-2018')
    sprintf('- Total time steps: %d', length(times_data))
    sprintf('- Variables: Hsig, Hswell, Direction, Period')
    ''
    'Ready for R analysis using exported CSV files'
};

% Write summary
fileID = fopen(fullfile(outputPath, 'export_summary.txt'), 'w');
for i = 1:length(export_summary)
    fprintf(fileID, '%s\n', export_summary{i});
end
fclose(fileID);

fprintf('\n=== Export Complete ===\n');
fprintf('Total files created: %d\n', 17); % Updated count
fprintf('All files saved to: %s\n', outputPath);