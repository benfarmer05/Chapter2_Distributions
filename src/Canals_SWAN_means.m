%% Wave Statistics Analysis and Plotting Script
% Script to analyze and plot yearly wave statistics data
% Adapted to use project structure and R export style

clear; clc; close all;

%% Get the project root directory and define paths
projectPath = matlab.project.rootProject().RootFolder;

% Define paths relative to the project root
dataPath = fullfile(projectPath, 'data', 'Canals_SWAN');
outputPath = fullfile(projectPath, 'output');
tempPath = fullfile(projectPath, 'temp');

% Create output directories if they don't exist
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end
if ~exist(tempPath, 'dir')
    mkdir(tempPath);
end

fprintf('=== PROJECT CONFIGURATION ===\n');
fprintf('Project root: %s\n', projectPath);
fprintf('Data directory: %s\n', dataPath);
fprintf('Output directory: %s\n', outputPath);
fprintf('==============================\n\n');

%% USER SETTINGS - Change these as needed
PLOT_YEAR = 2023;  % Change this to plot different years in Figures 1 & 2
                   % Available years: 2008-2023

% Multi-year analysis settings
START_YEAR = 2013;  % Start year for multi-year analysis
END_YEAR = 2018;    % End year for multi-year analysis

%% Load coordinates
fprintf('Loading coordinates...\n');
coords_file = fullfile(dataPath, 'swanv10_coordinates.mat');
if ~exist(coords_file, 'file')
    error('Coordinates file not found: %s', coords_file);
end

coords_data = load(coords_file);
lon = coords_data.swancoords.lon;
lat = coords_data.swancoords.lat;
tri = coords_data.swancoords.tri;

fprintf('Grid points: %d\n', length(lon));
fprintf('Coordinate range: Lon [%.2f, %.2f], Lat [%.2f, %.2f]\n', ...
    min(lon), max(lon), min(lat), max(lat));

%% Load yearly data
years = 2008:2023;
fprintf('Loading yearly data files...\n');

% Pre-allocate storage
hmean_data = zeros(length(lon), length(years));
hmax_data = zeros(length(lon), length(years));
tp_at_hmax_data = zeros(length(lon), length(years));
dir_at_hmax_data = zeros(length(lon), length(years));

valid_years = [];
for i = 1:length(years)
    filename = fullfile(dataPath, sprintf('yearlystats_%d.mat', years(i)));
    if exist(filename, 'file')
        fprintf('  Loading %d...\n', years(i));
        data = load(filename);
        
        % Extract the yearly statistics
        stats = data.yearlystats;
        hmean_data(:, i) = stats.hmean;
        hmax_data(:, i) = stats.hmax;
        tp_at_hmax_data(:, i) = stats.tp_at_hmax;
        dir_at_hmax_data(:, i) = stats.dir_at_hmax;
        
        valid_years = [valid_years, years(i)];
    else
        fprintf('  Warning: %s not found\n', filename);
        % Fill with NaN for missing years
        hmean_data(:, i) = NaN;
        hmax_data(:, i) = NaN;
        tp_at_hmax_data(:, i) = NaN;
        dir_at_hmax_data(:, i) = NaN;
    end
end

fprintf('Successfully loaded %d years of data\n', length(valid_years));

%% Find index for selected plot year
plot_idx = find(years == PLOT_YEAR);

%% Plot 1: Spatial distribution of mean wave height (selected year)
figure(1);
if ~isempty(plot_idx) && ~all(isnan(hmean_data(:, plot_idx)))
    % Create scatter plot
    scatter(lon, lat, 2, hmean_data(:, plot_idx), 'filled');
    colorbar;
    colormap(jet);
    
    xlabel('Longitude');
    ylabel('Latitude');
    title(sprintf('Mean Wave Height (m) - %d', PLOT_YEAR));
    
    % Set reasonable color limits
    valid_data = hmean_data(~isnan(hmean_data(:, plot_idx)), plot_idx);
    if ~isempty(valid_data)
        low_pct = prctile(valid_data, 5);
        high_pct = prctile(valid_data, 95);
        clim([low_pct, high_pct]);
    end
    
    axis equal;
    grid on;
else
    fprintf('Warning: No data available for year %d\n', PLOT_YEAR);
    text(0.5, 0.5, sprintf('No data for %d', PLOT_YEAR), 'HorizontalAlignment', 'center');
end

%% Plot 2: Spatial distribution of maximum wave height (selected year)
figure(2);
if ~isempty(plot_idx) && ~all(isnan(hmax_data(:, plot_idx)))
    scatter(lon, lat, 2, hmax_data(:, plot_idx), 'filled');
    colorbar;
    colormap(jet);
    
    xlabel('Longitude');
    ylabel('Latitude');
    title(sprintf('Maximum Wave Height (m) - %d', PLOT_YEAR));
    
    % Set reasonable color limits
    valid_data = hmax_data(~isnan(hmax_data(:, plot_idx)), plot_idx);
    if ~isempty(valid_data)
        low_pct = prctile(valid_data, 5);
        high_pct = prctile(valid_data, 95);
        clim([low_pct, high_pct]);
    end
    
    axis equal;
    grid on;
else
    fprintf('Warning: No data available for year %d\n', PLOT_YEAR);
    text(0.5, 0.5, sprintf('No data for %d', PLOT_YEAR), 'HorizontalAlignment', 'center');
end

%% Plot 3: Time series of spatially averaged wave statistics
figure(3);
subplot(2,2,1);
% Calculate spatial means for each year
hmean_yearly = mean(hmean_data, 1, 'omitnan');
plot(years, hmean_yearly, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Year');
ylabel('Mean Wave Height (m)');
title('Spatial Average: Mean Wave Height');
grid on;

% Add trend line
valid_idx = ~isnan(hmean_yearly);
if sum(valid_idx) > 1
    p = polyfit(years(valid_idx), hmean_yearly(valid_idx), 1);
    trend = polyval(p, years);
    hold on;
    plot(years, trend, '--r', 'LineWidth', 1.5);
    if p(1) > 0
        trend_text = sprintf('Trend: +%.3f m/year', p(1));
    else
        trend_text = sprintf('Trend: %.3f m/year', p(1));
    end
    legend('Data', trend_text, 'Location', 'best');
    hold off;
end

subplot(2,2,2);
hmax_yearly = mean(hmax_data, 1, 'omitnan');
plot(years, hmax_yearly, 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [1 0.5 0]);
xlabel('Year');
ylabel('Max Wave Height (m)');
title('Spatial Average: Maximum Wave Height');
grid on;

% Add trend line
valid_idx = ~isnan(hmax_yearly);
if sum(valid_idx) > 1
    p = polyfit(years(valid_idx), hmax_yearly(valid_idx), 1);
    trend = polyval(p, years);
    hold on;
    plot(years, trend, '--r', 'LineWidth', 1.5);
    if p(1) > 0
        trend_text = sprintf('Trend: +%.3f m/year', p(1));
    else
        trend_text = sprintf('Trend: %.3f m/year', p(1));
    end
    legend('Data', trend_text, 'Location', 'best');
    hold off;
end

subplot(2,2,3);
tp_yearly = mean(tp_at_hmax_data, 1, 'omitnan');
plot(years, tp_yearly, 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0 0.7 0]);
xlabel('Year');
ylabel('Peak Period (s)');
title('Spatial Average: Peak Period at Hmax');
grid on;

subplot(2,2,4);
% For direction, use circular mean
dir_yearly = zeros(size(years));
for i = 1:length(years)
    valid_dirs = dir_at_hmax_data(~isnan(dir_at_hmax_data(:,i)), i);
    if ~isempty(valid_dirs)
        % Convert to radians, compute circular mean, convert back
        dir_rad = deg2rad(valid_dirs);
        mean_dir = atan2(mean(sin(dir_rad)), mean(cos(dir_rad)));
        dir_yearly(i) = rad2deg(mean_dir);
        if dir_yearly(i) < 0
            dir_yearly(i) = dir_yearly(i) + 360;
        end
    else
        dir_yearly(i) = NaN;
    end
end

plot(years, dir_yearly, 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.7 0 0.7]);
xlabel('Year');
ylabel('Direction (degrees)');
title('Spatial Average: Direction at Hmax');
ylim([0 360]);
grid on;

%% Plot 4: Multi-year average analysis
figure(4);

% Find indices for the selected year range
start_idx = find(years == START_YEAR);
end_idx = find(years == END_YEAR);

if ~isempty(start_idx) && ~isempty(end_idx) && start_idx <= end_idx
    % Extract data for the selected years
    selected_years = years(start_idx:end_idx);
    hmean_selected = hmean_data(:, start_idx:end_idx);
    hmax_selected = hmax_data(:, start_idx:end_idx);
    tp_selected = tp_at_hmax_data(:, start_idx:end_idx);
    dir_selected = dir_at_hmax_data(:, start_idx:end_idx);
    
    % Calculate statistics across the selected years
    % Mean of means, max of maxes, etc.
    hmean_multiyear = mean(hmean_selected, 2, 'omitnan');  % Mean across years for each location
    hmax_multiyear = max(hmax_selected, [], 2, 'omitnan');  % Max across years for each location
    tp_multiyear = mean(tp_selected, 2, 'omitnan');        % Mean across years for each location
    
    % For direction, calculate circular mean across years for each location
    dir_multiyear = zeros(size(dir_selected, 1), 1);
    for loc = 1:size(dir_selected, 1)
        valid_dirs = dir_selected(loc, ~isnan(dir_selected(loc, :)));
        if ~isempty(valid_dirs)
            dir_rad = deg2rad(valid_dirs);
            mean_dir = atan2(mean(sin(dir_rad)), mean(cos(dir_rad)));
            dir_multiyear(loc) = rad2deg(mean_dir);
            if dir_multiyear(loc) < 0
                dir_multiyear(loc) = dir_multiyear(loc) + 360;
            end
        else
            dir_multiyear(loc) = NaN;
        end
    end
    
    % Create subplot for multi-year analysis
    subplot(2,2,1);
    scatter(lon, lat, 2, hmean_multiyear, 'filled');
    colorbar;
    colormap(jet);
    xlabel('Longitude');
    ylabel('Latitude');
    title(sprintf('Mean of Mean Wave Heights (%d-%d)', START_YEAR, END_YEAR));
    
    % Set reasonable color limits
    valid_data = hmean_multiyear(~isnan(hmean_multiyear));
    if ~isempty(valid_data)
        low_pct = prctile(valid_data, 5);
        high_pct = prctile(valid_data, 95);
        clim([low_pct, high_pct]);
    end
    axis equal;
    grid on;
    
    subplot(2,2,2);
    scatter(lon, lat, 2, hmax_multiyear, 'filled');
    colorbar;
    colormap(jet);
    xlabel('Longitude');
    ylabel('Latitude');
    title(sprintf('Max of Maximum Wave Heights (%d-%d)', START_YEAR, END_YEAR));
    
    % Set reasonable color limits
    valid_data = hmax_multiyear(~isnan(hmax_multiyear));
    if ~isempty(valid_data)
        low_pct = prctile(valid_data, 5);
        high_pct = prctile(valid_data, 95);
        clim([low_pct, high_pct]);
    end
    axis equal;
    grid on;
    
    subplot(2,2,3);
    scatter(lon, lat, 2, tp_multiyear, 'filled');
    colorbar;
    colormap(jet);
    xlabel('Longitude');
    ylabel('Latitude');
    title(sprintf('Mean Peak Period at Hmax (%d-%d)', START_YEAR, END_YEAR));
    
    % Set reasonable color limits
    valid_data = tp_multiyear(~isnan(tp_multiyear));
    if ~isempty(valid_data)
        low_pct = prctile(valid_data, 5);
        high_pct = prctile(valid_data, 95);
        clim([low_pct, high_pct]);
    end
    axis equal;
    grid on;
    
    subplot(2,2,4);
    scatter(lon, lat, 2, dir_multiyear, 'filled');
    colorbar;
    colormap(hsv); % HSV colormap works better for circular data like direction
    xlabel('Longitude');
    ylabel('Latitude');
    title(sprintf('Mean Direction at Hmax (%d-%d)', START_YEAR, END_YEAR));
    clim([0, 360]);
    axis equal;
    grid on;
    
    % Print summary statistics for multi-year analysis
    fprintf('\n=== MULTI-YEAR ANALYSIS (%d-%d) ===\n', START_YEAR, END_YEAR);
    fprintf('Years included: %s\n', mat2str(selected_years));
    fprintf('Number of years: %d\n', length(selected_years));
    
    valid_hmean_multi = hmean_multiyear(~isnan(hmean_multiyear));
    valid_hmax_multi = hmax_multiyear(~isnan(hmax_multiyear));
    valid_tp_multi = tp_multiyear(~isnan(tp_multiyear));
    
    if ~isempty(valid_hmean_multi)
        fprintf('Mean of mean wave heights: %.2f ± %.2f m\n', ...
            mean(valid_hmean_multi), std(valid_hmean_multi));
    end
    if ~isempty(valid_hmax_multi)
        fprintf('Max of maximum wave heights: %.2f ± %.2f m\n', ...
            mean(valid_hmax_multi), std(valid_hmax_multi));
    end
    if ~isempty(valid_tp_multi)
        fprintf('Mean of peak periods: %.2f ± %.2f s\n', ...
            mean(valid_tp_multi), std(valid_tp_multi));
    end
    
else
    fprintf('Warning: Invalid year range specified (%d-%d)\n', START_YEAR, END_YEAR);
    text(0.5, 0.5, sprintf('Invalid year range: %d-%d', START_YEAR, END_YEAR), ...
        'HorizontalAlignment', 'center', 'FontSize', 14);
end

%% Calculate Comprehensive Statistics Across All Years
fprintf('\n=== Calculating Comprehensive Statistics ===\n');

% Calculate statistics similar to ERDDAP script style
stats = struct();

% Store coordinate data
stats.longitude = lon;
stats.latitude = lat;
stats.times = years; % Store years as time dimension

% Calculate comprehensive statistics across all years
fprintf('Processing wave height statistics...\n');
stats.mean_hmean = mean(hmean_data, 2, 'omitnan');
stats.max_hmean = max(hmean_data, [], 2, 'omitnan');
stats.min_hmean = min(hmean_data, [], 2, 'omitnan');
stats.std_hmean = std(hmean_data, 0, 2, 'omitnan');

stats.mean_hmax = mean(hmax_data, 2, 'omitnan');
stats.max_hmax = max(hmax_data, [], 2, 'omitnan');
stats.min_hmax = min(hmax_data, [], 2, 'omitnan');
stats.std_hmax = std(hmax_data, 0, 2, 'omitnan');

fprintf('Processing wave period statistics...\n');
stats.mean_tp = mean(tp_at_hmax_data, 2, 'omitnan');
stats.max_tp = max(tp_at_hmax_data, [], 2, 'omitnan');
stats.min_tp = min(tp_at_hmax_data, [], 2, 'omitnan');
stats.std_tp = std(tp_at_hmax_data, 0, 2, 'omitnan');

% Wave direction (use circular statistics)
fprintf('Processing wave direction statistics...\n');
stats.mean_dir = zeros(size(dir_at_hmax_data, 1), 1);
for loc = 1:size(dir_at_hmax_data, 1)
    valid_dirs = dir_at_hmax_data(loc, ~isnan(dir_at_hmax_data(loc, :)));
    if ~isempty(valid_dirs)
        dir_rad = deg2rad(valid_dirs);
        mean_dir = atan2(mean(sin(dir_rad)), mean(cos(dir_rad)));
        stats.mean_dir(loc) = rad2deg(mean_dir);
        if stats.mean_dir(loc) < 0
            stats.mean_dir(loc) = stats.mean_dir(loc) + 360;
        end
    else
        stats.mean_dir(loc) = NaN;
    end
end

% Overall statistics (scalar values across entire domain and time)
stats.overall_mean_hmean = mean(hmean_data(:), 'omitnan');
stats.overall_max_hmean = max(hmean_data(:), [], 'omitnan');
stats.overall_min_hmean = min(hmean_data(:), [], 'omitnan');

stats.overall_mean_hmax = mean(hmax_data(:), 'omitnan');
stats.overall_max_hmax = max(hmax_data(:), [], 'omitnan');
stats.overall_min_hmax = min(hmax_data(:), [], 'omitnan');

stats.overall_mean_tp = mean(tp_at_hmax_data(:), 'omitnan');
stats.overall_max_tp = max(tp_at_hmax_data(:), [], 'omitnan');
stats.overall_min_tp = min(tp_at_hmax_data(:), [], 'omitnan');

% Print comprehensive summary
fprintf('\n=== COMPREHENSIVE STATISTICS SUMMARY ===\n');
fprintf('Data period: %d - %d (%d years)\n', min(valid_years), max(valid_years), length(valid_years));
fprintf('Total grid points: %,d\n', length(lon));
fprintf('Years with data: %d/%d\n', length(valid_years), length(years));

fprintf('\nOverall Domain Statistics (all years, all locations):\n');
fprintf('  Mean Wave Height: %.2f-%.2fm (overall mean: %.2fm)\n', ...
    stats.overall_min_hmean, stats.overall_max_hmean, stats.overall_mean_hmean);
fprintf('  Max Wave Height: %.2f-%.2fm (overall mean: %.2fm)\n', ...
    stats.overall_min_hmax, stats.overall_max_hmax, stats.overall_mean_hmax);
fprintf('  Peak Period: %.2f-%.2fs (overall mean: %.2fs)\n', ...
    stats.overall_min_tp, stats.overall_max_tp, stats.overall_mean_tp);

%% Export Data for R Analysis (ERDDAP Script Style)
fprintf('\n=== Exporting Data for R Analysis ===\n');

% Your coordinates are already vectors (unstructured grid), not a regular lat/lon grid
% So we can use them directly without meshgrid
lon_vec = stats.longitude(:);  % Ensure column vector
lat_vec = stats.latitude(:);   % Ensure column vector

fprintf('Coordinate vectors: %d longitude points, %d latitude points\n', ...
    length(lon_vec), length(lat_vec));

% Your statistics are already vectors matching the coordinate points
% Just ensure they're column vectors
mean_hmean_vec = stats.mean_hmean(:);
max_hmean_vec = stats.max_hmean(:);
min_hmean_vec = stats.min_hmean(:);
std_hmean_vec = stats.std_hmean(:);

mean_hmax_vec = stats.mean_hmax(:);
max_hmax_vec = stats.max_hmax(:);
min_hmax_vec = stats.min_hmax(:);
std_hmax_vec = stats.std_hmax(:);

mean_dir_vec = stats.mean_dir(:);

mean_tp_vec = stats.mean_tp(:);
max_tp_vec = stats.max_tp(:);
min_tp_vec = stats.min_tp(:);
std_tp_vec = stats.std_tp(:);

% Create comprehensive summary table (following ERDDAP script variable naming)
fprintf('Creating comprehensive summary table...\n');
summary_table = table(lon_vec, lat_vec, ...
    mean_hmean_vec, max_hmean_vec, min_hmean_vec, std_hmean_vec, ...
    mean_hmax_vec, max_hmax_vec, min_hmax_vec, std_hmax_vec, ...
    mean_dir_vec, mean_tp_vec, max_tp_vec, min_tp_vec, std_tp_vec, ...
    'VariableNames', {'longitude', 'latitude', ...
    'mean_hmean', 'max_hmean', 'min_hmean', 'std_hmean', ...
    'mean_hmax', 'max_hmax', 'min_hmax', 'std_hmax', ...
    'mean_dir', 'mean_tp', 'max_tp', 'min_tp', 'std_tp'});

% Export comprehensive CSV
csv_filename = 'swan_canals_comprehensive_summary_for_R.csv';
csv_fullpath = fullfile(outputPath, csv_filename);
writetable(summary_table, csv_fullpath);

fprintf('  Exported comprehensive statistics: %s\n', csv_filename);
fprintf('  File location: %s\n', csv_fullpath);
fprintf('  Variables exported: %s\n', strjoin(summary_table.Properties.VariableNames, ', '));
fprintf('  Total records: %,d\n', height(summary_table));

%% Export Individual Summary Matrices (for direct raster creation in R)
fprintf('\nExporting individual summary matrices for R raster creation...\n');

% Create a function to save matrices in R-readable format (like your old script)
save_matrix_for_R = @(matrix, filename) ...
    writematrix(matrix', fullfile(outputPath, [filename '.csv']));

% Export all summary statistics as individual matrices
save_matrix_for_R(stats.mean_hmean, 'mean_hmean_matrix');
save_matrix_for_R(stats.max_hmean, 'max_hmean_matrix');
save_matrix_for_R(stats.min_hmean, 'min_hmean_matrix');
save_matrix_for_R(stats.std_hmean, 'std_hmean_matrix');

save_matrix_for_R(stats.mean_hmax, 'mean_hmax_matrix');
save_matrix_for_R(stats.max_hmax, 'max_hmax_matrix');
save_matrix_for_R(stats.min_hmax, 'min_hmax_matrix');
save_matrix_for_R(stats.std_hmax, 'std_hmax_matrix');

save_matrix_for_R(stats.mean_dir, 'mean_dir_matrix');
save_matrix_for_R(stats.mean_tp, 'mean_tp_matrix');
save_matrix_for_R(stats.max_tp, 'max_tp_matrix');
save_matrix_for_R(stats.min_tp, 'min_tp_matrix');
save_matrix_for_R(stats.std_tp, 'std_tp_matrix');

fprintf('Individual matrices saved (13 files)\n');

%% Export Coordinate Information (matching your R workflow exactly)
fprintf('Creating coordinate reference files for R...\n');

% Note: For unstructured grids like SWAN, we need to handle this differently
% Your original script used regular grids, but SWAN uses triangular meshes

% Export coordinate vectors as your R code expects
lon_table = table(stats.longitude(:), 'VariableNames', {'longitude'});
lat_table = table(stats.latitude(:), 'VariableNames', {'latitude'});
writetable(lon_table, fullfile(outputPath, 'swan_longitude.csv'));
writetable(lat_table, fullfile(outputPath, 'swan_latitude.csv'));

% Create grid info - but note this is for unstructured grid
grid_info = table({'longitude'; 'latitude'; 'npoints'; 'npoints'}, ...
    {min(stats.longitude); min(stats.latitude); length(stats.longitude); length(stats.latitude)}, ...
    {max(stats.longitude); max(stats.latitude); length(stats.longitude); length(stats.latitude)}, ...
    'VariableNames', {'dimension', 'min_value', 'max_value'});
writetable(grid_info, fullfile(outputPath, 'swan_grid_info.csv'));

fprintf('Coordinate files saved: swan_longitude.csv, swan_latitude.csv, swan_grid_info.csv\n');

% WARNING MESSAGE for R users about grid structure
fprintf('\n*** IMPORTANT NOTE FOR R USERS ***\n');
fprintf('Your SWAN data is on an UNSTRUCTURED triangular grid, not a regular lat/lon grid.\n');
fprintf('The matrix exports will NOT work directly with your R raster code.\n');
fprintf('Consider using the comprehensive CSV instead for spatial analysis in R.\n');
fprintf('For raster creation, you may need to interpolate to a regular grid first.\n');
fprintf('***********************************\n');

%% Export Multi-Year Period Analysis for R
if ~isempty(start_idx) && ~isempty(end_idx) && start_idx <= end_idx
    fprintf('\nExporting multi-year period analysis for R...\n');
    
    % Create vectors for multi-year analysis (already vectors, just ensure column format)
    multiyear_hmean_vec = hmean_multiyear(:);
    multiyear_hmax_vec = hmax_multiyear(:);
    multiyear_tp_vec = tp_multiyear(:);
    multiyear_dir_vec = dir_multiyear(:);
    
    % Create multi-year table
    multiyear_table = table(lon_vec, lat_vec, ...
        multiyear_hmean_vec, multiyear_hmax_vec, multiyear_tp_vec, multiyear_dir_vec, ...
        'VariableNames', {'longitude', 'latitude', ...
        'mean_of_hmean', 'max_of_hmax', 'mean_tp', 'mean_dir'});
    
    % Export multi-year CSV
    multiyear_csv_filename = sprintf('swan_canals_multiyear_%d_%d_for_R.csv', START_YEAR, END_YEAR);
    multiyear_csv_fullpath = fullfile(outputPath, multiyear_csv_filename);
    writetable(multiyear_table, multiyear_csv_fullpath);
    
    fprintf('  Exported multi-year analysis: %s\n', multiyear_csv_filename);
    fprintf('  Period: %d-%d (%d years)\n', START_YEAR, END_YEAR, length(selected_years));
    fprintf('  Variables: %s\n', strjoin(multiyear_table.Properties.VariableNames, ', '));
end

%% Create Maps (ERDDAP Script Style)
fprintf('\n=== Creating Maps ===\n');

% Create subdirectory for maps
maps_output_dir = fullfile(outputPath, 'swan_canals_maps');
if ~exist(maps_output_dir, 'dir')
    mkdir(maps_output_dir);
end

fprintf('Creating comprehensive statistics maps...\n');

% Mean wave height maps
create_wave_map(stats.longitude, stats.latitude, stats.mean_hmean, ...
    'Mean of Mean Wave Heights (All Years)', ...
    fullfile(maps_output_dir, 'mean_hmean_map.png'));

create_wave_map(stats.longitude, stats.latitude, stats.max_hmean, ...
    'Maximum of Mean Wave Heights (All Years)', ...
    fullfile(maps_output_dir, 'max_hmean_map.png'));

create_wave_map(stats.longitude, stats.latitude, stats.std_hmean, ...
    'Std Dev of Mean Wave Heights (All Years)', ...
    fullfile(maps_output_dir, 'std_hmean_map.png'));

% Maximum wave height maps
create_wave_map(stats.longitude, stats.latitude, stats.mean_hmax, ...
    'Mean of Maximum Wave Heights (All Years)', ...
    fullfile(maps_output_dir, 'mean_hmax_map.png'));

create_wave_map(stats.longitude, stats.latitude, stats.max_hmax, ...
    'Maximum of Maximum Wave Heights (All Years)', ...
    fullfile(maps_output_dir, 'max_hmax_map.png'));

create_wave_map(stats.longitude, stats.latitude, stats.std_hmax, ...
    'Std Dev of Maximum Wave Heights (All Years)', ...
    fullfile(maps_output_dir, 'std_hmax_map.png'));

% Wave direction map
create_direction_map(stats.longitude, stats.latitude, stats.mean_dir, ...
    'Mean Wave Direction (All Years)', ...
    fullfile(maps_output_dir, 'mean_dir_map.png'));

% Wave period maps
create_wave_map(stats.longitude, stats.latitude, stats.mean_tp, ...
    'Mean Wave Period (All Years)', ...
    fullfile(maps_output_dir, 'mean_tp_map.png'));

create_wave_map(stats.longitude, stats.latitude, stats.max_tp, ...
    'Maximum Wave Period (All Years)', ...
    fullfile(maps_output_dir, 'max_tp_map.png'));

fprintf('  Maps saved to: %s\n', maps_output_dir);

%% Create Summary Report (ERDDAP Script Style)
fprintf('\n=== Creating Final Summary Report ===\n');

report_lines = {
    sprintf('SWAN Canals Wave Analysis Report - %s', datestr(now))
    '=================================================================='
    ''
    sprintf('Data Source: %s', dataPath)
    sprintf('Time Period Analyzed: %d to %d (%d years)', min(valid_years), max(valid_years), length(valid_years))
    sprintf('Grid Points: %,d', length(lon))
    ''
    'Dataset Statistics:'
    sprintf('  Mean Wave Height: %.2f-%.2fm (overall mean: %.2fm)', ...
        stats.overall_min_hmean, stats.overall_max_hmean, stats.overall_mean_hmean)
    sprintf('  Max Wave Height: %.2f-%.2fm (overall mean: %.2fm)', ...
        stats.overall_min_hmax, stats.overall_max_hmax, stats.overall_mean_hmax)
    sprintf('  Peak Period: %.2f-%.2fs (overall mean: %.2fs)', ...
        stats.overall_min_tp, stats.overall_max_tp, stats.overall_mean_tp)
    ''
    'Output Files Created:'
    sprintf('  - Comprehensive statistics CSV: %s', csv_filename)
};

if ~isempty(start_idx) && ~isempty(end_idx) && start_idx <= end_idx
    report_lines{end+1} = sprintf('  - Multi-year analysis CSV (%d-%d): %s', START_YEAR, END_YEAR, multiyear_csv_filename);
end

report_lines{end+1} = '  - Statistical maps in: swan_canals_maps/';
report_lines{end+1} = '  - This summary report';
report_lines{end+1} = '';
report_lines{end+1} = 'Configuration Used:';
report_lines{end+1} = sprintf('  - Single year plots: %d', PLOT_YEAR);
report_lines{end+1} = sprintf('  - Multi-year analysis: %d-%d', START_YEAR, END_YEAR);

% Write report
report_file = fullfile(outputPath, 'swan_canals_analysis_report.txt');
fileID = fopen(report_file, 'w');
for i = 1:length(report_lines)
    fprintf(fileID, '%s\n', report_lines{i});
end
fclose(fileID);

%% Summary statistics
fprintf('\n=== FINAL SUMMARY ===\n');
fprintf('Data period: %d - %d\n', min(valid_years), max(valid_years));
fprintf('Total grid points: %,d\n', length(lon));
fprintf('Years with data: %d/%d\n', length(valid_years), length(years));

if ~isempty(valid_years)
    fprintf('\nSelected year (%d) statistics:\n', PLOT_YEAR);
    if ~isempty(plot_idx) && ~all(isnan(hmean_data(:, plot_idx)))
        valid_hmean_plot = hmean_data(~isnan(hmean_data(:,plot_idx)), plot_idx);
        valid_hmax_plot = hmax_data(~isnan(hmax_data(:,plot_idx)), plot_idx);
        valid_tp_plot = tp_at_hmax_data(~isnan(tp_at_hmax_data(:,plot_idx)), plot_idx);
        
        fprintf('  Mean wave height: %.2f ± %.2f m (mean ± std)\n', ...
            mean(valid_hmean_plot, 'omitnan'), std(valid_hmean_plot, 'omitnan'));
        fprintf('  Max wave height: %.2f ± %.2f m (mean ± std)\n', ...
            mean(valid_hmax_plot, 'omitnan'), std(valid_hmax_plot, 'omitnan'));
        fprintf('  Peak period: %.2f ± %.2f s (mean ± std)\n', ...
            mean(valid_tp_plot, 'omitnan'), std(valid_tp_plot, 'omitnan'));
    else
        fprintf('  No data available for year %d\n', PLOT_YEAR);
    end
end

fprintf('\nPlots created:\n');
fprintf('  Figure 1: Spatial map of mean wave height (%d)\n', PLOT_YEAR);
fprintf('  Figure 2: Spatial map of maximum wave height (%d)\n', PLOT_YEAR);
fprintf('  Figure 3: Time series of spatially averaged statistics\n');
fprintf('  Figure 4: Multi-year composite analysis (%d-%d)\n', START_YEAR, END_YEAR);

fprintf('\nTo change years:\n');
fprintf('  - Modify PLOT_YEAR for single-year spatial plots (Figures 1 & 2)\n');
fprintf('  - Modify START_YEAR and END_YEAR for multi-year analysis (Figure 4)\n');
fprintf('Script completed successfully!\n');