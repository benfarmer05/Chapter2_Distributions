%% Streamlined CARICOOS SWAN Wave Data Download Script
%   Updated: October 2025
%   Enhanced with dynamic progress tracking and ETA

clear;clc

%% CONFIGURATION
USE_ALREADY_DOWNLOADED = true;  % Toggle: skip download and use existing files
USE_CUSTOM_DATES = true;
CUSTOM_START_DATE = '2024-01-03';
CUSTOM_END_DATE = '2024-12-31';
CHUNK_SIZE_DAYS = 1;
SAVE_RAW_DATA = true;  % Toggle: save downloaded NetCDF files locally

% Output paths
dataPath = 'D:\SWAN_ERDDAP\data';
outputPath = 'D:\SWAN_ERDDAP\output';

% Create directories if they don't exist
if ~exist(dataPath, 'dir'), mkdir(dataPath); end
if ~exist(outputPath, 'dir'), mkdir(outputPath); end

%% DATASET DEFINITIONS
% Initialize as struct array with explicit scalar fields
datasets = struct([]);

datasets(1).name = 'USVI_HR';
datasets(1).id = 'SWAN_HighRes_USVI';
datasets(1).server = 'http://dm3.caricoos.org/erddap';
datasets(1).variables = {'dir','tp'};
datasets(1).description = 'USVI High-Res (100m)';
datasets(1).resolution = 100;
datasets(1).lat_bounds = '(18.18):(18.8)';
datasets(1).lon_bounds = '(-65.2):(-64.0)';

datasets(2).name = 'PuertoRico_HR';
datasets(2).id = 'SWAN_HighRes_PR';
datasets(2).server = 'http://dm3.caricoos.org/erddap';
datasets(2).variables = {'dir','tp'};
datasets(2).description = 'PR High-Res (120m)';
datasets(2).resolution = 120;
datasets(2).lat_bounds = '(17.8):(18.6)';
datasets(2).lon_bounds = '(-67.5):(-65.1)';

datasets(3).name = 'StCroix_HR';
datasets(3).id = 'SWAN_HighRes_StCroix';
datasets(3).server = 'http://dm3.caricoos.org/erddap';
datasets(3).variables = {'dir','tp'};
datasets(3).description = 'StCroix High-Res (180m)';
datasets(3).resolution = 180;
datasets(3).lat_bounds = '(17.63):(17.825)';
datasets(3).lon_bounds = '(-65.0):(-64.48)';

datasets(4).name = 'Regional_PRVI';
datasets(4).id = 'swan_1km_caricoos';
datasets(4).server = 'http://dm3.caricoos.org/erddap';
datasets(4).variables = {'Dir','Per'};
datasets(4).description = 'Regional (1km)';
datasets(4).resolution = 1000;
datasets(4).lat_bounds = '(17.0):(19.5)';
datasets(4).lon_bounds = '(-68.0):(-64.0)';

%% GET TIME RANGES
fprintf('=== Checking Time Ranges ===\n');
for i = 1:length(datasets)
    info_url = sprintf('%s/griddap/%s.das', datasets(i).server, datasets(i).id);
    try
        info = webread(info_url);
        datasets(i).start = regexp(info, 'time_coverage_start\s+"([^"]+)"', 'tokens', 'once');
        datasets(i).end = regexp(info, 'time_coverage_end\s+"([^"]+)"', 'tokens', 'once');
        fprintf('%s: %s to %s\n', datasets(i).name, datasets(i).start{1}, datasets(i).end{1});
    catch
        datasets(i).start = {'2024-01-02T12:00:00Z'};
        datasets(i).end = {'2025-10-10T12:00:00Z'};
    end
end

%% SET DOWNLOAD PERIOD
if USE_CUSTOM_DATES
    % Get actual data availability
    actual_starts = cellfun(@(x) datetime(x{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z'''), {datasets.start});
    actual_ends = cellfun(@(x) datetime(x{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z'''), {datasets.end});
    server_start = max(actual_starts);
    server_end = min(actual_ends);
    
    % Parse custom dates
    custom_start = datetime(sprintf('%sT12:00:00Z', CUSTOM_START_DATE), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    custom_end = datetime(sprintf('%sT12:00:00Z', CUSTOM_END_DATE), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    
    % Constrain to available range
    overlap_start = max(custom_start, server_start);
    overlap_end = min(custom_end, server_end);
    
    if custom_start < server_start || custom_end > server_end
        fprintf('WARNING: Custom dates adjusted to available range\n');
        fprintf('  Requested: %s to %s\n', CUSTOM_START_DATE, CUSTOM_END_DATE);
        fprintf('  Available: %s to %s\n', char(server_start, 'yyyy-MM-dd'), char(server_end, 'yyyy-MM-dd'));
    end
    
    overlap_start_str = char(overlap_start, 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    overlap_end_str = char(overlap_end, 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    fprintf('Using dates: %s to %s\n', char(overlap_start, 'yyyy-MM-dd'), char(overlap_end, 'yyyy-MM-dd'));
else
    actual_starts = cellfun(@(x) datetime(x{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z'''), {datasets.start});
    actual_ends = cellfun(@(x) datetime(x{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z'''), {datasets.end});
    overlap_start = max(actual_starts);
    overlap_end = min(actual_ends);
    overlap_start_str = char(overlap_start, 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    overlap_end_str = char(overlap_end, 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    fprintf('Using server overlap: %s to %s\n', overlap_start_str, overlap_end_str);
end

chunk_starts = overlap_start:days(CHUNK_SIZE_DAYS):overlap_end;

%% DOWNLOAD DATA WITH PROGRESS TRACKING
all_data = struct();

if USE_ALREADY_DOWNLOADED
    fprintf('\n=== Loading Already Downloaded Data ===\n');
    
    for d = 1:length(datasets)
        fprintf('Loading %s... ', datasets(d).name);
        
        % Find all NetCDF files for this dataset
        pattern = fullfile(dataPath, sprintf('%s_*.nc', datasets(d).name));
        files = dir(pattern);
        
        if isempty(files)
            warning('No files found for %s. Skipping.', datasets(d).name);
            continue;
        end
        
        fprintf('found %d files\n', length(files));
        
        % Pre-scan to get dimensions and pre-allocate
        fprintf('  Pre-scanning dimensions... ');
        first_file = fullfile(files(1).folder, files(1).name);
        info = ncinfo(first_file);
        
        % Get spatial dimensions
        nlat = 0; nlon = 0;
        for i = 1:length(info.Dimensions)
            if strcmp(info.Dimensions(i).Name, 'latitude')
                nlat = info.Dimensions(i).Length;
            elseif strcmp(info.Dimensions(i).Name, 'longitude')
                nlon = info.Dimensions(i).Length;
            end
        end
        
        % Count total time steps across all files
        total_timesteps = 0;
        for f = 1:length(files)
            filepath = fullfile(files(f).folder, files(f).name);
            finfo = ncinfo(filepath);
            for i = 1:length(finfo.Dimensions)
                if strcmp(finfo.Dimensions(i).Name, 'time')
                    total_timesteps = total_timesteps + finfo.Dimensions(i).Length;
                    break;
                end
            end
        end
        fprintf('%d total timesteps\n', total_timesteps);
        
        % Pre-allocate arrays
        fprintf('  Allocating memory... ');
        data = struct();
        for v = 1:length(datasets(d).variables)
            data.(datasets(d).variables{v}) = NaN(nlon, nlat, total_timesteps);
        end
        data.times = NaN(total_timesteps, 1);
        fprintf('done\n');
        
        % Load data into pre-allocated arrays
        fprintf('  Reading files: ');
        time_idx = 1;
        for f = 1:length(files)
            if mod(f, 10) == 0, fprintf('%d ', f); end
            filepath = fullfile(files(f).folder, files(f).name);
            
            try
                % Get number of timesteps in this file
                finfo = ncinfo(filepath);
                file_timesteps = 0;
                for i = 1:length(finfo.Dimensions)
                    if strcmp(finfo.Dimensions(i).Name, 'time')
                        file_timesteps = finfo.Dimensions(i).Length;
                        break;
                    end
                end
                
                % Read all variables directly into position
                for v = 1:length(datasets(d).variables)
                    chunk_data = ncread(filepath, datasets(d).variables{v});
                    data.(datasets(d).variables{v})(:, :, time_idx:time_idx+file_timesteps-1) = chunk_data;
                end
                
                % Read time
                chunk_time = ncread(filepath, 'time');
                data.times(time_idx:time_idx+file_timesteps-1) = chunk_time(:);
                
                % Read lat/lon on first file only
                if f == 1
                    data.latitude = ncread(filepath, 'latitude');
                    data.longitude = ncread(filepath, 'longitude');
                end
                
                time_idx = time_idx + file_timesteps;
                
            catch ME
                warning('Failed to read %s: %s', files(f).name, ME.message);
            end
        end
        fprintf('done\n');
        
        all_data.(datasets(d).name) = data;
        fprintf('  Complete: %d time steps loaded\n', length(data.times));
    end
    
    fprintf('=== Load Complete ===\n');
    
else
    fprintf('\n=== Downloading Data ===\n');
    
    % Initialize overall progress tracker
    total_chunks = length(chunk_starts) * length(datasets);
    overall_progress = 0;
    overall_start_time = tic;
    chunk_times = [];

    for d = 1:length(datasets)
        fprintf('\n[%d/%d] %s\n', d, length(datasets), datasets(d).description);
        data = struct('times', [], 'first', true);
        for v = 1:length(datasets(d).variables), data.(datasets(d).variables{v}) = []; end
        
        h = waitbar(0, datasets(d).name);
        dataset_start_time = tic;
        
        for c = 1:length(chunk_starts)
            chunk_start_time = tic;
            chunk_start = chunk_starts(c);
            chunk_end = min(chunk_start + days(CHUNK_SIZE_DAYS-1), overlap_end);
            
            % Calculate progress percentages
            dataset_progress = (c-1) / length(chunk_starts) * 100;
            global_progress = (overall_progress + (c-1)) / total_chunks * 100;
            
            % Build URL
            start_str = sprintf('%sT00:00:00Z', char(chunk_start, 'yyyy-MM-dd'));
            end_str = sprintf('%sT23:59:59Z', char(chunk_end, 'yyyy-MM-dd'));
            var_list = strjoin(arrayfun(@(i) sprintf('%s[(%s):(%s)][%s][%s]', ...
                datasets(d).variables{i}, start_str, end_str, ...
                datasets(d).lat_bounds, datasets(d).lon_bounds), ...
                1:length(datasets(d).variables), 'UniformOutput', false), ',');
            url = sprintf('%s/griddap/%s.nc?%s', datasets(d).server, datasets(d).id, var_list);
            
            % Download
            if SAVE_RAW_DATA
                filename = fullfile(dataPath, sprintf('%s_%s.nc', datasets(d).name, char(chunk_start, 'yyyyMMdd')));
            else
                filename = sprintf('%s_temp.nc', datasets(d).name);
            end
            
            try
                websave(filename, url);
                
                % Read data
                for v = 1:length(datasets(d).variables)
                    chunk_data = ncread(filename, datasets(d).variables{v});
                    if data.first
                        data.(datasets(d).variables{v}) = chunk_data;
                    else
                        data.(datasets(d).variables{v}) = cat(3, data.(datasets(d).variables{v}), chunk_data);
                    end
                end
                
                chunk_time = ncread(filename, 'time');
                data.times = [data.times; chunk_time(:)];
                
                if data.first
                    data.latitude = ncread(filename, 'latitude');
                    data.longitude = ncread(filename, 'longitude');
                    data.first = false;
                end
                
                if ~SAVE_RAW_DATA && exist(filename, 'file'), delete(filename); end
                
                % Track chunk processing time
                chunk_elapsed = toc(chunk_start_time);
                chunk_times = [chunk_times; chunk_elapsed];
                
                % Calculate ETAs
                if length(chunk_times) >= 3
                    avg_chunk_time = mean(chunk_times(max(1,end-9):end));
                    
                    % Dataset ETA
                    chunks_remaining_dataset = length(chunk_starts) - c;
                    dataset_eta_sec = chunks_remaining_dataset * avg_chunk_time;
                    
                    % Overall ETA
                    chunks_remaining_overall = total_chunks - (overall_progress + c);
                    overall_eta_sec = chunks_remaining_overall * avg_chunk_time;
                    
                    % Format time strings
                    if dataset_eta_sec < 60
                        dataset_eta_str = sprintf('%.0fs', dataset_eta_sec);
                    elseif dataset_eta_sec < 3600
                        dataset_eta_str = sprintf('%dm %.0fs', floor(dataset_eta_sec/60), mod(dataset_eta_sec,60));
                    else
                        dataset_eta_str = sprintf('%dh %dm', floor(dataset_eta_sec/3600), floor(mod(dataset_eta_sec,3600)/60));
                    end
                    
                    if overall_eta_sec < 60
                        overall_eta_str = sprintf('%.0fs', overall_eta_sec);
                    elseif overall_eta_sec < 3600
                        overall_eta_str = sprintf('%dm %.0fs', floor(overall_eta_sec/60), mod(overall_eta_sec,60));
                    else
                        overall_eta_str = sprintf('%dh %dm', floor(overall_eta_sec/3600), floor(mod(overall_eta_sec,3600)/60));
                    end
                    
                    % Update waitbar with detailed info
                    waitbar_msg = sprintf('%s: %d/%d (%.1f%%)\nDataset ETA: %s | Overall ETA: %s\nSpeed: %.1fs/chunk | Global: %.1f%%', ...
                        datasets(d).name, c, length(chunk_starts), dataset_progress, ...
                        dataset_eta_str, overall_eta_str, avg_chunk_time, global_progress);
                    waitbar(c/length(chunk_starts), h, waitbar_msg);
                    
                    % Console output every 10 chunks or at completion
                    if mod(c, 10) == 0 || c == length(chunk_starts)
                        fprintf('  Progress: %d/%d | Dataset ETA: %s | Overall: %.1f%% (ETA: %s)\n', ...
                            c, length(chunk_starts), dataset_eta_str, global_progress, overall_eta_str);
                    end
                else
                    % Initial progress without ETA
                    waitbar(c/length(chunk_starts), h, sprintf('%s: %d/%d (%.1f%%) - Calculating ETA...', ...
                        datasets(d).name, c, length(chunk_starts), dataset_progress));
                end
                
            catch ME
                if contains(ME.message, 'InterruptException') || contains(ME.message, 'ctrl-c')
                    close(h);
                    error('User interrupted download at dataset %d, chunk %d', d, c);
                end
                warning('Failed chunk %d: %s', c, ME.message);
            end
        end
        
        close(h);
        overall_progress = overall_progress + length(chunk_starts);
        
        % Dataset completion summary
        dataset_elapsed = toc(dataset_start_time);
        all_data.(datasets(d).name) = data;
        
        if dataset_elapsed < 60
            elapsed_str = sprintf('%.0fs', dataset_elapsed);
        elseif dataset_elapsed < 3600
            elapsed_str = sprintf('%dm %.0fs', floor(dataset_elapsed/60), mod(dataset_elapsed,60));
        else
            elapsed_str = sprintf('%dh %dm', floor(dataset_elapsed/3600), floor(mod(dataset_elapsed,3600)/60));
        end
        
        fprintf('Complete: %d time steps in %s\n', length(data.times), elapsed_str);
    end

    % Overall completion summary
    total_elapsed = toc(overall_start_time);

    if total_elapsed < 60
        total_str = sprintf('%.0fs', total_elapsed);
    elseif total_elapsed < 3600
        total_str = sprintf('%dm %.0fs', floor(total_elapsed/60), mod(total_elapsed,60));
    else
        total_str = sprintf('%dh %dm', floor(total_elapsed/3600), floor(mod(total_elapsed,3600)/60));
    end

    fprintf('\n=== Download Complete ===\n');
    fprintf('Total time: %s\n', total_str);
    fprintf('Average speed: %.2f chunks/min\n', total_chunks / (total_elapsed/60));
end

%% CALCULATE STATISTICS
fprintf('\n=== Calculating Statistics ===\n');

for d = 1:length(datasets)
    name = datasets(d).name;
    if ~isfield(all_data, name), continue; end
    
    data = all_data.(name);
    
    % Map variables
    if isfield(data, 'Dir')
        dir = data.Dir; per = data.Per;
    else
        dir = data.dir; per = data.tp;
    end
    
    stats.longitude = data.longitude;
    stats.latitude = data.latitude;
    stats.times = data.times;
    stats.mean_dir = atan2d(mean(sind(dir), 3, 'omitnan'), mean(cosd(dir), 3, 'omitnan'));
    stats.mean_dir(stats.mean_dir < 0) = stats.mean_dir(stats.mean_dir < 0) + 360;
    stats.mean_per = mean(per, 3, 'omitnan');
    stats.max_per = max(per, [], 3, 'omitnan');
    
    all_data.(name).stats = stats;
    fprintf('%s: Per=%.1f-%.1fs\n', name, min(per(:), [], 'omitnan'), max(per(:), [], 'omitnan'));
end

%% CREATE PRIORITY-BASED COMPOSITE
fprintf('\n=== Creating Priority Composite ===\n');

% Sort datasets by resolution (highest first) - extract scalar values
n_datasets = length(datasets);
resolutions = zeros(1, n_datasets);
for i = 1:n_datasets
    resolutions(i) = datasets(i).resolution;
end

[~, sort_idx] = sort(resolutions, 'ascend');  % Ascending = smallest resolution first = highest priority

% Build priority order array
priority_order = cell(1, n_datasets);
for i = 1:n_datasets
    priority_order{i} = datasets(sort_idx(i)).name;
end

fprintf('Priority order (high to low res): %s\n', strjoin(priority_order, ', '));

% Collect all valid points - FAST VERSION
composite_data = [];

% Pre-compute bounds for all datasets
bounds = struct();
for p = 1:length(priority_order)
    name = priority_order{p};
    if ~isfield(all_data, name) || ~isfield(all_data.(name), 'stats'), continue; end
    stats = all_data.(name).stats;
    bounds.(name).lat_min = min(stats.latitude);
    bounds.(name).lat_max = max(stats.latitude);
    bounds.(name).lon_min = min(stats.longitude);
    bounds.(name).lon_max = max(stats.longitude);
end

for p = 1:length(priority_order)
    name = priority_order{p};
    if ~isfield(all_data, name) || ~isfield(all_data.(name), 'stats'), continue; end
    
    stats = all_data.(name).stats;
    [lon_grid, lat_grid] = meshgrid(stats.longitude, stats.latitude);
    
    fprintf('Processing %s (priority %d)... ', name, p);
    
    % Vectorize the coverage check
    mask = true(size(lon_grid));  % Start with all points valid
    
    % Check against higher priority datasets
    for h = 1:(p-1)
        higher_name = priority_order{h};
        if ~isfield(bounds, higher_name), continue; end
        
        % Vectorized bounds check
        covered = (lat_grid >= bounds.(higher_name).lat_min) & ...
                  (lat_grid <= bounds.(higher_name).lat_max) & ...
                  (lon_grid >= bounds.(higher_name).lon_min) & ...
                  (lon_grid <= bounds.(higher_name).lon_max);
        
        mask = mask & ~covered;  % Remove covered points
    end
    
    % Get uncovered points
    uncovered_idx = find(mask);
    n_points = length(uncovered_idx);
    
    if n_points > 0
        % Pre-allocate block
        new_data = zeros(n_points, 6);
        new_data(:, 1) = lon_grid(uncovered_idx);
        new_data(:, 2) = lat_grid(uncovered_idx);
        
        % Transpose stats for linear indexing
        temp_dir = stats.mean_dir';
        temp_per = stats.mean_per';
        temp_max = stats.max_per';
        
        new_data(:, 3) = temp_dir(uncovered_idx);
        new_data(:, 4) = temp_per(uncovered_idx);
        new_data(:, 5) = temp_max(uncovered_idx);
        new_data(:, 6) = p;
        
        composite_data = [composite_data; new_data];
    end
    
    fprintf('%d points added\n', n_points);
end

fprintf('Total composite points: %d\n', size(composite_data, 1));

%% CREATE MAPS
fprintf('\n=== Creating Maps ===\n');

% Create composite maps
if ~isempty(composite_data)
    fprintf('Creating composite maps...\n');
    
    % Get unique coordinates
    unique_lons = unique(composite_data(:,1));
    unique_lats = unique(composite_data(:,2));
    
    % Create grids
    [lon_grid, lat_grid] = meshgrid(unique_lons, unique_lats);
    mean_dir_grid = NaN(size(lon_grid));
    mean_per_grid = NaN(size(lon_grid));
    max_per_grid = NaN(size(lon_grid));
    
    % Fill grids from composite data
    for i = 1:size(composite_data, 1)
        lon_idx = find(unique_lons == composite_data(i,1), 1);
        lat_idx = find(unique_lats == composite_data(i,2), 1);
        if ~isempty(lon_idx) && ~isempty(lat_idx)
            mean_dir_grid(lat_idx, lon_idx) = composite_data(i,3);
            mean_per_grid(lat_idx, lon_idx) = composite_data(i,4);
            max_per_grid(lat_idx, lon_idx) = composite_data(i,5);
        end
    end
    
    % Create composite output directory
    comp_outdir = fullfile(outputPath, 'Composite');
    if ~exist(comp_outdir, 'dir'), mkdir(comp_outdir); end
    
    % Plot composite maps
    plot_dir_map(unique_lons, unique_lats, mean_dir_grid, ...
        'Mean Direction - Priority Composite', fullfile(comp_outdir, 'composite_mean_dir.png'));
    plot_map(unique_lons, unique_lats, mean_per_grid, ...
        'Mean Period - Priority Composite', fullfile(comp_outdir, 'composite_mean_per.png'));
    plot_map(unique_lons, unique_lats, max_per_grid, ...
        'Max Period - Priority Composite', fullfile(comp_outdir, 'composite_max_per.png'));
    
    fprintf('Composite maps saved to: %s\n', comp_outdir);
end

% Create individual dataset maps
for d = 1:length(datasets)
    name = datasets(d).name;
    if ~isfield(all_data, name) || ~isfield(all_data.(name), 'stats'), continue; end
    
    stats = all_data.(name).stats;
    outdir = fullfile(outputPath, name);
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    
    % Direction and period maps only
    plot_dir_map(stats.longitude, stats.latitude, stats.mean_dir, ...
        sprintf('Mean Direction - %s', datasets(d).description), fullfile(outdir, 'mean_dir.png'));
    plot_map(stats.longitude, stats.latitude, stats.mean_per, ...
        sprintf('Mean Period - %s', datasets(d).description), fullfile(outdir, 'mean_per.png'));
    plot_map(stats.longitude, stats.latitude, stats.max_per, ...
        sprintf('Max Period - %s', datasets(d).description), fullfile(outdir, 'max_per.png'));
end

%% EXPORT FOR R
fprintf('\n=== Exporting for R ===\n');

% Check if composite data exists

% Export individual datasets
for d = 1:length(datasets)
    name = datasets(d).name;
    if ~isfield(all_data, name) || ~isfield(all_data.(name), 'stats'), continue; end
    
    stats = all_data.(name).stats;
    [lon_grid, lat_grid] = meshgrid(stats.longitude, stats.latitude);
    
    lon_vec = lon_grid(:);
    lat_vec = lat_grid(:);
    temp = stats.mean_dir';
    mean_dir_vec = temp(:);
    temp = stats.mean_per';
    mean_per_vec = temp(:);
    temp = stats.max_per';
    max_per_vec = temp(:);
    
    T = table(lon_vec, lat_vec, mean_dir_vec, mean_per_vec, max_per_vec, ...
        'VariableNames', {'lon', 'lat', 'mean_dir', 'mean_per', 'max_per'});
    
    writetable(T, fullfile(outputPath, sprintf('swan_%s_summary.csv', name)));
    fprintf('Exported %s\n', name);
end

fprintf('\n=== Complete ===\n');
if SAVE_RAW_DATA
    fprintf('Data saved to: %s\n', dataPath);
else
    fprintf('Data saved to: Not saved (temp files only)\n');
end
fprintf('Output saved to: %s\n', outputPath);

%% HELPER FUNCTIONS
function plot_map(lon, lat, data, ttl, fname)
    figure('Position', [100 100 600 500], 'Visible', 'off');
    imagesc(lon, lat, data');
    set(gca, 'YDir', 'normal');
    colorbar; colormap('jet');
    title(ttl, 'FontSize', 12);
    xlabel('Longitude'); ylabel('Latitude');
    saveas(gcf, fname); close(gcf);
end

function plot_dir_map(lon, lat, data, ttl, fname)
    figure('Position', [100 100 600 500], 'Visible', 'off');
    imagesc(lon, lat, data');
    set(gca, 'YDir', 'normal');
    colorbar; colormap('hsv'); clim([0 360]);
    title(ttl, 'FontSize', 12);
    xlabel('Longitude'); ylabel('Latitude');
    saveas(gcf, fname); close(gcf);
end