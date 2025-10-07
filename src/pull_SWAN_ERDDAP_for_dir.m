%% Streamlined CARICOOS SWAN Wave Data Download Script
%   Updated: October 2025
%   Focus: Wave Direction Composite

clear;clc

%% CONFIGURATION
USE_ALREADY_DOWNLOADED = true;
USE_CUSTOM_DATES = false;
CUSTOM_START_DATE = '2024-01-03';
CUSTOM_END_DATE = '2024-12-31';
CHUNK_SIZE_DAYS = 1;
SAVE_RAW_DATA = true;
SKIP_EXISTING_FILES = true;
DOWNLOAD_TIMEOUT = 300;
MAX_RETRIES = 3;

dataPath = 'D:\SWAN_ERDDAP\data';
outputPath = 'D:\SWAN_ERDDAP\output';

if ~exist(dataPath, 'dir'), mkdir(dataPath); end
if ~exist(outputPath, 'dir'), mkdir(outputPath); end

%% DATASET DEFINITIONS
datasets = struct([]);

datasets(1).name = 'USVI_HR';
datasets(1).id = 'SWAN_HighRes_USVI';
datasets(1).server = 'http://dm3.caricoos.org/erddap';
datasets(1).variables = {'dir'};
datasets(1).description = 'USVI High-Res (100m)';
datasets(1).resolution = 100;
datasets(1).lat_bounds = '(18.18):(18.8)';
datasets(1).lon_bounds = '(-65.2):(-64.0)';

datasets(2).name = 'PuertoRico_HR';
datasets(2).id = 'SWAN_HighRes_PR';
datasets(2).server = 'http://dm3.caricoos.org/erddap';
datasets(2).variables = {'dir'};
datasets(2).description = 'PR High-Res (120m)';
datasets(2).resolution = 120;
datasets(2).lat_bounds = '(17.8):(18.6)';
datasets(2).lon_bounds = '(-67.5):(-65.1)';

datasets(3).name = 'StCroix_HR';
datasets(3).id = 'SWAN_HighRes_StCroix';
datasets(3).server = 'http://dm3.caricoos.org/erddap';
datasets(3).variables = {'dir'};
datasets(3).description = 'StCroix High-Res (180m)';
datasets(3).resolution = 180;
datasets(3).lat_bounds = '(17.63):(17.825)';
datasets(3).lon_bounds = '(-65.0):(-64.48)';

datasets(4).name = 'Regional_PRVI';
datasets(4).id = 'swan_1km_caricoos';
datasets(4).server = 'http://dm3.caricoos.org/erddap';
datasets(4).variables = {'Dir'};
datasets(4).description = 'Regional (1km)';
datasets(4).resolution = 1000;
datasets(4).lat_bounds = '(17.0):(19.5)';
datasets(4).lon_bounds = '(-68.0):(-64.0)';

%% KEEP PC AWAKE
try
    awake_file = fullfile(tempdir, 'matlab_keep_awake.txt');
    keep_awake_timer = timer('ExecutionMode', 'fixedRate', 'Period', 120, ...
        'TimerFcn', @(~,~) fopen(awake_file, 'w'));
    start(keep_awake_timer);
    fprintf('Keep-awake timer started (pings every 2 min)\n');
catch
    fprintf('Warning: Could not start keep-awake timer (not critical)\n');
    keep_awake_timer = [];
end

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
    actual_starts = cellfun(@(x) datetime(x{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z'''), {datasets.start});
    actual_ends = cellfun(@(x) datetime(x{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z'''), {datasets.end});
    server_start = max(actual_starts);
    server_end = min(actual_ends);
    
    custom_start = datetime(sprintf('%sT12:00:00Z', CUSTOM_START_DATE), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    custom_end = datetime(sprintf('%sT12:00:00Z', CUSTOM_END_DATE), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    
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

%% DOWNLOAD OR SKIP
if ~USE_ALREADY_DOWNLOADED
    fprintf('\n=== Downloading Data ===\n');
    
    total_chunks = length(chunk_starts) * length(datasets);
    overall_progress = 0;
    overall_start_time = tic;
    chunk_times = [];

    for d = 1:length(datasets)
        fprintf('\n[%d/%d] %s\n', d, length(datasets), datasets(d).description);
        
        h = waitbar(0, datasets(d).name);
        dataset_start_time = tic;
        
        for c = 1:length(chunk_starts)
            chunk_start_time = tic;
            chunk_start = chunk_starts(c);
            chunk_end = min(chunk_start + days(CHUNK_SIZE_DAYS-1), overlap_end);
            
            dataset_progress = (c-1) / length(chunk_starts) * 100;
            global_progress = (overall_progress + (c-1)) / total_chunks * 100;
            
            start_str = sprintf('%sT00:00:00Z', char(chunk_start, 'yyyy-MM-dd'));
            end_str = sprintf('%sT23:59:59Z', char(chunk_end, 'yyyy-MM-dd'));
            var_list = strjoin(arrayfun(@(i) sprintf('%s[(%s):(%s)][%s][%s]', ...
                datasets(d).variables{i}, start_str, end_str, ...
                datasets(d).lat_bounds, datasets(d).lon_bounds), ...
                1:length(datasets(d).variables), 'UniformOutput', false), ',');
            url = sprintf('%s/griddap/%s.nc?%s', datasets(d).server, datasets(d).id, var_list);
            
            filename = fullfile(dataPath, sprintf('%s_%s.nc', datasets(d).name, char(chunk_start, 'yyyyMMdd')));
            
            if SKIP_EXISTING_FILES && exist(filename, 'file')
                chunk_elapsed = toc(chunk_start_time);
                chunk_times = [chunk_times; chunk_elapsed];
                continue;
            end
            
            download_success = false;
            retry_count = 0;
            
            while ~download_success && retry_count < MAX_RETRIES
                try
                    options = weboptions('Timeout', DOWNLOAD_TIMEOUT);
                    websave(filename, url, options);
                    download_success = true;
                catch ME
                    retry_count = retry_count + 1;
                    if retry_count < MAX_RETRIES
                        fprintf('  Download failed (attempt %d/%d), retrying in 10s...\n', retry_count, MAX_RETRIES);
                        pause(10);
                    else
                        warning('Failed chunk %d after %d attempts: %s', c, MAX_RETRIES, ME.message);
                    end
                end
            end
            
            chunk_elapsed = toc(chunk_start_time);
            chunk_times = [chunk_times; chunk_elapsed];
            
            if length(chunk_times) >= 3
                avg_chunk_time = mean(chunk_times(max(1,end-9):end));
                chunks_remaining_dataset = length(chunk_starts) - c;
                dataset_eta_sec = chunks_remaining_dataset * avg_chunk_time;
                chunks_remaining_overall = total_chunks - (overall_progress + c);
                overall_eta_sec = chunks_remaining_overall * avg_chunk_time;
                
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
                
                waitbar_msg = sprintf('%s: %d/%d (%.1f%%)\nDataset ETA: %s | Overall ETA: %s\nSpeed: %.1fs/chunk | Global: %.1f%%', ...
                    datasets(d).name, c, length(chunk_starts), dataset_progress, ...
                    dataset_eta_str, overall_eta_str, avg_chunk_time, global_progress);
                waitbar(c/length(chunk_starts), h, waitbar_msg);
                
                if mod(c, 10) == 0 || c == length(chunk_starts)
                    fprintf('  Progress: %d/%d | Dataset ETA: %s | Overall: %.1f%% (ETA: %s)\n', ...
                        c, length(chunk_starts), dataset_eta_str, global_progress, overall_eta_str);
                end
            else
                waitbar(c/length(chunk_starts), h, sprintf('%s: %d/%d (%.1f%%) - Calculating ETA...', ...
                    datasets(d).name, c, length(chunk_starts), dataset_progress));
            end
        end
        
        close(h);
        overall_progress = overall_progress + length(chunk_starts);
        
        dataset_elapsed = toc(dataset_start_time);
        if dataset_elapsed < 60
            elapsed_str = sprintf('%.0fs', dataset_elapsed);
        elseif dataset_elapsed < 3600
            elapsed_str = sprintf('%dm %.0fs', floor(dataset_elapsed/60), mod(dataset_elapsed,60));
        else
            elapsed_str = sprintf('%dh %dm', floor(dataset_elapsed/3600), floor(mod(dataset_elapsed,3600)/60));
        end
        
        fprintf('Complete in %s\n', elapsed_str);
    end

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
    fprintf('Data saved to: %s\n', dataPath);
    
else
    fprintf('\n=== Using Already Downloaded Data ===\n');
    fprintf('Skipping download, proceeding to statistics...\n');
end

%% CALCULATE STATISTICS (only if using already downloaded data)
if USE_ALREADY_DOWNLOADED
    fprintf('\n=== Calculating Statistics ===\n');
    all_data = struct();

    for d = 1:length(datasets)
        fprintf('Processing %s... ', datasets(d).name);
        
        pattern = fullfile(dataPath, sprintf('%s_*.nc', datasets(d).name));
        files = dir(pattern);
        
        if isempty(files)
            warning('No files found for %s. Skipping.', datasets(d).name);
            continue;
        end
        
        fprintf('found %d files\n', length(files));
        
        % Read metadata from first file
        first_file = fullfile(files(1).folder, files(1).name);
        data.latitude = ncread(first_file, 'latitude');
        data.longitude = ncread(first_file, 'longitude');
        
        % Initialize accumulators
        nlat = length(data.latitude);
        nlon = length(data.longitude);
        sum_sin = zeros(nlon, nlat);
        sum_cos = zeros(nlon, nlat);
        data.times = [];
        
        var_name = datasets(d).variables{1};
        
        fprintf('  Computing mean: ');
        for f = 1:length(files)
            if mod(f, 50) == 0, fprintf('%d ', f); end
            filepath = fullfile(files(f).folder, files(f).name);
            
            dir_chunk = ncread(filepath, var_name);
            sum_sin = sum_sin + sum(sind(dir_chunk), 3, 'omitnan');
            sum_cos = sum_cos + sum(cosd(dir_chunk), 3, 'omitnan');
            
            chunk_time = ncread(filepath, 'time');
            data.times = [data.times; chunk_time(:)];
        end
        fprintf('done\n');
        
        % Calculate circular mean
        data.stats.mean_dir = atan2d(sum_sin, sum_cos);
        data.stats.mean_dir(data.stats.mean_dir < 0) = data.stats.mean_dir(data.stats.mean_dir < 0) + 360;
        data.stats.longitude = data.longitude;
        data.stats.latitude = data.latitude;
        data.stats.times = data.times;
        
        all_data.(datasets(d).name) = data;
        fprintf('  Complete: %d timesteps, mean calculated\n', length(data.times));
    end

    %% CREATE PRIORITY-BASED COMPOSITE
    fprintf('\n=== Creating Priority Composite ===\n');

    n_datasets = length(datasets);
    resolutions = zeros(1, n_datasets);
    for i = 1:n_datasets
        resolutions(i) = datasets(i).resolution;
    end

    [~, sort_idx] = sort(resolutions, 'ascend');

    priority_order = cell(1, n_datasets);
    for i = 1:n_datasets
        priority_order{i} = datasets(sort_idx(i)).name;
    end

    fprintf('Priority order (high to low res): %s\n', strjoin(priority_order, ', '));

    % Pre-compute bounding boxes
    bounds = struct();
    for p = 1:length(priority_order)
        name = priority_order{p};
        if ~isfield(all_data, name) || ~isfield(all_data.(name), 'stats'), continue; end
        
        stats = all_data.(name).stats;
        temp_dir = stats.mean_dir';
        valid_mask = ~isnan(temp_dir(:));
        
        if any(valid_mask)
            [lon_grid, lat_grid] = meshgrid(stats.longitude, stats.latitude);
            lon_valid = lon_grid(valid_mask);
            lat_valid = lat_grid(valid_mask);
            
            bounds.(name).lon_min = min(lon_valid);
            bounds.(name).lon_max = max(lon_valid);
            bounds.(name).lat_min = min(lat_valid);
            bounds.(name).lat_max = max(lat_valid);
        end
    end

    composite_data = [];

    for p = 1:length(priority_order)
        name = priority_order{p};
        if ~isfield(all_data, name) || ~isfield(all_data.(name), 'stats'), continue; end
        
        stats = all_data.(name).stats;
        [lon_grid, lat_grid] = meshgrid(stats.longitude, stats.latitude);
        
        fprintf('Processing %s (priority %d)... ', name, p);
        
        temp_dir = stats.mean_dir';
        valid_in_current = ~isnan(temp_dir(:));
        
        lon_vec = lon_grid(:);
        lat_vec = lat_grid(:);
        
        covered = false(size(lon_vec));
        for h = 1:(p-1)
            higher_name = priority_order{h};
            if ~isfield(bounds, higher_name), continue; end
            
            in_bounds = (lon_vec >= bounds.(higher_name).lon_min) & ...
                        (lon_vec <= bounds.(higher_name).lon_max) & ...
                        (lat_vec >= bounds.(higher_name).lat_min) & ...
                        (lat_vec <= bounds.(higher_name).lat_max);
            
            covered = covered | in_bounds;
        end
        
        keep_idx = find(valid_in_current & ~covered);
        n_points = length(keep_idx);
        
        if n_points > 0
            new_data = zeros(n_points, 4);
            new_data(:, 1) = lon_vec(keep_idx);
            new_data(:, 2) = lat_vec(keep_idx);
            new_data(:, 3) = temp_dir(keep_idx);
            new_data(:, 4) = p;
            
            composite_data = [composite_data; new_data];
        end
        
        fprintf('%d points added\n', n_points);
    end

    fprintf('Total composite points: %d\n', size(composite_data, 1));

    %% CREATE MAPS
    fprintf('\n=== Creating Maps ===\n');

    if ~isempty(composite_data)
        fprintf('Creating composite map...\n');
        
        comp_outdir = fullfile(outputPath, 'Composite');
        if ~exist(comp_outdir, 'dir'), mkdir(comp_outdir); end
        
        lon = composite_data(:,1);
        lat = composite_data(:,2);
        mean_dir = composite_data(:,3);
        
        figure('Position', [100 100 600 500], 'Visible', 'off');
        scatter(lon, lat, 2, mean_dir, 'filled');
        colorbar; colormap('hsv'); clim([0 360]);
        title('Mean Wave Direction - Priority Composite', 'FontSize', 12);
        xlabel('Longitude'); ylabel('Latitude');
        axis equal; grid on;
        saveas(gcf, fullfile(comp_outdir, 'composite_mean_dir.png')); 
        close(gcf);
        
        fprintf('Composite map saved to: %s\n', comp_outdir);
    end

    for d = 1:length(datasets)
        name = datasets(d).name;
        if ~isfield(all_data, name) || ~isfield(all_data.(name), 'stats'), continue; end
        
        stats = all_data.(name).stats;
        outdir = fullfile(outputPath, name);
        if ~exist(outdir, 'dir'), mkdir(outdir); end
        
        [lon_grid, lat_grid] = meshgrid(stats.longitude, stats.latitude);
        lon_vec = lon_grid(:);
        lat_vec = lat_grid(:);
        
        temp_dir = stats.mean_dir';
        mean_dir_vec = temp_dir(:);
        
        figure('Position', [100 100 600 500], 'Visible', 'off');
        scatter(lon_vec, lat_vec, 2, mean_dir_vec, 'filled');
        colorbar; colormap('hsv'); clim([0 360]);
        title(sprintf('Mean Wave Direction - %s', datasets(d).description), 'FontSize', 12);
        xlabel('Longitude'); ylabel('Latitude');
        axis equal; grid on;
        saveas(gcf, fullfile(outdir, 'mean_dir.png')); 
        close(gcf);
    end

    %% EXPORT FOR R
    fprintf('\n=== Exporting Composite for R ===\n');

    if ~isempty(composite_data)
        T = table(composite_data(:,1), composite_data(:,2), composite_data(:,3), ...
                  'VariableNames', {'lon', 'lat', 'mean_dir'});
        
        writetable(T, fullfile(outputPath, 'swan_composite_direction.csv'));
        fprintf('Exported composite CSV: swan_composite_direction.csv\n');
    else
        warning('No composite data to export');
    end

    %% INTERACTIVE FIGURES
    fprintf('\n=== Creating Interactive Figures ===\n');

    if ~isempty(composite_data)
        lon = composite_data(:,1);
        lat = composite_data(:,2);
        mean_dir = composite_data(:,3);
        
        figure('Position', [100 100 800 600], 'Name', 'Composite Wave Direction');
        scatter(lon, lat, 2, mean_dir, 'filled');
        colorbar; colormap('hsv'); clim([0 360]);
        title('Mean Wave Direction - Priority Composite', 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Longitude'); ylabel('Latitude');
        axis equal; grid on;
        
        fprintf('Composite figure created\n');
    end

    for d = 1:length(datasets)
        name = datasets(d).name;
        if ~isfield(all_data, name) || ~isfield(all_data.(name), 'stats'), continue; end
        
        stats = all_data.(name).stats;
        [lon_grid, lat_grid] = meshgrid(stats.longitude, stats.latitude);
        lon_vec = lon_grid(:);
        lat_vec = lat_grid(:);
        
        temp_dir = stats.mean_dir';
        mean_dir_vec = temp_dir(:);
        
        figure('Position', [100 100 800 600], 'Name', sprintf('%s Wave Direction', name));
        scatter(lon_vec, lat_vec, 2, mean_dir_vec, 'filled');
        colorbar; colormap('hsv'); clim([0 360]);
        title(sprintf('Mean Wave Direction - %s', datasets(d).description), 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Longitude'); ylabel('Latitude');
        axis equal; grid on;
        
        fprintf('%s figure created\n', name);
    end

    fprintf('Output saved to: %s\n', outputPath);
    fprintf('Interactive figures ready for inspection\n');
end

%% CLEANUP
fprintf('\n=== Complete ===\n');

if exist('keep_awake_timer', 'var') && ~isempty(keep_awake_timer) && isvalid(keep_awake_timer)
    stop(keep_awake_timer);
    delete(keep_awake_timer);
    fprintf('Keep-awake timer stopped\n');
end