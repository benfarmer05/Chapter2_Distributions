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
datasets = struct(...
    'name', {'USVI_HR', 'PuertoRico_HR', 'StCroix_HR', 'Regional_PRVI'}, ...
    'id', {'SWAN_HighRes_USVI', 'SWAN_HighRes_PR', 'SWAN_HighRes_StCroix', 'swan_1km_caricoos'}, ...
    'server', repmat({'http://dm3.caricoos.org/erddap'}, 1, 4), ...
    'variables', {{'hs','dir','tp'}, {'hs','dir','tp'}, {'hs','dir','tp'}, {'Hsig','Dir','Per'}}, ...
    'description', {'USVI High-Res (100m)', 'PR High-Res (120m)', 'StCroix High-Res (180m)', 'Regional (1km)'}, ...
    'lat_bounds', {'(18.18):(18.8)', '(17.8):(18.6)', '(17.63):(17.825)', '(17.0):(19.5)'}, ...
    'lon_bounds', {'(-65.2):(-64.0)', '(-67.5):(-65.1)', '(-65.0):(-64.48)', '(-68.0):(-64.0)'});

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
                    avg_chunk_time = mean(chunk_times(max(1,end-9):end)); % Use last 10 chunks for better accuracy
                    
                    % Dataset ETA
                    chunks_remaining_dataset = length(chunk_starts) - c;
                    dataset_eta_sec = chunks_remaining_dataset * avg_chunk_time;
                    
                    % Overall ETA
                    chunks_remaining_overall = total_chunks - (overall_progress + c);
                    overall_eta_sec = chunks_remaining_overall * avg_chunk_time;
                    
                    % Format time strings (inline to avoid scope issues)
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
                % Check for user interrupt (Ctrl-C)
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
        
        % Format elapsed time (inline)
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

    % Format total elapsed time (inline)
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
    if isfield(data, 'Hsig')
        hs = data.Hsig; dir = data.Dir; per = data.Per;
    else
        hs = data.hs; dir = data.dir; per = data.tp;
    end
    
    stats.longitude = data.longitude;
    stats.latitude = data.latitude;
    stats.times = data.times;
    stats.mean_hsig = mean(hs, 3, 'omitnan');
    stats.max_hsig = max(hs, [], 3, 'omitnan');
    stats.std_hsig = std(hs, 0, 3, 'omitnan');
    stats.mean_dir = atan2d(mean(sind(dir), 3, 'omitnan'), mean(cosd(dir), 3, 'omitnan'));
    stats.mean_dir(stats.mean_dir < 0) = stats.mean_dir(stats.mean_dir < 0) + 360;
    stats.mean_per = mean(per, 3, 'omitnan');
    stats.max_per = max(per, [], 3, 'omitnan');
    
    all_data.(name).stats = stats;
    fprintf('%s: Hsig=%.2f-%.2fm, Per=%.1f-%.1fs\n', name, ...
        min(hs(:), [], 'omitnan'), max(hs(:), [], 'omitnan'), ...
        min(per(:), [], 'omitnan'), max(per(:), [], 'omitnan'));
end

%% CREATE MAPS
fprintf('\n=== Creating Maps ===\n');

for d = 1:length(datasets)
    name = datasets(d).name;
    if ~isfield(all_data, name) || ~isfield(all_data.(name), 'stats'), continue; end
    
    stats = all_data.(name).stats;
    outdir = fullfile(outputPath, name);
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    
    % Wave height maps
    plot_map(stats.longitude, stats.latitude, stats.mean_hsig, ...
        sprintf('Mean Hsig - %s', datasets(d).description), fullfile(outdir, 'mean_hsig.png'));
    plot_map(stats.longitude, stats.latitude, stats.max_hsig, ...
        sprintf('Max Hsig - %s', datasets(d).description), fullfile(outdir, 'max_hsig.png'));
    
    % Direction and period maps
    plot_dir_map(stats.longitude, stats.latitude, stats.mean_dir, ...
        sprintf('Mean Direction - %s', datasets(d).description), fullfile(outdir, 'mean_dir.png'));
    plot_map(stats.longitude, stats.latitude, stats.mean_per, ...
        sprintf('Mean Period - %s', datasets(d).description), fullfile(outdir, 'mean_per.png'));
end

%% EXPORT FOR R
fprintf('\n=== Exporting for R ===\n');

for d = 1:length(datasets)
    name = datasets(d).name;
    if ~isfield(all_data, name) || ~isfield(all_data.(name), 'stats'), continue; end
    
    stats = all_data.(name).stats;
    [lon_grid, lat_grid] = meshgrid(stats.longitude, stats.latitude);
    
    lon_vec = lon_grid(:);
    lat_vec = lat_grid(:);
    temp = stats.mean_hsig';
    mean_hsig_vec = temp(:);
    temp = stats.max_hsig';
    max_hsig_vec = temp(:);
    temp = stats.std_hsig';
    std_hsig_vec = temp(:);
    temp = stats.mean_dir';
    mean_dir_vec = temp(:);
    temp = stats.mean_per';
    mean_per_vec = temp(:);
    temp = stats.max_per';
    max_per_vec = temp(:);
    
    T = table(lon_vec, lat_vec, mean_hsig_vec, max_hsig_vec, std_hsig_vec, ...
        mean_dir_vec, mean_per_vec, max_per_vec, ...
        'VariableNames', {'lon', 'lat', 'mean_hsig', 'max_hsig', 'std_hsig', ...
        'mean_dir', 'mean_per', 'max_per'});
    
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