% Create progress bar for this dataset
    h = waitbar(0, sprintf('Starting %s...', current_dataset.description), ...
        'Name', sprintf('Dataset %d/%d Progress', dataset_idx, total_datasets));%% Enhanced Script to pull CARICOOS SWAN wave data from multiple sources
%   Pulls both regional low-resolution and nearshore high-resolution data
%   Updated: August 2025 - DAILY downloads with proper variable names

clear;clc

%% USER CONFIGURATION - TESTING PARAMETERS
% Set USE_CUSTOM_DATES = true for testing, false for full production run
USE_CUSTOM_DATES = true;  % Toggle this for testing vs production

% Custom date range for testing (only used if USE_CUSTOM_DATES = true)
CUSTOM_START_DATE = '2024-01-06';  % Format: 'YYYY-MM-DD' 
CUSTOM_END_DATE = '2024-01-07';    % Format: 'YYYY-MM-DD'

% Flexible timestep configuration - specify chunk size in days
CHUNK_SIZE_DAYS = 1;  % Options: 1 (daily), 7 (weekly), 30 (monthly), or any number

% Quick test presets (uncomment one to use):
% CUSTOM_START_DATE = '2024-01-06'; CUSTOM_END_DATE = '2024-01-08'; CHUNK_SIZE_DAYS = 1;   % 3 days, daily chunks
% CUSTOM_START_DATE = '2024-01-06'; CUSTOM_END_DATE = '2024-01-19'; CHUNK_SIZE_DAYS = 7;   % 2 weeks, weekly chunks  
% CUSTOM_START_DATE = '2024-01-06'; CUSTOM_END_DATE = '2024-02-05'; CHUNK_SIZE_DAYS = 30;  % 1 month, monthly chunks
% CUSTOM_START_DATE = '2024-01-06'; CUSTOM_END_DATE = '2025-08-17'; CHUNK_SIZE_DAYS = 7;   % Full run, weekly chunks

fprintf('=== CONFIGURATION ===\n');
if USE_CUSTOM_DATES
    fprintf('Using CUSTOM date range: %s to %s\n', CUSTOM_START_DATE, CUSTOM_END_DATE);
    custom_start = datetime(CUSTOM_START_DATE, 'InputFormat', 'yyyy-MM-dd');
    custom_end = datetime(CUSTOM_END_DATE, 'InputFormat', 'yyyy-MM-dd');
    estimated_days = days(custom_end - custom_start) + 1;
    estimated_chunks = ceil(estimated_days / CHUNK_SIZE_DAYS);
    estimated_files = estimated_chunks * 4; % 4 datasets
    
    fprintf('Chunk size: %d day(s) per file\n', CHUNK_SIZE_DAYS);
    fprintf('Total period: %d days in %d chunks\n', estimated_days, estimated_chunks);
    fprintf('Estimated files: %d total (%d per dataset)\n', estimated_files, estimated_chunks);
    
    % Estimate file sizes and timing
    approx_mb_per_day = 6.5; % Average across datasets
    approx_mb_per_chunk = approx_mb_per_day * CHUNK_SIZE_DAYS;
    approx_seconds_per_chunk = max(2, approx_mb_per_chunk * 0.3); % Rough estimate
    estimated_minutes = (estimated_files * approx_seconds_per_chunk) / 60;
    
    fprintf('Estimated chunk size: ~%.1f MB each\n', approx_mb_per_chunk);
    fprintf('Estimated total time: ~%.0f minutes\n', estimated_minutes);
    
    if approx_mb_per_chunk > 200
        fprintf('WARNING: Large chunks (>200MB) may cause timeouts!\n');
    elseif approx_mb_per_chunk > 100
        fprintf('CAUTION: Medium chunks (>100MB) - monitor for timeouts\n');
    end
else
    fprintf('Using FULL server date range (will be determined from metadata)\n');
    fprintf('Chunk size: %d day(s) per file\n', CHUNK_SIZE_DAYS);
end
fprintf('======================\n\n');

% Get the project root directory
projectPath = matlab.project.rootProject().RootFolder;

% Define paths relative to the project root
dataPath = fullfile(projectPath, 'data');
srcPath = fullfile(projectPath, 'src');
outputPath = fullfile(projectPath, 'output');
tempPath = fullfile(projectPath, 'temp');

%% Dataset Configuration
% Define all datasets to download - HIGH-RESOLUTION FIRST
datasets = struct();

% HIGH-RESOLUTION datasets first (to test these work properly)
datasets(1).name = 'USVI_HR';
datasets(1).id = 'Historical_7719_73b6_fa7d';
datasets(1).server = 'http://dm3.caricoos.org/erddap';
datasets(1).variables = {'hs', 'dir', 'tp'};
datasets(1).description = 'US Virgin Islands High-Resolution (100m)';
datasets(1).lat_bounds = '(18.18):(18.8)';
datasets(1).lon_bounds = '(-65.2):(-64.0)';

datasets(2).name = 'PuertoRico_HR';
datasets(2).id = 'Historical_9810_6460_3129';
datasets(2).server = 'http://dm3.caricoos.org/erddap';
datasets(2).variables = {'hs', 'dir', 'tp'};
datasets(2).description = 'Puerto Rico High-Resolution (120m)';
datasets(2).lat_bounds = '(17.8):(18.6)';
datasets(2).lon_bounds = '(-67.5):(-65.1)';

datasets(3).name = 'StCroix_HR';
datasets(3).id = 'Historical_257a_2af8_9075';
datasets(3).server = 'http://dm3.caricoos.org/erddap';
datasets(3).variables = {'hs', 'dir', 'tp'};
datasets(3).description = 'St. Croix High-Resolution (180m)';
datasets(3).lat_bounds = '(17.63):(17.825)';
datasets(3).lon_bounds = '(-65.0):(-64.48)';

% Regional dataset last (broader domain, lower resolution)
datasets(4).name = 'Regional_PRVI';
datasets(4).id = 'caricoos_dm3_d848_78a5_baa6';
datasets(4).server = 'http://dm3.caricoos.org/erddap';
datasets(4).variables = {'Hsig', 'Hswell', 'Dir', 'Per'};
datasets(4).description = 'Regional PR-VI Grid (Low Resolution)';
datasets(4).lat_bounds = '(17.0):(19.5)';
datasets(4).lon_bounds = '(-68.0):(-64.0)';

%% Get available time ranges for each dataset and find overlap
fprintf('=== Checking Available Time Ranges ===\n');

% Initialize arrays to store parsed dates
start_dates = {};
end_dates = {};

for i = 1:length(datasets)
    fprintf('Checking %s (%s)...\n', datasets(i).description, datasets(i).name);
    
    % Get dataset info
    info_url = sprintf('%s/griddap/%s.das', datasets(i).server, datasets(i).id);
    try
        info_text = webread(info_url);
        
        % Extract time coverage start
        start_pattern = 'time_coverage_start\s+"([^"]+)"';
        start_match = regexp(info_text, start_pattern, 'tokens');
        
        % Extract time coverage end
        end_pattern = 'time_coverage_end\s+"([^"]+)"';
        end_match = regexp(info_text, end_pattern, 'tokens');
        
        if ~isempty(start_match) && ~isempty(end_match)
            start_date_str = start_match{1}{1};
            end_date_str = end_match{1}{1};
            
            % Store the date strings
            start_dates{i} = start_date_str;
            end_dates{i} = end_date_str;
            
            fprintf('  Available from: %s\n', start_date_str);
            fprintf('  Available to:   %s\n', end_date_str);
        else
            fprintf('  Time range information not found\n');
            % Set defaults if parsing fails
            start_dates{i} = '2024-01-01T00:00:00Z';
            end_dates{i} = '2025-08-17T12:00:00Z';
        end
        
    catch ME
        fprintf('  Error checking dataset info: %s\n', ME.message);
        % Set defaults if request fails
        start_dates{i} = '2024-01-01T00:00:00Z';
        end_dates{i} = '2025-08-17T12:00:00Z';
    end
    
    datasets(i).available = true; % Assume available unless error
end

%% Find overlapping time period
fprintf('\n=== Finding Overlapping Time Period ===\n');

if USE_CUSTOM_DATES
    % Use custom dates instead of server metadata
    fprintf('Using custom date range (overriding server metadata):\n');
    overlap_start_str = sprintf('%sT12:00:00Z', CUSTOM_START_DATE);
    overlap_end_str = sprintf('%sT12:00:00Z', CUSTOM_END_DATE);
    
    overlap_start = datetime(overlap_start_str, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    overlap_end = datetime(overlap_end_str, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    
    fprintf('  Custom Start: %s\n', overlap_start_str);
    fprintf('  Custom End:   %s\n', overlap_end_str);
    
else
    % Use server metadata to find overlap (original behavior)
    fprintf('Using server metadata to find overlap:\n');
    
    % Convert date strings to datetime objects for comparison
    start_datetimes = datetime(start_dates, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    end_datetimes = datetime(end_dates, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    
    % Find the latest start date and earliest end date
    overlap_start = max(start_datetimes);
    overlap_end = min(end_datetimes);
    
    % Convert back to the required format
    overlap_start_str = datestr(overlap_start, 'yyyy-mm-ddTHH:MM:SSZ');
    overlap_end_str = datestr(overlap_end, 'yyyy-mm-ddTHH:MM:SSZ');
    
    fprintf('  Server Start: %s\n', overlap_start_str);
    fprintf('  Server End:   %s\n', overlap_end_str);
end

% Extract timing info for display
start_year = year(overlap_start);
end_year = year(overlap_end);
start_month = month(overlap_start);
start_day = day(overlap_start);
end_month = month(overlap_end);
end_day = day(overlap_end);

fprintf('Download plan:\n');
fprintf('  Years: %d to %d\n', start_year, end_year);
fprintf('  Start month/day: %d/%d\n', start_month, start_day);
fprintf('  End month/day: %d/%d\n', end_month, end_day);

% Calculate and display total days and chunks
total_days_to_download = days(overlap_end - overlap_start) + 1;
total_chunks = ceil(total_days_to_download / CHUNK_SIZE_DAYS);
fprintf('  Total days: %d\n', total_days_to_download);
fprintf('  Chunk size: %d day(s) per file\n', CHUNK_SIZE_DAYS);
fprintf('  Total chunks: %d per dataset\n', total_chunks);
fprintf('  Total files: %d (4 datasets × %d chunks)\n', total_chunks * 4, total_chunks);

%% Initialize storage for all datasets
all_data = struct();

%% Main download loop for each dataset
total_datasets = length(datasets);
dataset_start_time = tic;

% Create overall progress tracking
fprintf('\n=== Starting Multi-Dataset Download ===\n');

for dataset_idx = 1:total_datasets
    current_dataset = datasets(dataset_idx);
    
    fprintf('\n--- Processing Dataset %d/%d: %s ---\n', ...
        dataset_idx, total_datasets, current_dataset.description);
    
    % Skip if dataset not available
    if ~current_dataset.available
        fprintf('Skipping unavailable dataset\n');
        continue;
    end
    
    % Initialize data arrays for current dataset
    dataset_data = struct();
    for var_idx = 1:length(current_dataset.variables)
        var_name = current_dataset.variables{var_idx};
        dataset_data.(var_name) = [];
    end
    dataset_data.times = [];
    dataset_data.first_chunk = true;
    
    % Calculate chunks to download based on overlap period and chunk size
    overlap_start_date = datetime(overlap_start_str, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    overlap_end_date = datetime(overlap_end_str, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z''');
    
    % Create array of all chunk start dates
    chunk_starts = overlap_start_date:days(CHUNK_SIZE_DAYS):overlap_end_date;
    total_chunks = length(chunk_starts);
    current_chunk_count = 0;
    
    % Create progress bar for this dataset
    h = waitbar(0, sprintf('Starting %s...', current_dataset.description), ...
        'Name', sprintf('Dataset %d/%d Progress', dataset_idx, total_datasets));
    
    % Download data CHUNK BY CHUNK within the overlap period
    try
        for chunk_idx = 1:total_chunks
            current_chunk_count = current_chunk_count + 1;
            chunk_start_time = tic;
            
            % Calculate chunk boundaries
            chunk_start_date = chunk_starts(chunk_idx);
            chunk_end_date = min(chunk_start_date + days(CHUNK_SIZE_DAYS - 1), overlap_end_date);
            actual_days_in_chunk = days(chunk_end_date - chunk_start_date) + 1;
            
            % Update progress bar
            progress = current_chunk_count / total_chunks;
            if CHUNK_SIZE_DAYS == 1
                status_msg = sprintf('%s: %s (%d/%d chunks)', ...
                    current_dataset.name, datestr(chunk_start_date, 'yyyy-mm-dd'), current_chunk_count, total_chunks);
            else
                status_msg = sprintf('%s: %s to %s (%d/%d chunks, %d days)', ...
                    current_dataset.name, datestr(chunk_start_date, 'yyyy-mm-dd'), ...
                    datestr(chunk_end_date, 'yyyy-mm-dd'), current_chunk_count, total_chunks, actual_days_in_chunk);
            end
            waitbar(progress, h, status_msg);
            
            % Construct chunk URL boundaries
            start_date = sprintf('%sT00:00:00Z', datestr(chunk_start_date, 'yyyy-mm-dd'));
            end_date = sprintf('%sT23:59:59Z', datestr(chunk_end_date, 'yyyy-mm-dd'));
            
            % Build variable list for URL using CORRECT coordinate bounds for each dataset
            var_list = '';
            for var_idx = 1:length(current_dataset.variables)
                var_name = current_dataset.variables{var_idx};
                if var_idx > 1
                    var_list = [var_list ','];
                end
                
                % Use dataset-specific coordinate bounds
                coord_bounds = sprintf('[%s][%s]', current_dataset.lat_bounds, current_dataset.lon_bounds);
                
                var_list = sprintf('%s%s[(%s):(%s)]%s', var_list, var_name, ...
                    start_date, end_date, coord_bounds);
            end
            
            % Construct full URL
            url = sprintf('%s/griddap/%s.nc?%s', ...
                current_dataset.server, current_dataset.id, var_list);
            
            % DETAILED LOGGING - Show exactly what we're requesting
            fprintf('\n  === CHUNK DOWNLOAD DETAILS ===\n');
            fprintf('  Dataset: %s\n', current_dataset.description);
            fprintf('  Variables: %s\n', strjoin(current_dataset.variables, ', '));
            if CHUNK_SIZE_DAYS == 1
                fprintf('  Date: %s (1 day)\n', datestr(chunk_start_date, 'yyyy-mm-dd'));
            else
                fprintf('  Period: %s to %s (%d days)\n', ...
                    datestr(chunk_start_date, 'yyyy-mm-dd'), datestr(chunk_end_date, 'yyyy-mm-dd'), actual_days_in_chunk);
            end
            fprintf('  Temporal: %s to %s\n', start_date, end_date);
            fprintf('  Spatial: Lat %s, Lon %s\n', current_dataset.lat_bounds, current_dataset.lon_bounds);
            
            % Calculate approximate grid size for estimation
            lat_range = sscanf(current_dataset.lat_bounds, '(%f):(%f)');
            lon_range = sscanf(current_dataset.lon_bounds, '(%f):(%f)');
            if length(lat_range) == 2 && length(lon_range) == 2
                lat_span = lat_range(2) - lat_range(1);
                lon_span = lon_range(2) - lon_range(1);
                fprintf('  Spatial span: %.3f° lat x %.3f° lon\n', lat_span, lon_span);
                
                % Rough grid size estimates (these are approximate)
                if strcmp(current_dataset.name, 'USVI_HR')
                    approx_res = 0.0009; % ~100m in degrees
                elseif strcmp(current_dataset.name, 'PuertoRico_HR')
                    approx_res = 0.0011; % ~120m in degrees
                elseif strcmp(current_dataset.name, 'StCroix_HR')
                    approx_res = 0.0016; % ~180m in degrees
                else
                    approx_res = 0.01; % ~1km in degrees
                end
                
                approx_nlat = round(lat_span / approx_res);
                approx_nlon = round(lon_span / approx_res);
                fprintf('  Estimated grid: %d x %d points\n', approx_nlon, approx_nlat);
                
                % Chunk time steps (3-hourly data = 8 per day)
                approx_time_steps = actual_days_in_chunk * 8;
                fprintf('  Estimated time steps: %d (%d days × 8 per day)\n', approx_time_steps, actual_days_in_chunk);
                
                % Rough file size estimate (4 bytes per float, 3-4 variables)
                num_vars = length(current_dataset.variables);
                approx_size_mb = (approx_nlat * approx_nlon * approx_time_steps * num_vars * 4) / (1024^2);
                fprintf('  Estimated file size: %.1f MB\n', approx_size_mb);
                
                if approx_size_mb > 200
                    fprintf('  WARNING: Large chunk - high timeout risk!\n');
                elseif approx_size_mb > 100
                    fprintf('  CAUTION: Medium chunk - monitor for timeouts\n');
                elseif approx_size_mb > 50
                    fprintf('  NOTE: Larger chunk - may take longer\n');
                end
            end
            fprintf('  ===============================\n\n');
            
            % Download with retry logic
            if CHUNK_SIZE_DAYS == 1
                filename = sprintf('%s_%s.nc', current_dataset.name, datestr(chunk_start_date, 'yyyymmdd'));
            else
                filename = sprintf('%s_%s_to_%s.nc', current_dataset.name, ...
                    datestr(chunk_start_date, 'yyyymmdd'), datestr(chunk_end_date, 'yyyymmdd'));
            end
            
            download_success = false;
            retry_count = 0;
            max_retries = 3; % Reduced retries for large files
            
            fprintf('  Attempting download: %s\n', filename);
            download_start_time = tic;
            
            while ~download_success && retry_count < max_retries
                try
                    websave(filename, url);
                    download_time = toc(download_start_time);
                    
                    % Check actual file size
                    file_info = dir(filename);
                    actual_size_mb = file_info.bytes / (1024^2);
                    
                    fprintf('  ✓ Downloaded successfully: %s (%.1f MB in %.1f sec)\n', ...
                        filename, actual_size_mb, download_time);
                    download_success = true;
                    
                catch ME
                    retry_count = retry_count + 1;
                    download_time = toc(download_start_time);
                    
                    if retry_count < max_retries
                        fprintf('  ✗ Download failed (attempt %d/%d, %.1f sec): %s\n', ...
                            retry_count, max_retries, download_time, ME.message);
                        fprintf('    Retrying in 5 seconds...\n');
                        pause(5);
                        download_start_time = tic; % Reset timer for retry
                    else
                        fprintf('  ✗ Download FAILED after %d attempts (%.1f sec): %s\n', ...
                            max_retries, download_time, ME.message);
                        fprintf('    Error details: %s\n', ME.message);
                        
                        % Suggest potential solutions
                        if contains(ME.message, 'timeout') || contains(ME.message, 'time')
                            fprintf('    → Likely cause: Chunk too large for server timeout\n');
                            fprintf('    → Suggestion: Reduce CHUNK_SIZE_DAYS (try %d days)\n', max(1, CHUNK_SIZE_DAYS/2));
                        elseif contains(ME.message, '500') || contains(ME.message, 'server')
                            fprintf('    → Likely cause: Server error or overload\n');
                            fprintf('    → Suggestion: Try again later or contact CARICOOS\n');
                        elseif contains(ME.message, 'coordinate') || contains(ME.message, 'bound')
                            fprintf('    → Likely cause: Invalid coordinate bounds\n');
                            fprintf('    → Suggestion: Check spatial domain specifications\n');
                        end
                    end
                end
            end
            
            if ~download_success
                continue; % Skip to next chunk
            end
            
            % Read chunk data
            try
                for var_idx = 1:length(current_dataset.variables)
                    var_name = current_dataset.variables{var_idx};
                    chunk_data = ncread(filename, var_name);
                    
                    if dataset_data.first_chunk
                        dataset_data.(var_name) = chunk_data;
                    else
                        % Concatenate along time dimension (3rd dimension)
                        dataset_data.(var_name) = cat(3, dataset_data.(var_name), chunk_data);
                    end
                end
                
                % Read time data
                chunk_time = ncread(filename, 'time');
                if dataset_data.first_chunk
                    dataset_data.times = chunk_time(:);
                    
                    % Read coordinates on first iteration
                    dataset_data.latitude = ncread(filename, 'latitude');
                    dataset_data.longitude = ncread(filename, 'longitude');
                    
                    fprintf('  Grid dimensions: %d x %d\n', ...
                        length(dataset_data.longitude), length(dataset_data.latitude));
                    
                    dataset_data.first_chunk = false;
                else
                    dataset_data.times = [dataset_data.times; chunk_time(:)];
                end
                
                if CHUNK_SIZE_DAYS == 1
                    fprintf('  Day %s: Added %d time steps (Total so far: %d)\n', ...
                        datestr(chunk_start_date, 'yyyy-mm-dd'), length(chunk_time), length(dataset_data.times));
                else
                    fprintf('  Chunk %d: Added %d time steps from %d days (Total so far: %d)\n', ...
                        chunk_idx, length(chunk_time), actual_days_in_chunk, length(dataset_data.times));
                end
                
            catch ME
                warning('Failed to read data from %s: %s', filename, ME.message);
            end
            
            % Clean up file
            if exist(filename, 'file')
                delete(filename);
            end
            
            % Show progress
            chunk_time = toc(chunk_start_time);
            fprintf('  Progress: %d/%d chunks (%.1f%%) | Chunk time: %.1f sec\n', ...
                current_chunk_count, total_chunks, progress*100, chunk_time);
        end
        
        % Close progress bar
        close(h);
        
    catch ME
        if ishandle(h)
            close(h);
        end
        warning('Error processing dataset %s: %s', current_dataset.name, ME.message);
        continue;
    end
    
    % Store dataset in main structure
    all_data.(current_dataset.name) = dataset_data;
    
    fprintf('Dataset %s complete: %d time steps total\n', ...
        current_dataset.name, length(dataset_data.times));
end

total_download_time = toc(dataset_start_time);

%% Calculate Statistics Across All Datasets
fprintf('\n=== Calculating Multi-Dataset Statistics ===\n');

% Process each dataset and calculate statistics
for dataset_idx = 1:length(datasets)
    dataset_name = datasets(dataset_idx).name;
    
    if ~isfield(all_data, dataset_name)
        fprintf('Skipping missing dataset: %s\n', dataset_name);
        continue;
    end
    
    data = all_data.(dataset_name);
    fprintf('\nProcessing %s (%s)...\n', datasets(dataset_idx).description, dataset_name);
    
    % Map variable names to standard names for processing
    if isfield(data, 'Hsig')
        % Regional dataset
        hsig_data = data.Hsig;
        dir_data = data.Dir;
        per_data = data.Per;
        has_swell = isfield(data, 'Hswell');
        if has_swell
            hswell_data = data.Hswell;
        end
    else
        % High-resolution dataset
        hsig_data = data.hs;  % hs maps to significant wave height
        dir_data = data.dir;
        per_data = data.tp;   % tp maps to peak period
        has_swell = false;    % High-res datasets don't have swell
    end
    
    % Calculate statistics
    stats = struct();
    stats.longitude = data.longitude;
    stats.latitude = data.latitude;
    stats.times = data.times;
    
    % Significant wave height statistics
    stats.mean_hsig = mean(hsig_data, 3, 'omitnan');
    stats.max_hsig = max(hsig_data, [], 3, 'omitnan');
    stats.min_hsig = min(hsig_data, [], 3, 'omitnan');
    stats.std_hsig = std(hsig_data, 0, 3, 'omitnan');
    
    % Wave direction (use circular statistics)
    stats.mean_dir = atan2d(mean(sind(dir_data), 3, 'omitnan'), ...
                           mean(cosd(dir_data), 3, 'omitnan'));
    stats.mean_dir(stats.mean_dir < 0) = stats.mean_dir(stats.mean_dir < 0) + 360;
    
    % Wave period statistics
    stats.mean_per = mean(per_data, 3, 'omitnan');
    stats.max_per = max(per_data, [], 3, 'omitnan');
    stats.min_per = min(per_data, [], 3, 'omitnan');
    stats.std_per = std(per_data, 0, 3, 'omitnan');
    
    % Swell statistics (if available)
    if has_swell
        stats.mean_hswell = mean(hswell_data, 3, 'omitnan');
        stats.max_hswell = max(hswell_data, [], 3, 'omitnan');
        stats.min_hswell = min(hswell_data, [], 3, 'omitnan');
        stats.std_hswell = std(hswell_data, 0, 3, 'omitnan');
    end
    
    % Overall statistics
    stats.overall_mean_hsig = mean(hsig_data(:), 'omitnan');
    stats.overall_max_hsig = max(hsig_data(:), [], 'omitnan');
    stats.overall_min_hsig = min(hsig_data(:), [], 'omitnan');
    stats.overall_mean_per = mean(per_data(:), 'omitnan');
    stats.overall_max_per = max(per_data(:), [], 'omitnan');
    stats.overall_min_per = min(per_data(:), [], 'omitnan');
    
    if has_swell
        stats.overall_mean_hswell = mean(hswell_data(:), 'omitnan');
        stats.overall_max_hswell = max(hswell_data(:), [], 'omitnan');
        stats.overall_min_hswell = min(hswell_data(:), [], 'omitnan');
    end
    
    % Store statistics
    all_data.(dataset_name).stats = stats;
    
    % Print summary
    fprintf('  Significant Wave Height: Mean=%.2fm, Max=%.2fm, Min=%.2fm\n', ...
        stats.overall_mean_hsig, stats.overall_max_hsig, stats.overall_min_hsig);
    fprintf('  Wave Period: Mean=%.2fs, Max=%.2fs, Min=%.2fs\n', ...
        stats.overall_mean_per, stats.overall_max_per, stats.overall_min_per);
    if has_swell
        fprintf('  Swell Height: Mean=%.2fm, Max=%.2fm, Min=%.2fm\n', ...
            stats.overall_mean_hswell, stats.overall_max_hswell, stats.overall_min_hswell);
    end
end

%% Create Maps for All Datasets
fprintf('\n=== Creating Maps for All Datasets ===\n');

for dataset_idx = 1:length(datasets)
    dataset_name = datasets(dataset_idx).name;
    
    if ~isfield(all_data, dataset_name) || ~isfield(all_data.(dataset_name), 'stats')
        continue;
    end
    
    stats = all_data.(dataset_name).stats;
    dataset_desc = datasets(dataset_idx).description;
    
    fprintf('Creating maps for %s...\n', dataset_desc);
    
    % Create subdirectory for this dataset's maps
    dataset_output_dir = fullfile(outputPath, dataset_name);
    if ~exist(dataset_output_dir, 'dir')
        mkdir(dataset_output_dir);
    end
    
    % Significant Wave Height Maps
    create_wave_map(stats.longitude, stats.latitude, stats.mean_hsig, ...
        sprintf('Mean Significant Wave Height - %s', dataset_desc), ...
        fullfile(dataset_output_dir, 'mean_hsig_map.png'));
    
    create_wave_map(stats.longitude, stats.latitude, stats.max_hsig, ...
        sprintf('Maximum Significant Wave Height - %s', dataset_desc), ...
        fullfile(dataset_output_dir, 'max_hsig_map.png'));
    
    create_wave_map(stats.longitude, stats.latitude, stats.std_hsig, ...
        sprintf('Std Dev Significant Wave Height - %s', dataset_desc), ...
        fullfile(dataset_output_dir, 'std_hsig_map.png'));
    
    % Wave Direction Map
    create_direction_map(stats.longitude, stats.latitude, stats.mean_dir, ...
        sprintf('Mean Wave Direction - %s', dataset_desc), ...
        fullfile(dataset_output_dir, 'mean_dir_map.png'));
    
    % Wave Period Maps
    create_wave_map(stats.longitude, stats.latitude, stats.mean_per, ...
        sprintf('Mean Wave Period - %s', dataset_desc), ...
        fullfile(dataset_output_dir, 'mean_per_map.png'));
    
    create_wave_map(stats.longitude, stats.latitude, stats.max_per, ...
        sprintf('Maximum Wave Period - %s', dataset_desc), ...
        fullfile(dataset_output_dir, 'max_per_map.png'));
    
    % Swell maps (if available)
    if isfield(stats, 'mean_hswell')
        create_wave_map(stats.longitude, stats.latitude, stats.mean_hswell, ...
            sprintf('Mean Swell Height - %s', dataset_desc), ...
            fullfile(dataset_output_dir, 'mean_hswell_map.png'));
        
        create_wave_map(stats.longitude, stats.latitude, stats.max_hswell, ...
            sprintf('Maximum Swell Height - %s', dataset_desc), ...
            fullfile(dataset_output_dir, 'max_hswell_map.png'));
    end
    
    fprintf('  Maps saved to: %s\n', dataset_output_dir);
end

%% Export Data for R Analysis (All Datasets)
fprintf('\n=== Exporting All Data for R Analysis ===\n');

for dataset_idx = 1:length(datasets)
    dataset_name = datasets(dataset_idx).name;
    
    if ~isfield(all_data, dataset_name) || ~isfield(all_data.(dataset_name), 'stats')
        continue;
    end
    
    stats = all_data.(dataset_name).stats;
    fprintf('Exporting %s for R...\n', datasets(dataset_idx).description);
    
    % Create coordinate grids
    [lon_grid, lat_grid] = meshgrid(stats.longitude, stats.latitude);
    lon_vec = lon_grid(:);
    lat_vec = lat_grid(:);
    
    % Flatten summary statistics
    mean_hsig_vec = (stats.mean_hsig')(:);
    max_hsig_vec = (stats.max_hsig')(:);
    min_hsig_vec = (stats.min_hsig')(:);
    std_hsig_vec = (stats.std_hsig')(:);
    mean_dir_vec = (stats.mean_dir')(:);
    mean_per_vec = (stats.mean_per')(:);
    max_per_vec = (stats.max_per')(:);
    min_per_vec = (stats.min_per')(:);
    std_per_vec = (stats.std_per')(:);
    
    % Create summary table
    if isfield(stats, 'mean_hswell')
        % Include swell data
        mean_hswell_vec = (stats.mean_hswell')(:);
        max_hswell_vec = (stats.max_hswell')(:);
        min_hswell_vec = (stats.min_hswell')(:);
        std_hswell_vec = (stats.std_hswell')(:);
        
        summary_table = table(lon_vec, lat_vec, ...
            mean_hsig_vec, max_hsig_vec, min_hsig_vec, std_hsig_vec, ...
            mean_hswell_vec, max_hswell_vec, min_hswell_vec, std_hswell_vec, ...
            mean_dir_vec, mean_per_vec, max_per_vec, min_per_vec, std_per_vec, ...
            'VariableNames', {'longitude', 'latitude', ...
            'mean_hsig', 'max_hsig', 'min_hsig', 'std_hsig', ...
            'mean_hswell', 'max_hswell', 'min_hswell', 'std_hswell', ...
            'mean_dir', 'mean_per', 'max_per', 'min_per', 'std_per'});
    else
        % No swell data
        summary_table = table(lon_vec, lat_vec, ...
            mean_hsig_vec, max_hsig_vec, min_hsig_vec, std_hsig_vec, ...
            mean_dir_vec, mean_per_vec, max_per_vec, min_per_vec, std_per_vec, ...
            'VariableNames', {'longitude', 'latitude', ...
            'mean_hsig', 'max_hsig', 'min_hsig', 'std_hsig', ...
            'mean_dir', 'mean_per', 'max_per', 'min_per', 'std_per'});
    end
    
    % Export dataset-specific CSV
    csv_filename = sprintf('swan_%s_summary_for_R.csv', dataset_name);
    writetable(summary_table, fullfile(outputPath, csv_filename));
    
    fprintf('  Exported: %s\n', csv_filename);
end

%% Create Summary Report
fprintf('\n=== Creating Final Summary Report ===\n');

report_lines = {
    sprintf('CARICOOS Multi-Dataset SWAN Wave Analysis Report - %s', datestr(now))
    '=================================================================='
    ''
    sprintf('Total Processing Time: %.2f minutes', total_download_time/60)
    sprintf('Time Period Analyzed: %s to %s', overlap_start_str, overlap_end_str)
    ''
    'Datasets Processed:'
};

for dataset_idx = 1:length(datasets)
    dataset_name = datasets(dataset_idx).name;
    if isfield(all_data, dataset_name)
        stats = all_data.(dataset_name).stats;
        report_lines{end+1} = sprintf('%d. %s (%s)', dataset_idx, ...
            datasets(dataset_idx).description, dataset_name);
        report_lines{end+1} = sprintf('   Grid: %d x %d, Time steps: %d', ...
            length(stats.longitude), length(stats.latitude), length(stats.times));
        report_lines{end+1} = sprintf('   Hsig: %.2f-%.2fm (mean: %.2fm)', ...
            stats.overall_min_hsig, stats.overall_max_hsig, stats.overall_mean_hsig);
        report_lines{end+1} = '';
    end
end

report_lines{end+1} = 'Output Files Created:';
report_lines{end+1} = '- Individual maps for each dataset in subdirectories';
report_lines{end+1} = '- R-compatible CSV files for each dataset';
report_lines{end+1} = '- This summary report';

% Write report
report_file = fullfile(outputPath, 'multi_dataset_analysis_report.txt');
fileID = fopen(report_file, 'w');
for i = 1:length(report_lines)
    fprintf(fileID, '%s\n', report_lines{i});
end
fclose(fileID);

fprintf('\n=== Analysis Complete ===\n');
fprintf('Total datasets processed: %d\n', sum(isfield(all_data, {datasets.name})));
fprintf('Summary report: %s\n', report_file);
fprintf('Output directory: %s\n', outputPath);

%% Helper Functions

function create_wave_map(longitude, latitude, data, title_str, filename)
    figure('Position', [100, 100, 600, 500], 'Visible', 'off');
    imagesc(longitude, latitude, data');
    set(gca, 'YDir', 'normal');
    colorbar;
    colormap('jet');
    title(title_str, 'FontSize', 14);
    xlabel('Longitude (°W)', 'FontSize', 12);
    ylabel('Latitude (°N)', 'FontSize', 12);
    clim([0 max(data(:))]);
    grid on;
    saveas(gcf, filename);
    close(gcf);
end

function create_direction_map(longitude, latitude, data, title_str, filename)
    figure('Position', [100, 100, 600, 500], 'Visible', 'off');
    imagesc(longitude, latitude, data');
    set(gca, 'YDir', 'normal');
    colorbar;
    colormap('hsv'); % Better for directional data
    title(title_str, 'FontSize', 14);
    xlabel('Longitude (°W)', 'FontSize', 12);
    ylabel('Latitude (°N)', 'FontSize', 12);
    clim([0 360]);
    grid on;
    saveas(gcf, filename);
    close(gcf);
end