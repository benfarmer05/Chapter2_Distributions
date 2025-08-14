% Simple MATLAB script to download Ocean Color data from ERDDAP
% Downloads NOAA VIIRS/OLCI DINEOF gap-filled data for one timepoint
% Variables: Chlorophyll-a and Kd490

clear; clc;

%% Define parameters
base_url = 'https://coastwatch.noaa.gov/erddap/griddap/';

% Time: Use 'last' for most recent available data, or specify a date like '(2024-08-01T12:00:00Z)'
time_str = '(last)';

% Geographic bounds (matching your other datasets)
lat_min = 17.0;
lat_max = 19.5;
lon_min = -68.0;
lon_max = -64.0;

% Depth (surface level)
depth_str = '(0.0)';

% Define the three datasets
datasets = {
    'noaacwNPPN20S3ASCIDINEOF2kmDaily', 'chlor_a', 'Chlorophyll-a Concentration';
    'noaacwNPPN20S3AkdSCIDINEOF2kmDaily', 'kd_490', 'Diffuse Attenuation Coefficient at 490nm';
    'noaacwNPPN20S3AspmSCIDINEOF2kmDaily', 'spm', 'Suspended Particulate Matter'
};

fprintf('Triple Ocean Color Data Download\n');
fprintf('=================================\n');
fprintf('Variables: Chlorophyll-a, Kd490, and Suspended Particulate Matter\n');
fprintf('Geographic bounds: %.1f°N to %.1f°N, %.1f°W to %.1f°W\n\n', lat_min, lat_max, abs(lon_max), abs(lon_min));

% Storage for all three datasets
chlor_a_data = [];
kd490_data = [];
spm_data = [];
longitude = [];
latitude = [];
time_data = [];
altitude = [];

%% Download and process each dataset
for i = 1:size(datasets, 1)
    dataset_id = datasets{i, 1};
    variable = datasets{i, 2};
    var_name = datasets{i, 3};
    
    fprintf('Processing %s (%s)...\n', var_name, variable);
    
    %% Construct download URL for NetCDF format
    query_str = sprintf('%s[%s][%s][(%f):(%f)][(%f):(%f)]', ...
        variable, time_str, depth_str, lat_min, lat_max, lon_min, lon_max);
    
    download_url = sprintf('%s%s.nc?%s', base_url, dataset_id, query_str);
    
    %% Download the data with robust retry logic
    filename = sprintf('%s_data.nc', variable);
    fprintf('Downloading data from ERDDAP...\n');
    fprintf('URL: %s\n', download_url);
    
    % Download with robust retry logic for ERDDAP server issues
    download_success = false;
    retry_count = 0;
    max_retries = 15;  % Increased for problematic datasets
    
    % Set longer timeout for high-resolution data requests
    options = weboptions('Timeout', 60, 'RequestMethod', 'get');
    
    fprintf('Attempting download with retry logic...\n');
    
    try
        while ~download_success && retry_count < max_retries
            retry_count = retry_count + 1;
            
            try
                fprintf('Download attempt %d/%d...\n', retry_count, max_retries);
                websave(filename, download_url, options);
                download_success = true;
                fprintf('Data successfully downloaded to: %s\n', filename);
                
            catch ME
                if retry_count < max_retries
                    fprintf('Download failed (attempt %d/%d): %s\n', retry_count, max_retries, ME.message);
                    
                    if contains(ME.message, '503') || contains(ME.message, 'Service Unavailable')
                        fprintf('Server temporarily unavailable (503 error)\n');
                    elseif contains(ME.message, 'timeout') || contains(ME.message, 'Timeout')
                        fprintf('Request timed out\n');
                    else
                        fprintf('Network error: %s\n', ME.message);
                    end
                    
                    fprintf('Waiting 15 seconds before retry...\n');
                    pause(15);  % Simple 15-second wait for all retries
                else
                    fprintf('\nFailed to download %s after %d attempts.\n', var_name, max_retries);
                    fprintf('Final error: %s\n', ME.message);
                    fprintf('Skipping this variable and continuing...\n\n');
                    break;
                end
            end
        end
        
        if ~download_success
            fprintf('Download failed for %s after all retry attempts. Skipping...\n\n', var_name);
            continue;
        end
        
        %% Read the NetCDF file
        fprintf('Reading NetCDF file...\n');
        
        % Read variables
        temp_longitude = ncread(filename, 'longitude');
        temp_latitude = ncread(filename, 'latitude');
        temp_time = ncread(filename, 'time');
        temp_altitude = ncread(filename, 'altitude');
        temp_data = ncread(filename, variable);
        
        % Store coordinates from first successful download
        if i == 1
            longitude = temp_longitude;
            latitude = temp_latitude;
            time_data = temp_time;
            altitude = temp_altitude;
        end
        
        % Store the specific variable data
        if strcmp(variable, 'chlor_a')
            chlor_a_data = temp_data;
        elseif strcmp(variable, 'kd_490')
            kd490_data = temp_data;
        elseif strcmp(variable, 'spm')
            spm_data = temp_data;
        end
        
        % Get time units and convert to readable date
        try
            time_units = ncreadatt(filename, 'time', 'units');
            fprintf('Time units: %s\n', time_units);
        catch
            fprintf('Time units: not available\n');
        end
        
        % Display basic information
        fprintf('\n--- Dataset Information ---\n');
        fprintf('Dataset: %s\n', var_name);
        fprintf('Time: %s\n', time_str);
        fprintf('Longitude range: %.4f to %.4f (degrees East)\n', min(temp_longitude), max(temp_longitude));
        fprintf('Latitude range: %.4f to %.4f (degrees North)\n', min(temp_latitude), max(temp_latitude));
        fprintf('Altitude: %.1f m (surface)\n', temp_altitude);
        fprintf('Data dimensions: %d x %d x %d x %d (lon x lat x alt x time)\n', size(temp_data));
        
        % Determine units
        if strcmp(variable, 'chlor_a')
            var_units = 'mg m-3';
        elseif strcmp(variable, 'kd_490')
            var_units = 'm-1';
        elseif strcmp(variable, 'spm')
            var_units = 'g m-3';
        end
        fprintf('Units: %s\n', var_units);
        
        % Display data statistics (excluding NaN values)
        valid_data = temp_data(~isnan(temp_data));
        if ~isempty(valid_data)
            fprintf('\n--- Data Statistics ---\n');
            fprintf('Valid data points: %d\n', length(valid_data));
            fprintf('Min %s: %.4f %s\n', variable, min(valid_data), var_units);
            fprintf('Max %s: %.4f %s\n', variable, max(valid_data), var_units);
            fprintf('Mean %s: %.4f %s\n', variable, mean(valid_data), var_units);
            fprintf('Median %s: %.4f %s\n', variable, median(valid_data), var_units);
        else
            fprintf('No valid data found (all NaN values)\n');
        end
        
        % Display additional metadata
        try
            dataset_title = ncreadatt(filename, '/', 'title');
            fprintf('\nDataset title: %s\n', dataset_title);
        catch
        end
        
        try
            institution = ncreadatt(filename, '/', 'institution');
            fprintf('Institution: %s\n', institution);
        catch
        end
        
        try
            source = ncreadatt(filename, '/', 'platform');
            fprintf('Platform/Source: %s\n', source);
        catch
        end
        
        % Clean up file
        if exist(filename, 'file')
            delete(filename);
        end
        
    catch ME
        fprintf('Error processing %s after successful download:\n', var_name);
        fprintf('%s\n', ME.message);
        fprintf('The file was downloaded but there may be an issue with the data format.\n');
    end
    
    fprintf('\n' + string(repmat('=', 1, 40)) + '\n\n');
end

%% Create plots for successfully downloaded variables
fprintf('Creating plots...\n');

% Debug: Check what data we have stored
fprintf('Debug info:\n');
fprintf('- chlor_a_data is empty: %s\n', mat2str(isempty(chlor_a_data)));
fprintf('- kd490_data is empty: %s\n', mat2str(isempty(kd490_data)));
fprintf('- spm_data is empty: %s\n', mat2str(isempty(spm_data)));

% Count successful downloads
plot_count = 0;
if ~isempty(chlor_a_data)
    plot_count = plot_count + 1;
    fprintf('- Chlor_a available for plotting\n');
end
if ~isempty(kd490_data)
    plot_count = plot_count + 1;
    fprintf('- Kd490 available for plotting\n');
end
if ~isempty(spm_data)
    plot_count = plot_count + 1;
    fprintf('- SPM available for plotting\n');
end

fprintf('Total plots to create: %d\n', plot_count);

if plot_count > 0 && ~isempty(longitude)
    % Create meshgrid for plotting
    [lon_grid, lat_grid] = meshgrid(longitude, latitude');
    
    % Determine subplot layout
    if plot_count == 1
        figure('Position', [100, 100, 600, 500]);
    elseif plot_count == 2
        figure('Position', [100, 100, 1200, 500]);
    else
        figure('Position', [100, 100, 1800, 500]);  % Wide layout for 3 plots
    end
    
    current_plot = 0;
    
    % Plot Chlorophyll-a if available
    if ~isempty(chlor_a_data)
        current_plot = current_plot + 1;
        if plot_count > 1
            subplot(1, plot_count, current_plot);
        end
        
        % Plot the data (squeeze to remove singleton dimensions and transpose)
        chlor_a_plot = squeeze(chlor_a_data)';
        pcolor(lon_grid, lat_grid, chlor_a_plot);
        shading interp;
        colormap(jet);
        
        title('Chlorophyll-a Concentration (DINEOF gap-filled)');
        xlabel('Longitude (degrees East)');
        ylabel('Latitude (degrees North)');
        
        % Set color limits with logarithmic scale
        valid_chlor = chlor_a_data(~isnan(chlor_a_data));
        if ~isempty(valid_chlor)
            caxis([0.03, 30]);  % Use the dataset's suggested colorbar range
            set(gca, 'ColorScale', 'log');  % Set logarithmic color scale
        end
        
        % Add colorbar
        c = colorbar;
        c.Label.String = 'Chlorophyll-a (mg m^{-3})';
        
        % Improve plot appearance
        grid on;
        set(gca, 'GridAlpha', 0.3);
    end
    
    % Plot Kd490 if available
    if ~isempty(kd490_data)
        current_plot = current_plot + 1;
        if plot_count > 1
            subplot(1, plot_count, current_plot);
        end
        
        % Plot the data (squeeze to remove singleton dimensions and transpose)
        kd490_plot = squeeze(kd490_data)';
        pcolor(lon_grid, lat_grid, kd490_plot);
        shading interp;
        colormap(jet);
        
        title('Diffuse Attenuation Coefficient at 490nm (DINEOF gap-filled)');
        xlabel('Longitude (degrees East)');
        ylabel('Latitude (degrees North)');
        
        % Set color limits with logarithmic scale
        valid_kd = kd490_data(~isnan(kd490_data));
        if ~isempty(valid_kd)
            caxis([0.01, 10]);  % Appropriate range for Kd490
            set(gca, 'ColorScale', 'log');  % Set logarithmic color scale
        end
        
        % Add colorbar
        c = colorbar;
        c.Label.String = 'Kd490 (m^{-1})';  % Fixed units for Kd490
        
        % Improve plot appearance
        grid on;
        set(gca, 'GridAlpha', 0.3);
    end
    
    % Plot SPM if available
    if ~isempty(spm_data)
        current_plot = current_plot + 1;
        if plot_count > 1
            subplot(1, plot_count, current_plot);
        end
        
        % Plot the data (squeeze to remove singleton dimensions and transpose)
        spm_plot = squeeze(spm_data)';
        pcolor(lon_grid, lat_grid, spm_plot);
        shading interp;
        colormap(jet);
        
        title('Suspended Particulate Matter (DINEOF gap-filled)');
        xlabel('Longitude (degrees East)');
        ylabel('Latitude (degrees North)');
        
        % Set color limits with logarithmic scale
        valid_spm = spm_data(~isnan(spm_data));
        if ~isempty(valid_spm)
            caxis([0.01, 100]);  % Appropriate range for SPM
            set(gca, 'ColorScale', 'log');  % Set logarithmic color scale
        end
        
        % Add colorbar
        c = colorbar;
        c.Label.String = 'SPM (g m^{-3})';
        
        % Improve plot appearance
        grid on;
        set(gca, 'GridAlpha', 0.3);
    end
    
    % Add overall title if multiple variables plotted
    if plot_count > 1
        sgtitle('Ocean Color Variables - NOAA VIIRS/OLCI DINEOF Data', 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    fprintf('Plot(s) created successfully.\n');
else
    fprintf('No data was successfully downloaded for plotting.\n');
end

fprintf('\nScript completed.\n');