clear;clc


matFileDir = '/Users/benja/Downloads/';
matFiles = dir(fullfile(matFileDir, '*.mat'));

if length(matFiles) ~= 1
    error('There should be exactly one .mat file in the specified directory.');
end

matFilePath = fullfile(matFileDir, matFiles(1).name);

loadedData = load(matFilePath);
structName = fieldnames(loadedData);
loadedData = loadedData.(structName{1});
time_values = loadedData.time;


% Initialize an empty datetime array to store the converted dates
converted_dates = datetime.empty;

% Convert the Unix timestamps to datetime values
for i = 1:length(time_values)
    converted_dates(i) = datetime(time_values(i), 'ConvertFrom', 'posixtime');
end

% Display the converted dates and times
disp('Converted date and time values:');
disp(converted_dates);

% 1. Calculate the size of the .mat file
fileInfo = dir(matFilePath);
fileSizeMB = fileInfo.bytes / (1024 * 1024); % Convert bytes to megabytes

% Display the size of the .mat file
disp(['Size of the .mat file: ', num2str(fileSizeMB), ' MB']);

% 2. Calculate the mode periodicity of the time values
time_diffs = diff(time_values); % Calculate the differences between consecutive time values
mode_periodicity_seconds = mode(time_diffs); % Calculate the mode of these differences
mode_periodicity_hours = mode_periodicity_seconds / 3600; % Convert seconds to hours

% Display the mode periodicity
disp(['Mode periodicity of the time values: ', num2str(mode_periodicity_hours), ' hours']);

% 3. Calculate the average MB required per entry
total_entries = length(time_values);
average_mb_per_entry = fileSizeMB / total_entries;

% Display the average MB required per entry
disp(['Average MB required per entry: ', num2str(average_mb_per_entry), ' MB']);


% 5. Calculate the storage requirements for a given date range
% User-specified date range
start_date = datetime(2018, 1, 1);
end_date = datetime(2018, 6, 30);

% Calculate the number of periods within the date range
duration_hours = hours(end_date - start_date);
num_periods = ceil(duration_hours / mode_periodicity_hours);

% Calculate the total storage required for the specified date range
total_storage_mb = num_periods * average_mb_per_entry;

% Display the total storage requirements for the specified date range
disp(['Total storage required from ', datestr(start_date), ' to ', datestr(end_date), ': ', num2str(total_storage_mb), ' MB']);

%%

% Assuming loadedData contains the structure as before with Dir field and total_entries calculated

% Extract dimensions from 'Dir' field
dim2 = size(loadedData.Dir, 2);
dim3 = size(loadedData.Dir, 3);

% Use the total_entries as the size of the 1st dimension
dim1 = num_periods;

% Create a matrix of zeros with the specified dimensions
matrix = zeros(dim1, dim2, dim3);

% Calculate the size of the matrix in bytes
matrix_size_bytes = numel(matrix) * 8;  % Assuming double precision (8 bytes per element)

% Convert size to megabytes
matrix_size_mb = matrix_size_bytes / (1024 * 1024);

% Display the size of the matrix in megabytes
disp(['Size of the matrix: ', num2str(matrix_size_mb), ' MB']);

disp(['Size of all matrices: ', num2str(matrix_size_mb*4), ' MB']);

%%

% Assuming 'matrix' has been defined and initialized with zeros previously

% Get the total number of elements in the matrix
total_elements = numel(matrix);

% Generate random values
random_values = rand(total_elements, 1);  % Generate random values between 0 and 1

tic
% Assign random values to the matrix
matrix(:) = random_values;  % Assign randomly generated values to all elements in 'matrix'
toc

% % Display a sample of the updated matrix
% disp('Sample of the updated matrix:');
% disp(matrix(:,:,1));  % Display the first page of the updated matrix

%%

% Estimate maximum number of elements based on available memory
% For example, assuming 8 GB of available memory
available_memory_gb = 12;  % Adjust this based on your system
bytes_per_element = 4;    % Double precision (8 bytes per element)

max_elements_estimation = (available_memory_gb * 1024^3) / bytes_per_element;

% Display the estimation results
disp(['Total number of elements in the matrix: ', num2str(total_elements)]);
disp(['Estimated maximum number of elements MATLAB can handle: ', num2str(max_elements_estimation)]);
disp(['Percentage of estimated limit used: ', num2str(100 * total_elements / max_elements_estimation), '%']);




%%
clear;clc

% Define the coordinates of the ocean extent
lonRange = [17, 19.5];
latRange = [-68, -64];

% Load coastline data (you need the Mapping Toolbox for this)
load coastlines

% Create figure
figure
hold on

% Plot coastline data
geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'k')

% Set axis properties
axis equal
axis([latRange, lonRange])

% Add labels and title
xlabel('Longitude')
ylabel('Latitude')
title('Map of the Ocean - Caribbean Sea Region')

% Create a web map using ESRI World Imagery basemap and zoom to the specified limits
webmap('oceanbasemap');
