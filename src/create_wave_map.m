function create_wave_map(longitude, latitude, data, title_str, filename)
    figure('Position', [100, 100, 600, 500], 'Visible', 'off');
    scatter(longitude, latitude, 4, data, 'filled');
    colorbar;
    colormap('jet');
    title(title_str, 'FontSize', 14);
    xlabel('Longitude (°W)', 'FontSize', 12);
    ylabel('Latitude (°N)', 'FontSize', 12);
    
    % Better color limits
    valid_data = data(~isnan(data));
    if ~isempty(valid_data) && max(valid_data) > 0
        clim([0 max(valid_data)]);
    end
    
    grid on;
    axis equal;
    saveas(gcf, filename);
    close(gcf);
end