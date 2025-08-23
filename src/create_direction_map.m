function create_direction_map(longitude, latitude, data, title_str, filename)
    figure('Position', [100, 100, 600, 500], 'Visible', 'off');
    scatter(longitude, latitude, 4, data, 'filled');
    colorbar;
    colormap('hsv');
    title(title_str, 'FontSize', 14);
    xlabel('Longitude (°W)', 'FontSize', 12);
    ylabel('Latitude (°N)', 'FontSize', 12);
    clim([0 360]);
    grid on;
    axis equal;
    saveas(gcf, filename);
    close(gcf);
end