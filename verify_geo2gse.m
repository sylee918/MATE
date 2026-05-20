%% Coordinate Transformation Verification Script (GSE to GEO)
% This script verifies the correctness of the coordinate transformations in 'GEO2GSE_2020001.txt'.
% It compares the values from the file against an independent, high-precision
% analytical astronomical coordinate conversion model (Hapgood 1992, 1995 / IAU 1982).
%
% Created for: NASA MATE project
% Language: MATLAB

clear; clc; close all;

%% 1. Configuration & File Loading
filename = 'C:\Users\slee122\OneDrive - NASA\Desktop/Work/git/MATE\GEO2GSE_2020001.txt';
gif_filename = 'gse_to_geo_rotation_animation.gif';

fprintf('========================================================================\n');
fprintf('       GEO2GSE Coordinate Transformation Verification Tool             \n');
fprintf('========================================================================\n');

if ~exist(filename, 'file')
    error('오류: %s 파일이 현재 디렉토리에 존재하지 않습니다.', filename);
end

fprintf('Loading "%s"... (이 과정은 수 초가 소요될 수 있습니다.)\n', filename);
tic;
% Optimized reading using detectImportOptions
opts = detectImportOptions(filename);
opts.VariableNamingRule = 'preserve';
data = readtable(filename, opts);
loading_time = toc;
fprintf('Loaded %d data rows successfully in %.2f seconds.\n\n', height(data), loading_time);

% Normalize variable names based on number of detected columns (handles whitespace-split dates)
num_cols = width(data);
fprintf('Detected %d columns in the imported table.\n', num_cols);

if num_cols == 6
    % Date and Time columns were split due to whitespace delimiter!
    % Column 1: Date, Column 2: Time, Column 3: GSE.long, Column 4: GSE.lati, Column 5: GEO.long, Column 6: GEO.lati
    date_col = data{:, 1};
    time_col = data{:, 2};
    
    % Safe conversion to string arrays and concatenate with a space (avoids strcat space-stripping gotcha)
    date_str = string(date_col);
    time_str = string(time_col);
    combined_time = cellstr(date_str + " " + time_str);
    
    % Reconstruct table with standard names
    data = table(combined_time, data{:, 3}, data{:, 4}, data{:, 5}, data{:, 6}, ...
        'VariableNames', {'Time', 'GSE_long', 'GSE_lati', 'GEO_long', 'GEO_lati'});
    
elseif num_cols == 5
    % Columns were read with correct 5-column layout
    cols = data.Properties.VariableNames;
    cols{1} = 'Time';
    cols{2} = 'GSE_long';
    cols{3} = 'GSE_lati';
    cols{4} = 'GEO_long';
    cols{5} = 'GEO_lati';
    data.Properties.VariableNames = cols;
else
    disp('Variable names in imported table:');
    disp(data.Properties.VariableNames);
    error('오류: 데이터 파일의 열(column) 개수가 예상과 다릅니다 (%d개). 파일 포맷을 확인하세요.', num_cols);
end


%% 2. Grid & Time Analysis
unique_times = unique(data.Time);
num_steps = length(unique_times);

% Determine the grid size of the first time step
time1 = unique_times{1};
idx_time1 = strcmp(data.Time, time1);
grid_data = data(idx_time1, :);
num_points = sum(idx_time1);

% GSE long and lat unique values to determine dimensions
gse_longs = unique(grid_data.GSE_long);
gse_latis = unique(grid_data.GSE_lati);
n_lon = length(gse_longs);
n_lat = length(gse_latis);

fprintf('--- Grid Information ---\n');
fprintf('총 시간 스텝 수 (Time Steps): %d steps\n', num_steps);
fprintf('시간 간격 (Time Step Interval): %s ~ %s\n', unique_times{1}, unique_times{min(2, num_steps)});
fprintf('그리드 크기 (Grid Size per Step): %d points (%d Longitude x %d Latitude)\n', num_points, n_lon, n_lat);
fprintf('GSE Longitude Range: [%.1f, %.1f] (Step: %.1f deg)\n', min(gse_longs), max(gse_longs), diff(gse_longs(1:2)));
fprintf('GSE Latitude Range: [%.1f, %.1f] (Step: %.1f deg)\n', min(gse_latis), max(gse_latis), diff(gse_latis(1:2)));
fprintf('------------------------\n\n');

%% 3. Perform Analytical Verification
fprintf('Performing analytical coordinate transformation verification...\n');
tic;

% Initialize error arrays
max_errors = zeros(num_steps, 1);
mean_errors = zeros(num_steps, 1);
all_mean_lon_err = zeros(num_steps, 1);
all_mean_lat_err = zeros(num_steps, 1);

% Prepare coordinate containers for visualization
% We will sample a few representative time steps for plotting (e.g., start, 1/4, 1/2, 3/4)
plot_steps = round(linspace(1, num_steps, 4));
stored_grids = cell(length(plot_steps), 1);

for t = 1:num_steps
    t_str = unique_times{t};
    t_idx = strcmp(data.Time, t_str);
    
    file_gse_lon = data.GSE_long(t_idx);
    file_gse_lat = data.GSE_lati(t_idx);
    file_geo_lon = data.GEO_long(t_idx);
    file_geo_lat = data.GEO_lati(t_idx);
    
    % Compute our high-precision analytical values
    [calc_geo_lon, calc_geo_lat] = gse2geo_analytical(file_gse_lon, file_gse_lat, t_str);
    
    % Handle longitude wrapping (-180 to 180 degrees)
    lon_diff = mod(file_geo_lon - calc_geo_lon + 180, 360) - 180;
    lat_diff = file_geo_lat - calc_geo_lat;
    
    % Great circle distance error (true angle error in degrees)
    % dist = sqrt( (d_lon * cos(lat))^2 + (d_lat)^2 )
    dist_error = sqrt( (lon_diff .* cosd(file_geo_lat)).^2 + lat_diff.^2 );
    
    max_errors(t) = max(dist_error);
    mean_errors(t) = mean(dist_error);
    all_mean_lon_err(t) = mean(abs(lon_diff));
    all_mean_lat_err(t) = mean(abs(lat_diff));
    
    % Store a few grids for the 2D static overview plot
    p_idx = find(plot_steps == t);
    if ~isempty(p_idx)
        stored_grids{p_idx}.time = t_str;
        stored_grids{p_idx}.geo_lon = reshape(file_geo_lon, [n_lon, n_lat]);
        stored_grids{p_idx}.geo_lat = reshape(file_geo_lat, [n_lon, n_lat]);
    end
end
verification_time = toc;
fprintf('Verification completed in %.2f seconds.\n', verification_time);
fprintf('전체 평균 오차 (Overall Mean Great Circle Error): %.6f degrees\n', mean(mean_errors));
fprintf('최대 오차 (Overall Max Great Circle Error): %.6f degrees\n', max(max_errors));
fprintf('Longitude 평균 오차 (Mean Lon Absolute Error): %.6f degrees\n', mean(all_mean_lon_err));
fprintf('Latitude 평균 오차 (Mean Lat Absolute Error): %.6f degrees\n', mean(all_mean_lat_err));
fprintf('========================================================================\n\n');

%% 4. Graphical Visualizations

% Set custom aesthetics (Clean theme)
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

% Load standard coastlines for geospatial reference
try
    coast = load('coastlines.mat');
catch
    % Fallback if coastlines is not in the path (though it should be)
    coast.coastlon = [];
    coast.coastlat = [];
end

%% Figure 1: Static Overview of GSE Grid projected on GEO at 4 Time Steps
fig1 = figure('Name', 'GSE Grid projected on GEO Coordinates', 'Position', [100, 100, 1200, 800], 'Color', 'w');
sgtitle('GSE Grid Projected onto GEO Coordinates (시간에 따른 변화)', 'FontSize', 16, 'FontWeight', 'bold');

colors = cell(4,1);
colors{1} = [0, 0.4470, 0.7410];  % Blue
colors{2} = [0.8500, 0.3250, 0.0980]; % Orange
colors{3} = [0.9290, 0.6940, 0.1250]; % Yellow
colors{4} = [0.4940, 0.1840, 0.5560]; % Purple

for p = 1:4
    subplot(2, 2, p);
    hold on;
    grid on;
    box on;
    
    % Draw background world coastlines
    if ~isempty(coast.coastlon)
        plot(coast.coastlon, coast.coastlat, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    end
    
    grid_data_p = stored_grids{p};
    lon_m = grid_data_p.geo_lon;
    lat_m = grid_data_p.geo_lat;
    
    % Plot lines of constant GSE latitude
    for lat_idx = 1:size(lon_m, 2)
        plot_map_path(lon_m(:, lat_idx), lat_m(:, lat_idx), 'Color', [colors{p}, 0.25], 'LineWidth', 0.5);
    end
    
    % Plot lines of constant GSE longitude
    for lon_idx = 1:size(lon_m, 1)
        plot_map_path(lon_m(lon_idx, :), lat_m(lon_idx, :), 'Color', [colors{p}, 0.25], 'LineWidth', 0.5);
    end
    
    % Plot GSE Equator (GSE latitude closest to 0)
    eq_idx = round(size(lon_m, 2) / 2);
    plot_map_path(lon_m(:, eq_idx), lat_m(:, eq_idx), 'Color', colors{p}, 'LineWidth', 1.8, 'DisplayName', 'GSE Equator');
    
    % Plot GSE Meridian (GSE longitude closest to 0) in Orange
    [~, lon0_idx] = min(abs(gse_longs));
    plot_map_path(lon_m(lon0_idx, :), lat_m(lon0_idx, :), 'Color', [1, 0.4, 0], 'LineWidth', 1.8, 'DisplayName', 'GSE Meridian (Lon=0)');
    
    title(sprintf('Time: %s UTC', grid_data_p.time), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('GEO Longitude (deg)', 'FontSize', 10);
    ylabel('GEO Latitude (deg)', 'FontSize', 10);
    xlim([-180, 180]);
    ylim([-90, 90]);
    xticks(-180:60:180);
    yticks(-90:30:90);
end

%% Figure 2: Error Analysis & Time Series Comparison
fig2 = figure('Name', 'Coordinate Error Analysis', 'Position', [200, 150, 1100, 700], 'Color', 'w');
sgtitle('GEO2GSE 오차 분석 및 시계열 비교', 'FontSize', 16, 'FontWeight', 'bold');

% Subplot 1: Maximum & Mean Error over Time
subplot(2, 1, 1);
hold on;
grid on;
box on;
time_idx = 1:num_steps;
plot(time_idx, max_errors, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Max Great Circle Error');
plot(time_idx, mean_errors, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Mean Great Circle Error');
xlabel('Time Step Index (5-min interval)', 'FontSize', 11);
ylabel('Error (degrees)', 'FontSize', 11);
title('시간에 따른 분석 모델과의 각도 오차 변화', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10, 'AutoUpdate', 'off');
set(gca, 'FontSize', 10);

% Subplot 2: Single Point Trace Comparison
% Let's trace the GEO position corresponding to the center of the GSE grid (GSE.long = 0, GSE.lati = 0) over time.'
target_lon = 0;
target_lat = 0;
trace_lon_file = zeros(num_steps, 1);
trace_lat_file = zeros(num_steps, 1);
trace_lon_calc = zeros(num_steps, 1);
trace_lat_calc = zeros(num_steps, 1);

for t = 1:num_steps
    t_str = unique_times{t};
    t_idx = strcmp(data.Time, t_str);
    
    file_gse_lon = data.GSE_long(t_idx);
    file_gse_lat = data.GSE_lati(t_idx);
    file_geo_lon = data.GEO_long(t_idx);
    file_geo_lat = data.GEO_lati(t_idx);
    
    % Find point closest to GSE(0,0)
    [~, min_idx] = min(abs(file_gse_lon - target_lon) + abs(file_gse_lat - target_lat));
    
    trace_lon_file(t) = file_geo_lon(min_idx);
    trace_lat_file(t) = file_geo_lat(min_idx);
    
    [calc_lon, calc_lat] = gse2geo_analytical(target_lon, target_lat, t_str);
    trace_lon_calc(t) = calc_lon;
    trace_lat_calc(t) = calc_lat;
end

subplot(2, 2, 3);
hold on;
grid on;
box on;
plot(time_idx, trace_lon_file, 'ro', 'MarkerSize', 4, 'DisplayName', 'File (GEO)');
plot(time_idx, trace_lon_calc, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Analytical Model');
xlabel('Time Step Index', 'FontSize', 11);
ylabel('GEO Longitude (degrees)', 'FontSize', 11);
title('GSE(0,0) 지점의 GEO 경도 변화 비교', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'AutoUpdate', 'off');
set(gca, 'FontSize', 10);

subplot(2, 2, 4);
hold on;
grid on;
box on;
plot(time_idx, trace_lat_file, 'ro', 'MarkerSize', 4, 'DisplayName', 'File (GEO)');
plot(time_idx, trace_lat_calc, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Analytical Model');
xlabel('Time Step Index', 'FontSize', 11);
ylabel('GEO Latitude (degrees)', 'FontSize', 11);
title('GSE(0,0) 지점의 GEO 위도 변화 비교', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'AutoUpdate', 'off');
set(gca, 'FontSize', 10);

%% 5. Generating 2D Animation & Saving to GIF
fprintf('Generating interactive map animation and saving to GIF: %s...\n', gif_filename);

fig3 = figure('Name', 'Coordinate Transformation Animation', 'Position', [150, 100, 950, 600], 'Color', 'w');
hold on;
grid on;
box on;

% Draw background world coastlines
if ~isempty(coast.coastlon)
    plot(coast.coastlon, coast.coastlat, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
end

xlim([-180, 180]);
ylim([-90, 90]);
xticks(-180:60:180);
yticks(-90:30:90);
xlabel('GEO Longitude (deg)', 'FontSize', 11);
ylabel('GEO Latitude (deg)', 'FontSize', 11);
set(gca, 'FontSize', 10);

% Pre-allocate graphic handles for animation speedup
h_grid = [];
h_eq = [];
h_mer = [];
h_title = title('', 'FontSize', 14, 'FontWeight', 'bold');

% Loop to generate animation and write to GIF
% (To avoid slow execution, we sample every 4th step for the GIF)
animation_step = 4;
first_frame = true;

for t = 1:animation_step:num_steps
    t_str = unique_times{t};
    t_idx = strcmp(data.Time, t_str);
    
    file_geo_lon = data.GEO_long(t_idx);
    file_geo_lat = data.GEO_lati(t_idx);
    
    lon_m = reshape(file_geo_lon, [n_lon, n_lat]);
    lat_m = reshape(file_geo_lat, [n_lon, n_lat]);
    
    % Self-healing: Ensure fig3 is a valid figure handle
    if ~ishandle(fig3)
        fig3 = figure('Name', 'Coordinate Transformation Animation', 'Position', [150, 100, 950, 600], 'Color', 'w');
        hold on;
        grid on;
        box on;
        if ~isempty(coast.coastlon)
            plot(coast.coastlon, coast.coastlat, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        end
        xlim([-180, 180]);
        ylim([-90, 90]);
        xticks(-180:60:180);
        yticks(-90:30:90);
        xlabel('GEO Longitude (deg)', 'FontSize', 11);
        ylabel('GEO Latitude (deg)', 'FontSize', 11);
        set(gca, 'FontSize', 10);
        h_title = title('', 'FontSize', 14, 'FontWeight', 'bold');
        h_grid = [];
        h_eq = [];
        h_mer = [];
    else
        % Clear previous step graphics if they exist
        if ~isempty(h_grid)
            delete(h_grid);
        end
        if ~isempty(h_eq)
            delete(h_eq);
        end
        if ~isempty(h_mer)
            delete(h_mer);
        end
    end
    
    h_grid = [];
    % Draw constant GSE latitude lines
    for lat_idx = 1:size(lon_m, 2)
        h = plot_map_path_anim(lon_m(:, lat_idx), lat_m(:, lat_idx), 'Color', [0, 0.4470, 0.7410, 0.15], 'LineWidth', 0.5);
        h_grid = [h_grid; h];
    end
    
    % Draw constant GSE longitude lines
    for lon_idx = 1:size(lon_m, 1)
        h = plot_map_path_anim(lon_m(lon_idx, :), lat_m(lon_idx, :), 'Color', [0, 0.4470, 0.7410, 0.15], 'LineWidth', 0.5);
        h_grid = [h_grid; h];
    end
    
    % Draw GSE equator
    eq_idx = round(size(lon_m, 2) / 2);
    h_eq = plot_map_path_anim(lon_m(:, eq_idx), lat_m(:, eq_idx), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.8);
    
    % Draw GSE meridian (longitude closest to 0) in Orange
    [~, lon0_idx] = min(abs(gse_longs));
    h_mer = plot_map_path_anim(lon_m(lon0_idx, :), lat_m(lon0_idx, :), 'Color', [1, 0.4, 0], 'LineWidth', 1.8);
    
    set(h_title, 'String', sprintf('GSE Grid in GEO Coordinates - Time: %s UTC', t_str));
    drawnow;
    
    % Write to GIF
    frame = getframe(fig3);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if first_frame
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
        first_frame = false;
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end
fprintf('Animation GIF saved successfully as "%s"!\n\n', gif_filename);

%% 6. Final Report Summary
fprintf('==================== VERIFICATION SUMMARY ====================\n');
avg_err = mean(mean_errors);
max_err = max(max_errors);

if avg_err < 0.05
    fprintf('결과: 성공 (SUCCESS)!\n');
    fprintf('파일의 GSE-to-GEO 변환은 분석 모델과 거의 완벽하게 일치합니다.\n');
    fprintf('평균 Great Circle 각도 오차가 %.6f도(<0.05도)로 매우 작습니다.\n', avg_err);
    fprintf('이는 지구 자전 속도 및 태양 황경 계산이 올바르게 반영되었음을 증명합니다.\n');
elseif avg_err < 1.0
    fprintf('결과: 부분 성공 (PARTIAL SUCCESS - LOW ACCURACY MODEL)!\n');
    fprintf('파일의 GSE-to-GEO 변환이 분석 모델과 대체로 일치하지만, 평균 %.4f도의 미세한 오차가 존재합니다.\n', avg_err);
    fprintf('이는 분석 모델 간의 Greenwich Sidereal Time(GMST) 공식 차이(예: Apparent vs Mean Sidereal Time)나\n');
    fprintf('태양 황도 경도 계산 시 고차 섭동항의 생략 여부에 따른 차이일 가능성이 큽니다.\n');
else
    fprintf('결과: 검토 필요 (WARNING - SIGNIFICANT DISCREPANCY)!\n');
    fprintf('평균 오차가 %.4f도로, 분석 좌표계 변환 모델과 큰 차이가 발생했습니다.\n', avg_err);
    fprintf('시간에 따른 지구의 자전(Z축 회전)이나 태양-지구 황도면의 기울기(X축 회전) 계산 코드를 재검토하세요.\n');
end
fprintf('==============================================================\n');

%% Helper Functions

function [geo_lon_calc, geo_lat_calc] = gse2geo_analytical(gse_lon, gse_lat, t_str)
    % Perform high-precision analytical coordinate transformation from GSE to GEO
    % Input angles are in DEGREES.
    
    % Parse string to datetime
    dt = datetime(t_str, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    [yr, mo, dy, hr, mn, sc] = datevec(dt);
    
    % Manual JD calculation for maximum robustness and precision
    % Handles month vectors and Gregorian calendar transition details
    jd = zeros(size(yr));
    for i = 1:length(yr)
        y = yr(i); m = mo(i); d = dy(i);
        h = hr(i); min_val = mn(i); s = sc(i);
        if m <= 2
            y = y - 1;
            m = m + 12;
        end
        A = floor(y / 100);
        B = 2 - A + floor(A / 4);
        jd(i) = floor(365.25 * (y + 4716)) + floor(30.6001 * (m + 1)) + d + B - 1524.5 + (h + min_val/60 + s/3600)/24;
    end
    
    % T: Julian centuries since J2000.0
    T = (jd - 2451545.0) / 36525;
    
    % Mean longitude of the Sun (degrees)
    L = 280.460 + 36000.770 * T;
    
    % Mean anomaly of the Sun (degrees)
    M = 357.528 + 35999.050 * T;
    
    % Ecliptic longitude of the Sun (degrees)
    lambda_sun = L + 1.915 * sind(M) + 0.020 * sind(2*M);
    lambda_sun = mod(lambda_sun, 360);
    
    % Obliquity of the ecliptic (degrees)
    epsilon = 23.439 - 0.013 * T;
    
    % Greenwich Mean Sidereal Time (degrees)
    % IAU 1982 formula for GMST in degrees at any time:
    gmst = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T.^2 - T.^3 / 38710000;
    gmst = mod(gmst, 360);
    
    % Convert spherical GSE to Cartesian unit vectors
    cos_lat = cosd(gse_lat);
    x_gse = cos_lat .* cosd(gse_lon);
    y_gse = cos_lat .* sind(gse_lon);
    z_gse = sind(gse_lat);
    
    % Step 1: Rotate about Z by -lambda_sun
    cos_lam = cosd(lambda_sun);
    sin_lam = sind(lambda_sun);
    x1 = x_gse .* cos_lam - y_gse .* sin_lam;
    y1 = x_gse .* sin_lam + y_gse .* cos_lam;
    z1 = z_gse;
    
    % Step 2: Rotate about X by -epsilon
    cos_eps = cosd(epsilon);
    sin_eps = sind(epsilon);
    x2 = x1;
    y2 = y1 .* cos_eps - z1 .* sin_eps;
    z2 = y1 .* sin_eps + z1 .* cos_eps;
    
    % Step 3: Rotate about Z by +gmst
    cos_gmst = cosd(gmst);
    sin_gmst = sind(gmst);
    x_geo = x2 .* cos_gmst + y2 .* sin_gmst;
    y_geo = -x2 .* sin_gmst + y2 .* cos_gmst;
    z_geo = z2;
    
    % Convert Cartesian GEO back to spherical coordinates
    geo_lat_calc = asind(z_geo);
    geo_lon_calc = atan2d(y_geo, x_geo);
end

function plot_map_path(lon, lat, varargin)
    % Split path at -180/180 crossing to prevent visual stretch streaks
    diff_lon = diff(lon);
    jump_idx = find(abs(diff_lon) > 180);
    start_idx = 1;
    for k = 1:length(jump_idx)
        idx_range = start_idx:jump_idx(k);
        plot(lon(idx_range), lat(idx_range), varargin{:});
        start_idx = jump_idx(k) + 1;
    end
    plot(lon(start_idx:end), lat(start_idx:end), varargin{:});
end

function h = plot_map_path_anim(lon, lat, varargin)
    % Same as plot_map_path but returns handles of all plotted segments
    diff_lon = diff(lon);
    jump_idx = find(abs(diff_lon) > 180);
    start_idx = 1;
    h = [];
    for k = 1:length(jump_idx)
        idx_range = start_idx:jump_idx(k);
        h_part = plot(lon(idx_range), lat(idx_range), varargin{:});
        h = [h; h_part];
        start_idx = jump_idx(k) + 1;
    end
    h_part = plot(lon(start_idx:end), lat(start_idx:end), varargin{:});
    h = [h; h_part];
end

