%% plot_MATE_1D_nH_10Re.m
% =========================================================================
% Purpose: Read and compare MATE exospheric hydrogen density (nH) simulations
%          at 3 Re and 10 Re across DOY 164-174 (2008 storm event).
%
% Compared Runs:
%   1. With Charge Exchange (With CX):    MATE_nH_GRCPX1_test_2008*.data
%   2. Without Charge Exchange (No CX):   MATE_nH_GRCP_RCCX2_2008*.data
%
% Panels:
%   1. Top:    Geomagnetic Dst index [nT]
%   2. Middle: nH [cm^-3] at 3 Re (With CX vs No CX)
%   3. Bottom: nH [cm^-3] at 10 Re (With CX vs No CX)
% =========================================================================

clear; close all; clc;

%% 1. Simulation & Comparison Settings
year = 2008;
doys = 164:174; % 11 days (2008164 to 2008174)

runs = [ ...
    struct('id', 'With_CX',    'name', 'With CX (GRCPX1\_RCCX2)',   'prefix', 'MATE_nH_GRCPX1_RCCX2_', ...
           'color', [0.00, 0.45, 0.74], 'linestyle', '-',  'marker', 'none', 'linewidth', 2.0), ...
    struct('id', 'Without_CX', 'name', 'Without CX (GRCP\_RCCX2)', 'prefix', 'MATE_nH_GRCP_RCCX2_',  ...
           'color', [0.85, 0.33, 0.10], 'linestyle', '--', 'marker', 'none', 'linewidth', 2.0)  ...
];

% Candidate directories to search
candidate_dirs = { ...
    'C:\Users\slee122\OneDrive - NASA\Desktop\Work\git\MATE', ...
    '\\wsl.localhost\Ubuntu-22.04\home\sylee\exospherecode\MATE\output\0728', ...
    '\\wsl$\Ubuntu-22.04\home\sylee\exospherecode\MATE\output\0728', ...
    './' ...
};

%% 2. Load Data for Each Simulation Run
sim_data = struct();

for r = 1:length(runs)
    prefix = runs(r).prefix;
    fprintf('\n--- Processing Run: %s (%s) ---\n', runs(r).name, prefix);
    
    % Locate data directory for this prefix
    data_dir = '';
    for k = 1:length(candidate_dirs)
        cur_dir = candidate_dirs{k};
        if isempty(cur_dir) || ~exist(cur_dir, 'dir'), continue; end
        test_file = fullfile(cur_dir, sprintf('%s%04d%03d.data', prefix, year, doys(1)));
        if exist(test_file, 'file')
            data_dir = cur_dir;
            break;
        end
    end
    
    % Fallback if prefix was previously named _test_
    if isempty(data_dir) && strcmp(prefix, 'MATE_nH_GRCPX1_RCCX2_')
        fallback_prefix = 'MATE_nH_GRCPX1_test_';
        for k = 1:length(candidate_dirs)
            cur_dir = candidate_dirs{k};
            if isempty(cur_dir) || ~exist(cur_dir, 'dir'), continue; end
            test_file = fullfile(cur_dir, sprintf('%s%04d%03d.data', fallback_prefix, year, doys(1)));
            if exist(test_file, 'file')
                data_dir = cur_dir;
                prefix = fallback_prefix;
                fprintf('  (Using fallback prefix: %s)\n', prefix);
                break;
            end
        end
    end
    
    if isempty(data_dir)
        error('Could not locate files for prefix "%s" in candidate directories.', prefix);
    end
    fprintf('  Data directory: %s\n', data_dir);
    
    all_t   = [];
    all_3R  = [];
    all_10R = [];
    
    for idoy = doys
        fn = sprintf('%s%04d%03d.data', prefix, year, idoy);
        fp = fullfile(data_dir, fn);
        
        if ~exist(fp, 'file')
            warning('File not found: %s', fn);
            continue;
        end
        
        fid = fopen(fp, 'rb');
        raw = fread(fid, 'float32');
        fclose(fid);
        
        n_floats = length(raw);
        
        if n_floats == 1086912 % 17 * 72 * 37 * 24
            nR = 17; nLon = 72; nLat = 37; ntpd = 24;
            iLat_eq = fix(nLat / 2) + 1; % 19
            data4D = reshape(raw, [nR, nLon, nLat, ntpd]);
            
            % Index 3 is 3.0 Re, Index 17 is 10.0 Re
            nH_3Re  = squeeze(data4D(3,  1, iLat_eq, :));
            nH_10Re = squeeze(data4D(17, 1, iLat_eq, :));
            t_vec   = idoy + (0:ntpd-1)' / 24.0;
            
        elseif n_floats == (17 * 24) % Pure 1D [17, 24]
            nR = 17; ntpd = 24;
            data2D  = reshape(raw, [nR, ntpd]);
            nH_3Re  = data2D(3,  :)';
            nH_10Re = data2D(17, :)';
            t_vec   = idoy + (0:ntpd-1)' / 24.0;
            
        elseif mod(n_floats, 17) == 0
            nR = 17;
            cur_nt  = n_floats / nR;
            data2D  = reshape(raw, [nR, cur_nt]);
            nH_3Re  = data2D(3,  :)';
            nH_10Re = data2D(17, :)';
            t_vec   = idoy + (0:cur_nt-1)' / 24.0;
        else
            warning('Unrecognized length (%d) in %s. Skipping.', n_floats, fn);
            continue;
        end
        
        all_t   = [all_t;   t_vec(:)];
        all_3R  = [all_3R;  nH_3Re(:)];
        all_10R = [all_10R; nH_10Re(:)];
    end
    
    % Exclude 1st time step (initial startup transient)
    if length(all_t) > 1
        all_t   = all_t(2:end);
        all_3R  = all_3R(2:end);
        all_10R = all_10R(2:end);
        fprintf('  Excluded 1st time step (t=0) as initial transient.\n');
    end
    
    sim_data(r).time    = all_t;
    sim_data(r).nH_3Re  = all_3R;
    sim_data(r).nH_10Re = all_10R;
    sim_data(r).name    = runs(r).name;
    sim_data(r).color   = runs(r).color;
    sim_data(r).ls      = runs(r).linestyle;
    sim_data(r).lw      = runs(r).linewidth;
    sim_data(r).marker  = runs(r).marker;
    
    fprintf('  Total time steps: %d (DOY %.2f to %.2f)\n', length(all_t), all_t(1), all_t(end));
    fprintf('   3 Re nH: Min = %.2f, Max = %.2f, Mean = %.2f cm^-3\n', ...
        min(all_3R), max(all_3R), mean(all_3R));
    fprintf('  10 Re nH: Min = %.2f, Max = %.2f, Mean = %.2f cm^-3\n', ...
        min(all_10R), max(all_10R), mean(all_10R));
end

%% 3. Quantitative Comparison Summary
if length(sim_data) >= 2
    fprintf('\n==========================================================\n');
    fprintf('           CHARGE EXCHANGE (CX) IMPACT COMPARISON          \n');
    fprintf('==========================================================\n');
    
    mean_3R_with   = mean(sim_data(1).nH_3Re);
    mean_3R_without = mean(sim_data(2).nH_3Re);
    depletion_3R   = (1 - mean_3R_with / mean_3R_without) * 100;
    
    mean_10R_with   = mean(sim_data(1).nH_10Re);
    mean_10R_without = mean(sim_data(2).nH_10Re);
    diff_10R        = (mean_10R_with - mean_10R_without) / mean_10R_without * 100;
    
    fprintf('  3 Re Mean Density:\n');
    fprintf('    Without CX: %7.2f cm^-3\n', mean_3R_without);
    fprintf('    With CX:    %7.2f cm^-3\n', mean_3R_with);
    fprintf('    Depletion:  %7.2f %% (CX loss effect)\n', depletion_3R);
    
    fprintf(' 10 Re Mean Density:\n');
    fprintf('    Without CX: %7.2f cm^-3\n', mean_10R_without);
    fprintf('    With CX:    %7.2f cm^-3\n', mean_10R_with);
    fprintf('    Difference: %+7.2f %%\n', diff_10R);
    fprintf('==========================================================\n');
end

%% 4. Load Geomagnetic Dst Index
has_dst  = false;
dst_time = [];
dst_val  = [];

dst_candidates = { ...
    fullfile(pwd, '2008164_Dst.txt'), ...
    'C:\Users\slee122\OneDrive - NASA\Documents\MATLAB\2008164_Dst.txt', ...
    'C:\Users\slee122\OneDrive - NASA\Desktop\Work\git\MATE\2008164_Dst.txt', ...
    '\\wsl.localhost\Ubuntu-22.04\home\sylee\exospherecode\MATE\output\0728\2008164_Dst.txt' ...
};

for k = 1:length(dst_candidates)
    if exist(dst_candidates{k}, 'file')
        fid_dst = fopen(dst_candidates{k}, 'r');
        dst_raw = textscan(fid_dst, '%f %f %f %f', 'HeaderLines', 1, 'MultipleDelimsAsOne', true);
        fclose(fid_dst);
        dst_time = dst_raw{2} + dst_raw{3} / 24.0;
        dst_val  = dst_raw{4};
        has_dst  = true;
        fprintf('Loaded Dst index from: %s\n', dst_candidates{k});
        break;
    end
end

%% 5. Plotting (3-Panel Comparison Layout)
fig = figure('Name', 'MATE 1D nH CX Comparison (3 Re & 10 Re)', ...
             'Color', 'w', 'Position', [120, 60, 1000, 780]);

clr_dst = [0.15, 0.15, 0.15]; % Charcoal Black for neutral Dst
x_limits = [floor(min(sim_data(1).time)), max(sim_data(1).time)];

if has_dst
    tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Panel 1 (Top): Dst Index
    ax1 = nexttile;
    plot(ax1, dst_time, dst_val, 'LineWidth', 2.0, 'Color', clr_dst);
    grid on; box on;
    xlim(ax1, x_limits);
    ylabel(ax1, 'Dst [nT]', 'FontSize', 12, 'FontWeight', 'bold');
    title(ax1, sprintf('Geomagnetic Activity (Dst) & MATE 1D Exospheric H Density: CX Effect (%d)', year), ...
          'FontSize', 13, 'FontWeight', 'bold');
    yline(ax1, 0, '--k', 'Alpha', 0.4);
    set(ax1, 'FontSize', 11, 'LineWidth', 1.0);
    
    % Panel 2 (Middle): nH at 3 Re Comparison
    ax2 = nexttile;
    hold(ax2, 'on');
    h_plots_3R = gobjects(length(sim_data), 1);
    for r = 1:length(sim_data)
        h_plots_3R(r) = plot(ax2, sim_data(r).time, sim_data(r).nH_3Re, ...
            'LineStyle', sim_data(r).ls, 'LineWidth', sim_data(r).lw, ...
            'Color', sim_data(r).color, 'DisplayName', sim_data(r).name);
    end
    hold(ax2, 'off');
    grid on; box on;
    xlim(ax2, x_limits);
    ylabel(ax2, 'n_H [cm^{-3}] at 3 R_E', 'FontSize', 12, 'FontWeight', 'bold');
    lgd2 = legend(ax2, h_plots_3R, 'Location', 'northwest', 'Box', 'off', 'FontSize', 11);
    set(ax2, 'FontSize', 11, 'LineWidth', 1.0);
    
    % Panel 3 (Bottom): nH at 10 Re Comparison
    ax3 = nexttile;
    hold(ax3, 'on');
    h_plots_10R = gobjects(length(sim_data), 1);
    for r = 1:length(sim_data)
        h_plots_10R(r) = plot(ax3, sim_data(r).time, sim_data(r).nH_10Re, ...
            'LineStyle', sim_data(r).ls, 'LineWidth', sim_data(r).lw, ...
            'Color', sim_data(r).color, 'DisplayName', sim_data(r).name);
    end
    hold(ax3, 'off');
    grid on; box on;
    xlim(ax3, x_limits);
    xlabel(ax3, sprintf('Day of Year (DOY in %d)', year), 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax3, 'n_H [cm^{-3}] at 10 R_E', 'FontSize', 12, 'FontWeight', 'bold');
    lgd3 = legend(ax3, h_plots_10R, 'Location', 'northwest', 'Box', 'off', 'FontSize', 11);
    set(ax3, 'FontSize', 11, 'LineWidth', 1.0);
    
    linkaxes([ax1, ax2, ax3], 'x');
else
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Panel 1 (Top): nH at 3 Re Comparison
    ax1 = nexttile;
    hold(ax1, 'on');
    h_plots_3R = gobjects(length(sim_data), 1);
    for r = 1:length(sim_data)
        h_plots_3R(r) = plot(ax1, sim_data(r).time, sim_data(r).nH_3Re, ...
            'LineStyle', sim_data(r).ls, 'LineWidth', sim_data(r).lw, ...
            'Color', sim_data(r).color, 'DisplayName', sim_data(r).name);
    end
    hold(ax1, 'off');
    grid on; box on;
    xlim(ax1, x_limits);
    ylabel(ax1, 'n_H [cm^{-3}] at 3 R_E', 'FontSize', 12, 'FontWeight', 'bold');
    title(ax1, sprintf('MATE 1D Exospheric H Density: CX Effect (3 R_E & 10 R_E, %d)', year), ...
          'FontSize', 13, 'FontWeight', 'bold');
    legend(ax1, h_plots_3R, 'Location', 'northwest', 'Box', 'off', 'FontSize', 11);
    set(ax1, 'FontSize', 11, 'LineWidth', 1.0);
    
    % Panel 2 (Bottom): nH at 10 Re Comparison
    ax2 = nexttile;
    hold(ax2, 'on');
    h_plots_10R = gobjects(length(sim_data), 1);
    for r = 1:length(sim_data)
        h_plots_10R(r) = plot(ax2, sim_data(r).time, sim_data(r).nH_10Re, ...
            'LineStyle', sim_data(r).ls, 'LineWidth', sim_data(r).lw, ...
            'Color', sim_data(r).color, 'DisplayName', sim_data(r).name);
    end
    hold(ax2, 'off');
    grid on; box on;
    xlim(ax2, x_limits);
    xlabel(ax2, sprintf('Day of Year (DOY in %d)', year), 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax2, 'n_H [cm^{-3}] at 10 R_E', 'FontSize', 12, 'FontWeight', 'bold');
    legend(ax2, h_plots_10R, 'Location', 'northwest', 'Box', 'off', 'FontSize', 11);
    set(ax2, 'FontSize', 11, 'LineWidth', 1.0);
    
    linkaxes([ax1, ax2], 'x');
end

%% 6. Save Figure
output_png = 'MATE_1D_nH_CX_comparison_3Re_10Re.png';
exportgraphics(fig, output_png, 'Resolution', 300);
fprintf('\nFigure saved to: %s\n', fullfile(pwd, output_png));

% Also update default timeseries filename for convenience
exportgraphics(fig, 'MATE_1D_nH_3Re_10Re_timeseries.png', 'Resolution', 300);
