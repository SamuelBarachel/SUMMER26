%% ============================================================================
%  analyze_TOPAS_sim.m  —  UITF ND3 Target Simulation Results Analysis (MATLAB)
%% ============================================================================
%  Place this file in the same directory as your TOPAS CSV output files and run:
%
%      analyze_TOPAS_sim
%
%  No arguments needed. All outputs are written to the same directory.
%
%  This trimmed version produces only:
%    fig01_rz_heatmap.png    - 2D Dose map (R-Z plane, phi-averaged)
%    fig02_rz_contours.png   - 2D Dose map (R-Z plane) with isodose contours
%    fig03_radial_profile.png - radial dose profile (actual line graph)
%    fig04_axial_profile.png  - axial dose profile (actual line graph)
%
%  Each profile graph (03/04) reports the Non-Uniformity Index
%  NUI = sigma / xbar  (std dev / mean) computed over the plotted profile.
%
%  EXPECTED INPUT FILES (3-D voxel files, 25 R x 72 Phi x 50 Z = 90000 rows,
%  cols: r,phi,z,value):
%        Dose_ND3_Binned3D.csv   (preferred primary dose scorer)
%        Dose_RZ.csv / Dose_R.csv / Dose_Z.csv  (alternates, same format)
%% ============================================================================
clear; clc; close all;
HERE = fileparts(mfilename('fullpath'));
if isempty(HERE)
    HERE = pwd;   % called from command window
end
fprintf('\n%s\n', repmat('=',1,60));
fprintf('  TOPAS UITF ND3  —  ANALYSIS\n');
fprintf('  Directory: %s\n', HERE);
fprintf('%s\n\n', repmat('=',1,60));
%% ============================================================================
%  1.  GEOMETRY CONSTANTS  (from UITF_ND3.txt)
%% ============================================================================
ND3_RMIN_CM   = 1.25;        % cm   inner radius
ND3_RMAX_CM   = 2.50;        % cm   outer radius
ND3_HL_CM     = 1.25;        % cm   half-length  (full length 2.50 cm)
ND3_RASTER_X  = 2.5;         % cm   beam raster full width
ND3_RASTER_Y  = 2.0;         % cm   beam raster full height
%% ============================================================================
%  2.  COLOUR PALETTE  (dark theme)
%% ============================================================================
BG      = [0.043 0.059 0.102];
PANEL   = [0.075 0.098 0.161];
BORDER  = [0.165 0.208 0.314];
TEXT_C  = [0.910 0.929 0.961];
DIM     = [0.353 0.416 0.541];
CYAN    = [0.000 0.898 1.000];
AMBER   = [1.000 0.702 0.278];
GREEN_C = [0.298 0.686 0.314];
RED_C   = [1.000 0.420 0.420];
BLUE_C  = [0.259 0.647 0.961];
%% ============================================================================
%  3.  FILE LOOKUP TABLE  (3-D dose only)
%% ============================================================================
DOSE_3D_NAMES = {'Dose_ND3_Binned3D.csv', 'Dose_RZ.csv', 'Dose_R.csv', 'Dose_Z.csv', ...
                 'DoseAtSample.csv', 'DoseND3Binned3D.csv'};
%% ============================================================================
%  4.  LOAD PRIMARY 3-D DOSE
%% ============================================================================
d3 = []; hd = struct(); r_ax = []; phi_ax = []; z_ax = [];
p = first_existing(DOSE_3D_NAMES, HERE);
if isempty(p)
    error('No 3-D dose file found (looked for: %s)', strjoin(DOSE_3D_NAMES, ', '));
end
fprintf('[dose 3D]  %s\n', p);
[hd, data] = read_topas_csv(p);
if isempty(data) || size(data,2) ~= 4
    error('Could not parse 3-D dose data from %s', p);
end
[d3, r_ax, phi_ax, z_ax] = to_3d(data, hd);
fprintf('\n');
%% ============================================================================
%  5.  FIGURES
%% ============================================================================
fprintf('Generating figures...\n');
TABGROUP = make_tab_group(BG, PANEL, TEXT_C, BORDER);

fig01_rz_heatmap(d3, r_ax, z_ax, hd, TABGROUP, HERE, BG, PANEL, BORDER, TEXT_C, DIM, AMBER, CYAN, ...
    ND3_RMIN_CM, ND3_RMAX_CM, ND3_HL_CM, ND3_RASTER_X, ND3_RASTER_Y);

fig02_rz_contours(d3, r_ax, z_ax, hd, TABGROUP, HERE, BG, PANEL, BORDER, TEXT_C, DIM, AMBER, CYAN, BLUE_C, GREEN_C, RED_C, ...
    ND3_RMIN_CM, ND3_RMAX_CM, ND3_HL_CM);

fig03_radial(d3, r_ax, hd, TABGROUP, HERE, BG, PANEL, BORDER, TEXT_C, DIM, AMBER, CYAN, GREEN_C, ...
    ND3_RMIN_CM, ND3_RMAX_CM);

fig04_axial(d3, z_ax, hd, TABGROUP, HERE, BG, PANEL, BORDER, TEXT_C, DIM, AMBER, GREEN_C, ND3_HL_CM);

fprintf('\n%s\n', repmat('=',1,60));
fprintf('  Done.  4 figures loaded into layout panel.\n');
fprintf('  All outputs saved to -> %s\n', HERE);
fprintf('%s\n\n', repmat('=',1,60));

%% ============================================================================
%  LOCAL FUNCTIONS
%% ============================================================================
function p = first_existing(names, base)
    p = [];
    for k = 1:numel(names)
        candidate = fullfile(base, names{k});
        if exist(candidate, 'file')
            p = candidate;
            return;
        end
    end
end

function h = parse_header(filepath)
    h.version   = ''; h.param_file = ''; h.scorer = ''; h.component = '';
    h.quantity  = ''; h.unit = '';       h.stat = '';
    h.r_bins    = []; h.r_size = [];
    h.phi_bins  = []; h.phi_size = [];
    h.z_bins    = []; h.z_size  = [];
    h.warnings  = {};
    h.n_hdr     = 0;
    fid = fopen(filepath, 'r');
    if fid < 0
        warning('parse_header: cannot open %s', filepath);
        return;
    end
    while true
        raw = fgetl(fid);
        if ~ischar(raw); break; end
        raw = strtrim(raw);
        if isempty(raw) || raw(1) ~= '#'; break; end
        h.n_hdr = h.n_hdr + 1;
        text = strtrim(raw(2:end));
        lo   = lower(text);
        if contains(lo,'topas version')
            parts = strsplit(text,':'); if numel(parts)>1; h.version = strtrim(parts{2}); end
        elseif contains(lo,'parameter file')
            parts = strsplit(text,':'); if numel(parts)>1; h.param_file = strtrim(parts{2}); end
        elseif contains(lo,'results for scorer')
            parts = strsplit(text,':'); if numel(parts)>1; h.scorer = strtrim(parts{2}); end
        elseif contains(lo,'scored in component')
            parts = strsplit(text,':'); if numel(parts)>1; h.component = strtrim(parts{2}); end
        elseif contains(lo,'warning')
            h.warnings{end+1} = text; 
        elseif startsWith(lo,'r in')
            tok = regexp(text, '([\d.]+)\s*bins?\s*of\s*([\d.eE+-]+)', 'tokens', 'once');
            if ~isempty(tok)
                h.r_bins = round(str2double(tok{1}));
                h.r_size = str2double(tok{2});
            end
        elseif startsWith(lo,'phi in')
            tok = regexp(text, '([\d.]+)\s*bins?\s*of\s*([\d.eE+-]+)', 'tokens', 'once');
            if ~isempty(tok)
                h.phi_bins = round(str2double(tok{1}));
                h.phi_size = str2double(tok{2});
            end
        elseif startsWith(lo,'z in')
            tok = regexp(text, '([\d.]+)\s*bins?\s*of\s*([\d.eE+-]+)', 'tokens', 'once');
            if ~isempty(tok)
                h.z_bins = round(str2double(tok{1}));
                h.z_size = str2double(tok{2});
            end
        elseif contains(lo,'scored quantity')
            parts = strsplit(text,':'); if numel(parts)>1; h.quantity = strtrim(parts{2}); end
        elseif contains(lo,'in units of') || contains(lo,'units of')
            tok = regexp(text, 'units of\s*(\S+)', 'tokens', 'once');
            if ~isempty(tok); h.unit = tok{1}; end
        end
    end
    fclose(fid);
    if isempty(h.quantity)
        h.quantity = h.scorer;
        if isempty(h.quantity); h.quantity = 'Dose'; end
    end
    if isempty(h.unit); h.unit = 'Gy'; end
    if isempty(h.r_bins);   h.r_bins   = 25; h.r_size   = 0.05; end
    if isempty(h.phi_bins); h.phi_bins = 72; h.phi_size = 5.0;  end
    if isempty(h.z_bins);   h.z_bins   = 50; h.z_size   = 0.05; end
end

function [h, data] = read_topas_csv(filepath)
    h    = parse_header(filepath);
    data = [];
    fid = fopen(filepath, 'r');
    if fid < 0
        warning('read_topas_csv: cannot open %s', filepath);
        return;
    end
    for k = 1:h.n_hdr
        fgetl(fid);
    end
    rows = {};
    while true
        raw = fgetl(fid);
        if ~ischar(raw); break; end
        raw = strtrim(raw);
        while ~isempty(raw) && (raw(end) == ',' || raw(end) == ' ')
            raw(end) = [];
        end
        if isempty(raw) || raw(1) == '#'; continue; end
        parts = strsplit(raw, ',');
        nums = str2double(parts);
        if all(~isnan(nums))
            rows{end+1} = nums; 
        end
    end
    fclose(fid);
    if isempty(rows); return; end
    ncols = max(cellfun(@numel, rows));
    mat   = NaN(numel(rows), ncols);
    for k = 1:numel(rows)
        mat(k, 1:numel(rows{k})) = rows{k};
    end
    data = mat;
    if size(data,2) == 1
        data = data(:,1);
    end
end

function [arr3d, r, phi, z] = to_3d(data, h)
    R   = h.r_bins;   Phi = h.phi_bins;  Z = h.z_bins;
    arr3d = zeros(R, Phi, Z);
    ri   = min(max(round(data(:,1)), 0), R-1);
    phii = min(max(round(data(:,2)), 0), Phi-1);
    zi   = min(max(round(data(:,3)), 0), Z-1);
    vals = data(:,4);
    lin = sub2ind([R, Phi, Z], ri+1, phii+1, zi+1);
    for k = 1:numel(lin)
        arr3d(lin(k)) = arr3d(lin(k)) + vals(k);
    end
    r   = (0:R-1)'   * h.r_size   + h.r_size/2;
    phi = (0:Phi-1)' * h.phi_size + h.phi_size/2;
    z   = (0:Z-1)'   * h.z_size   + h.z_size/2;
end

function out = rz_mean(arr3d)
    out = mean(arr3d, 2);
    out = squeeze(out); 
end

function out = rad_prof(arr3d)
    tmp = mean(arr3d, 3);
    out = mean(tmp, 2);
end

function out = axl_prof(arr3d)
    tmp = mean(arr3d, 1);
    tmp = mean(tmp, 2);
    out = squeeze(tmp);
end

function nui = compute_nui(profile)
    profile = profile(:);
    xbar = mean(profile);
    if xbar == 0
        nui = 0;
    else
        nui = std(profile) / xbar;
    end
end

function out = gauss_smooth2(data2d, sigma)
    if sigma <= 0; out = data2d; return; end
    hw   = ceil(3*sigma);
    x    = -hw:hw;
    kern = exp(-x.^2 / (2*sigma^2));
    kern = kern / sum(kern);
    out  = conv2(kern, kern, data2d, 'same');
end

% ── Tab Framework Helpers ───────────────────────────────────────────────────
function tg = make_tab_group(BG, PANEL, TEXT_C, BORDER)
    fig = figure('Color', BG, 'Units','normalized', 'Position',[0.08 0.08 0.84 0.84], ...
                  'Name','TOPAS UITF ND3 -- Dose Analysis', 'NumberTitle','off');
    tg = uitabgroup('Parent', fig);
end

function ax = make_tab_axes(tg, tab_title, w, h, bg, panel)
    tab = uitab('Parent', tg, 'Title', tab_title, 'BackgroundColor', bg);
    ax  = axes('Parent', tab, 'Units','normalized', 'Position',[0.09 0.12 0.85 0.78]);
    ax.UserData.aspect = [w h];
    ax.Color = panel;
end

function apply_dark_axes(ax, bg, panel, border, text_c, dim)
    set(ax, 'Color',panel, 'XColor',dim, 'YColor',dim, ...
            'GridColor',[0.118 0.165 0.259], 'GridAlpha',1, ...
            'FontName','Courier New', 'FontSize',9);
    set(ax.Title, 'Color',text_c, 'FontSize',11, 'FontWeight','bold');
    set(ax.XLabel,'Color',text_c,'FontSize',10);
    set(ax.YLabel,'Color',text_c,'FontSize',10);
    grid(ax,'on');
    ax.GridLineStyle = '--';
end

function h_cb = add_colorbar(ax, label, text_c, dim)
    h_cb = colorbar(ax);
    h_cb.Label.String = label;
    h_cb.Label.Color  = text_c;
    h_cb.Color = dim;
    h_cb.FontName = 'Courier New';
end

function add_nui_text(ax, nui, text_c)
    xl = xlim(ax); yl = ylim(ax);
    text(ax, xl(1) + 0.02*diff(xl), yl(2) - 0.05*diff(yl), ...
        sprintf('NUI = \\sigma/\\bar{x} = %.4f', nui), ...
        'Color', text_c, 'FontSize', 9, 'FontName','Courier New', ...
        'VerticalAlignment','top', 'BackgroundColor',[0 0 0 0.5], ...
        'Margin', 4);
end

function add_stamp(ax, dim)
    xl = xlim(ax); yl = ylim(ax);
    text(ax, xl(2) - 0.02*diff(xl), yl(1) + 0.04*diff(yl), ...
        'UITF ND3 Target Sim', 'Color', dim, 'FontSize', 8, ...
        'FontName', 'Courier New', 'HorizontalAlignment', 'right');
end

function save_tab_fig(ax, filename, folder)
    exportgraphics(ax, fullfile(folder, filename), 'Resolution', 150);
end

% ── FIGURE 01: R-Z 2D dose map (heatmap) ─────────────────────────────────────
function fig01_rz_heatmap(d3, r_ax, z_ax, hd, TABGROUP, HERE, BG, PANEL, BORDER, TEXT_C, DIM, AMBER, CYAN, ...
    ND3_RMIN_CM, ND3_RMAX_CM, ND3_HL_CM, ND3_RASTER_X, ND3_RASTER_Y)
    rz   = gauss_smooth2(rz_mean(d3), 0.8);
    nui  = compute_nui(rz(:));
    
    ax = make_tab_axes(TABGROUP, '2D Heatmap', 10, 6, BG, PANEL);
    pcolor(ax, z_ax(:)', r_ax(:), rz);
    shading(ax,'flat');
    colormap(ax, 'hot');
    apply_dark_axes(ax, BG, PANEL, BORDER, TEXT_C, DIM);
    add_colorbar(ax, sprintf('%s (%s)', hd.quantity, hd.unit), TEXT_C, DIM);
    hold(ax,'on');
    
    for r_b = [ND3_RMIN_CM, ND3_RMAX_CM]
        plot(ax, [z_ax(1) z_ax(end)], [r_b r_b], '--', 'Color', AMBER, 'LineWidth', 0.9);
    end
    z_full = 2*ND3_HL_CM;
    z_edge = (z_full - ND3_RASTER_X)/2;
    for ze = [max(z_edge,z_ax(1)), min(z_full-z_edge, z_ax(end))]
        plot(ax, [ze ze], [r_ax(1) r_ax(end)], ':', 'Color', CYAN, 'LineWidth', 0.7);
    end
    xlabel(ax,'Z (cm)'); ylabel(ax,'R (cm)');
    title(ax,'2D Dose Map -- R-Z Plane (phi-averaged, Gaussian smooth \sigma=0.8)');
    add_nui_text(ax, nui, TEXT_C);
    add_stamp(ax, DIM);
    save_tab_fig(ax, 'fig01_rz_heatmap.png', HERE);
end

% ── FIGURE 02: R-Z 2D dose map with isodose contours ─────────────────────────
function fig02_rz_contours(d3, r_ax, z_ax, hd, TABGROUP, HERE, BG, PANEL, BORDER, TEXT_C, DIM, AMBER, CYAN, BLUE_C, GREEN_C, RED_C, ...
    ND3_RMIN_CM, ND3_RMAX_CM, ND3_HL_CM)
    rz     = gauss_smooth2(rz_mean(d3), 1.0);
    nui    = compute_nui(rz(:));
    peak   = max(rz(:));
    r_prof = mean(rz, 2);   
    iso_fracs  = [0.20, 0.50, 0.80, 0.90, 0.95];
    iso_colors = {BLUE_C; CYAN; GREEN_C; AMBER; RED_C};
    
    ax = make_tab_axes(TABGROUP, 'Isodose Contours', 11, 6, BG, PANEL);
    pcolor(ax, z_ax(:)', r_ax(:), rz);
    shading(ax,'flat');
    hold(ax,'on');
    colormap(ax,'hot');
    apply_dark_axes(ax, BG, PANEL, BORDER, TEXT_C, DIM);
    add_colorbar(ax, sprintf('%s (%s)', hd.quantity, hd.unit), TEXT_C, DIM);
    
    legend_entries = {};
    for k = 1:numel(iso_fracs)
        frac = iso_fracs(k);
        threshold = peak * frac;
        idxs = find(r_prof >= threshold);
        if isempty(idxs); continue; end
        r_iso = r_ax(idxs(end));
        plot(ax, [z_ax(1) z_ax(end)], [r_iso r_iso], '-', ...
            'Color', iso_colors{k}, 'LineWidth', 1.5);
        legend_entries{end+1} = sprintf('%d%%  (R=%.3f cm)', round(frac*100), r_iso); 
    end
    for r_b = [ND3_RMIN_CM, ND3_RMAX_CM]
        plot(ax, [z_ax(1) z_ax(end)], [r_b r_b], '--', 'Color', TEXT_C, 'LineWidth', 1.0, 'HandleVisibility','off');
    end
    if ~isempty(legend_entries)
        legend(ax, legend_entries, 'Location','northeast','TextColor',TEXT_C, ...
               'EdgeColor',BORDER,'Color',PANEL,'FontSize',8,'FontName','Courier New');
    end
    xlabel(ax,'Z (cm)'); ylabel(ax,'R (cm)');
    title(ax, sprintf('2D Dose Map with Isodose Contours  |  Peak: %.4e %s', peak, hd.unit));
    add_nui_text(ax, nui, TEXT_C);
    add_stamp(ax, DIM);
    save_tab_fig(ax, 'fig02_rz_contours.png', HERE);
end

% ── FIGURE 03: radial profile (graph) ────────────────────────────────────────
function fig03_radial(d3, r_ax, hd, TABGROUP, HERE, BG, PANEL, BORDER, TEXT_C, DIM, AMBER, CYAN, GREEN_C, ...
    ND3_RMIN_CM, ND3_RMAX_CM)
    prof = rad_prof(d3);
    nui  = compute_nui(prof);
    
    ax = make_tab_axes(TABGROUP, 'Radial Profile', 9, 5, BG, PANEL);
    hold(ax,'on');
    fill(ax, [r_ax; flipud(r_ax)]', [prof; zeros(size(prof))]', CYAN, 'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');
    plot(ax, r_ax, prof, 'Color', CYAN, 'LineWidth', 2, 'DisplayName','Radial profile (phi-Z avg)');
    
    xpatch = [ND3_RMIN_CM ND3_RMAX_CM ND3_RMAX_CM ND3_RMIN_CM];
    yl = [0 0 max(prof)*1.1 max(prof)*1.1];
    fill(ax, xpatch, yl, GREEN_C, 'FaceAlpha',0.07,'EdgeColor','none','HandleVisibility','off');
    plot(ax,[ND3_RMIN_CM ND3_RMIN_CM],[0 max(prof)],'--','Color',GREEN_C,'LineWidth',0.9,'HandleVisibility','off');
    plot(ax,[ND3_RMAX_CM ND3_RMAX_CM],[0 max(prof)],'--','Color',GREEN_C,'LineWidth',0.9,'HandleVisibility','off');
    
    [~, ip] = max(prof);
    text(ax, r_ax(ip)+0.05, prof(ip)*1.05, sprintf('%.3e', prof(ip)), ...
        'Color',AMBER,'FontSize',8,'FontName','Courier New');
    legend(ax, 'Location','best','TextColor',TEXT_C,'EdgeColor',BORDER,'Color',PANEL,'FontSize',8,'FontName','Courier New');
    apply_dark_axes(ax, BG, PANEL, BORDER, TEXT_C, DIM);
    xlabel(ax,'R (cm)'); ylabel(ax, sprintf('%s (%s)', hd.quantity, hd.unit));
    title(ax, sprintf('Radial Dose Profile  (averaged over all phi and Z)  |  NUI = %.4f', nui));
    add_stamp(ax, DIM);
    save_tab_fig(ax, 'fig03_radial_profile.png', HERE);
end

% ── FIGURE 04: axial profile (graph) ─────────────────────────────────────────
function fig04_axial(d3, z_ax, hd, TABGROUP, HERE, BG, PANEL, BORDER, TEXT_C, DIM, AMBER, GREEN_C, ND3_HL_CM) %#ok<INUSD>
    prof = axl_prof(d3);
    nui  = compute_nui(prof);
    
    ax = make_tab_axes(TABGROUP, 'Axial Profile', 9, 5, BG, PANEL);
    hold(ax,'on');
    fill(ax, [z_ax; flipud(z_ax)]', [prof; zeros(size(prof))]', GREEN_C, 'FaceAlpha',0.12,'EdgeColor','none', 'HandleVisibility','off');
    plot(ax, z_ax, prof, 'Color', GREEN_C, 'LineWidth', 2, 'DisplayName','Axial profile (phi-R avg)');
    
    [~, ip] = max(prof);
    text(ax, z_ax(ip)+0.05, prof(ip)*1.05, sprintf('%.3e', prof(ip)), ...
        'Color',AMBER,'FontSize',8,'FontName','Courier New');
    legend(ax,'Location','best','TextColor',TEXT_C,'EdgeColor',BORDER,'Color',PANEL,'FontSize',8,'FontName','Courier New');
    apply_dark_axes(ax, BG, PANEL, BORDER, TEXT_C, DIM);
    xlabel(ax,'Z (cm)'); ylabel(ax, sprintf('%s (%s)', hd.quantity, hd.unit));
    title(ax, sprintf('Axial Dose Profile  (averaged over all phi and R)  |  NUI = %.4f', nui));
    add_stamp(ax, DIM);
    save_tab_fig(ax, 'fig04_axial_profile.png', HERE);
end