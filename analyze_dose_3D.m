%% ============================================================================
%  analyze_dose_profiles.m  —  UITF ND3  Dose-to-Medium  Radial & Axial Profiles
%% ============================================================================
%  Reads Dose_R.csv and Dose_Z.csv (or Dose_ND3_Binned3D.csv / DoseAtSample.csv
%  as fallback) produced by TOPAS scorers Sc/Dose_R and Sc/Dose_Z, both scored
%  on ND3Target (25 R × 72 Phi × 50 Z bins).
%
%  For each scorer the script:
%    1.  Parses the TOPAS CSV header for bin geometry.
%    2.  Reconstructs the full 3-D array (R × Phi × Z).
%    3.  Collapses it to a 1-D profile by averaging over the other two axes.
%    4.  Calculates the Non-Uniformity Index  NUI = sigma / x_bar
%        (coefficient of variation, dimensionless) over the active bins.
%    5.  Saves a publication-quality PNG figure (dark theme, matching the
%        existing analysis suite).
%    6.  Saves a plain-background PNG figure ("raw" version, white background,
%        suitable for reports / quick inspection).
%    7.  Exports a two-column CSV  (axis_cm, dose_Gy).
%
%  TOPAS CSV FORMAT ASSUMED
%  ────────────────────────
%  Header:  comment lines beginning with '#'.  Relevant tags parsed:
%    "R in N bins of X cm"
%    "Phi in N bins of X deg"
%    "Z in N bins of X cm"
%    "DoseToMedium (Gy) : Sum"   (or Mean – handled)
%  Data rows:  r_bin_idx , phi_bin_idx , z_bin_idx , value [, variance_col]
%    Indices are 0-based integers.
%    A trailing comma or variance column is silently ignored.
%
%  USAGE
%  ─────
%    Place this file in the same directory as the TOPAS output CSVs and run:
%        analyze_dose_profiles
%    from the MATLAB command window or R-click → Run.
%
%  OUTPUT FILES (written to same directory)
%  ────────────────────────────────────────
%    radial_profile_dark.png      dark-theme figure
%    radial_profile_raw.png       white-background figure
%    radial_profile.csv           R_cm, Dose_Gy
%    axial_profile_dark.png
%    axial_profile_raw.png
%    axial_profile.csv            Z_cm, Dose_Gy
%    dose_profiles_report.txt     plain-text summary of NUI and statistics
%
%  Dependencies: base MATLAB only (no toolboxes required).
%% ============================================================================

clear; clc; close all;

%% ── Working directory ───────────────────────────────────────────────────────
HERE = fileparts(mfilename('fullpath'));
if isempty(HERE)
    HERE = pwd;
end

fprintf('\n%s\n', repmat('=',1,60));
fprintf('  UITF ND3  —  Dose Profile Analysis\n');
fprintf('  Directory: %s\n', HERE);
fprintf('%s\n\n', repmat('=',1,60));

%% ── Geometry constants (locked to TOPAS parameter file) ────────────────────
ND3_RMIN   = 1.25;   % cm  inner radius of ND3 annulus
ND3_RMAX   = 2.50;   % cm  outer radius
ND3_HL     = 1.25;   % cm  half-length (full length = 2.50 cm)
ND3_RASTER_X = 2.5;  % cm  beam raster full width in Z direction
ND3_RASTER_Y = 2.0;  % cm  beam raster full height

%% ── Colour palettes ─────────────────────────────────────────────────────────
% Dark theme (matches existing suite)
DK.bg     = [0.043 0.059 0.102];
DK.panel  = [0.075 0.098 0.161];
DK.border = [0.165 0.208 0.314];
DK.text   = [0.910 0.929 0.961];
DK.dim    = [0.353 0.416 0.541];
DK.cyan   = [0.000 0.898 1.000];
DK.amber  = [1.000 0.702 0.278];
DK.green  = [0.298 0.686 0.314];
DK.red    = [1.000 0.420 0.420];

% Raw / white theme
WH.bg     = [1 1 1];
WH.panel  = [0.97 0.97 0.97];
WH.border = [0.6  0.6  0.6];
WH.text   = [0.1  0.1  0.1];
WH.dim    = [0.4  0.4  0.4];
WH.cyan   = [0.0  0.45 0.70];  % muted blue
WH.amber  = [0.80 0.40 0.00];  % burnt orange
WH.green  = [0.10 0.50 0.20];
WH.red    = [0.80 0.10 0.10];

%% ── Candidate file lists ────────────────────────────────────────────────────
% 1-D scorer files (preferred: two columns or four columns)
DOSE_R_CANDIDATES = {'Dose_R.csv', 'DoseR.csv'};
DOSE_Z_CANDIDATES = {'Dose_Z.csv', 'DoseZ.csv'};
% Full 3-D file as fallback
DOSE_3D_CANDIDATES = {'Dose_ND3_Binned3D.csv', 'DoseAtSample.csv', ...
                       'DoseND3Binned3D.csv'};

%% ── Load data ───────────────────────────────────────────────────────────────
fprintf('Loading CSV files...\n');

[hR, dataR, r_ax, phi_ax, z_ax, d3_from_R] = load_scorer( ...
    DOSE_R_CANDIDATES, DOSE_3D_CANDIDATES, HERE);

% If R scorer loaded its own 3-D array, reuse it for Z as well when Dose_Z
% is absent — avoids reading the 3-D file twice.
[hZ, dataZ, ~, ~, ~, d3_from_Z] = load_scorer( ...
    DOSE_Z_CANDIDATES, DOSE_3D_CANDIDATES, HERE);

% Reconcile 3-D arrays
d3 = [];
if ~isempty(d3_from_R); d3 = d3_from_R; end
if ~isempty(d3_from_Z) && isempty(d3); d3 = d3_from_Z; end

% Build axis vectors if not yet populated
if isempty(r_ax) && ~isempty(d3) && ~isempty(hR) && ~isempty(hR.r_bins)
    r_ax   = make_axis(hR.r_bins,   hR.r_size);
    phi_ax = make_axis(hR.phi_bins, hR.phi_size);
    z_ax   = make_axis(hR.z_bins,   hR.z_size);
end
if isempty(r_ax) && ~isempty(d3) && ~isempty(hZ) && ~isempty(hZ.r_bins)
    r_ax   = make_axis(hZ.r_bins,   hZ.r_size);
    phi_ax = make_axis(hZ.phi_bins, hZ.phi_size);
    z_ax   = make_axis(hZ.z_bins,   hZ.z_size);
end

if isempty(d3) && isempty(dataR) && isempty(dataZ)
    error(['No recognised dose CSV files found in:\n  %s\n\n' ...
           'Expected one or more of:\n  %s\n  %s\n  %s'], ...
        HERE, DOSE_R_CANDIDATES{1}, DOSE_Z_CANDIDATES{1}, DOSE_3D_CANDIDATES{1});
end

%% ── Derive profiles ─────────────────────────────────────────────────────────

% --- Radial profile ---
if ~isempty(dataR) && size(dataR,2) == 4 && isempty(d3_from_R)
    % 4-column scorer: same bin structure, reduce to radial by averaging Phi&Z
    [r_prof, r_ax] = profile_from_4col(dataR, hR, 'R');
elseif ~isempty(d3_from_R)
    r_prof = phi_z_mean(d3_from_R);      % (R,)
elseif ~isempty(d3)
    r_prof = phi_z_mean(d3);
else
    r_prof = [];
    fprintf('[WARNING]  Cannot derive radial profile — no usable data\n');
end

% --- Axial profile ---
if ~isempty(dataZ) && size(dataZ,2) == 4 && isempty(d3_from_Z)
    [z_prof, z_ax] = profile_from_4col(dataZ, hZ, 'Z');
elseif ~isempty(d3_from_Z)
    z_prof = phi_r_mean(d3_from_Z);      % (Z,)
elseif ~isempty(d3)
    z_prof = phi_r_mean(d3);
else
    z_prof = [];
    fprintf('[WARNING]  Cannot derive axial profile — no usable data\n');
end

% ── Unit / quantity labels ───────────────────────────────────────────────────
qty_lbl  = 'DoseToMedium';
unit_lbl = 'Gy';
if ~isempty(hR) && ~isempty(hR.quantity); qty_lbl  = hR.quantity; end
if ~isempty(hR) && ~isempty(hR.unit);     unit_lbl = hR.unit;     end

%% ── Statistics and NUI ──────────────────────────────────────────────────────

radial_stats = struct('nui',NaN,'sigma',NaN,'mean',NaN,'max',NaN,'min',NaN, ...
                      'peak_r',NaN,'n_active',0);
axial_stats  = struct('nui',NaN,'sigma',NaN,'mean',NaN,'max',NaN,'min',NaN, ...
                      'peak_z',NaN,'n_active',0);

if ~isempty(r_prof) && ~isempty(r_ax)
    radial_stats = calc_stats(r_prof, r_ax);
    fprintf('Radial  NUI (sigma/mean) = %.4f   peak R = %.3f cm   N_active = %d\n', ...
        radial_stats.nui, radial_stats.peak_r, radial_stats.n_active);
end

if ~isempty(z_prof) && ~isempty(z_ax)
    axial_stats = calc_stats(z_prof, z_ax);
    fprintf('Axial   NUI (sigma/mean) = %.4f   peak Z = %.3f cm   N_active = %d\n', ...
        axial_stats.nui, axial_stats.peak_z, axial_stats.n_active);
end

%% ── Figures ─────────────────────────────────────────────────────────────────
fprintf('\nGenerating figures...\n');

if ~isempty(r_prof) && ~isempty(r_ax)
    make_radial_figure(r_ax, r_prof, radial_stats, qty_lbl, unit_lbl, ...
        ND3_RMIN, ND3_RMAX, DK, 'radial_profile_dark.png', HERE, true);
    make_radial_figure(r_ax, r_prof, radial_stats, qty_lbl, unit_lbl, ...
        ND3_RMIN, ND3_RMAX, WH, 'radial_profile_raw.png',  HERE, false);
    export_profile_csv(r_ax, r_prof, 'R_cm', ['Dose_' unit_lbl], ...
        fullfile(HERE,'radial_profile.csv'));
end

if ~isempty(z_prof) && ~isempty(z_ax)
    make_axial_figure(z_ax, z_prof, axial_stats, qty_lbl, unit_lbl, ...
        ND3_HL, ND3_RASTER_X, DK, 'axial_profile_dark.png', HERE, true);
    make_axial_figure(z_ax, z_prof, axial_stats, qty_lbl, unit_lbl, ...
        ND3_HL, ND3_RASTER_X, WH, 'axial_profile_raw.png',  HERE, false);
    export_profile_csv(z_ax, z_prof, 'Z_cm', ['Dose_' unit_lbl], ...
        fullfile(HERE,'axial_profile.csv'));
end

%% ── Plain-text report ────────────────────────────────────────────────────────
write_report(radial_stats, axial_stats, qty_lbl, unit_lbl, ...
    ND3_RMIN, ND3_RMAX, ND3_HL, HERE);

fprintf('\n%s\n', repmat('=',1,60));
fprintf('  Done.\n');
fprintf('%s\n\n', repmat('=',1,60));


%% ============================================================================
%  LOCAL FUNCTIONS
%% ============================================================================

% ── File finder ──────────────────────────────────────────────────────────────

function p = first_existing(names, base)
    p = [];
    for k = 1:numel(names)
        c = fullfile(base, names{k});
        if exist(c,'file')
            p = c;
            return;
        end
    end
end

% ── TOPAS header parser ───────────────────────────────────────────────────────

function h = parse_header(filepath)
    h.version = ''; h.scorer = ''; h.component = '';
    h.quantity = ''; h.unit = ''; h.stat = '';
    h.r_bins = []; h.r_size = [];
    h.phi_bins = []; h.phi_size = [];
    h.z_bins = []; h.z_size = [];
    h.n_hdr = 0;

    fid = fopen(filepath,'r');
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
            parts = strsplit(text,':');
            if numel(parts)>1; h.version = strtrim(parts{2}); end
        elseif contains(lo,'results for scorer')
            parts = strsplit(text,':');
            if numel(parts)>1; h.scorer = strtrim(parts{2}); end
        elseif contains(lo,'scored in component')
            parts = strsplit(text,':');
            if numel(parts)>1; h.component = strtrim(parts{2}); end
        end

        % "R in N bins of X cm"
        tok = regexp(text, ...
            '([A-Za-z]+)\s+in\s+(\d+)\s+bins?\s+of\s+([\d.eE+\-]+)', ...
            'tokens','ignorecase');
        if ~isempty(tok)
            aname = upper(tok{1}{1});
            n  = str2double(tok{1}{2});
            sz = str2double(tok{1}{3});
            if     strcmp(aname,'R');                h.r_bins=n;   h.r_size=sz;
            elseif contains(aname,'PH');             h.phi_bins=n; h.phi_size=sz;
            elseif strcmp(aname,'Z');                h.z_bins=n;   h.z_size=sz;
            end
        end

        % "DoseToMedium (Gy) : Sum"
        tok2 = regexp(text,'(\w+)\s*\(\s*([^)]+)\s*\)\s*:\s*(\w+)','tokens');
        if ~isempty(tok2)
            h.quantity = strtrim(tok2{1}{1});
            h.unit     = strtrim(tok2{1}{2});
            h.stat     = strtrim(tok2{1}{3});
        end
    end
    fclose(fid);
end

% ── CSV data reader ───────────────────────────────────────────────────────────

function [h, data] = read_topas_csv(filepath)
    h    = parse_header(filepath);
    data = [];

    fid = fopen(filepath,'r');
    if fid < 0
        warning('read_topas_csv: cannot open %s', filepath);
        return;
    end
    for k = 1:h.n_hdr; fgetl(fid); end

    rows = {};
    while true
        raw = fgetl(fid);
        if ~ischar(raw); break; end
        raw = strtrim(raw);
        % Strip trailing comma
        while ~isempty(raw) && (raw(end)==',' || raw(end)==' ')
            raw(end) = [];
        end
        if isempty(raw) || raw(1)=='#'; continue; end
        parts = strsplit(raw,',');
        nums  = str2double(parts);
        if all(~isnan(nums))
            rows{end+1} = nums; %#ok<AGROW>
        end
    end
    fclose(fid);

    if isempty(rows); return; end
    ncols = max(cellfun(@numel,rows));
    mat   = NaN(numel(rows),ncols);
    for k = 1:numel(rows)
        mat(k,1:numel(rows{k})) = rows{k};
    end
    data = mat;
    if size(data,2)==1; data=data(:,1); end
end

% ── 3-D reconstruction ────────────────────────────────────────────────────────

function [arr3d, r, phi, z] = to_3d(data, h)
    R=h.r_bins; Phi=h.phi_bins; Z=h.z_bins;
    arr3d = zeros(R,Phi,Z);

    ri   = min(max(round(data(:,1)),0),R-1);
    phii = min(max(round(data(:,2)),0),Phi-1);
    zi   = min(max(round(data(:,3)),0),Z-1);
    vals = data(:,4);

    lin  = sub2ind([R,Phi,Z], ri+1, phii+1, zi+1);
    for k = 1:numel(lin)
        arr3d(lin(k)) = arr3d(lin(k)) + vals(k);
    end

    r   = make_axis(R,   h.r_size);
    phi = make_axis(Phi, h.phi_size);
    z   = make_axis(Z,   h.z_size);
end

function ax = make_axis(n_bins, bin_size)
    ax = ((0:n_bins-1)' + 0.5) * bin_size;
end

% ── Profile extractors ────────────────────────────────────────────────────────

function out = phi_z_mean(arr3d)
    % Average over Phi and Z -> (R,)
    tmp = mean(arr3d,3);   % (R,Phi)
    out = mean(tmp,2);     % (R,)
    out = out(:);
end

function out = phi_r_mean(arr3d)
    % Average over Phi and R -> (Z,)
    tmp = mean(arr3d,1);   % (1,Phi,Z)
    tmp = mean(tmp,2);     % (1,1,Z)
    out = squeeze(tmp);    % (Z,)
    out = out(:);
end

% ── Profile from 4-column data with header ───────────────────────────────────

function [prof, ax] = profile_from_4col(data, h, axis_type)
    % Reconstruct 3-D then project
    [arr3d, r_a, ~, z_a] = to_3d(data, h);
    if strcmp(axis_type,'R')
        prof = phi_z_mean(arr3d);
        ax   = r_a;
    else
        prof = phi_r_mean(arr3d);
        ax   = z_a;
    end
end

% ── Master load function ──────────────────────────────────────────────────────

function [h, data, r_ax, phi_ax, z_ax, d3] = load_scorer(primary_names, fallback_names, HERE)
    h=[]; data=[]; r_ax=[]; phi_ax=[]; z_ax=[]; d3=[];

    p = first_existing(primary_names, HERE);
    if ~isempty(p)
        fprintf('  Loading %s\n', p);
        [h, data] = read_topas_csv(p);
        if isempty(data)
            fprintf('  [WARNING]  File empty or unreadable: %s\n', p);
            data = [];
        end
        % If 4-column, build 3-D now so we can derive both profiles from one file
        if ~isempty(data) && size(data,2)>=4 && ...
                ~isempty(h.r_bins) && ~isempty(h.phi_bins) && ~isempty(h.z_bins)
            [d3, r_ax, phi_ax, z_ax] = to_3d(data, h);
        end
        return;
    end

    % Fallback: full 3-D file
    p = first_existing(fallback_names, HERE);
    if ~isempty(p)
        fprintf('  Loading (fallback) %s\n', p);
        [h, data] = read_topas_csv(p);
        if ~isempty(data) && size(data,2)>=4 && ...
                ~isempty(h.r_bins) && ~isempty(h.phi_bins) && ~isempty(h.z_bins)
            [d3, r_ax, phi_ax, z_ax] = to_3d(data, h);
        end
        return;
    end

    fprintf('  [INFO]  None of the candidate files found (primary + fallback)\n');
end

% ── Statistics / NUI ─────────────────────────────────────────────────────────

function s = calc_stats(prof, ax)
    v    = prof(:);
    mask = v > 0;
    vact = v(mask);
    aact = ax(mask);

    s.n_active = sum(mask);
    s.mean     = 0;
    s.sigma    = 0;
    s.nui      = 0;
    s.max      = max(v);
    s.min      = min(v(mask));
    s.peak_r   = NaN;
    s.peak_z   = NaN;   % one struct covers both axes; caller uses the right field

    if s.n_active < 2; return; end

    s.mean    = mean(vact);
    s.sigma   = std(vact);
    s.nui     = s.sigma / s.mean;   % coefficient of variation

    [~, ip]   = max(prof);
    s.peak_r  = ax(ip);   % works for either axis
    s.peak_z  = ax(ip);
end

% ── Figure helpers ────────────────────────────────────────────────────────────

function fig = new_fig(w_in, h_in, bg)
    fig = figure('Visible','off','Color',bg, ...
                 'Units','inches','Position',[1 1 w_in h_in]);
end

function style_axes(ax, pal, dark_mode)
    set(ax,'Color',pal.panel, ...
           'XColor',pal.dim,'YColor',pal.dim, ...
           'GridColor',pal.border,'GridAlpha',1.0, ...
           'GridLineStyle','--', ...
           'FontName','Courier New','FontSize',10);
    set(ax.Title, 'Color',pal.text,'FontSize',12,'FontWeight','bold');
    set(ax.XLabel,'Color',pal.text,'FontSize',11);
    set(ax.YLabel,'Color',pal.text,'FontSize',11);
    grid(ax,'on');
    if ~dark_mode
        ax.GridAlpha = 0.3;
    end
end

function stamp(fig, pal)
    annotation(fig,'textbox',[0.55 0.002 0.44 0.030], ...
        'String',sprintf('TOPAS UITF ND3  •  %s',datestr(now,'yyyy-mm-dd')), ...
        'EdgeColor','none','Color',pal.dim, ...
        'FontSize',7,'FontName','Courier New', ...
        'HorizontalAlignment','right','FitBoxToText','off');
end

function save_fig(fig, fname, HERE)
    fpath = fullfile(HERE, fname);
    print(fig, fpath, '-dpng', '-r200');
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

% ── RADIAL PROFILE FIGURE ─────────────────────────────────────────────────────

function make_radial_figure(r_ax, r_prof, stats, qty, unit, ...
                             r_min, r_max, pal, fname, HERE, dark_mode)

    fig = new_fig(10, 6, pal.bg);
    ax  = axes('Parent',fig,'Color',pal.panel, ...
               'Position',[0.10 0.13 0.82 0.77]);
    hold(ax,'on');

    % Filled area under curve
    fill(ax, [r_ax; flipud(r_ax)], [r_prof; zeros(size(r_prof))], ...
         pal.cyan, 'FaceAlpha',0.18, 'EdgeColor','none');

    % Profile line
    h_line = plot(ax, r_ax, r_prof, '-', 'Color',pal.cyan, ...
                  'LineWidth',2.2, 'DisplayName','Dose (phi-Z avg)');

    % ND3 annulus boundaries
    y_top = max(r_prof) * 1.15;
    xpatch = [r_min r_max r_max r_min];
    ypatch = [0 0 y_top y_top];
    fill(ax, xpatch, ypatch, pal.green, 'FaceAlpha',0.07, 'EdgeColor','none', ...
         'HandleVisibility','off');
    plot(ax,[r_min r_min],[0 y_top],'--','Color',pal.amber,'LineWidth',1.0, ...
         'DisplayName',sprintf('R_{min} = %.2f cm', r_min));
    plot(ax,[r_max r_max],[0 y_top],'--','Color',pal.red,  'LineWidth',1.0, ...
         'DisplayName',sprintf('R_{max} = %.2f cm', r_max));

    % Peak marker
    [~,ip] = max(r_prof);
    plot(ax, r_ax(ip), r_prof(ip), 'o', 'Color',pal.amber, ...
         'MarkerFaceColor',pal.amber,'MarkerSize',7, ...
         'DisplayName',sprintf('Peak = %.3e %s  @  R=%.3f cm', ...
             r_prof(ip), unit, r_ax(ip)));

    % NUI legend entry (invisible line for legend text)
    nui_str = sprintf('NUI  \\sigma/\\bar{x} = %.4f', stats.nui);
    if stats.n_active > 1
        nui_str = sprintf('%s   (\\sigma=%.3e,  \\bar{x}=%.3e %s)', ...
                          nui_str, stats.sigma, stats.mean, unit);
    end
    plot(ax,NaN,NaN,'s','Color',pal.dim,'MarkerFaceColor',pal.dim, ...
         'MarkerSize',0,'DisplayName',nui_str);

    % Axis limits and labels
    xlim(ax,[0 r_ax(end)*1.05]);
    ylim(ax,[0 y_top]);
    xlabel(ax,'R  (cm)');
    ylabel(ax,sprintf('%s  (%s)',qty,unit));
    title(ax,'Radial Dose Profile  —  DoseToMedium  (averaged over all \phi and Z)');

    style_axes(ax, pal, dark_mode);

    lg = legend(ax,'Location','best','FontSize',9,'FontName','Courier New');
    set(lg,'TextColor',pal.text,'EdgeColor',pal.border,'Color',pal.panel);

    stamp(fig,pal);
    save_fig(fig, fname, HERE);
end

% ── AXIAL PROFILE FIGURE ──────────────────────────────────────────────────────

function make_axial_figure(z_ax, z_prof, stats, qty, unit, ...
                            nd3_hl, raster_x, pal, fname, HERE, dark_mode)

    fig = new_fig(10, 6, pal.bg);
    ax  = axes('Parent',fig,'Color',pal.panel, ...
               'Position',[0.10 0.13 0.82 0.77]);
    hold(ax,'on');

    z_full       = 2 * nd3_hl;
    raster_start = (z_full - raster_x) / 2;
    raster_end   = raster_start + raster_x;
    y_top        = max(z_prof) * 1.15;

    % Raster footprint shading
    rs = max(raster_start, z_ax(1));
    re = min(raster_end,   z_ax(end));
    fill(ax,[rs re re rs],[0 0 y_top y_top], pal.amber, ...
         'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');

    % Filled area under curve
    fill(ax,[z_ax; flipud(z_ax)],[z_prof; zeros(size(z_prof))], ...
         pal.green,'FaceAlpha',0.18,'EdgeColor','none');

    % Profile line
    plot(ax, z_ax, z_prof, '-', 'Color',pal.green, ...
         'LineWidth',2.2,'DisplayName','Dose (phi-R avg)');

    % Target length boundaries
    plot(ax,[z_ax(1) z_ax(1)],[0 y_top],   '--','Color',pal.dim, 'LineWidth',0.8, ...
         'HandleVisibility','off');
    plot(ax,[z_ax(end) z_ax(end)],[0 y_top],'--','Color',pal.dim,'LineWidth',0.8, ...
         'DisplayName',sprintf('Target length = %.2f cm', z_full));

    % Raster extent markers
    plot(ax,[rs rs],[0 y_top],':','Color',pal.amber,'LineWidth',1.0, ...
         'DisplayName',sprintf('Raster start = %.2f cm', rs));
    plot(ax,[re re],[0 y_top],':','Color',pal.red,  'LineWidth',1.0, ...
         'DisplayName',sprintf('Raster end   = %.2f cm', re));

    % Peak marker
    [~,ip] = max(z_prof);
    plot(ax, z_ax(ip), z_prof(ip), 'o','Color',pal.amber, ...
         'MarkerFaceColor',pal.amber,'MarkerSize',7, ...
         'DisplayName',sprintf('Peak = %.3e %s  @  Z=%.3f cm', ...
             z_prof(ip), unit, z_ax(ip)));

    % NUI legend entry
    nui_str = sprintf('NUI  \\sigma/\\bar{x} = %.4f', stats.nui);
    if stats.n_active > 1
        nui_str = sprintf('%s   (\\sigma=%.3e,  \\bar{x}=%.3e %s)', ...
                          nui_str, stats.sigma, stats.mean, unit);
    end
    plot(ax,NaN,NaN,'s','Color',pal.dim,'MarkerFaceColor',pal.dim, ...
         'MarkerSize',0,'DisplayName',nui_str);

    xlim(ax,[0 z_ax(end)*1.05]);
    ylim(ax,[0 y_top]);
    xlabel(ax,'Z  (cm)');
    ylabel(ax,sprintf('%s  (%s)',qty,unit));
    title(ax,'Axial Dose Profile  —  DoseToMedium  (averaged over all \phi and R)');

    style_axes(ax, pal, dark_mode);

    lg = legend(ax,'Location','best','FontSize',9,'FontName','Courier New');
    set(lg,'TextColor',pal.text,'EdgeColor',pal.border,'Color',pal.panel);

    stamp(fig,pal);
    save_fig(fig, fname, HERE);
end

% ── CSV export ────────────────────────────────────────────────────────────────

function export_profile_csv(ax, prof, col1, col2, fpath)
    fid = fopen(fpath,'w');
    if fid < 0
        warning('export_profile_csv: cannot write %s', fpath);
        return;
    end
    fprintf(fid,'%s,%s\n',col1,col2);
    for k = 1:numel(ax)
        fprintf(fid,'%.6f,%.6e\n', ax(k), prof(k));
    end
    fclose(fid);
    [~,fn,ext] = fileparts(fpath);
    fprintf('  Saved: %s%s\n', fn, ext);
end

% ── Plain-text report ─────────────────────────────────────────────────────────

function write_report(rs, as, qty, unit, r_min, r_max, nd3_hl, HERE)
    fpath = fullfile(HERE,'dose_profiles_report.txt');
    fid   = fopen(fpath,'w');
    if fid < 0
        warning('write_report: cannot write report');
        return;
    end
    SEP = repmat('-',1,60);
    fprintf(fid,'%s\n',repmat('=',1,60));
    fprintf(fid,'  UITF ND3  —  Dose Profile Report\n');
    fprintf(fid,'  %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid,'%s\n\n',repmat('=',1,60));

    fprintf(fid,'GEOMETRY\n%s\n', SEP);
    fprintf(fid,'  ND3 annulus  : R = %.2f - %.2f cm\n',r_min,r_max);
    fprintf(fid,'  ND3 length   : %.2f cm (half-length %.2f cm)\n\n', ...
            2*nd3_hl, nd3_hl);

    fprintf(fid,'QUANTITY  : %s  (%s)\n\n', qty, unit);

    fprintf(fid,'RADIAL PROFILE (phi-Z averaged)\n%s\n', SEP);
    fprintf_stats(fid, rs, unit, 'R');

    fprintf(fid,'AXIAL PROFILE (phi-R averaged)\n%s\n', SEP);
    fprintf_stats(fid, as, unit, 'Z');

    fprintf(fid,'OUTPUT FILES\n%s\n', SEP);
    fprintf(fid,'  radial_profile_dark.png  (dark theme)\n');
    fprintf(fid,'  radial_profile_raw.png   (white background)\n');
    fprintf(fid,'  radial_profile.csv\n');
    fprintf(fid,'  axial_profile_dark.png\n');
    fprintf(fid,'  axial_profile_raw.png\n');
    fprintf(fid,'  axial_profile.csv\n');
    fprintf(fid,'%s\n',repmat('=',1,60));
    fclose(fid);
    fprintf('  Saved: dose_profiles_report.txt\n');
end

function fprintf_stats(fid, s, unit, axis_lbl)
    if isnan(s.nui)
        fprintf(fid,'  [no data]\n\n');
        return;
    end
    fprintf(fid,'  N active bins : %d\n', s.n_active);
    fprintf(fid,'  Max           : %.6e %s\n', s.max, unit);
    fprintf(fid,'  Min (nonzero) : %.6e %s\n', s.min, unit);
    fprintf(fid,'  Mean          : %.6e %s\n', s.mean, unit);
    fprintf(fid,'  Sigma         : %.6e %s\n', s.sigma, unit);
    fprintf(fid,'  NUI = sigma/mean : %.6f  (%.2f%%)\n', s.nui, s.nui*100);
    fprintf(fid,'  Peak %s       : %.3f cm\n\n', axis_lbl, s.peak_r);
end
