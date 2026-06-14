%% analyze_nui_profiles.m
%  Companion TOPAS post-processing for ND3 1D profiles and azimuthal R-Z map.
%  Run from the directory containing Dose_Z.csv, Dose_R.csv, and Dose_RZ.csv.

clear; clc; close all;
set(groot, 'defaultFigureVisible', 'on');

HERE = fileparts(mfilename('fullpath'));
if isempty(HERE)
    HERE = pwd;
end

TOPAS_CONFIG = fullfile(HERE, 'UITF_ND3.txt');
BEAM_ENERGY_MEV = read_beam_energy_mev(TOPAS_CONFIG, 7.0);

% Native ND3Target geometry from UITF_ND3.txt.
ND3_RMIN_CM = 1.25;
ND3_Z_HL_CM = 1.25;

axial = load_projected_profile(fullfile(HERE, 'Dose_Z.csv'), 'axial', ND3_RMIN_CM, ND3_Z_HL_CM);
radial = load_projected_profile(fullfile(HERE, 'Dose_R.csv'), 'radial', ND3_RMIN_CM, ND3_Z_HL_CM);
rz = load_projected_profile(fullfile(HERE, 'Dose_RZ.csv'), 'rz', ND3_RMIN_CM, ND3_Z_HL_CM);
energySummary = load_target_energy_summary(HERE, TOPAS_CONFIG, BEAM_ENERGY_MEV);
energyComponents = load_energy_component_summaries(HERE, TOPAS_CONFIG, BEAM_ENERGY_MEV);

axialNorm = normalize_to_peak(axial.dose);
radialNorm = normalize_to_peak(radial.dose);
rzNorm = normalize_to_peak(rz.dose);
brems = load_bremsstrahlung_spectrum(HERE, BEAM_ENERGY_MEV);

plotItems = build_plot_items(axial, radial, rz, axialNorm, radialNorm, rzNorm, ...
    energySummary, energyComponents, brems, BEAM_ENERGY_MEV);
open_analysis_panel(plotItems, axial, radial, energyComponents);

%% Local functions

function style_axes_text(ax, accentColor)
% Recolor every text element on an axes (axis lines, ticks, tick labels,
% title, x/y/z labels) to a vivid, non-grey color so nothing on the plot
% defaults to MATLAB's standard grey/near-black text.
    set(ax, 'XColor', accentColor, 'YColor', accentColor, 'ZColor', accentColor, ...
        'GridColor', accentColor, 'GridAlpha', 0.25, 'FontWeight', 'bold');
    if ~isempty(ax.Title)
        set(ax.Title, 'Color', accentColor);
    end
    if ~isempty(ax.XLabel)
        set(ax.XLabel, 'Color', accentColor);
    end
    if ~isempty(ax.YLabel)
        set(ax.YLabel, 'Color', accentColor);
    end
end

function set_profile_ylim_zero(ax, values)
    finiteValues = values(isfinite(values));
    if isempty(finiteValues)
        ylim(ax, [0 1]);
        return;
    end

    ymax = max(finiteValues);
    if ymax <= 0
        ylim(ax, [0 1]);
    else
        ylim(ax, [0 ymax * 1.10]);
    end
end

function out = normalize_to_peak(values)
    out = values;
    finiteValues = values(isfinite(values));
    if isempty(finiteValues)
        return;
    end

    peak = max(finiteValues);
    if peak > 0
        out = values ./ peak;
    end
end

function add_energy_legend(ax, dataHandles, dataLabels, energySummary, location)
    hold(ax, 'on');
    handles = dataHandles;
    labels = dataLabels;

    if energySummary.available
        hTarget = plot(ax, NaN, NaN, 'LineStyle', 'none', 'Marker', 'none');
        hPercent = plot(ax, NaN, NaN, 'LineStyle', 'none', 'Marker', 'none');
        handles = [handles(:).' hTarget hPercent];
        labels = [labels(:).' {energySummary.targetLabel, energySummary.percentLabel}];
    end

    if ~isempty(handles)
        legend(ax, handles, labels, 'Location', location);
    end
end

function components = load_energy_component_summaries(baseDir, configPath, beamEnergyMeV)
    files = dir(fullfile(baseDir, 'EnergyDep_*.csv'));
    components = repmat(default_energy_component(), 0, 1);

    primaryElectrons = read_primary_electron_count(configPath);
    incidentMeV = beamEnergyMeV * primaryElectrons;

    for k = 1:numel(files)
        filePath = fullfile(files(k).folder, files(k).name);
        item = default_energy_component();
        item.file = files(k).name;
        item.component = component_name_from_file(files(k).name);
        item.unit = 'MeV';
        item.incidentMeV = incidentMeV;

        if files(k).bytes == 0
            item.status = 'empty file';
            components(end+1) = item; %#ok<AGROW>
            continue;
        end

        try
            [h, data] = read_topas_csv(filePath);
            [energyValues, ~] = select_energy_values(data, h);
            item.component = h.component;
            item.unit = h.unit;
            item.totalMeV = sum(energyValues, 'omitnan');
            item.available = isfinite(item.totalMeV);
            item.status = 'ok';
            if isfinite(incidentMeV) && incidentMeV > 0
                item.percentIncident = 100.0 * item.totalMeV / incidentMeV;
            end
        catch ME
            item.status = ME.message;
        end

        components(end+1) = item; %#ok<AGROW>
    end

    if isempty(components)
        return;
    end

    totalScoredMeV = sum([components([components.available]).totalMeV], 'omitnan');
    if isfinite(totalScoredMeV) && totalScoredMeV > 0
        for k = 1:numel(components)
            if components(k).available
                components(k).percentScored = 100.0 * components(k).totalMeV / totalScoredMeV;
            end
        end
    end
end

function item = default_energy_component()
    item = struct( ...
        'file', '', ...
        'component', '', ...
        'unit', 'MeV', ...
        'totalMeV', NaN, ...
        'incidentMeV', NaN, ...
        'percentIncident', NaN, ...
        'percentScored', NaN, ...
        'available', false, ...
        'status', 'not loaded');
end

function [values, valueCol] = select_energy_values(data, h)
    nIndexCols = sum(~cellfun(@isempty, {h.rBins, h.phiBins, h.zBins}));
    nIndexCols = min(nIndexCols, size(data, 2) - 1);
    nIndexCols = max(0, nIndexCols);
    valueCol = nIndexCols + 1;
    values = data(:, valueCol);
end

function name = component_name_from_file(fileName)
    [~, stem, ~] = fileparts(fileName);
    name = regexprep(stem, '^EnergyDep_', '');
end

function label = pretty_component_name(name)
    label = regexprep(name, '(?<!^)([A-Z])', ' $1');
    label = strrep(label, 'N D3', 'ND3');
end

function safe = sanitize_filename(name)
    safe = regexprep(name, '[^A-Za-z0-9]+', '_');
    safe = regexprep(safe, '^_+|_+$', '');
end

function summary = load_target_energy_summary(baseDir, configPath, beamEnergyMeV)
    summary.available = false;
    summary.targetMeV = NaN;
    summary.unit = 'MeV';
    summary.primaryElectrons = NaN;
    summary.incidentMeV = NaN;
    summary.percent = NaN;
    summary.targetLabel = 'Target Edep: unavailable';
    summary.percentLabel = 'Target fraction: unavailable';

    energyPath = fullfile(baseDir, 'EnergyDep_ND3.csv');
    if ~isfile(energyPath)
        return;
    end
    info = dir(energyPath);
    if isempty(info) || info.bytes == 0
        return;
    end

    [h, data] = read_topas_csv(energyPath);
    [energyValues, ~] = select_energy_values(data, h);
    totalMeV = sum(energyValues, 'omitnan');
    primaryElectrons = read_primary_electron_count(configPath);
    incidentMeV = beamEnergyMeV * primaryElectrons;

    summary.available = isfinite(totalMeV);
    summary.targetMeV = totalMeV;
    summary.unit = h.unit;
    summary.primaryElectrons = primaryElectrons;
    summary.incidentMeV = incidentMeV;
    if isfinite(incidentMeV) && incidentMeV > 0
        summary.percent = 100.0 * totalMeV / incidentMeV;
    end

    summary.targetLabel = sprintf('Target Edep: %.4g %s', totalMeV, h.unit);
    if isfinite(summary.percent)
        summary.percentLabel = sprintf('Target share: %.3f%% of incident', summary.percent);
    end
end

function primaryElectrons = read_primary_electron_count(configPath)
    primaryElectrons = read_integer_parameter(configPath, 'MC/NumberOfPrimaryElectrons');
    if ~isfinite(primaryElectrons)
        primaryElectrons = read_integer_parameter(configPath, 'MC/NumberOfHistories');
    end
end

function value = read_integer_parameter(configPath, parameterName)
    value = NaN;
    if ~isfile(configPath)
        return;
    end

    fid = fopen(configPath, 'r');
    if fid < 0
        return;
    end
    cleanup = onCleanup(@() fclose(fid));

    escapedName = regexptranslate('escape', parameterName);
    pattern = ['^\s*i:' escapedName '\s*=\s*([0-9.eE+-]+)'];
    while true
        raw = fgetl(fid);
        if ~ischar(raw)
            break;
        end
        tok = regexp(raw, pattern, 'tokens', 'once');
        if ~isempty(tok)
            value = str2double(tok{1});
            return;
        end
    end
end

function beamEnergyMeV = read_beam_energy_mev(configPath, fallbackMeV)
    beamEnergyMeV = fallbackMeV;
    if ~isfile(configPath)
        return;
    end

    fid = fopen(configPath, 'r');
    if fid < 0
        return;
    end
    cleanup = onCleanup(@() fclose(fid));

    while true
        raw = fgetl(fid);
        if ~ischar(raw)
            break;
        end
        tok = regexp(raw, '^\s*d:So/BeamSource/BeamEnergy\s*=\s*([0-9.eE+-]+)\s*MeV', 'tokens', 'once');
        if ~isempty(tok)
            beamEnergyMeV = str2double(tok{1});
            return;
        end
    end
end

function brems = load_bremsstrahlung_spectrum(baseDir, beamEnergyMeV)
    files = dir(fullfile(baseDir, 'PhaseSpace_Escape_*.phsp'));
    nBins = 140;
    edges = linspace(0, beamEnergyMeV * 1.05, nBins + 1);

    brems.available = ~isempty(files);
    brems.edges = edges;
    brems.binCenters = edges(1:end-1) + diff(edges) / 2;
    brems.totalCounts = zeros(1, nBins);
    brems.forwardCounts = zeros(1, nBins);
    brems.endpointMeV = NaN;
    brems.gammaCount = 0;
    brems.filesUsed = 0;
    brems.duaneHuntLambdaPm = 1.239841984 / beamEnergyMeV;

    if ~brems.available
        return;
    end

    maxEnergy = -Inf;
    for k = 1:numel(files)
        path = fullfile(files(k).folder, files(k).name);
        isForward = contains(files(k).name, 'Escape_ZPlus');
        [counts, forwardCounts, fileMaxEnergy, gammaCount] = phase_space_gamma_histogram(path, edges, isForward);
        brems.totalCounts = brems.totalCounts + counts;
        brems.forwardCounts = brems.forwardCounts + forwardCounts;
        brems.gammaCount = brems.gammaCount + gammaCount;
        brems.filesUsed = brems.filesUsed + 1;
        if isfinite(fileMaxEnergy)
            maxEnergy = max(maxEnergy, fileMaxEnergy);
        end
    end

    if isfinite(maxEnergy)
        brems.endpointMeV = maxEnergy;
    end
end

function [counts, forwardCounts, maxEnergy, gammaCount] = phase_space_gamma_histogram(path, edges, isForward)
    nBins = numel(edges) - 1;
    counts = zeros(1, nBins);
    forwardCounts = zeros(1, nBins);
    maxEnergy = NaN;
    gammaCount = 0;

    fid = fopen(path, 'r');
    if fid < 0
        warning('Could not open phase-space file: %s', path);
        return;
    end
    cleanup = onCleanup(@() fclose(fid));

    seek_to_first_numeric_line(fid);
    chunkRows = 250000;
    format = '%*f %*f %*f %*f %*f %f %f %f %*f %*f %*f %*f %*f %*f %*[^\n]';

    while ~feof(fid)
        block = textscan(fid, format, chunkRows, ...
            'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
        if isempty(block) || isempty(block{1})
            break;
        end

        energy = block{1};
        weight = block{2};
        pdg = block{3};

        isGamma = pdg == 22 & isfinite(energy) & energy > 0 & isfinite(weight);
        if ~any(isGamma)
            continue;
        end

        gammaEnergy = energy(isGamma);
        gammaWeight = weight(isGamma);
        gammaCount = gammaCount + numel(gammaEnergy);
        maxEnergy = max([maxEnergy; gammaEnergy(:)], [], 'omitnan');

        h = weighted_histcounts(gammaEnergy, edges, gammaWeight);
        counts = counts + h;
        if isForward
            forwardCounts = forwardCounts + h;
        end
    end
end

function counts = weighted_histcounts(values, edges, weights)
    nBins = numel(edges) - 1;
    [~, ~, bin] = histcounts(values, edges);
    valid = bin >= 1 & bin <= nBins & isfinite(weights);
    if any(valid)
        counts = accumarray(bin(valid), weights(valid), [nBins, 1], @sum, 0).';
    else
        counts = zeros(1, nBins);
    end
end

function seek_to_first_numeric_line(fid)
    while true
        pos = ftell(fid);
        raw = fgetl(fid);
        if ~ischar(raw)
            return;
        end
        if ~isempty(regexp(raw, '^\s*[-+]?\d', 'once'))
            fseek(fid, pos, 'bof');
            return;
        end
    end
end

function out = load_projected_profile(filepath, mode, rMinCm, zHalfLengthCm)
    [h, data] = read_topas_csv(filepath);
    [doseValues, varianceValues, nIndexCols, selected] = select_value_columns(data, h);
    axesInfo = make_axes(h, rMinCm, zHalfLengthCm);

    doseGrid = values_to_grid(data(:, 1:nIndexCols), doseValues, axesInfo);
    varianceGrid = values_to_grid(data(:, 1:nIndexCols), varianceValues, axesInfo);

    switch lower(mode)
        case 'axial'
            [profileDose, profileVariance] = project_axial(doseGrid, varianceGrid, axesInfo.r);
            axisValues = axesInfo.z;
            [nui, method] = compute_profile_nui(profileDose, profileVariance);
            out = base_output(h, selected, profileDose, nui, method);
            out.axis = axisValues;
        case 'radial'
            [profileDose, profileVariance] = project_radial(doseGrid, varianceGrid);
            axisValues = axesInfo.r;
            [nui, method] = compute_profile_nui(profileDose, profileVariance);
            out = base_output(h, selected, profileDose, nui, method);
            out.axis = axisValues;
        case 'rz'
            [mapDose, mapVariance] = project_rz(doseGrid, varianceGrid); %#ok<ASGLU>
            out = base_output(h, selected, mapDose(:), NaN, 'not computed for map');
            out.dose = mapDose;
            out.r_axis = axesInfo.r;
            out.z_axis = axesInfo.z;
        otherwise
            error('Unknown projection mode: %s', mode);
    end
end

function out = base_output(h, selected, dose, nui, method)
    out.dose = dose(:);
    out.meanDose = mean(out.dose, 'omitnan');
    out.nui = nui;
    out.nuiMethod = method;
    out.unit = h.unit;
    out.doseLabel = selected.doseLabel;
    out.varianceLabel = selected.varianceLabel;
end

function [h, data] = read_topas_csv(filepath)
    if ~isfile(filepath)
        error('Missing TOPAS CSV file: %s', filepath);
    end

    h = default_header();
    rows = {};

    fid = fopen(filepath, 'r');
    if fid < 0
        error('Could not open TOPAS CSV file: %s', filepath);
    end

    cleanup = onCleanup(@() fclose(fid));
    while true
        raw = fgetl(fid);
        if ~ischar(raw)
            break;
        end
        raw = strtrim(raw);
        if isempty(raw)
            continue;
        end

        if startsWith(raw, '#')
            h = parse_header_line(h, raw);
            continue;
        end

        nums = parse_numeric_row(raw);
        if ~isempty(nums)
            rows{end+1} = nums; %#ok<AGROW>
        elseif any(isletter(raw)) && isempty(h.dataColumnNames)
            h.dataColumnNames = split_header_names(raw);
        end
    end

    if isempty(rows)
        error('No numeric TOPAS data rows found in %s', filepath);
    end

    nCols = max(cellfun(@numel, rows));
    data = NaN(numel(rows), nCols);
    for k = 1:numel(rows)
        data(k, 1:numel(rows{k})) = rows{k};
    end
end

function h = default_header()
    h.scorer = '';
    h.component = '';
    h.quantity = 'Dose';
    h.unit = 'Gy';
    h.valueLabels = {};
    h.dataColumnNames = {};
    h.rBins = [];
    h.phiBins = [];
    h.zBins = [];
    h.rSize = [];
    h.phiSize = [];
    h.zSize = [];
end

function h = parse_header_line(h, raw)
    text = strtrim(regexprep(raw, '^#\s*', ''));
    lo = lower(text);

    if contains(lo, 'results for scorer')
        h.scorer = value_after_colon(text);
    elseif contains(lo, 'scored in component')
        h.component = value_after_colon(text);
    elseif startsWith(lo, 'r in')
        [h.rBins, h.rSize] = parse_axis_bins(text);
    elseif startsWith(lo, 'phi in')
        [h.phiBins, h.phiSize] = parse_axis_bins(text);
    elseif startsWith(lo, 'z in')
        [h.zBins, h.zSize] = parse_axis_bins(text);
    else
        tok = regexp(text, '^(.+?)\s*\(\s*([^)]+)\s*\)\s*:\s*(.*)$', 'tokens', 'once');
        if ~isempty(tok)
            h.quantity = strtrim(tok{1});
            h.unit = strtrim(tok{2});
            h.valueLabels = split_header_names(tok{3});
        end
    end
end

function val = value_after_colon(text)
    parts = strsplit(text, ':');
    if numel(parts) >= 2
        val = strtrim(strjoin(parts(2:end), ':'));
    else
        val = '';
    end
end

function [bins, binSize] = parse_axis_bins(text)
    tok = regexp(text, '(\d+)\s+bins?\s+of\s+([0-9.eE+-]+)', 'tokens', 'once');
    if isempty(tok)
        bins = [];
        binSize = [];
    else
        bins = str2double(tok{1});
        binSize = str2double(tok{2});
    end
end

function names = split_header_names(raw)
    parts = regexp(strtrim(raw), '[,\s]+', 'split');
    names = parts(~cellfun('isempty', parts));
end

function nums = parse_numeric_row(raw)
    raw = regexprep(raw, '[,\s]+$', '');
    parts = regexp(raw, '[,\s]+', 'split');
    parts = parts(~cellfun('isempty', parts));
    nums = str2double(parts);
    if isempty(nums) || any(isnan(nums))
        nums = [];
    end
end

function [doseValues, varianceValues, nIndexCols, selected] = select_value_columns(data, h)
    nIndexCols = guess_index_columns(data, h);
    nValueCols = size(data, 2) - nIndexCols;
    if nValueCols < 1
        error('TOPAS table has no value columns after index columns.');
    end

    labels = h.valueLabels;
    if numel(h.dataColumnNames) == size(data, 2)
        labels = h.dataColumnNames(nIndexCols + 1:end);
    end
    while numel(labels) < nValueCols
        labels{end+1} = sprintf('Value_%d', numel(labels) + 1); %#ok<AGROW>
    end

    labelsLower = lower(labels);
    doseLocal = find(contains(labelsLower, 'dose') | contains(labelsLower, 'value_sum') | ...
        strcmp(labelsLower, 'sum') | contains(labelsLower, 'mean'), 1, 'first');
    if isempty(doseLocal)
        doseLocal = 1;
    end

    varianceLocal = find(contains(labelsLower, 'variance') | contains(labelsLower, 'value_variance') | ...
        strcmp(labelsLower, 'var'), 1, 'first');

    doseValues = data(:, nIndexCols + doseLocal);
    if isempty(varianceLocal)
        varianceValues = NaN(size(doseValues));
        varianceLabel = 'not present; using spatial profile std/mean fallback';
    else
        varianceValues = data(:, nIndexCols + varianceLocal);
        varianceLabel = labels{varianceLocal};
    end

    selected.doseLabel = labels{doseLocal};
    selected.varianceLabel = varianceLabel;
end

function nIndexCols = guess_index_columns(data, h)
    nByHeader = sum(~cellfun(@isempty, {h.rBins, h.phiBins, h.zBins}));
    if nByHeader > 0
        nIndexCols = nByHeader;
    else
        nIndexCols = min(3, size(data, 2) - 1);
    end
    nIndexCols = max(1, min(nIndexCols, size(data, 2) - 1));
end

function axesInfo = make_axes(h, rMinCm, zHalfLengthCm)
    if isempty(h.rBins); h.rBins = 1; end
    if isempty(h.phiBins); h.phiBins = 1; end
    if isempty(h.zBins); h.zBins = 1; end
    if isempty(h.rSize); h.rSize = 1; end
    if isempty(h.zSize); h.zSize = 1; end

    axesInfo.dims = [h.rBins, h.phiBins, h.zBins];
    axesInfo.r = rMinCm + ((0:h.rBins-1)' + 0.5) * h.rSize;
    axesInfo.phi = ((0:h.phiBins-1)' + 0.5);
    axesInfo.z = -zHalfLengthCm + ((0:h.zBins-1)' + 0.5) * h.zSize;
end

function grid = values_to_grid(indexData, values, axesInfo)
    dims = axesInfo.dims;
    sums = NaN(dims);
    counts = zeros(dims);

    for row = 1:size(indexData, 1)
        idx = infer_grid_index(indexData(row, :), dims);
        if any(idx < 1) || any(idx > dims) || isnan(values(row))
            continue;
        end

        lin = sub2ind(dims, idx(1), idx(2), idx(3));
        if isnan(sums(lin))
            sums(lin) = 0;
        end
        sums(lin) = sums(lin) + values(row);
        counts(lin) = counts(lin) + 1;
    end

    grid = sums;
    filled = counts > 0;
    grid(filled) = sums(filled) ./ counts(filled);
end

function idx = infer_grid_index(rowIndexData, dims)
    idx = [1, 1, 1];
    n = numel(rowIndexData);

    if n >= 3
        idx = round(rowIndexData(1:3)) + 1;
        return;
    end

    activeAxes = find(dims > 1);
    if isempty(activeAxes)
        activeAxes = 1;
    end
    for k = 1:min(n, numel(activeAxes))
        idx(activeAxes(k)) = round(rowIndexData(k)) + 1;
    end
end

function [profile, profileVariance] = project_axial(doseGrid, varianceGrid, rAxis)
    dims = size(doseGrid);
    if numel(dims) < 3; dims(3) = 1; end
    profile = NaN(dims(3), 1);
    profileVariance = NaN(dims(3), 1);
    radialWeights = repmat(rAxis(:), 1, dims(2));

    for z = 1:dims(3)
        values = doseGrid(:, :, z);
        weights = radialWeights;
        valid = isfinite(values) & isfinite(weights);
        if any(valid(:))
            w = weights(valid);
            v = values(valid);
            wsum = sum(w);
            profile(z) = sum(v .* w) / wsum;

            varSlice = varianceGrid(:, :, z);
            validVar = valid & isfinite(varSlice) & varSlice >= 0;
            if any(validVar(:))
                wv = weights(validVar) / sum(weights(validVar));
                vv = varSlice(validVar);
                profileVariance(z) = sum((wv .^ 2) .* vv);
            end
        end
    end
end

function [profile, profileVariance] = project_radial(doseGrid, varianceGrid)
    dims = size(doseGrid);
    if numel(dims) < 3; dims(3) = 1; end
    profile = NaN(dims(1), 1);
    profileVariance = NaN(dims(1), 1);

    for r = 1:dims(1)
        values = squeeze(doseGrid(r, :, :));
        valid = isfinite(values);
        if any(valid(:))
            profile(r) = mean(values(valid));
        end

        vars = squeeze(varianceGrid(r, :, :));
        validVar = isfinite(vars) & vars >= 0;
        n = nnz(validVar);
        if n > 0
            profileVariance(r) = sum(vars(validVar)) / (n ^ 2);
        end
    end
end

function [mapDose, mapVariance] = project_rz(doseGrid, varianceGrid)
    dims = size(doseGrid);
    if numel(dims) < 3; dims(3) = 1; end
    mapDose = NaN(dims(1), dims(3));
    mapVariance = NaN(dims(1), dims(3));

    for r = 1:dims(1)
        for z = 1:dims(3)
            values = squeeze(doseGrid(r, :, z));
            valid = isfinite(values);
            if any(valid(:))
                mapDose(r, z) = mean(values(valid));
            end

            vars = squeeze(varianceGrid(r, :, z));
            validVar = isfinite(vars) & vars >= 0;
            n = nnz(validVar);
            if n > 0
                mapVariance(r, z) = sum(vars(validVar)) / (n ^ 2);
            end
        end
    end
end

function [nui, method] = compute_profile_nui(profileDose, profileVariance)
    validDose = isfinite(profileDose);
    mu = mean(profileDose(validDose));
    if isempty(mu) || mu == 0 || ~isfinite(mu)
        nui = NaN;
        method = 'undefined; mean dose is zero or unavailable';
        return;
    end

    validVar = validDose & isfinite(profileVariance) & profileVariance >= 0;
    if any(validVar)
        sigma = sqrt(mean(profileVariance(validVar)));
        nui = sigma / abs(mu);
        method = 'sqrt(mean projected variance) / mean projected dose';
    else
        nui = std(profileDose(validDose), 0) / abs(mu);
        method = 'spatial std(projected profile) / mean projected dose';
    end
end

%% --- Plot catalog + interactive panel -----------------------------------

function items = build_plot_items(axial, radial, rz, axialNorm, radialNorm, rzNorm, ...
    energySummary, energyComponents, brems, beamEnergyMeV)
% Build a struct array describing every available plot. Each entry has a
% short 'name' (shown in the list box) and a 'draw' function handle that
% renders the plot into the currently active axes when called.

    items = struct('name', {}, 'draw', {});

    items(end+1) = struct('name', 'Axial Dose Profile', ...
        'draw', @(ax) draw_axial_profile(ax, axial));

    items(end+1) = struct('name', 'Radial Dose Profile', ...
        'draw', @(ax) draw_radial_profile(ax, radial));

    items(end+1) = struct('name', 'Axial Dose Profile (normalized)', ...
        'draw', @(ax) draw_normalized_profile(ax, axial.axis, axialNorm, ...
            'Z position (cm)', 'Axial Dose Profile (normalized to peak)'));

    items(end+1) = struct('name', 'Radial Dose Profile (normalized)', ...
        'draw', @(ax) draw_normalized_profile(ax, radial.axis, radialNorm, ...
            'R position (cm)', 'Radial Dose Profile (normalized to peak)'));

    items(end+1) = struct('name', 'R-Z Dose Map', ...
        'draw', @(ax) draw_rz_map(ax, rz));

    items(end+1) = struct('name', 'R-Z Dose Map (normalized)', ...
        'draw', @(ax) draw_rz_map_norm(ax, rz, rzNorm));

    if ~isempty(energyComponents) && any([energyComponents.available])
        items(end+1) = struct('name', 'Energy Deposition by Component', ...
            'draw', @(ax) draw_energy_components(ax, energyComponents));

        available = energyComponents([energyComponents.available]);
        for k = 1:numel(available)
            c = available(k);
            items(end+1) = struct('name', sprintf('Energy Dep: %s', pretty_component_name(c.component)), ...
                'draw', @(ax) draw_single_energy_component(ax, c)); %#ok<AGROW>
        end
    end

    if brems.available
        items(end+1) = struct('name', 'Bremsstrahlung Spectrum', ...
            'draw', @(ax) draw_brems_spectrum(ax, brems, beamEnergyMeV));
    end

    if energySummary.available
        items(end+1) = struct('name', 'Target Energy Summary', ...
            'draw', @(ax) draw_target_energy_summary(ax, energySummary));
    end
end

function draw_axial_profile(ax, axial)
    h1 = plot(ax, axial.axis, axial.dose, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.00 0.45 0.74]);
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'Z position (cm)');
    ylabel(ax, sprintf('%s (%s)', axial.doseLabel, axial.unit));
    title(ax, sprintf('Axial Dose Profile  |  NUI = %.4g (%s)', axial.nui, axial.nuiMethod));
    set_profile_ylim_zero(ax, axial.dose);
    legend(ax, h1, sprintf('Mean = %.4g %s', axial.meanDose, axial.unit), 'Location', 'best', 'TextColor', [0.00 0.45 0.74]);
    style_axes_text(ax, [0.00 0.45 0.74]);
end

function draw_radial_profile(ax, radial)
    h1 = plot(ax, radial.axis, radial.dose, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.85 0.33 0.10]);
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'R position (cm)');
    ylabel(ax, sprintf('%s (%s)', radial.doseLabel, radial.unit));
    title(ax, sprintf('Radial Dose Profile  |  NUI = %.4g (%s)', radial.nui, radial.nuiMethod));
    set_profile_ylim_zero(ax, radial.dose);
    legend(ax, h1, sprintf('Mean = %.4g %s', radial.meanDose, radial.unit), 'Location', 'best', 'TextColor', [0.85 0.33 0.10]);
    style_axes_text(ax, [0.85 0.33 0.10]);
end

function draw_normalized_profile(ax, axisValues, normValues, xlabelText, titleText)
    plot(ax, axisValues, normValues, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.47 0.67 0.19]);
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, xlabelText);
    ylabel(ax, 'Normalized dose (peak = 1)');
    title(ax, titleText);
    ylim(ax, [0 1.10]);
    style_axes_text(ax, [0.47 0.67 0.19]);
end

function draw_rz_map(ax, rz)
    imagesc(ax, rz.z_axis, rz.r_axis, rz.dose);
    set(ax, 'YDir', 'normal');
    set(ax, 'Color', [0.05 0.05 0.20]);
    axis(ax, 'tight');
    colormap(ax, 'turbo');
    cb = colorbar(ax);
    cb.Color = [0.00 0.45 0.74];
    xlabel(ax, 'Z position (cm)');
    ylabel(ax, 'R position (cm)');
    title(ax, sprintf('R-Z Dose Map (%s, %s)', rz.doseLabel, rz.unit));
    style_axes_text(ax, [0.00 0.45 0.74]);
end

function draw_rz_map_norm(ax, rz, rzNorm)
    rzNormGrid = reshape(rzNorm, size(rz.dose));
    imagesc(ax, rz.z_axis, rz.r_axis, rzNormGrid);
    set(ax, 'YDir', 'normal');
    set(ax, 'Color', [0.05 0.05 0.20]);
    axis(ax, 'tight');
    colormap(ax, 'turbo');
    cb = colorbar(ax);
    cb.Label.String = 'Normalized dose (peak = 1)';
    cb.Color = [0.00 0.45 0.74];
    xlabel(ax, 'Z position (cm)');
    ylabel(ax, 'R position (cm)');
    title(ax, 'R-Z Dose Map (normalized to peak)');
    caxis(ax, [0 1]); %#ok<CAXIS>
    style_axes_text(ax, [0.00 0.45 0.74]);
end

function draw_energy_components(ax, components)
    available = components([components.available]);
    [~, order] = sort([available.percentIncident], 'descend');
    available = available(order);
    labels = arrayfun(@(c) pretty_component_name(c.component), available, 'UniformOutput', false);
    values = [available.percentIncident];

    bar(ax, values, 'FaceColor', [0.46 0.05 0.68], 'EdgeColor', 'none');
    grid(ax, 'on'); box(ax, 'on');
    ylabel(ax, 'Energy deposited (% of incident)');
    title(ax, 'Energy Deposition by Component');
    set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels);
    xtickangle(ax, 35);
    set_profile_ylim_zero(ax, values);
    for k = 1:numel(available)
        text(ax, k, values(k), sprintf('%.4g%%\n%.4g %s', ...
            available(k).percentIncident, available(k).totalMeV, available(k).unit), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'Color', [0.46 0.05 0.68], 'FontWeight', 'bold');
    end
    style_axes_text(ax, [0.46 0.05 0.68]);
end

function draw_single_energy_component(ax, c)
    label = pretty_component_name(c.component);
    bar(ax, 1, c.percentIncident, 'FaceColor', [0.46 0.05 0.68], 'EdgeColor', 'none');
    grid(ax, 'on'); box(ax, 'on');
    ylabel(ax, 'Energy deposited (% of incident)');
    title(ax, sprintf('%s Energy Deposition', label));
    set(ax, 'XTick', 1, 'XTickLabel', {label});
    set_profile_ylim_zero(ax, c.percentIncident);
    yText = max(c.percentIncident * 1.05, 0.05);
    text(ax, 1, yText, sprintf('Total: %.6g %s\nIncident: %.6g%%\nScored deposits: %.6g%%', ...
        c.totalMeV, c.unit, c.percentIncident, c.percentScored), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Color', [0.46 0.05 0.68], 'FontWeight', 'bold');
    style_axes_text(ax, [0.46 0.05 0.68]);
end

function draw_brems_spectrum(ax, brems, beamEnergyMeV)
    bar(ax, brems.binCenters, brems.totalCounts, 1.0, 'FaceColor', [0.30 0.45 0.69], 'EdgeColor', 'none');
    hold(ax, 'on');
    if any(brems.forwardCounts > 0)
        bar(ax, brems.binCenters, brems.forwardCounts, 1.0, 'FaceColor', [0.85 0.33 0.10], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        legend(ax, {'All escaping \gamma', 'Forward (+Z) escaping \gamma'}, 'Location', 'best', 'TextColor', [0.30 0.45 0.69]);
    end
    xline(ax, beamEnergyMeV, '--', sprintf('Beam energy = %.3g MeV', beamEnergyMeV), ...
        'Color', [0.00 0.60 0.30], 'LineWidth', 1.5, 'LabelColor', [0.00 0.60 0.30], ...
        'LabelHorizontalAlignment', 'left', 'LabelOrientation', 'horizontal');
    if isfinite(brems.endpointMeV)
        xline(ax, brems.endpointMeV, ':', sprintf('Max observed = %.3g MeV', brems.endpointMeV), ...
            'Color', [0.85 0.10 0.10], 'LineWidth', 1.5, 'LabelColor', [0.85 0.10 0.10], ...
            'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal');
    end
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'Photon energy (MeV)');
    ylabel(ax, 'Weighted counts');
    title(ax, sprintf('Bremsstrahlung Escape Spectrum (%d file(s), %d \\gamma)', brems.filesUsed, brems.gammaCount));
    hold(ax, 'off');
    style_axes_text(ax, [0.30 0.45 0.69]);
end

function draw_target_energy_summary(ax, energySummary)
    axis(ax, 'off');
    lines = {
        sprintf('Target energy deposition: %.6g %s', energySummary.targetMeV, energySummary.unit)
        sprintf('Incident beam energy: %.6g MeV per primary', energySummary.incidentMeV / max(energySummary.primaryElectrons, 1))
        sprintf('Primary electrons simulated: %.6g', energySummary.primaryElectrons)
        sprintf('Total incident energy: %.6g MeV', energySummary.incidentMeV)
        sprintf('Target share of incident: %.4f %%', energySummary.percent)
    };
    text(ax, 0.05, 0.8, lines, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontSize', 12, 'FontName', 'FixedWidth', 'Color', [0.00 0.30 0.60], 'FontWeight', 'bold');
    title(ax, 'Target Energy Summary', 'Color', [0.00 0.30 0.60]);
end

function open_analysis_panel(plotItems, axial, radial, energyComponents)
% Open a single figure containing a list box of available plots and an
% axes that updates to show whichever plot is selected, instead of
% generating and saving a separate figure file per plot.

    if isempty(plotItems)
        warning('No plots available to display.');
        return;
    end

    fig = figure('Name', 'ND3 Profile Analysis', 'NumberTitle', 'off', ...
        'Units', 'normalized', 'Position', [0.08 0.08 0.84 0.80], ...
        'Color', 'w');

    listPanel = uipanel('Parent', fig, 'Units', 'normalized', ...
        'Position', [0.00 0.00 0.22 1.00], 'Title', 'Available plots');

    plotNames = {plotItems.name};
    listBox = uicontrol('Parent', listPanel, 'Style', 'listbox', ...
        'Units', 'normalized', 'Position', [0.04 0.06 0.92 0.90], ...
        'String', plotNames, 'FontSize', 10, 'Value', 1);

    ax = axes('Parent', fig, 'Units', 'normalized', ...
        'Position', [0.30 0.10 0.66 0.82]);

    statusText = uicontrol('Parent', fig, 'Style', 'text', ...
        'Units', 'normalized', 'Position', [0.24 0.005 0.74 0.05], ...
        'HorizontalAlignment', 'left', 'FontSize', 10, 'FontWeight', 'bold', ...
        'ForegroundColor', [0.00 0.30 0.60], 'BackgroundColor', 'w', ...
        'String', summary_status_string(axial, radial, energyComponents));

    listBox.Callback = @(src, evt) render_selected_plot(ax, plotItems, src.Value);

    % Render the first plot immediately.
    render_selected_plot(ax, plotItems, 1);

    drawnow;
end

function render_selected_plot(ax, plotItems, idx)
    idx = max(1, min(idx, numel(plotItems)));
    cla(ax, 'reset');
    try
        plotItems(idx).draw(ax);
    catch ME
        cla(ax, 'reset');
        axis(ax, 'off');
        text(ax, 0.05, 0.5, sprintf('Could not render "%s":\n%s', ...
            plotItems(idx).name, ME.message), 'Units', 'normalized', ...
            'Color', 'r', 'FontSize', 11);
    end
end

function str = summary_status_string(axial, radial, energyComponents)
    nComponents = 0;
    if ~isempty(energyComponents)
        nComponents = sum([energyComponents.available]);
    end

    str = sprintf('Axial NUI = %.4g (%s)   |   Radial NUI = %.4g (%s)   |   Energy components loaded: %d', ...
        axial.nui, axial.nuiMethod, radial.nui, radial.nuiMethod, nComponents);
end