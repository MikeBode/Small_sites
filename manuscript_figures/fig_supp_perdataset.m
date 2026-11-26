function fig_supp_perdataset()
% FIG_SUPP_PERDATASET  Supplementary per-dataset figures for the 13 non-tier-1 studies.
%
%   Two-row layout matching the main-text figure. Up to 5 datasets per page.
%
%   Output: manuscript_figures/output/fig_supp_perdataset_page{1,2,3}.png

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
warning('off', 'MATLAB:print:ContentTypeImageSuggested');

here = fileparts(mfilename('fullpath'));
data_dir = fullfile(fileparts(here), 'data');
outdir = fullfile(here, 'output');
if ~exist(outdir, 'dir'); mkdir(outdir); end

tier1_ids = {'Bovendorp_2017', 'Arroyo_Rodriguez_2009', 'Pinho_2019a'};

files = dir(fullfile(data_dir, '*.csv'));
files = files(~strcmp({files.name}, 'study_metadata.csv'));

studies = struct('study_id', {}, 'taxon', {}, 'region', {}, ...
                 'patch_area', {}, 'abundance_matrix', {}, 'fragment_id', {});
for f = 1:numel(files)
    s = load_study(fullfile(data_dir, files(f).name));
    if any(strcmp(s.study_id, tier1_ids)); continue; end
    studies(end+1) = s; %#ok<AGROW>
end

n_species_each = arrayfun(@(x) size(x.abundance_matrix, 2), studies);
[~, order] = sort(n_species_each, 'descend');
studies = studies(order);

n_studies = numel(studies);
n_per_page = 5;
n_pages = ceil(n_studies / n_per_page);

for page = 1:n_pages
    fig = figure('Position', [50 50 1500 800], 'Color', 'w');
    set(fig, 'PaperPositionMode', 'auto');
    ix_start = (page - 1) * n_per_page + 1;
    ix_end = min(page * n_per_page, n_studies);
    n_cols = ix_end - ix_start + 1;

    for c = 1:n_cols
        s = studies(ix_start + c - 1);
        draw_column(fig, s, c, n_cols, 2);
    end

    outpath = fullfile(outdir, sprintf('fig_supp_perdataset_page%d.png', page));
    print(fig, outpath, '-dpng', '-r200');
    close(fig);
    fprintf('  wrote %s\n', outpath);
end
end


function draw_column(fig, study, col_idx, n_cols, n_rows)
area = study.patch_area;
abundance_matrix = study.abundance_matrix;
mean_area = mean(area);
log_area = log(area);

ax = subplot(n_rows, n_cols, col_idx);
plot_combined_overlay(ax, area, abundance_matrix, log_area);

ax = subplot(n_rows, n_cols, n_cols + col_idx);
plot_cumulative_curve(ax, area, abundance_matrix, mean_area);
end


function plot_combined_overlay(ax, area, abundance_matrix, log_area)
occupancy = sum(abundance_matrix > 0, 1);
rare_idx = find(occupancy == 1);
rare_log_areas = zeros(numel(rare_idx), 1);
for k = 1:numel(rare_idx)
    p = find(abundance_matrix(:, rare_idx(k)) > 0, 1);
    rare_log_areas(k) = log_area(p);
end

selected_min_area = greedy_setcover_area(abundance_matrix, area);
log_area_min_area = log_area(selected_min_area);

n_bins = 15;
edges = linspace(min(log_area), max(log_area), n_bins + 1);
centres = (edges(1:end-1) + edges(2:end)) / 2;
all_counts = histc(log_area, edges); all_counts = all_counts(1:end-1);
rare_counts = histc(rare_log_areas, edges); rare_counts = rare_counts(1:end-1);
min_area_counts = histc(log_area_min_area, edges); min_area_counts = min_area_counts(1:end-1);

all_density = all_counts / max(sum(all_counts), 1);
rare_density = rare_counts / max(sum(rare_counts), 1);
min_area_density = min_area_counts / max(sum(min_area_counts), 1);

hold(ax, 'on'); box(ax, 'on');
bar(ax, centres, all_density, 1.0, 'FaceColor', 'none', ...
    'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.0);
bar(ax, centres, rare_density, 0.75, 'FaceColor', [0.84 0.15 0.15], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6);
bar(ax, centres, min_area_density, 0.45, 'FaceColor', [0.10 0.30 0.65], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.7);
set(ax, 'FontSize', 8);
xlabel(ax, '$\ln(\mathrm{patch\ area})$', 'FontSize', 10, 'Interpreter', 'latex');
ylabel(ax, 'Relative frequency', 'FontSize', 10, 'Interpreter', 'latex');
end


function plot_cumulative_curve(ax, area, abundance_matrix, mean_area, show_legend)
% Row 2: cumulative fraction of single-patch rare species accumulated as the
% patch-size threshold sweeps from 0 to 10 percent of the landscape mean.
% Nested null intervals show the 99%, 95%, and 67% central ranges.

if nargin < 5
    show_legend = false;
end

cutoffs_pct = 0:0.2:10;
cutoff_areas = (cutoffs_pct / 100) * mean_area;
n_cutoffs = numel(cutoffs_pct);

occupancy = sum(abundance_matrix > 0, 1);
rare_idx = find(occupancy == 1);
n_rare = numel(rare_idx);
if n_rare == 0
    axis(ax, 'off');
    return;
end

% Empirical cumulative curve
min_area_rare = zeros(1, n_rare);
for k = 1:n_rare
    p = find(abundance_matrix(:, rare_idx(k)) > 0, 1);
    min_area_rare(k) = area(p);
end

emp_curve = zeros(1, n_cutoffs);
for c = 1:n_cutoffs
    emp_curve(c) = 100 * sum(min_area_rare <= cutoff_areas(c)) / n_rare;
end

% Null model: independently permute each species' abundance vector across fragments
n_iter = 2000;
rng(42, 'twister');
null_curves = zeros(n_iter, n_cutoffs);
n_fragments = size(abundance_matrix, 1);

for it = 1:n_iter
    permuted = zeros(size(abundance_matrix));
    for sp = 1:size(abundance_matrix, 2)
        permuted(:, sp) = randsample(abundance_matrix(:, sp), n_fragments);
    end

    perm_occupancy = sum(permuted > 0, 1);
    perm_rare_idx = find(perm_occupancy == 1);

    if isempty(perm_rare_idx)
        continue;
    end

    perm_min = zeros(1, numel(perm_rare_idx));
    for k = 1:numel(perm_rare_idx)
        p = find(permuted(:, perm_rare_idx(k)) > 0, 1);
        perm_min(k) = area(p);
    end

    for c = 1:n_cutoffs
        null_curves(it, c) = 100 * sum(perm_min <= cutoff_areas(c)) / numel(perm_rare_idx);
    end
end

% Nested central intervals
null_mean = mean(null_curves, 1);

null_99_lo = prctile(null_curves, 0.5, 1);
null_99_hi = prctile(null_curves, 99.5, 1);

null_95_lo = prctile(null_curves, 2.5, 1);
null_95_hi = prctile(null_curves, 97.5, 1);

% Central 67% interval, approximately +/- 1 SD for a Gaussian distribution
null_67_lo = prctile(null_curves, 16.5, 1);
null_67_hi = prctile(null_curves, 83.5, 1);

hold(ax, 'on');
box(ax, 'on');

% Plot widest interval first, narrowest last
h99 = patch(ax, [cutoffs_pct, fliplr(cutoffs_pct)], ...
    [null_99_lo, fliplr(null_99_hi)], ...
    [0.88 0.88 0.88], ...
    'FaceAlpha', 1.0, 'EdgeColor', 'none');

h95 = patch(ax, [cutoffs_pct, fliplr(cutoffs_pct)], ...
    [null_95_lo, fliplr(null_95_hi)], ...
    [0.72 0.72 0.72], ...
    'FaceAlpha', 1.0, 'EdgeColor', 'none');

h67 = patch(ax, [cutoffs_pct, fliplr(cutoffs_pct)], ...
    [null_67_lo, fliplr(null_67_hi)], ...
    [0.52 0.52 0.52], ...
    'FaceAlpha', 1.0, 'EdgeColor', 'none');

hmean = plot(ax, cutoffs_pct, null_mean, '--', ...
    'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);

hemp = plot(ax, cutoffs_pct, emp_curve, 'r-', 'LineWidth', 2.2);

set(ax, 'FontSize', 9);
xlabel(ax, 'Patch size threshold (\% of landscape mean)', ...
    'FontSize', 11, 'Interpreter', 'latex');
ylabel(ax, '\% of rare species accumulated', ...
    'FontSize', 11, 'Interpreter', 'latex');
xlim(ax, [0 10]);
ylim(ax, [0 100]);

if show_legend
    legend(ax, [hemp hmean h67 h95 h99], ...
        {'observed', 'null mean', '67\% null interval', ...
         '95\% null interval', '99\% null interval'}, ...
        'Location', 'southeast', ...
        'FontSize', 8, ...
        'Box', 'off', ...
        'Interpreter', 'latex');
end
end



function selected = greedy_setcover_area(abundance_matrix, area)
n_fragments = size(abundance_matrix, 1);
n_species = size(abundance_matrix, 2);
uncovered = true(1, n_species);
selected = false(1, n_fragments);
while any(uncovered)
    marginal = zeros(1, n_fragments);
    for p = 1:n_fragments
        if selected(p); continue; end
        marginal(p) = sum(abundance_matrix(p, :) > 0 & uncovered);
    end
    efficiency = marginal ./ max(area(:)', eps);
    [best_score, best_p] = max(efficiency);
    if best_score == 0; break; end
    selected(best_p) = true;
    uncovered = uncovered & ~(abundance_matrix(best_p, :) > 0);
end
selected = find(selected);
end


function study = load_study(filepath)
meta = struct();
fid = fopen(filepath, 'r');
while true
    line = fgetl(fid);
    if ~ischar(line) || isempty(line) || line(1) ~= '#'; break; end
    body = strtrim(line(2:end));
    colon = strfind(body, ':');
    if ~isempty(colon)
        key = strrep(strtrim(body(1:colon(1)-1)), ' ', '_');
        meta.(key) = strtrim(body(colon(1)+1:end));
    end
end
fclose(fid);

is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if is_octave
    [study.fragment_id, study.patch_area, study.abundance_matrix] = read_table_manual(filepath);
else
    t = readtable(filepath, 'CommentStyle', '#', 'VariableNamingRule', 'preserve');
    study.fragment_id = cellstr(string(t{:, 1}));
    study.patch_area = t{:, 2};
    study.abundance_matrix = t{:, 3:end};
end
study.study_id = get_field_or(meta, 'study_id', '');
study.taxon = get_field_or(meta, 'taxon', '');
study.region = get_field_or(meta, 'region', '');
end


function [fragment_id, patch_area, abundance_matrix] = read_table_manual(filepath)
fid = fopen(filepath, 'r');
header_seen = false; fragment_id = {}; areas = []; rows = {};
while true
    line = fgetl(fid);
    if ~ischar(line); break; end
    if isempty(line) || line(1) == '#'; continue; end
    if ~header_seen; header_seen = true; continue; end
    parts = strsplit(line, ',');
    fragment_id{end+1, 1} = parts{1}; %#ok<AGROW>
    nums = str2double(parts(2:end));
    areas(end+1, 1) = nums(1); %#ok<AGROW>
    rows{end+1, 1} = nums(2:end); %#ok<AGROW>
end
fclose(fid);
patch_area = areas;
abundance_matrix = cell2mat(rows);
end


function v = get_field_or(s, field, default)
if isfield(s, field); v = s.(field); else; v = default; end
end
