function fig_main_perdataset()

% Produces the main-text figure: 2-panel column figure for the three tier-1 datasets.
%
%   Layout is 2 rows x 3 columns. Each column is one tier-1 dataset.
%
%       Row 1: log(patch area) distribution of all patches (grey open bars) with
%              rare-species patches (red fill) and minimum-area set-cover patches
%              (blue fill) both superimposed.
%       Row 2: cumulative fraction of single-patch rare species accumulated as the
%              small-patch threshold sweeps from 0 to 10 percent of the landscape
%              mean, with empirical curve and null 95% CI band.
%
%   Output: manuscript_figures/output/fig_main_perdataset.png

% Identify the best datasets
tier1_ids = {'Bovendorp_2017', 'Arroyo_Rodriguez_2009', 'Pinho_2019a'};

% Panel label settings
panel_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
label_cfg.x_frac = 0.96; label_cfg.y_frac = 0.96;
label_cfg.font_size = 14; label_cfg.font_weight = 'bold'; label_cfg.interpreter = 'latex';
label_cfg.h_align = 'right'; label_cfg.v_align = 'top';

% Where to save the outputs
here = fileparts(mfilename('fullpath'));
data_dir = fullfile(fileparts(here), 'data');
outdir = fullfile(here, 'output');
if ~exist(outdir, 'dir'); mkdir(outdir); end

% Prep the figure
fig = figure(10); clf

for col = 1:numel(tier1_ids)
    study = load_study(fullfile(data_dir, [tier1_ids{col} '.csv']));
    draw_column(fig, study, col, numel(tier1_ids), 2, panel_labels, label_cfg);
end

Make_TIFF('output/fig_main_perdataset.tiff',[0 0 30 18]*1,'-r450')
end

function draw_column(fig, study, col_idx, n_cols, n_rows, panel_labels, label_cfg)
% Draw both panels for one study in one column of the figure grid.

area = study.patch_area;
abundance_matrix = study.abundance_matrix;
mean_area = mean(area);
log_area = log(area);

% Row 1: rare species + min-area set cover overlaid on all patches
ax = subplot(n_rows, n_cols, col_idx);
plot_combined_overlay(ax, area, abundance_matrix, log_area, col_idx == 2);
add_panel_label(ax, panel_labels{col_idx}, label_cfg);

% Row 2: cumulative rare-in-small curve with null CI
ax = subplot(n_rows, n_cols, n_cols + col_idx);
show_legend = (col_idx == 2);
plot_cumulative_curve(ax, area, abundance_matrix, mean_area, show_legend);
add_panel_label(ax, panel_labels{n_cols + col_idx}, label_cfg);
end


function plot_combined_overlay(ax, area, abundance_matrix, log_area, show_legend)
% Row 1: distribution of log(patch area) for all patches (grey open bars),
% with rare-species patches (red) and minimum-area set-cover patches (blue)
% both superimposed on the same underlying patch distribution.

% Find the rare species
occupancy = sum(abundance_matrix > 0, 1);
rare_idx = find(occupancy == 1);

% Find the size of their patches
rare_log_areas = zeros(numel(rare_idx), 1);
for k = 1:numel(rare_idx)
    p = find(abundance_matrix(:, rare_idx(k)) > 0, 1);
    rare_log_areas(k) = log_area(p);
end

% How much land is needed to protect them
selected_min_area = greedy_setcover_area(abundance_matrix, area);
log_area_min_area = log_area(selected_min_area);

n_bins = 15;
edges = linspace(min(log_area), max(log_area), n_bins + 1);
centres = (edges(1:end-1) + edges(2:end)) / 2;
width = centres(2)-centres(1);

% Compute normalized histogram counts for each overlaid patch category.
all_counts = histc(log_area, edges); all_counts = all_counts(1:end-1);
rare_counts = histc(rare_log_areas, edges); rare_counts = rare_counts(1:end-1);
min_area_counts = histc(log_area_min_area, edges); min_area_counts = min_area_counts(1:end-1);

all_density = all_counts / max(sum(all_counts), 1);
rare_density = rare_counts / max(sum(rare_counts), 1);
min_area_density = min_area_counts / max(sum(min_area_counts), 1);

hold(ax, 'on'); box(ax, 'on');
h1 = bar(ax, centres, all_density, 1.0, 'FaceColor', 0.65.*ones(1,3), ...
    'EdgeColor', 0.45.*ones(1,3));
h2 = bar(ax, centres-width*0.15, rare_density, 0.33, 'FaceColor', [0.84 0.15 0.15], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6);
h3 = bar(ax, centres+width*0.15, min_area_density, 0.33, 'FaceColor', [0.10 0.30 0.65], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.7);

xlabel(ax, '$\ln(\mathrm{patch\ area})$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel(ax, 'Relative frequency', 'FontSize', 12, 'Interpreter', 'latex');
ylim([0 0.25])

if show_legend
    legend(ax, [h1 h2 h3], ...
        {'all patches', 'rare species patches', 'protected patches'}, ...
        'Location', 'northwest', 'FontSize', 10, 'Box', 'off', 'Interpreter', 'latex');
end
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

% What proportion of (observed) rare species are found in the smallest % of patches
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

    % Randomly re-allocate each species across patches
    permuted = zeros(size(abundance_matrix));
    for sp = 1:size(abundance_matrix, 2)
        permuted(:, sp) = randsample(abundance_matrix(:, sp), n_fragments);
    end

    % Work out occupancy and rare species
    perm_occupancy = sum(permuted > 0, 1);
    perm_rare_idx = find(perm_occupancy == 1);

    perm_min = zeros(1, numel(perm_rare_idx));
    for k = 1:numel(perm_rare_idx)
        p = find(permuted(:, perm_rare_idx(k)) > 0, 1);
        perm_min(k) = area(p);
    end

    % What proportion of (randomised) rare species are found in the smallest % of patches
    for c = 1:n_cutoffs
        null_curves(it, c) = 100 * sum(perm_min <= cutoff_areas(c)) / numel(perm_rare_idx);
    end
end

% Construct nested central intervals
null_mean = mean(null_curves, 1);
null_99_lo = prctile(null_curves, 0.5, 1); null_99_hi = prctile(null_curves, 99.5, 1);
null_95_lo = prctile(null_curves, 2.5, 1); null_95_hi = prctile(null_curves, 97.5, 1);
null_67_lo = prctile(null_curves, 16.5, 1); null_67_hi = prctile(null_curves, 83.5, 1);
hold(ax, 'on'); box(ax, 'on');

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
         '95\%', '99\%'}, ...
        'Location', 'northwest', ...
        'FontSize', 10, ...
        'Box', 'off', ...
        'Interpreter', 'latex');
end
end



function selected = greedy_setcover_area(abundance_matrix, area)

% Greedy min-area set cover: at each step add the patch offering the most
% currently-uncovered species per hectare, until all species are covered.

n_fragments = size(abundance_matrix, 1);
n_species = size(abundance_matrix, 2);
uncovered = true(1, n_species);
selected = false(1, n_fragments);
while any(uncovered)
    marginal = zeros(1, n_fragments);
    for p = 1:n_fragments
        if selected(p); continue; end % Can't protect a site twice
        marginal(p) = sum(abundance_matrix(p, :) > 0 & uncovered); % How many new species per site?
    end

    % Pro-rata by area
    efficiency = marginal ./ max(area(:)', eps);
    [best_score, best_p] = max(efficiency);
    if best_score == 0; break; end
    selected(best_p(1)) = true;
    uncovered = uncovered & ~(abundance_matrix(best_p, :) > 0); % Any species in that site are now considered protected
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

t = readtable(filepath, 'CommentStyle', '#', 'VariableNamingRule', 'preserve');
study.fragment_id = cellstr(string(t{:, 1}));
study.patch_area = t{:, 2};
study.abundance_matrix = t{:, 3:end};
study.study_id = get_field_or(meta, 'study_id', '');
study.taxon = get_field_or(meta, 'taxon', '');
study.region = get_field_or(meta, 'region', '');
end


function v = get_field_or(s, field, default)
if isfield(s, field); v = s.(field); else; v = default; end
end

function add_panel_label(ax, label_text, cfg)

% ADD_PANEL_LABEL  Add a panel label at a fractional position within the axes.
%
%   cfg.x_frac and cfg.y_frac set the position in axes data coordinates:
%       x_frac = 0 left edge, 1 right edge
%       y_frac = 0 bottom edge, 1 top edge

xl = xlim(ax);
yl = ylim(ax);

if strcmp(get(ax, 'XScale'), 'log')
    xpos = 10^(log10(xl(1)) + cfg.x_frac * (log10(xl(2)) - log10(xl(1))));
else
    xpos = xl(1) + cfg.x_frac * (xl(2) - xl(1));
end

if strcmp(get(ax, 'YScale'), 'log')
    ypos = 10^(log10(yl(1)) + cfg.y_frac * (log10(yl(2)) - log10(yl(1))));
else
    ypos = yl(1) + cfg.y_frac * (yl(2) - yl(1));
end

text(ax, xpos, ypos, label_text, ...
    'FontSize', cfg.font_size, ...
    'FontWeight', cfg.font_weight, ...
    'Interpreter', cfg.interpreter, ...
    'HorizontalAlignment', cfg.h_align, ...
    'VerticalAlignment', cfg.v_align);
end

