function run_setcover_and_table()

% RUN_SETCOVER_AND_TABLE  Produce numeric summary tables for the manuscript.
%
%   Writes three CSVs to output/, one per rarity definition:
%       table_rarity_occ1.csv                 rare = species occupying exactly 1 patch
%       table_rarity_occ2.csv                 rare = species occupying 1 or 2 patches
%       table_rarity_bottom25pct_species.csv  rare = the bottom 25% of species ranked by occupancy; ties at the threshold are included

here = fileparts(mfilename('fullpath'));
data_dir = fullfile(fileparts(here), 'data');
outdir = fullfile(here, 'output');
if ~exist(outdir, 'dir'); mkdir(outdir); end

% Compile a list of the studies we're going to analyse
files = dir(fullfile(data_dir, '*.csv'));
files = files(~strcmp({files.name}, 'study_metadata.csv'));
n = numel(files);
studies = struct([]);
for f = 1:n
    studies(end+1).s = load_study(fullfile(data_dir, files(f).name));
end

% Set cover doesn't depend on rarity definition, compute once per study
setcov = struct([]);
for i = 1:numel(studies)
    s = studies(i).s;
    setcov(i).sel_area = greedy_setcover_area(s.abundance_matrix, s.patch_area);
end

definitions = { ...
    struct('name', 'occ1', 'fn', @(occ) find(occ == 1)), ...
    struct('name', 'occ2', 'fn', @(occ) find(occ <= 2)), ...
    struct('name', 'bottom25pct_species', 'fn', @rare_by_bottom_species_pct) };

for d = 1:numel(definitions)
    def = definitions{d};
    fprintf('Building table for rarity definition: %s\n', def.name);
    rows = cell(numel(studies), 1);
    for i = 1:numel(studies)
        s = studies(i).s;
        occ = sum(s.abundance_matrix > 0, 1);
        rare_idx = def.fn(occ);
        rows{i} = build_row(s, rare_idx, setcov(i));
    end
    outpath = fullfile(outdir, sprintf('table_rarity_%s.csv', def.name));
    write_csv(outpath, table_header(), rows);
    fprintf('  wrote %s\n', outpath);
end

fprintf('Done.\n');
end


function rare_idx = rare_by_bottom_species_pct(occupancy)
% RARE_BY_BOTTOM_SPECIES_PCT  Classify the bottom X% of species by occupancy as rare.
%
% Species are ranked by occupancy (number of fragments occupied). The N
% species with the lowest occupancy values are treated as rare, where N is
% ceil(X/100 * total number of species). Ties at the occupancy threshold are
% included so that no rare species is arbitrarily excluded.

n_species = numel(occupancy);
if n_species == 0
    rare_idx = [];
    return;
end

target_n = max(1, ceil(0.25 * n_species));

% Sort species occupancies in ascending order and take the target_n lowest.
% Because ties may straddle the target_n boundary, include all species whose
% occupancy is <= the occupancy of the target_n-th species in the sorted list.
sorted_occ = sort(occupancy(:), 'ascend');
threshold = sorted_occ(target_n);
rare_idx = find(occupancy <= threshold);

end



function row = build_row(s, rare_idx, setcov)
% One row per study for the current rarity definition.

area = s.patch_area;
log_area = log(area);
mean_area = mean(area);
median_area = median(area);
total_area = sum(area);
pct_patches_below_mean = 100 * mean(area <= mean_area);
n_rare = numel(rare_idx);

if n_rare > 0
    species_mean_area = zeros(1, n_rare);
    rare_log_areas_pooled = zeros(1, 0);
    for k = 1:n_rare
        occupied = find(s.abundance_matrix(:, rare_idx(k)) > 0);
        species_mean_area(k) = mean(area(occupied));
        rare_log_areas_pooled = [rare_log_areas_pooled, log_area(occupied)']; %#ok<AGROW>
    end
    pct_rare_below_mean = 100 * mean(species_mean_area <= mean_area);
    pct_rare_below_median = 100 * mean(species_mean_area <= median_area);

    try
        [~, ks_p, ks_d] = kstest2(log_area, rare_log_areas_pooled);
    catch
        ks_p = NaN; ks_d = NaN;
    end
else
    pct_rare_below_mean = NaN;
    pct_rare_below_median = NaN;
    ks_d = NaN; ks_p = NaN;
end

sel_a = setcov.sel_area;
areas_a = area(sel_a);
n_patches_minarea = numel(sel_a);
total_area_minarea = sum(areas_a);
pct_minarea_below_mean = 100 * mean(areas_a <= mean_area);
pct_minarea_below_median = 100 * mean(areas_a <= median_area);

row = { ...
    s.study_id, ...
    csv_safe(s.taxon), ...
    csv_safe(s.region), ...
    numel(area), ...
    total_area, ...
    size(s.abundance_matrix, 2), ...
    n_rare, ...
    mean_area, ...
    median_area, ...
    pct_patches_below_mean, ...
    pct_rare_below_mean, ...
    pct_rare_below_median, ...
    ks_d, ...
    ks_p, ...
    n_patches_minarea, ...
    total_area_minarea, ...
    pct_minarea_below_mean, ...
    pct_minarea_below_median };
end


function h = table_header()
h = { ...
    'study_id', ...
    'taxon', ...
    'region', ...
    'n_fragments', ...
    'total_patch_area_ha', ...
    'n_species', ...
    'n_rare', ...
    'mean_patch_area_ha', ...
    'median_patch_area_ha', ...
    'pct_patches_below_mean', ...
    'pct_rare_below_mean', ...
    'pct_rare_below_median', ...
    'ks_D', ...
    'ks_p', ...
    'minarea_n_patches', ...
    'minarea_total_area_ha', ...
    'minarea_pct_below_mean', ...
    'minarea_pct_below_median' };
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
% Load a particular study into the struct "study"

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

function out = csv_safe(s)
if ~ischar(s); out = s; return; end
out = strrep(s, ',', ' ');
end

function write_csv(path, header, rows)
fid = fopen(path, 'w');
fprintf(fid, '%s', header{1});
for j = 2:numel(header); fprintf(fid, ',%s', header{j}); end
fprintf(fid, '\n');
for i = 1:numel(rows)
    r = rows{i};
    for j = 1:numel(r)
        if j > 1; fprintf(fid, ','); end
        v = r{j};
        if ischar(v)
            fprintf(fid, '%s', v);
        elseif isnan(v)
            fprintf(fid, 'NA');
        elseif v == round(v) && abs(v) < 1e12
            fprintf(fid, '%d', v);
        else
            fprintf(fid, '%.4g', v);
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
end
