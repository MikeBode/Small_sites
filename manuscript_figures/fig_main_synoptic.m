function fig_main_synoptic(rarity_name)
% If no argument is given, default to the single-fragment definition.
if nargin < 1 || isempty(rarity_name)
    rarity_name = 'occ1';
end

% Map of the location of all datasets.
%
%   Reads the study metadata then plots:
%       x = % of all patches below the landscape mean patch area
%       y = % of single-patch rare species whose patch is below the landscape mean
%   Points coloured by taxon, tier-1 studies labelled, 1:1 line for reference.

here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'output');
if ~exist(outdir, 'dir'); mkdir(outdir); end

repo_root = fileparts(here);
csv_name = sprintf('table_rarity_%s.csv', rarity_name);
claim1 = read_csv(fullfile(here, 'output', csv_name));
tier1_ids = {'Bovendorp_2017', 'Arroyo_Rodriguez_2009', 'Pinho_2019a'};

% Extract the two columns we need
study_ids = claim1.study_id;
pct_patches = str2doublevec(claim1.pct_patches_below_mean);
pct_rare = str2doublevec(claim1.pct_rare_below_mean);

% Lookup taxon per study from the claim1 table (it already has it)
taxa = claim1.taxon;

% Group taxa: mammals include bats
taxa_grouped = taxa;
for i = 1:numel(taxa_grouped)
    t = lower(taxa_grouped{i});
    if strncmp(t, 'chiroptera', 10) || strncmp(t, 'mammalia', 8)
        taxa_grouped{i} = 'Mammalia';
    elseif strncmp(t, 'insecta', 7)
        taxa_grouped{i} = 'Insecta';
    elseif strncmp(t, 'aves', 4)
        taxa_grouped{i} = 'Aves';
    elseif strncmp(t, 'plantae', 7) || strncmp(t, 'plants', 6)
        taxa_grouped{i} = 'Plantae';
    elseif strncmp(t, 'mollusca', 8)
        taxa_grouped{i} = 'Mollusca';
    end
end

unique_taxa = unique(taxa_grouped);
n_taxa = numel(unique_taxa);
taxon_colours = get_taxon_palette(n_taxa);

fig = figure(1); clf
ax = axes('position',[0.07 0.1 0.9 0.58]);
hold(ax, 'on'); box(ax, 'on');

% 1:1 reference line
plot(ax, [0 110], [0 110], '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2);

% Plot points by taxon
handles = zeros(n_taxa, 1);
for t = 1:n_taxa
    sel = strcmp(taxa_grouped, unique_taxa{t});
    if ~any(sel); continue; end

    % Separate tier 1 (diamond) from tier 2 (circle)
    tier1_mask = false(sum(sel), 1);
    sel_idx = find(sel);
    for k = 1:numel(sel_idx)
        if any(strcmp(study_ids{sel_idx(k)}, tier1_ids))
            tier1_mask(k) = true;
        end
    end

    x_all = pct_patches(sel_idx);
    y_all = pct_rare(sel_idx);

    % Tier 2 (smaller circles)
    if any(~tier1_mask)
        scatter(ax, x_all(~tier1_mask), y_all(~tier1_mask), 35, ...
            taxon_colours(t, :), 'filled', 'MarkerEdgeColor', taxon_colours(t, :).^2, 'LineWidth', 1.2);
    end
    % Tier 1 (larger diamonds), plotted last so they sit on top
    if any(tier1_mask)
        scatter(ax, x_all(tier1_mask), y_all(tier1_mask), 35, ...
            taxon_colours(t, :), 'filled', 'MarkerEdgeColor', taxon_colours(t, :).^2, 'LineWidth', 1.2);
    end

    % Handle for the legend (invisible off-plot marker)
    handles(t) = scatter(ax, -10, -10, 110, taxon_colours(t, :), 'filled', 'MarkerEdgeColor', taxon_colours(t, :).^2, 'LineWidth', 1.2);
end

% Label tier-1 studies
for k = 1:numel(tier1_ids)
    ix = find(strcmp(study_ids, tier1_ids{k}), 1);
    if ~isempty(ix)
        label = short_label(tier1_ids{k});
        text(ax, pct_patches(ix) + 0.75, pct_rare(ix), label, ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
end

set(ax, 'FontSize', 11); axis square
xlabel(ax, 'Percent small patches', 'FontSize', 13);
ylabel(ax, 'Percent rare species in small patches', 'FontSize', 13);
xlim(ax, [50 101]); ylim(ax, [50 101]);

legend(handles(handles > 0), unique_taxa(handles > 0), 'Location', 'northwest', ...
    'FontSize', 11, 'Box', 'on','color',0.95.*ones(1,3),'edgecolor',0.95.*ones(1,3));

ax = axes('position',[0.05 0.69 0.9 0.275]);
map_study_locations
axis off

outpath = fullfile(outdir, sprintf('fig_main_synoptic_%s.tiff', rarity_name));
Make_TIFF(outpath,[0 0 30 18]*0.75,'-r550')

end


function palette = get_taxon_palette(n)
% A colourblind-friendly palette with strong contrast, up to 8 categories.
base = [
    0.12 0.47 0.71;   % blue
    0.89 0.10 0.11;   % red
    0.20 0.63 0.17;   % green
    0.99 0.50 0.00;   % orange
    0.42 0.24 0.60;   % purple
    0.55 0.34 0.16;   % brown
    0.65 0.65 0.10;   % olive
    0.10 0.65 0.60;   % teal
];
palette = base(1:min(n, size(base, 1)), :);
end


function tbl = read_csv(filepath)
% Read a plain CSV (no # comment rows) into a struct of column arrays.
fid = fopen(filepath, 'r');
header = strsplit(fgetl(fid), ',');
n_cols = numel(header);
cols = cell(n_cols, 1);
for i = 1:n_cols; cols{i} = {}; end
while true
    line = fgetl(fid);
    if ~ischar(line) || isempty(line); break; end
    parts = strsplit(line, ',', 'CollapseDelimiters', false);
    for i = 1:min(numel(parts), n_cols)
        cols{i}{end+1, 1} = parts{i};
    end
end
fclose(fid);
tbl = struct();
for i = 1:n_cols
    key = strrep(strtrim(header{i}), ' ', '_');
    tbl.(key) = cols{i};
end
end


function v = str2doublevec(c)
v = zeros(numel(c), 1);
for i = 1:numel(c)
    x = str2double(c{i});
    if isnan(x)
        v(i) = NaN;
    else
        v(i) = x;
    end
end
end


function out = short_label(study_id)
% Produce a short, human-friendly label from a study_id like "Bovendorp_2017".
parts = strsplit(study_id, '_');
out = strrep(parts{1}, '\_', ''); % just the author surname
end
