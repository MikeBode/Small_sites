function run_all_figures()
% RUN_ALL_FIGURES  Rebuild every manuscript figure and table from the occurrence data.
%
%   Clears the output/ directory, then produces:
%       output/fig_main_perdataset.png
%       output/fig_main_synoptic.png
%       output/fig_supp_perdataset_page{1,2,3}.png
%       output/table_rarity_occ1.csv
%       output/table_rarity_occ2.csv
%       out put/table_rarity_bottom5pct_occupancy.csv

here = fileparts(mfilename('fullpath'));
addpath(here);
outdir = fullfile(here, 'output');

% Clear the output directory
if exist(outdir, 'dir')
    old_files = dir(fullfile(outdir, '*'));
    for k = 1:numel(old_files)
        if old_files(k).isdir; continue; end
        delete(fullfile(outdir, old_files(k).name));
    end
    fprintf('Cleared existing outputs in %s\n', outdir);
else, mkdir(outdir);
end

run_setcover_and_table();
fig_main_perdataset();
fig_main_synoptic('occ1');
fig_main_synoptic('occ2');
fig_main_synoptic('bottom25pct_species');
fig_supp_perdataset();
end
