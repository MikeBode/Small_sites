function map_study_locations()

% Plot approximate locations of the 16 study landscapes.
%
% Study coordinates are approximate landscape centroids, intended for a
% small global overview map. Bovendorp_2017 is a regional Atlantic Forest
% compilation rather than a single landscape.

% CONFIGURATION
world_data_file = "/Users/bodem@qut.edu.au/Documents/MATLAB/PlotWorldData.mat";
output_dir = "/Users/bodem@qut.edu.au/Work/Simon_Mahood/Datasets/" + "LSS_final/manuscript_figures/output";
output_file = fullfile(output_dir, "figure_study_locations.png");

colour_by_taxon = true;
show_legend = 0;
marker_size = 25;
marker_edge_width = 1;

% STUDY LOCATIONS

study_id = { ...
    'Arroyo_Rodriguez_2009'
    'Banks_Leite_2009'
    'Bovendorp_2017'
    'Brosi_2009'
    'Giladi_2011_b'
    'Kappes_2009'
    'Manu_2007'
    'Pinho_2019a'
    'Pinho_2019c'
    'Raheem_2009'
    'Silva_2016'
    'Soto_2022'
    'Stireman_2014'
    'Struebig_2008'
    'Vasconcelos_2006'};

taxon = { ...
    'Plantae'
    'Aves'
    'Mammalia'
    'Insecta'
    'Plantae'
    'Mollusca'
    'Aves'
    'Plantae'
    'Plantae'
    'Mollusca'
    'Insecta'
    'Mammalia'
    'Insecta'
    'Mammalia'
    'Insecta'};

location = { ...
    'Los Tuxtlas Mexico'
    'Atlantic Forest Sao Paulo Brazil'
    'Atlantic Forest Brazil'
    'Southern Costa Rica'
    'Southern Judean Lowlands Israel'
    'Lower Rhine Embayment Germany'
    'Southwestern Nigeria'
    'Balbina Reservoir Brazilian Amazon'
    'Northeastern Atlantic Forest Brazil'
    'Southwestern Sri Lanka'
    'Central Brazil'
    'Central Brazil'
    'Southwestern Ohio USA'
    'Peninsular Malaysia'
    'Alter do Chao Brazilian Amazon'};

% Decimal degrees: longitude then latitude.
longitude = [ ...
    -95.15
    -47.10
    -45.00
    -82.95
     34.80
      6.90
      5.20
    -59.50
    -35.00
     80.35
    -47.40
    -48.60
    -84.20
    102.17
    -54.95];

latitude = [ ...
     18.60
    -23.80
    -22.50
      8.78
     31.55
     50.90
      6.40
     -1.90
     -8.50
      6.25
    -15.10
    -15.90
     39.70
      3.67
     -2.50];


% =========================================================================
% LOAD WORLD POLYGONS
% =========================================================================
load PlotWorldData

% =========================================================================
% DRAW MAP
% =========================================================================

cla
hold on
ax = gca;
for i = 1:numel(S)
    XX = S(i).X;
    YY = S(i).Y;
    F = [0 find(isnan(XX))];
    for f = 1:length(F)-1
        xx = XX(F(f)+1:F(f+1)-1);
        yy = YY(F(f)+1:F(f+1)-1);
        if length(xx) > 1e2
            patch(ax, xx, yy, 0.8.*ones(1,3), ...
                'EdgeColor', 'none', ...
                'LineWidth', 0.5);
        end
    end
end

set(ax, ...
    'XLim', [-180 180], ...
    'YLim', [-55 85], ...
    'DataAspectRatio', [1 1 1], ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'FontSize', 8, ...
    'XTick', -180:60:180, ...
    'YTick', -60:30:60, ...
    'XTickLabel', {}, ...
    'YTickLabel', {});

% =========================================================================
% ADD STUDY LOCATIONS
% =========================================================================

taxon_levels = {'Aves', 'Insecta', 'Mammalia', 'Mollusca', 'Plantae'};

taxon_colours = [ ...
    0.12 0.47 0.71   % Aves: blue
    0.89 0.10 0.11   % Insecta: red
    0.20 0.63 0.17   % Mammalia: green
    0.99 0.50 0.00   % Mollusca: orange
    0.42 0.24 0.60]; % Plantae: purple

if colour_by_taxon

    legend_handles = gobjects(numel(taxon_levels), 1);

    for t = 1:numel(taxon_levels)

        selected = strcmp(taxon, taxon_levels{t});

        legend_handles(t) = scatter(ax, ...
            longitude(selected), latitude(selected), marker_size, ...
            taxon_colours(t, :), ...
            'o', ...
            'filled', ...
            'MarkerEdgeColor', taxon_colours(t, :).^2, ...
            'LineWidth', marker_edge_width, ...
            'DisplayName', taxon_levels{t});
    end

    if show_legend
        legend(ax, legend_handles, taxon_levels, ...
            'Location', 'southoutside', ...
            'Orientation', 'horizontal', ...
            'Box', 'off', ...
            'FontSize', 8);
    end

else

    scatter(ax, longitude, latitude, marker_size, ...
        [0.12 0.30 0.55], ...
        's', ...
        'filled', ...
        'MarkerEdgeColor', taxon_colours(t, :).^2, ...
        'LineWidth', marker_edge_width);
end
% axis off

end
