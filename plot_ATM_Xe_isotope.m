%% plot_ATM_Xe_isotope.m
% Reproduces the atmospheric Xe isotope evolution figure
% (analogous to Extended Data Fig. 1 in Avice et al. 2025)
%
% Data source: ATM_Xe_isotope_table.xlsx
% Columns used: age_lower, age_upper (Ga), δXe_air (‰·u⁻¹), ± δXe_air (1σ)
%
% Usage: run this script directly in MATLAB (no toolboxes required)

clear; clc; close all;

%% ── 1. DATA ──────────────────────────────────────────────────────────────
% Each row: [age_lower  age_upper  dXe  dXe_1sigma  group_id]
% group_id maps to colour/marker (see legend below)
% age is plotted as midpoint; horizontal error bar spans [age_lower, age_upper]

% Group IDs
%  1 = Srinivasan (1976)
%  2 = Pujol et al. (2009)
%  3 = Pujol et al. (2013)
%  4 = Pujol et al. (2011)
%  5 = Avice et al. (2017/2018)
%  6 = Holland et al. (2013)   → box (age range × dXe range)
%  7 = Avice et al. (2025)

data = [
%  age_lo  age_hi   dXe    1sigma  grp
   3.390   3.570   21.00   3.00     2;  % Pujol 2009
   3.500   3.500   13.70   0.00     1;  % Srinivasan 1976 (no error)
   2.700   3.500   15.00   5.00     3;  % Pujol 2013
   2.800   3.200   10.00   5.00     4;  % Pujol 2011
   3.250   3.350   12.90   1.20     5;  % Avice 2018 Barberton
   2.700   2.700   13.00   1.20     5;  % Avice 2018 Fortescue
   2.400   2.700    3.80   2.50     5;  % Avice 2018 Quetico
   2.450   2.450    6.60   1.50     5;  % Avice 2018 Vetreny
   2.000   2.600    5.80   1.50     5;  % Avice 2018 Isua
   1.900   2.100    1.80   2.20     5;  % Avice 2018 Carnaiba
   2.030   2.170    2.60   2.10     5;  % Avice 2018 Gaoua
   1.700   1.700    0.32   0.78     5;  % Avice 2018 Caramal
   0.520   0.540    1.50   1.60     5;  % Avice 2018 Avranches
   0.403   0.405    0.10   1.90     5;  % Avice 2018 Rhynie
   % Avice 2025 new data
   1.802   2.428   -0.30   1.20     8;  % 2025-ref16 Ongeluk fm: Ardoin, 2022;
   2.439   2.503   -2.00   0.90     8;  % 2025-ref16 FD1A: Ardoin, 2022
   2.427   2.443    0.70   1.50     8;  % 2025-ref16 FD3A: Ardoin, 2022
   2.600   2.700    3.40   0.50     7;  % 2025-ref15 DM4: Almayrac 2021
   2.600   2.900   11.00   2.10     7;  % 2025-ref15 Z10/6: Almayrac 2021
   3.300   3.300   10.30   0.50     9;  % 2025-ref17 Barberton re-analysis A: Broadley, 2022
   3.300   3.300   11.30   0.65     9;  % 2025-ref17 Barberton re-analysis B: Broadley, 2022
   3.390   3.570   19.10   1.80     10;  % Avice 2025 North Pole reanalysis
];

% initial atmosphere composition
U_Xe = [4.6, 40.8];

% Holland et al. (2013) – fluid inclusions: shown as a shaded rectangle
% age range: ~1.5–2.64 Ga  |  dXe range: 0–5 ‰·u⁻¹  (box as in original)
holland_age  = [1.50, 2.64];   % Ga
holland_dXe  = [0.00, 5.00];   % ‰·u⁻¹

%% ── 2. COLOURS & MARKERS ─────────────────────────────────────────────────
% Matching the original figure palette
clr = [ ...
    0.55  0.10  0.10;   % 1 Srinivasan  – dark red diamond
    0.85  0.15  0.15;   % 2 Pujol 2009  – red circle
    0.55  0.60  0.10;   % 3 Pujol 2013  – olive circle
    0.25  0.75  0.85;   % 4 Pujol 2011  – cyan circle
    0.25  0.80  0.50;   % 5 Avice 2017/2018 – green circle
    0.95  0.70  0.40;   % 6 Holland 2013 – orange (box only)
    0.20  0.40  0.75;   % 7 Almayrac 2021  – blue circle
    0.60  0.30  0.70;   % 8 Ardoin, 2022 – Deep Purple
    0.90  0.40  0.70;   % 9 Broadley, 2022 – Magenta/Pink
    0.40  0.40  0.40;   % 10 Avice 2025 – Dark Grey
];

mkr = {'d','o','o','o','o','s','o','o','o','o'};  % marker shapes

%% ── 3. MODEL LINE ────────────────────────────────────────────────────────
% Linear decrease from 39 ‰·u⁻¹ at 4.0 Ga → 0 at 2.0 Ga, then 0
t_model = [4.6, 2.0, 0.0];
d_model = [39,  0,   0 ];


%% ── 4. PLOT ──────────────────────────────────────────────────────────────
% fig = figure('Color','w','Units','centimeters','Position',[3 3 16 13]);
fig = figure('color','white','Units','centimeters','Position',[2 5 17.8 12]);

ax  = axes('Parent',fig);
hold(ax,'on'); box(ax,'on');


% ── Holland box (behind everything) ──────────────────────────────────────
fill(ax, ...
    [holland_age(1), holland_age(2), holland_age(2), holland_age(1)], ...
    [holland_dXe(1), holland_dXe(1), holland_dXe(2), holland_dXe(2)], ...
    clr(6,:), 'FaceAlpha', 0.55, 'EdgeColor', 'none');

% ── Model line ────────────────────────────────────────────────────────────
hMod = plot(ax, t_model, d_model, '-', ...
    'Color', [0.10 0.20 0.75], 'LineWidth', 3.0);

% ── Initial composition ──
hini = scatter (ax, U_Xe(1),U_Xe(2),'red','filled');


% ── Data points with error bars ───────────────────────────────────────────
hLeg = gobjects(10,1);   % one handle per group for legend

for g = 1:10
    idx = data(:,5) == g;
    if ~any(idx), continue; end
    d_sub = data(idx,:);

    age_mid = 0.5*(d_sub(:,1) + d_sub(:,2));
    age_lo  = age_mid - d_sub(:,1);          % left  error
    age_hi  = d_sub(:,2) - age_mid;          % right error
    dXe     = d_sub(:,3);
    dXe_err = d_sub(:,4);

    % error bars (2σ shown, matching figure style)
    errorbar(ax, age_mid, dXe, 2*dXe_err, 2*dXe_err, age_lo, age_hi, ...
        'LineStyle','none', ...
        'Color', clr(g,:), ...
        'CapSize', 4, 'LineWidth', 1.0);

    % markers
    hLeg(g) = plot(ax, age_mid, dXe, mkr{g}, ...
        'MarkerFaceColor', clr(g,:), ...
        'MarkerEdgeColor', clr(g,:)*0.7, ...
        'MarkerSize', 8);
end

% Holland placeholder handle for legend
hHolland = fill(ax, nan, nan, clr(6,:), 'FaceAlpha', 0.55, 'EdgeColor','none');

%% ── 5. AXES FORMATTING ───────────────────────────────────────────────────
set(ax, ...
    'XDir',  'reverse', ...        % time axis: past on left, present on right
    'XLim',  [0, 4.6], ...
    'YLim',  [-5, 42], ...
    'XTick', 0:1:4, ...
    'YTick', 0:10:40, ...
    'FontName', 'Helvetica', ...
    'FontSize', 11, ...
    'TickDir', 'out', ...
    'LineWidth', 0.8);

xlabel(ax, 'time (Gyr ago)', 'FontSize', 12, 'FontName', 'Helvetica');
ylabel(ax, '\deltaXe_{atm} (per mil AMU^{-1})', 'FontSize', 12, 'FontName', 'Helvetica');

% "present" arrow annotation
% annotation('textarrow', ...
%     ax2fig(ax, 0,  -3.5, 'x'), ...
%     ax2fig(ax, 0,   0,   'y'), ...
%     'String', 'present', ...
%     'FontSize', 9, 'FontName', 'Helvetica', ...
%     'HeadLength', 7, 'HeadWidth', 5, 'LineWidth', 1.2);

%% ── 6. LEGEND ────────────────────────────────────────────────────────────
labels = { ...
    'U-Xe', ...
    'Srinivasan (1976)', ...
    'Pujol et al. (2009)', ...
    'Pujol et al. (2011)', ...
    'Pujol et al. (2013)', ...
    'Holland et al. (2013)', ...
    'Avice et al. (2017/2018)', ...
    'Almayrac et al. (2021)', ... 
    'Ardoin et al. (2022)', ...
    'Broadley et al. (2022)', ...
    'Avice et al. (2025)', ...
    'model atm evolution'};

% Build handle array (replace empty gobjects with NaN placeholders)
h_all = [hini; hLeg(1); hLeg(2); hLeg(4); hLeg(3); hHolland; hLeg(5); hLeg(7);hLeg(8);hLeg(9);hLeg(10); hMod];
valid = isgraphics(h_all);
legend(ax, h_all(valid), labels(valid), ...
    'Location', 'northeast', ...
    'FontSize',  9, ...
    'FontName', 'Helvetica', ...
    'Box', 'on', ...
    'EdgeColor', [0.6 0.6 0.6]);

hold(ax,'off');

% set(gca,'FontUnits','points','FontSize',6,'LabelFontSizeMultiplier',1.0);
box on;

%% ── 7. EXPORT ────────────────────────────────────────────────────────────
print(fig, 'ATM_Xe_isotope_evolution.pdf', '-dpdf', '-r300', '-bestfit');
% print(fig, 'ATM_Xe_isotope_evolution.png', '-dpng', '-r300');
fprintf('Saved: ATM_Xe_isotope_evolution.pdf / .png\n');


