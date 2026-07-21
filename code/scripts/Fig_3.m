clear; clc; close all;

%% ============================================================
% User settings
%% ============================================================
Delta_min = 0.5;
Delta_max = 5.0;

gamma_min = 0.05;
gamma_max = 10;

gamma_range = [linspace(gamma_min, 0.999, 140), 1, linspace(1.001, gamma_max, 180)];
Delta_range = linspace(Delta_min, Delta_max, 220);

% x-axis share of the under-parameterized regime (gamma < 1)
left_share = 0.30;     % 30% for gamma < 1, 70% for gamma > 1

alpha_left  = 0.45;    % gamma < 1 (Theorem 1)
alpha_right = 0.50;    % gamma > 1 (Theorem 2, isotropic)

% gamma ticks to display (labels show true gamma values)
gamma_ticks = [0.05, 0.5, 1, 5, 10];

% --- Colormap: deep blue -> deep red, mildly lightened for translucency
num_colors  = 256;
c_low       = [0.10, 0.15, 0.45];   % deep blue  (low error)
c_high      = [0.55, 0.05, 0.10];   % deep red   (high error)
white_blend = 0.20;                 % 0 = fully saturated, 1 = white
cmap_dark = [linspace(c_low(1), c_high(1), num_colors)', ...
             linspace(c_low(2), c_high(2), num_colors)', ...
             linspace(c_low(3), c_high(3), num_colors)'];
cmap_dark = (1 - white_blend) * cmap_dark + white_blend * ones(num_colors, 3);

% --- View: leave free for interactive rotation. After rotating,
%     run   [az, el] = view(gca)   in the command window, then set
%     fix_view = true and paste the values below and re-run.
fix_view = true;
view_az  = 6.1072;
view_el  = 47.1083;

save_figure = true;   % set true only after fixing the view
output_name = 'Figure_3D.pdf';

%% ============================================================
% Piecewise-linear x-axis rescaling:
%   [gamma_min, 1] -> [0, left_share]   (left part of the axis)
%   [1, gamma_max] -> [left_share, 1]   (right part of the axis)
%% ============================================================
u = @(g) (g <= 1) .* (left_share .* (g - gamma_min) ./ (1 - gamma_min)) ...
       + (g >  1) .* (left_share + (1 - left_share) .* (g - 1) ./ (gamma_max - 1));

%% ============================================================
% Meshgrid: x = rescaled gamma, y = Delta
%% ============================================================
[GAM, DEL] = meshgrid(gamma_range, Delta_range);
GAMPLOT = u(GAM);

idx_one = find(gamma_range == 1, 1, 'first');

ERR = limiting_lda_error(GAM, DEL);

%% ============================================================
% Figure
%% ============================================================
figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 5, 5]);
set(gcf, 'Renderer', 'opengl');

ax = axes;
hold(ax, 'on');

% gamma <= 1 : Theorem 1 (general Sigma), joined at gamma = 1
surf(ax, GAMPLOT(:, 1:idx_one), DEL(:, 1:idx_one), ...
    ERR(:, 1:idx_one), ERR(:, 1:idx_one), ...
    'FaceColor', 'interp', ...
    'FaceAlpha', alpha_left, ...
    'EdgeColor', 'none');

% gamma >= 1 : Theorem 2 (isotropic Sigma), joined at gamma = 1
surf(ax, GAMPLOT(:, idx_one:end), DEL(:, idx_one:end), ...
    ERR(:, idx_one:end), ERR(:, idx_one:end), ...
    'FaceColor', 'interp', ...
    'FaceAlpha', alpha_right, ...
    'EdgeColor', 'none');

colormap(ax, cmap_dark);
caxis(ax, [0, 0.5]);

xlabel(ax, '$\gamma$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax, '$\Delta$', 'Interpreter', 'latex', 'FontSize', 12);
zlabel(ax, 'Limiting $R_{\mathrm{LDA}}$', 'Interpreter', 'latex', 'FontSize', 12);

grid(ax, 'on');
box(ax, 'off');

xlim(ax, [0, 1]);
ylim(ax, [Delta_min, Delta_max]);
zlim(ax, [0, 0.5]);

xticks(ax, u(gamma_ticks));
xticklabels(ax, arrayfun(@(g) sprintf('%g', g), gamma_ticks, ...
    'UniformOutput', false));
zticks(ax, [0, 0.1, 0.2, 0.3, 0.4, 0.5]);

set(ax, ...
    'LineWidth', 0.8, ...
    'FontName', 'Times New Roman', ...
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 12);

if fix_view
    view(ax, view_az, view_el);
else
    rotate3d(gcf, 'on');
    disp('Rotate freely; then run  [az, el] = view(gca)  to record the angle.');
end

hold(ax, 'off');

%% ============================================================
% Export (only after fixing the view)
%% ============================================================
if save_figure && fix_view
    exportgraphics(gcf, output_name, ...
        'ContentType', 'image', ...
        'Resolution', 600);
end

%% ============================================================
% Local function: limiting error from Theorems 1 and 2
%   gamma < 1 : Phi( -sqrt(1-gamma) * D^2 / (2 sqrt(D^2 + 4 gamma)) )
%   gamma > 1 : Phi( -sqrt(gamma-1) * D^2 / (2 gamma sqrt(D^2 + 4 gamma)) )
%   gamma = 1 : 0.5, the common one-sided limit (continuous extension)
%% ============================================================
function ERR = limiting_lda_error(GAM, DEL)

    ERR = zeros(size(GAM));
    D2  = DEL.^2;

    % gamma < 1 (Theorem 1)
    left = (GAM < 1);
    arg_left = sqrt(1 - GAM(left)) .* D2(left) ...
        ./ (2 .* sqrt(D2(left) + 4 .* GAM(left)));
    ERR(left) = normcdf(-arg_left);

    % gamma > 1 (Theorem 2, isotropic)
    right = (GAM > 1);
    arg_right = sqrt(GAM(right) - 1) .* D2(right) ...
        ./ (2 .* GAM(right) .* sqrt(D2(right) + 4 .* GAM(right)));
    ERR(right) = normcdf(-arg_right);

    % interpolation threshold: continuous extension by the common limit
    ERR(abs(GAM - 1) < 1e-12) = 0.5;

end