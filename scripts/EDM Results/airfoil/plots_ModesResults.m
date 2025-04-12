clear all; close all; clc
load("AirfoilCalculations_ModeResults_normalizedMeanRef.mat")
% Figure Singular values of first eigenmode deformation
f = figure;
f.Position = [500 500 300 3*300/4];
semilogy(Re_vec, singular_values(1,:), 'o', ...
    'MarkerFaceColor', 'red', 'MarkerSize', 8, 'MarkerEdgeColor', 'black')
grid
yticks([10^(-15), 10^(-10), 10^(-5), 1])
set(gca, 'TickLabelInterpreter', 'latex')
xlim([Re_vec(1) Re_vec(end)])
xlabel('$\mu$', 'Interpreter', 'latex')
set(gca, "FontSize", 17)
%% Plot Deformation Modes
% First Deformation Mode
clim_abs = [0 0];
clim_abs(2) = max(abs(def_modes(n/2+1:end, 1, 1)));

clim_re = [0 0];
clim_re(2) = 0.4*max(real(def_modes(n/2+1:end, 1, 1)));
clim_re(1) = -clim_re(2);

plotField_nek(x, y, real(def_modes(n/2+1:end, 1, 1)), clim_re)
exportgraphics(gcf,'firstDeformationMode.png','Resolution', 500)
% Second Deformation Mode
clim_abs = [0 0];
clim_abs(2) = max(abs(def_modes(n/2+1:end, 2, 1)));

clim_re = [0 0];
clim_re(2) = 0.4*max(real(def_modes(n/2+1:end, 2, 1)));
clim_re(1) = -clim_re(2);
plotField_nek(x, y, real(def_modes(n/2+1:end, 2, 1)), clim_re)
exportgraphics(gcf,'secondDeformationMode.png','Resolution', 500)

% First Mean Mode
clim_abs = [0 0];
clim_abs(2) = max(abs( mean_modes(n/2+1:end,1) ));

clim_re = [0 0];
clim_re(2) = 0.4*max(real( mean_modes(n/2+1:end,1) ));
clim_re(1) = -clim_re(2);

plotField_nek(x, y, real(mean_modes(n/2+1:end,1)), clim_re)
exportgraphics(gcf,'firstMeanMode.png','Resolution', 500)
%% Plot coefficients
f4 = figure;
plot(Re_vec, abs(coeffs(1,:,1)), 'o-','MarkerSize', 7, 'Color', 'blue' ,'LineWidth', 1.3)
f4.Position = [500 500 300 7*300/16];
grid
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$\mu$', 'Interpreter', 'latex')
set(gca, "FontSize", 17)
xlim([Re_vec(1) Re_vec(end)])

f5 = figure;
plot(Re_vec, abs(coeffs(2,:,1)), 'o-', 'MarkerSize', 7, 'Color', 'blue', 'LineWidth', 1.3)
f5.Position = [500 500 300 7*300/16];
grid
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$\mu$', 'Interpreter', 'latex')
set(gca, "FontSize", 17)
xlim([Re_vec(1) Re_vec(end)])
%% Save Plots
%exportgraphics(f,'firstMode_singularValues.png','Resolution', 500)
%exportgraphics(f4,'firstDeformationCoefficient.png','Resolution',500)
%exportgraphics(f5,'secondDeformationCoefficient.png','Resolution',500)