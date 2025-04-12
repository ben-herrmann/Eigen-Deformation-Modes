clc; clear all; close all
load('Calculations_EDMsResults.mat')
%% Figure Singular values of first eigenmode deformation
f = figure;
f.Position = [500 500 300 3*300/4];
semilogy(h_vec, singular_values(1,:), 'o', ...
    'MarkerFaceColor', 'red', 'MarkerSize', 8, 'MarkerEdgeColor', 'black')
grid
yticks([10^(-15), 10^(-10), 10^(-5), 1])
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$\mu$', 'Interpreter', 'latex')
set(gca, "FontSize", 17)
%% EDMs
f2 = figure;
pdeplot3D(model.Mesh, "ColorMapData", def_modes(:,1,1))
colorbar("off")
f3 = figure;
pdeplot3D(model.Mesh, "ColorMapData", def_modes(:,2,1))
colorbar("off")
%% EDM Coefficients
f4 = figure;
plot(h_vec, coeffs(1,:,1), 'o-', 'MarkerSize', 7, 'Color', 'blue', 'LineWidth', 1.3)
f4.Position = [500 500 300 7*300/16];
grid
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$\mu$', 'Interpreter', 'latex')
set(gca, "FontSize", 17)
f5 = figure;
plot(h_vec, coeffs(2,:,1), 'o-', 'MarkerSize', 7, 'Color', 'blue', 'LineWidth', 1.3)
f5.Position = [500 500 300 7*300/16];
grid
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$\mu$', 'Interpreter', 'latex')
set(gca, "FontSize", 17)
%% Plot first mean eigenmode
f6 = figure;
pdeplot3D(model.Mesh, "ColorMapData", mean_modes(:,1))
colorbar("off")
%% Save plots
%exportgraphics(f,'firstMode_singularValues.png','Resolution',300)
%exportgraphics(f2,'firstDeformationMode.png','Resolution',300)
%exportgraphics(f3,'secondDeformationMode.png','Resolution',300)
%exportgraphics(f4,'firstDeformationCoefficient.png','Resolution',300)
%exportgraphics(f5,'secondDeformationCoefficient.png','Resolution',300)
%exportgraphics(f6,'firstMeanMode.png','Resolution',300)