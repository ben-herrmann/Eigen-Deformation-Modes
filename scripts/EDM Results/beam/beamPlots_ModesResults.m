clc; close all; clear all
load('beamCalculationsE2_ModesResults.mat')
%% Figure Singular values of first eigenmode deformation
f = figure;
f.Position = [500 500 300 3*300/4];
semilogy(xc_data, singular_values(1,:), 'o', ...
    'MarkerFaceColor', 'red', 'MarkerSize', 8, 'MarkerEdgeColor', 'black')
grid
yticks([10^(-15), 10^(-10), 10^(-5), 1])
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$\mu$', 'Interpreter', 'latex')
set(gca, "FontSize", 17)
%% Modal Deformation Modes of first Ux Mode
f2 = figure;
pdeplot3D(sModel.Mesh, "ColorMapData", def_modes(1:n,1,1))
colorbar("off")
f3 = figure;
pdeplot3D(sModel.Mesh, "ColorMapData", def_modes(1:n,2,1))
colorbar("off")
%% Modal Deformation Coefficients
f4 = figure;
plot(xc_data, coeffs(1,:,1), 'o-', 'MarkerSize', 7, 'Color', 'blue', 'LineWidth', 1.3)
f4.Position = [500 500 300 7*300/16];
grid
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$\mu$', 'Interpreter', 'latex')
set(gca, "FontSize", 17)
xlim([xc_data(1) xc_data(end)])

f5 = figure;
plot(xc_data, coeffs(2,:,1), 'o-', 'MarkerSize', 7, 'Color', 'blue', 'LineWidth', 1.3)
f5.Position = [500 500 300 7*300/16];
grid
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$\mu$', 'Interpreter', 'latex')
set(gca, "FontSize", 17)
xlim([xc_data(1) xc_data(end)])
%% Save Plots
%exportgraphics(f,'firstMode_singularValues.png','Resolution',300)
%exportgraphics(f2,'firstDeformationMode.png','Resolution',300)
%exportgraphics(f3,'secondDeformationMode.png','Resolution',300)
%exportgraphics(f4,'firstDeformationCoefficient.png','Resolution',300)
%exportgraphics(f5,'secondDeformationCoefficient.png','Resolution',300)