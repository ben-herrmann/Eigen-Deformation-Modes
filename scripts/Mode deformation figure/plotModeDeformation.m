clear all; close all; clc
load("AirfoilCalculations_ModeResults_normalizedMeanRef.mat")
%%
first_modes = permute(modes_matrix(:,1,:), [1 3 2]);

nx = 500;
ny = 500;
xu = linspace(min(x),max(x),nx);
yu = linspace(min(y),max(y),ny);
[X,Y] = meshgrid(xu, yu);
U_r = griddata(x,y,real(first_modes(n/2+1:end,1)), X,Y);
U_abs = griddata(x,y,abs(first_modes(n/2+1:end,1)), X,Y);

[M,c] = contour(X,Y, U_r, 6, 'LineWidth', 0.7, 'DisplayName', 'Re = 1850');

hold on

U_r = griddata(x,y,real(first_modes(n/2+1:end,p)), X,Y);
U_abs = griddata(x,y,abs(first_modes(n/2+1:end,p)), X,Y);
[M,c] = contour(X,Y, U_r, 'LevelList', c.LevelList, 'LineWidth', 0.7, 'LineStyle', '--', 'DisplayName', 'Re=2200');
%hold on
%[M,c] = contour(X,Y, U_abs, 'LevelList', c_abs.LevelList, 'LineStyle', '--');
colormap(parula)
colorbar();
plot_naca0012_aoa(1, 100, -5); % Plots the airfoil with a 5-degree angle of attack
xlim([-2 6.75])
ylim([-1 1])
lbl = strtrim(cellstr(num2str((1:5)', 'data%d')));
% Label only the desired lines
%legend('show', lbl(1:2))
set(gca, 'TickLabelInterpreter', 'latex')
exportgraphics(gcf,'ModeDeformationContours.png','Resolution', 1000) %1000 dpi para hacerle zoom en KeyNote