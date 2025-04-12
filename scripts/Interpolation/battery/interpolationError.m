clc; clear all; close all
load("BatteryData_EDMsResults.mat")
%%
FEM_M = assembleFEMatrices(model, 'M');
M = FEM_M.M;
clear FEM_M
[R,flag] = chol(M);
m=6; %number of eigenmodes used from the database per parameter 
[edms, coeffs] = themethod_meanRef(modes_matrix(:,1:m,:), R, 2);
mean_modes = mean(modes_matrix(:,1:m,:), 3);
%% First Eigenmode
[n,~,p] = size(modes_matrix);
first_modes_train = permute(modes_matrix(:,1,:), [1 3 2]);
%directly interpolated modes
di_first_modes = interp1(h_vec, first_modes_train', par_space)';
%reduced representation of modes
coeffs_int = interp1(h_vec, coeffs(:,:,1)', par_space)'; %edm coeffs interpolation
red_first_modes_r1 = mean_modes(:,1) + edms(:,1,1)*coeffs_int(1,:); %edm and coeffs multiplication
red_first_modes_r2 = mean_modes(:,1) + edms(:,:,1)*coeffs_int; %edm and coeffs multiplication

error_mean = vecnorm(R*(first_modes_gt - mean_modes(:,1)))./vecnorm(R*first_modes_gt);
error_di = vecnorm(R*(first_modes_gt - di_first_modes))./vecnorm(R*first_modes_gt);
error_red_r1 = vecnorm(R*(first_modes_gt - red_first_modes_r1))./vecnorm(R*first_modes_gt);
error_red_r2 = vecnorm(R*(first_modes_gt - red_first_modes_r2))./vecnorm(R*first_modes_gt);
%% Plot error comparison First Eigenmode
f = figure;
f.Position = [500 500 300 3*300/4];
semilogy(par_space, error_di+1e-15, 'b-', ...
    'LineWidth', 1.3, 'DisplayName', 'direct interpolation')
grid
hold on
semilogy(par_space, error_red_r1, 'r', ...
    'LineWidth', 1.3, 'DisplayName', 'r=1 EDMs')
hold on
semilogy(par_space, error_red_r2, '--', 'LineWidth', 1.3, ...
    'DisplayName', 'r=2 EDMs')
hold on
semilogy(par_space, error_mean, 'LineWidth', 1.3, ...
    'DisplayName', 'constant mean mode')
xlabel('$\mu$', 'Interpreter', 'latex')
ylabel('error', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
title('Battery heat transfer', 'Interpreter', 'latex')
xlim([h_vec(1) h_vec(end)])
ylim([1e-6 1])
%% export figure
exportgraphics(f,'Battery_FirstEigenmodeInterpolationError.png','Resolution', 500)
%% third eigenmode
third_modes_train = permute(modes_matrix(:,3,:), [1 3 2]);
%directly interpolated modes
di_third_modes = interp1(h_vec, third_modes_train', par_space)';
%reduced representation of modes
coeffs_int = interp1(h_vec, coeffs(:,:,3)', par_space)'; %edm coeffs interpolation
red_third_modes_r1 = mean_modes(:,3) + edms(:,1,3)*coeffs_int(1,:); %edm and coeffs multiplication
red_third_modes_r2 = mean_modes(:,3) + edms(:,:,3)*coeffs_int; %edm and coeffs multiplication

third_error_mean = vecnorm(R*(third_modes_gt - mean_modes(:,3)))./vecnorm(R*third_modes_gt);
third_error_di = vecnorm(R*(third_modes_gt - di_third_modes))./vecnorm(R*third_modes_gt);
third_error_red_r1 = vecnorm(R*(third_modes_gt - red_third_modes_r1))./vecnorm(R*third_modes_gt);
third_error_red_r2 = vecnorm(R*(third_modes_gt - red_third_modes_r2))./vecnorm(R*third_modes_gt);
%% Plot error comparison Third Eigenmode
f2 = figure;
f2.Position = [500 500 300 3*300/4];
semilogy(par_space, third_error_di+1e-15, 'b-', ...
    'LineWidth', 1.3, 'DisplayName', 'direct interpolation')
grid
hold on
semilogy(par_space, third_error_red_r1, 'r', ...
    'LineWidth', 1.3, 'DisplayName', 'r=1 EDMs')
hold on
semilogy(par_space, third_error_red_r2, '--', 'LineWidth', 1.3, ...
    'DisplayName', 'r=2 EDMs')
hold on
semilogy(par_space, third_error_mean, 'LineWidth', 1.3, ...
    'DisplayName', 'constant mean mode')
xlabel('$\mu$', 'Interpreter', 'latex')
ylabel('error', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
title('Battery heat transfer', 'Interpreter', 'latex')
xlim([h_vec(1) h_vec(end)])
ylim([1e-6 1])
%% export figure
exportgraphics(f2,'Battery_ThirdEigenmodeInterpolationError.png','Resolution', 500)