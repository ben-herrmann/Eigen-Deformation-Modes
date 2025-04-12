clc; close all; clear all
load("beamData_ModeResults.mat")
%%
dataset = modes_matrix(:,:,1:8:end); %we use 9 parameters from the original data set
xc_train = xc_vec(1:8:end);
xc_gt = xc_vec; % avalailable ground truth parameters
FEM_M = assembleFEMatrices(sModel);
M = FEM_M.M;
clear FEM_M
[R, flag] = chol(M);
[edms, coeffs] = themethod_meanRef(dataset, R, 5); %5
mean_modes = mean(dataset, 3);
%% First Eigenmode
first_modes_gt = permute(modes_matrix(:,1,:), [1 3 2]);
first_modes_train = permute(dataset(:,1,:), [1 3 2]);
[n,~,p] = size(modes_matrix);
%directly interpolated modes
di_first_modes = interp1(xc_train, first_modes_train', xc_gt)';
%reduced representation of modes
coeffs_int = interp1(xc_train, coeffs(:,:,1)', xc_gt)'; %edm coeffs interpolation
red_first_modes_r1 = mean_modes(:,1) + edms(:,1,1)*coeffs_int(1,:); %r=1 edms and coeffs multiplication
red_first_modes_r2 = mean_modes(:,1) + edms(:,1:2,1)*coeffs_int(1:2,:); %r=2 edms and coeffs multiplication
red_first_modes_r3 = mean_modes(:,1) + edms(:,1:3,1)*coeffs_int(1:3,:); %r=3 edms and coeffs multiplication
red_first_modes_r4 = mean_modes(:,1) + edms(:,1:4,1)*coeffs_int(1:4,:); %r=4 edms and coeffs multiplication
red_first_modes_r5 = mean_modes(:,1) + edms(:,:,1)*coeffs_int; %r=5 edms and coeffs multiplication

error_mean = vecnorm(R*(first_modes_gt - mean_modes(:,1)))./vecnorm(R*first_modes_gt);
error_di = vecnorm(R*(first_modes_gt - di_first_modes))./vecnorm(R*first_modes_gt);
error_red_r1 = vecnorm(R*(first_modes_gt - red_first_modes_r1))./vecnorm(R*first_modes_gt);
error_red_r2 = vecnorm(R*(first_modes_gt - red_first_modes_r2))./vecnorm(R*first_modes_gt);
error_red_r3 = vecnorm(R*(first_modes_gt - red_first_modes_r3))./vecnorm(R*first_modes_gt);
error_red_r4 = vecnorm(R*(first_modes_gt - red_first_modes_r4))./vecnorm(R*first_modes_gt);
error_red_r5 = vecnorm(R*(first_modes_gt - red_first_modes_r5))./vecnorm(R*first_modes_gt);
%% Plot error comparison First Eigenmode
f = figure;
f.Position = [500 500 300 3*300/4];
semilogy(xc_gt, error_di+1e-10, 'b', ...
    'LineWidth', 1.3, 'DisplayName', 'direct interpolation')
grid
hold on
semilogy(xc_gt, error_red_r1, 'r', ...
    'LineWidth', 1.3, 'DisplayName', 'r=1 EDMs')
hold on
semilogy(xc_gt, error_red_r2, 'color', [0.9290 0.6940 0.1250], 'LineWidth', 1.3, ...
    'DisplayName', 'r=2 EDMs')
semilogy(xc_gt, error_red_r3, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 1.3, ...
    'DisplayName', 'r=3 EDMs')
hold on
semilogy(xc_gt, error_red_r5, '-', 'color',[1 0.5529 0.1608] , 'LineWidth', 1.3, ...
    'DisplayName', 'r=5 EDMs')
hold on
semilogy(xc_gt, error_mean, 'color',[0.4940 0.1840 0.5560] ,'LineWidth', 1.3, ...
    'DisplayName', 'constant mean mode')
xlabel('$\mu$', 'Interpreter', 'latex')
ylabel('error', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
title('Vibrations on cantilever beam', 'Interpreter', 'latex')
xlim([xc_gt(1) xc_gt(end)])
%yscale(gca,'linear')
ylim([10^(-4) 10^(-1)])
%ylim([0 0.03])
%% export figure
exportgraphics(f,'Beam_FirstEigenmodeInterpolationError.png','Resolution', 500)
%% Third Eigenmode
third_modes_gt = permute(modes_matrix(:,3,:), [1 3 2]);
third_modes_train = permute(dataset(:,3,:), [1 3 2]);
[n,~,p] = size(modes_matrix);
%directly interpolated modes
di_third_modes = interp1(xc_train, third_modes_train', xc_gt)';
%reduced representation of modes
coeffs_int = interp1(xc_train, coeffs(:,:,3)', xc_gt)'; %edm coeffs interpolation
red_third_modes_r1 = mean_modes(:,3) + edms(:,1,3)*coeffs_int(1,:); %r=1 edms and coeffs multiplication
red_third_modes_r2 = mean_modes(:,3) + edms(:,1:2,3)*coeffs_int(1:2,:); %r=2 edms and coeffs multiplication
red_third_modes_r3 = mean_modes(:,3) + edms(:,1:3,3)*coeffs_int(1:3,:); %r=3 edms and coeffs multiplication
red_third_modes_r4 = mean_modes(:,3) + edms(:,1:4,3)*coeffs_int(1:4,:); %r=4 edms and coeffs multiplication
red_third_modes_r5 = mean_modes(:,3) + edms(:,:,3)*coeffs_int; %r=2 edms and coeffs multiplication

third_error_mean = vecnorm(R*(third_modes_gt - mean_modes(:,3)))./vecnorm(R*third_modes_gt);
third_error_di = vecnorm(R*(third_modes_gt - di_third_modes))./vecnorm(R*third_modes_gt);
third_error_red_r1 = vecnorm(R*(third_modes_gt - red_third_modes_r1))./vecnorm(R*third_modes_gt);
third_error_red_r2 = vecnorm(R*(third_modes_gt - red_third_modes_r2))./vecnorm(R*third_modes_gt);
third_error_red_r3 = vecnorm(R*(third_modes_gt - red_third_modes_r3))./vecnorm(R*third_modes_gt);
third_error_red_r4 = vecnorm(R*(third_modes_gt - red_third_modes_r4))./vecnorm(R*third_modes_gt);
third_error_red_r5 = vecnorm(R*(third_modes_gt - red_third_modes_r5))./vecnorm(R*third_modes_gt);
%% Plot error comparison Third Eigenmode
f2 = figure;
f2.Position = [500 500 300 3*300/4];
semilogy(xc_gt, third_error_di+1e-10, 'b', ...
    'LineWidth', 1.3, 'DisplayName', 'direct interpolation')
grid
hold on
semilogy(xc_gt, third_error_red_r1, 'r', ...
    'LineWidth', 1.3, 'DisplayName', 'r=1 EDMs')
hold on
semilogy(xc_gt, third_error_red_r2, 'color', [0.9290 0.6940 0.1250], 'LineWidth', 1.3, ...
    'DisplayName', 'r=2 EDMs')
semilogy(xc_gt, third_error_red_r3, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 1.3, ...
    'DisplayName', 'r=3 EDMs')
hold on
semilogy(xc_gt, third_error_red_r5, '-', 'color',[1 0.5529 0.1608] , 'LineWidth', 1.3, ...
    'DisplayName', 'r=5 EDMs')
hold on
semilogy(xc_gt, third_error_mean, 'color',[0.4940 0.1840 0.5560] ,'LineWidth', 1.3, ...
    'DisplayName', 'constant mean mode')
xlabel('$\mu$', 'Interpreter', 'latex')
ylabel('error', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
title('Vibrations on cantilever beam', 'Interpreter', 'latex')
xlim([xc_gt(1) xc_gt(end)])
%yscale(gca,'linear')
ylim([10^(-3) 5*10^(-1)])
%ylim([0 0.03])
%% export figure
exportgraphics(f2,'Beam_ThirdEigenmodeInterpolationError.png','Resolution', 500)
%% Plot Modes
% comp_idx = 1; %Ux,Uy,Uz
% mode_idx = 3;
% Nnodes = n/3;
% for i=1:4:32
%     figure
%     pdeplot3D(sModel, 'ColorMapData', modes_matrix((comp_idx-1)*Nnodes+1:comp_idx*Nnodes, mode_idx, i) )
%     %clim([min_min max_max])
%     title(['Position $x_c= $', num2str(xc_vec(i))] , 'Interpreter', 'latex')
% end