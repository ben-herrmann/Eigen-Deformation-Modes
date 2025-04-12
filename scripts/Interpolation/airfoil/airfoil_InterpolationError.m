clear all; close all; clc
load("mode_variation.mat");
load('lambdas.mat');
Re_vec = 1850:50:2200;
%%
% arrange modes acording to each eigenvalue
first_mode = cat(3, u_tensor(:,2,1:6), u_tensor(:,2,7:end));
second_mode = cat(3, u_tensor(:,4,1:5), u_tensor(:,6,6:end));
%third_mode = u_tensor(:, 5,:);
third_mode = cat(3, u_tensor(:,7,1), u_tensor(:,5,2:end)); %!!!! CORRECT MODE PAIRING
%fourth_mode = u_tensor(:,7,:);

modes_matrix = cat(2, first_mode, second_mode, third_mode);

[n,m,p] = size(modes_matrix);

% Normalization using Q = Mass^0.5
for k=1:p
    modes_matrix(:,:,k) = modes_matrix(:,:,k)./vecnorm(Q.*modes_matrix(:,:,k));
end
%% MODE SORTING, IMPORTANT PART
% We align the phase of each Complex Mode in modes_matrix
f = @(theta, ref_mode, rot_mode)norm( Q.*(rot_mode*exp(j*theta) - ref_mode) );
for i=[1 3]
    ref_mode = modes_matrix(:, i, 1); % Align to first parameter
    theta0 = 0; %initial phase guess
    for k=2:p
        rot_mode = modes_matrix(:,i,k);
        fun = @(theta)f(theta, ref_mode, rot_mode);
        thOpt = fminsearch(fun, theta0);
        modes_matrix(:,i,k) = modes_matrix(:,i,k)*exp(j*thOpt);
    end
end

% Sort real mode (second one) based in dot product
for i=1:p-1
    dot_prod_i = modes_matrix(:,2,i)'*modes_matrix(:,2, i+1);
    if dot_prod_i<0
        modes_matrix(:,2,i+1)=-modes_matrix(:,2,i+1);
    end
end
%% Method
[edms, coeffs] = themethod_meanRef(modes_matrix(:,:,1:2:end), Q, 3);
%% Error calculation First Eigenmode
mean_modes = mean(modes_matrix(:,:,1:2:end), 3);
first_modes = permute(modes_matrix(:,1,:), [1 3 2]);
first_modes_train = first_modes(:,1:2:end);
first_modes_gt = first_modes(:,1:end-1);
par_space = Re_vec(1:end-1);
Re_train = Re_vec(1:2:end);

di_first_modes = interp1(Re_train, first_modes_train', par_space)';
coeffs_int = interp1(Re_train, coeffs(:,:,1)', par_space)';
red_first_modes_r1 = mean_modes(:,1) + edms(:,1,1)*coeffs_int(1,:);
red_first_modes_r2 = mean_modes(:,1) + edms(:,1:2,1)*coeffs_int(1:2,:);
red_first_modes_r3 = mean_modes(:,1) + edms(:,:,1)*coeffs_int;

error_mean = vecnorm( first_modes_gt - mean_modes(:,1) )./vecnorm(first_modes_gt);
error_di = vecnorm( first_modes_gt - di_first_modes )./vecnorm(first_modes_gt);
error_red_r1 = vecnorm( first_modes_gt - red_first_modes_r1 )./vecnorm(first_modes_gt);
error_red_r2 = vecnorm( first_modes_gt - red_first_modes_r2 )./vecnorm(first_modes_gt);
error_red_r3 = vecnorm( first_modes_gt - red_first_modes_r3 )./vecnorm(first_modes_gt);
%% Plot error comparison First Eigenmode
f = figure;
f.Position = [500 500 300 3*300/4];
semilogy(par_space, error_di+1e-15, 'b', ...
    'LineWidth', 1.3, 'DisplayName', 'direct interpolation')
grid
hold on
semilogy(par_space, error_red_r1, 'r', ...
    'LineWidth', 1.3, 'DisplayName', 'r=1 EDMs')
hold on
semilogy(par_space, error_red_r2, 'color', [0.9290 0.6940 0.1250], 'LineWidth', 1.3, ...
    'DisplayName', 'r=2 EDMs')
semilogy(par_space, error_red_r3, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 1.3, ...
    'DisplayName', 'r=3 EDMs', 'LineStyle', '-')
hold on
semilogy(par_space, error_mean, 'color',[0.4940 0.1840 0.5560] ,'LineWidth', 1.3, ...
    'DisplayName', 'constant mean mode')
xlabel('$\mu$', 'Interpreter', 'latex')
ylabel('error', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
title('Flow around airfoil', 'Interpreter', 'latex')
xlim([par_space(1) par_space(end)])
%yscale(gca,'linear')
ylim([1e-6 1])
%ylim([0 0.03])
%% export figure
exportgraphics(f,'Airfoil_FirstEigenmodeInterpolationError.png','Resolution', 500)
%% Error calculation Third Eigenmode
third_modes = permute(modes_matrix(:,3,:), [1 3 2]);
third_modes_train = third_modes(:,1:2:end);
third_modes_gt = third_modes(:,1:end-1);

di_third_modes = interp1(Re_train, third_modes_train', par_space)';
coeffs_int = interp1(Re_train, coeffs(:,:,3)', par_space)';
red_third_modes_r1 = mean_modes(:,3) + edms(:,1,3)*coeffs_int(1,:);
red_third_modes_r2 = mean_modes(:,3) + edms(:,1:2,3)*coeffs_int(1:2,:);
red_third_modes_r3 = mean_modes(:,3) + edms(:,:,3)*coeffs_int;

third_error_mean = vecnorm( third_modes_gt - mean_modes(:,3) )./vecnorm(third_modes_gt);
third_error_di = vecnorm( third_modes_gt - di_third_modes )./vecnorm(third_modes_gt);
third_error_red_r1 = vecnorm( third_modes_gt - red_third_modes_r1 )./vecnorm(third_modes_gt);
third_error_red_r2 = vecnorm( third_modes_gt - red_third_modes_r2 )./vecnorm(third_modes_gt);
third_error_red_r3 = vecnorm( third_modes_gt - red_third_modes_r3 )./vecnorm(third_modes_gt);
%% Plot error comparison Third Eigenmode
f2 = figure;
f2.Position = [500 500 300 3*300/4];
semilogy(par_space, third_error_di+1e-15, 'b', ...
    'LineWidth', 1.3, 'DisplayName', 'direct interpolation')
grid
hold on
semilogy(par_space, third_error_red_r1, 'r', ...
    'LineWidth', 1.3, 'DisplayName', 'r=1 EDMs')
hold on
semilogy(par_space, third_error_red_r2, 'color', [0.9290 0.6940 0.1250], 'LineWidth', 1.3, ...
    'DisplayName', 'r=2 EDMs')
semilogy(par_space, third_error_red_r3, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 1.3, ...
    'DisplayName', 'r=3 EDMs', 'LineStyle', '-')
hold on
semilogy(par_space, third_error_mean, 'color',[0.4940 0.1840 0.5560] ,'LineWidth', 1.3, ...
    'DisplayName', 'constant mean mode')
xlabel('$\mu$', 'Interpreter', 'latex')
ylabel('error', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
title('Flow around airfoil', 'Interpreter', 'latex')
xlim([par_space(1) par_space(end)])
%yscale(gca,'linear')
ylim([1e-6 1])
%ylim([0 0.03])
%% export figure
exportgraphics(f2,'Airfoil_ThirdEigenmodeInterpolationError.png','Resolution', 500)