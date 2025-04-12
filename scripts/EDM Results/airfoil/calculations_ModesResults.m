clear all; close all; clc
load("mode_variation.mat");
load('lambdas.mat');
Re_vec = 1850:50:2200;
%%
% arrange modes acording to corresponding eigenvalue
first_mode = cat(3, u_tensor(:,2,1:6), u_tensor(:,2,7:end));
second_mode = cat(3, u_tensor(:,4,1:5), u_tensor(:,6,6:end));
third_mode = cat(3, u_tensor(:,7,1), u_tensor(:,5,2:end)); % CORRECT MODE PAIRING

modes_matrix = cat(2, first_mode, second_mode, third_mode);

[n,m,p] = size(modes_matrix);

% Normalization using Q = Mass^0.5
for k=1:p
    modes_matrix(:,:,k) = modes_matrix(:,:,k)./vecnorm(Q.*modes_matrix(:,:,k));
end
%% MODE SORT, IMPORTANT PART
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
%% Mean modes
mean_modes = mean(modes_matrix, 3);
%% Method
[def_modes, coeffs, singular_values] = themethod_meanRef(modes_matrix, Q, 2);

save AirfoilCalculations_ModeResults_normalizedMeanRef.mat modes_matrix def_modes coeffs singular_values Re_vec x y Q n m p mean_modes