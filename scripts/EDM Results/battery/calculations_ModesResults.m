clc; clear all; close all
load("BatteryData_EDMsResults.mat")
%% Mass Matrix, Cholesky decomposition and method application
FEM_M = assembleFEMatrices(model, 'M');
[R,flag] =  chol(FEM_M.M);
% Method
m=6; % Number of eigenmodes used for model reduction
modes_matrix_m = modes_matrix(:,1:m,:);
[def_modes, coeffs, singular_values] = themethod_meanRef(modes_matrix_m, R, 2);
%% Mean modes
mean_modes = mean(modes_matrix_m, 3);
%%
save Calculations_EDMsResults.mat model h_vec modes_matrix def_modes coeffs singular_values mean_modes