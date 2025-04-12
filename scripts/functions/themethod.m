function [def_modes, coefficients] = themethod(data_set, R_matrix, retained_modes)

%THEMETHOD Extracts parametric modal deformation modes and coefficientes of
%data set
%
%
% themethod(data_set, R_matrix, retained_modes)
%
% THEMETHOD Extracts parametric modal deformation modes and coefficientes
% of a data set of eigenmodes computed for different parameter values. The
% data set is a set of reduced basis {V_k} computed for mu_k values. Each
% one of the V_k basis is orthogonal using a non-standard


[n,m,p] = size(data_set);
Ri = pseudoinverse(R_matrix);
disp('pseudoinverse() function used')
%Ri = R_matrix\speye(n); %This is the one
%initialize resultsdef
def_modes = zeros(n,retained_modes,m);
coefficients = zeros(retained_modes,p,m);
% sort modes by position
sorted_data_set = permute(data_set, [1 3 2]);
for i=1:m
    % first parameter as reference mode
    ref_mode = sorted_data_set(:,1,i);
    % restamos el modo de referencia a cada columna
    phi_i = sorted_data_set(:,:,i) - ref_mode;
    % multiplicamos por matriz R
    phi_i = R_matrix*phi_i;
    % SVD
    [U,S,V] = svd(phi_i,"econ");
    Ur = U(:,1:retained_modes);
    Sr = S(1:retained_modes,1:retained_modes);
    Vr = V(:,1:retained_modes);
    def_modes(:,:,i) = Ri*Ur;
    %def_modes(:,:,i) = Ur;
    coefficients(:,:,i) = Sr*Vr';
end