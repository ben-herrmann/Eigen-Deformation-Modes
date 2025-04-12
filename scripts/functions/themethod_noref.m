function [def_modes, coefficients, singular_values] = themethod_noref(data_set, R, retained_modes)
%%%%%%
% Extracts parametric modal deformation modes and coefficientes
% of a data set of eigenmodes computed for different parameter values. The
% data set is a set of reduced basis {V_k} computed for mu_k values. Each
% one of the V_k basis is orthogonal using a non-standard inner product
% defined by M=R'*R.
%%%%%%

[n,m,p] = size(data_set);
[~,R_dim2] = size(R);
if R_dim2 == n
    Ri = pseudoinverse(R);
    disp('pseudoinverse() function used')
end
%Ri = R_matrix\speye(n); %This is the one
%initialize results
def_modes = zeros(n,retained_modes,m);
coefficients = zeros(retained_modes,p,m);
singular_values = zeros(m, p);
% sort modes by position
sorted_data_set = permute(data_set, [1 3 2]);
for j=1:m
    phi_j = sorted_data_set(:,:,j);
    % multiplicamos por matriz R, sea vector o matriz
    if R_dim2 == 1
        phi_j = R.*phi_j;
    else
        phi_j = R*phi_j;
    end
    % SVD
    [U,S,V] = svd(phi_j,"econ");
    Ur = U(:,1:retained_modes);
    Sr = S(1:retained_modes,1:retained_modes);
    Vr = V(:,1:retained_modes);
    if R_dim2 == 1
        def_modes(:,:,j) = (1./R).*Ur;
    else
        def_modes(:,:,j) = Ri*Ur;
    %def_modes(:,:,i) = Ur;
    coefficients(:,:,j) = Sr*Vr';

    %store singular values
    singular_values(j,:) = diag(S);
end