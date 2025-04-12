function [def_modes, coefficients, singular_values] = themethod_meanRef(data_set, R, retained_modes)
%%%%%%
% Extracts EDMs and EDM coefficients
% of a data set of eigenmodes computed for different parameter values. The
% data set is a set of reduced basis {V_k} computed for mu_k values. Each
% one of the V_k basis is orthogonal using a non-standard inner product
% defined by M=R'*R. The 
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
for i=1:m
    phi_i = sorted_data_set(:,:,i);
    phi_i = phi_i - mean(phi_i, 2); % Re
    % multiplicamos por matriz R, sea vector o matriz
    if R_dim2 == 1
        phi_i = R.*phi_i;
    else
        phi_i = R*phi_i;
    end
    % SVD
    [U,S,V] = svd(phi_i,"econ");
    Ur = U(:,1:retained_modes);
    Sr = S(1:retained_modes,1:retained_modes);
    Vr = V(:,1:retained_modes);
    if R_dim2 == 1
        def_modes(:,:,i) = (1./R).*Ur;
    else
        def_modes(:,:,i) = Ri*Ur;
    end
    %def_modes(:,:,i) = Ur;
    coefficients(:,:,i) = Sr*Vr';

    %store singular values
    singular_values(i,:) = diag(S);
end