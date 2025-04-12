function interp_basis=interpolate_whole(data_base, parameter_grid, mu_value)
[n,m,~] = size(data_base);
interp_basis=zeros(n,m);
for j=1:m
    modos_j = permute(data_base(:,j,:), [1 3 2]);
    interpolated_mode_j = interp1(parameter_grid, modos_j', mu_value, 'spline', 'extrap')';
    interp_basis(:,j) = interpolated_mode_j;
end
end