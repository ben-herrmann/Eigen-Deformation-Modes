function interp_basis=interpolate_model(deformation_modes, coefficients, parameter_grid, ref_basis, mu_value)
[n,r,m] = size(deformation_modes);
interp_basis=zeros(n,m);
for j=1:m
    def_modesj = deformation_modes(:,:,j);
    Aj = coefficients(:,:,j);
    interpolated_coefficientsj = interp1(parameter_grid, Aj', mu_value, 'spline', 'extrap')';
    interp_basis(:,j) = def_modesj*interpolated_coefficientsj;
end
interp_basis = interp_basis + ref_basis;
end