function rom_simulation = my_simulate_ROM(initial_condition, eigenvalues, t_list)
[m,~] = size(initial_condition);
rom_simulation = zeros(m, length(t_list));
% Propagate initial condition using eigenvalues
for i=1:length(t_list)
    rom_simulation(:,i) = spdiags( exp(-eigenvalues.*t_list(i)), 0, m, m )*initial_condition;
end
end