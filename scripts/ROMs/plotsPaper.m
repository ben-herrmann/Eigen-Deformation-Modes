clc; close all; clear all;
load('BatteryCalculations.mat')
%% error comparison plot
j = 1;
f = figure;
f.Position = [500 500 300 3*300/4];
semilogy(tlist, error_di(j,:), 'b-', ...
    'LineWidth', 1.3, 'DisplayName', 'direct interpolation')
grid
hold on
semilogy(tlist, error_ri(j,:), 'r--', ...
    'LineWidth', 1.3, 'DisplayName', 'method')
hold on
semilogy(tlist, error_rom_int(j,:), '-', 'LineWidth', 1.3, ...
    'DisplayName', 'solutions interpolation')
xlabel('t (s)', 'Interpreter', 'latex')
ylabel('error', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
xlim([tlist(1) tlist(end)])
yticks([10^(-6) 10^(-5) 10^(-4) 10^(-3) 10^(-2)])
%exportgraphics(f,'errorTimeEvolution_h20.png','Resolution', 500)
%% Plot: Simulation time reduction
f2 = figure;
f2.Position = [500 500 300 3*300/4];
semilogy(par_space, time_DI./time_FOM, "square", 'Color', 'b', 'DisplayName', 'direct interpolation');
hold on
semilogy(par_space, time_RI./time_FOM, "x", 'Color', 'r', 'DisplayName', 'method');
hold on
semilogy(par_space, time_rom_int./time_FOM, "v", 'DisplayName', 'solutions interpolation');
grid on
xlabel('$\mu$', 'Interpreter', 'latex')
ylabel('time reduction', 'Interpreter', 'latex')
ylim([10^(-4) 1])
xticks([0 40 80 120])
yticks([10^(-4) 10^(-3) 10^(-2) 10^(-1) 1])
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
%exportgraphics(f2,'NormalizedComputationalTime.png','Resolution', 500)
%% Plot: Normalized time integrated error over parameter space
f3 = figure;
f3.Position = [500 500 300 3*300/4];
plot(par_space, integral_DI./integral_fom, "square", 'Color', 'b', 'DisplayName', 'Basis direct interpolation');
hold on
plot(par_space, integral_RI./integral_fom, "x", 'Color', 'r', 'DisplayName', 'Basis reduced Interpolation');
hold on
plot(par_space, integral_rom_int./integral_fom, "v", 'DisplayName', 'ROM interpolation');
grid on
xlabel('$\mu$', 'Interpreter', 'latex')
ylabel('time integrated error', 'Interpreter', 'latex')
ylim([0 0.25])
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
xticks([0 40 80 120])
%exportgraphics(f3,'NormalizedIntegratedError.png','Resolution', 500)