function plot_naca0012_aoa(chord_length, num_points, angle_of_attack_deg)
    % PLOT_NACA0012_AOA Plots the NACA 0012 airfoil at a specified angle of attack.
    % 
    % Inputs:
    %   chord_length - Length of the airfoil chord (default: 1)
    %   num_points - Number of points to define the airfoil surface (default: 100)
    %   angle_of_attack_deg - Angle of attack in degrees (default: 0)
    %
    % Example usage:
    %   plot_naca0012_aoa(1, 100, 10);

    if nargin < 1, chord_length = 1; end
    if nargin < 2, num_points = 100; end
    if nargin < 3, angle_of_attack_deg = 0; end
    
    % Convert angle of attack from degrees to radians
    angle_of_attack_rad = deg2rad(angle_of_attack_deg);
    
    % Thickness distribution for NACA 0012
    t = 0.12; % Thickness as a fraction of the chord
    x = linspace(0, chord_length, num_points);
    y_t = 5 * t * (0.2969 * sqrt(x/chord_length) ...
                   - 0.1260 * (x/chord_length) ...
                   - 0.3516 * (x/chord_length).^2 ...
                   + 0.2843 * (x/chord_length).^3 ...
                   - 0.1015 * (x/chord_length).^4);
    
    % Upper and lower surfaces
    y_upper = y_t;
    y_lower = -y_t;

    % Translate the airfoil if necessary (example translation)
    translation = [-0.25, 0]; % Adjust as necessary
    x_translated = x + translation(1);
    y_upper_translated = y_upper + translation(2);
    y_lower_translated = y_lower + translation(2);

    % Rotation matrix
    R = [cos(angle_of_attack_rad) -sin(angle_of_attack_rad);
         sin(angle_of_attack_rad) cos(angle_of_attack_rad)];
    
    % Apply rotation to airfoil coordinates (vectorized)
    upper_rotated = R * [x_translated; y_upper_translated];
    lower_rotated = R * [x_translated; y_lower_translated];

    % Plotting the airfoil
    % Combine the upper and lower surfaces
    airfoil_x = [upper_rotated(1, :) lower_rotated(1, end:-1:1)];
    airfoil_y = [upper_rotated(2, :) lower_rotated(2, end:-1:1)];
    %color_airfoil = [1.0 1.0 1.0]; %white color
    color_airfoil = [.5 .5 .5]; %greycolor

    fill(airfoil_x, airfoil_y, color_airfoil, 'LineStyle', 'none'); % Fill with white color
    hold on;
    plot(upper_rotated(1, :), upper_rotated(2, :), 'w', 'LineWidth', 0.5);
    plot(lower_rotated(1, :), lower_rotated(2, :), 'w', 'LineWidth', 0.5);
    axis equal;
    hold off;
end