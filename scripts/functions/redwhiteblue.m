function c = redwhiteblue(m)
    if nargin < 1, m = size(get(gcf, 'colormap'), 1); end
    
    % If m is even, divide the colormap into two equal parts (red to white, white to blue)
    if mod(m, 2) == 0
        m1 = m / 2;
    else
        % If m is odd, adjust so the white color is centered
        m1 = floor(m / 2);
    end
    
    % Part 1: Red to white
    r1 = linspace(1, 1, m1); % Red stays 1
    g1 = linspace(0, 1, m1); % Green goes from 0 to 1
    b1 = linspace(0, 1, m1); % Blue goes from 0 to 1
    
    % Part 2: White to blue
    r2 = linspace(1, 0, m - m1); % Red goes from 1 to 0
    g2 = linspace(1, 0, m - m1); % Green goes from 1 to 0
    b2 = linspace(1, 1, m - m1); % Blue stays 1
    
    % Concatenate the color channels
    r = [r1 r2];
    g = [g1 g2];
    b = [b1 b2];
    
    % Create the colormap by combining r, g, and b
    c = [r' g' b'];
end
