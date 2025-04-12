function fhandle = plotField_nek(x,y,u, clims)
%
% This function is used to plot fields extracted from Nek5000 files.
%
%   INPUT
%   - x: x-coordinate of the state
%   - y: y-coordinate of the state
%   - u: nek field over the x,y domain
%   - clims: tuple of colormap min and max values
%   OUTPUT
%   - fhandle: figure with plot of the field, figure is adapted to scale
%   the colormap from min(u(:)) to max(u(:))
%
% Last edit: 20231115 Erick Kracht (erick.kracht@ing.uchile.cl)
%
    figure
    fhandle = gcf;


%-------------------------------------------------------------------------
% Create regular grid across data space
%-------------------------------------------------------------------------

    nx = 500;
    ny = 500;
    xu = linspace(min(x),max(x),nx);
    yu = linspace(min(y),max(y),ny);
    [X,Y] = meshgrid(xu, yu);
    U = griddata(x,y,u,X,Y);

%-------------------------------------------------------------------------
% Set custom colormap (red and blue) and black background
%-------------------------------------------------------------------------
    % Set the background color to black
    %set(gcf, 'Color', 'k');

    % Create contour plot
    %imagesc(xu,yu,U,[min(U(:)), max(U(:))])
    %imagesc(xu,yu,U,[clims(1), clims(2)])
    pcolor(xu,yu,U)
    %contour(xu, yu,U, 20)
    shading interp
    clim(clims)
%     max_val = max(abs(min(U(:))), abs(max(U(:))));
%     max_val = 1.7;
%     imagesc(xu,yu,U,[-max_val, max_val ])
%colormap('sky')
    %colormap(flipud(bone))
    %colormap(redblackblue(1000))
    colormap(redwhiteblue(1000))
    %colormap(redwhiteblue(1000))
    %colormap("jet")
    %colormap(redwhiteblue())
    %colormap(greenblackblue(1000))
    %colormap(summer)
%     colormap(hsv)

    % Set the tick labels and axis labels to white
    %set(gca, 'XColor', 'w', 'YColor', 'w');
    %set(gca, 'Color', 'k');

    hold on

%-------------------------------------------------------------------------
% Adds a cylinder or a naca0012 with aoa, comment one or the other before
% plotting
%-------------------------------------------------------------------------

%---------------------------
%       Cylinder
%---------------------------
%     t = (1:100)/100'*2*pi;
%     xc = 0.5*sin(t);
%     yc = 0.5*cos(t);
%     fill(xc,yc,[.8 .8 .8])
%     plot(xc,yc,'w','LineWidth',1.2)

%---------------------------
%       Naca0012 Airfoil
%---------------------------
    plot_naca0012_aoa(1, 100, -5); % Plots the airfoil with a 5-degree angle of attack
    xlim([-2 6.75])
    ylim([-2.5 2.5])

    %axes('position',[.65 .175 .25 .25])
    %box on

end