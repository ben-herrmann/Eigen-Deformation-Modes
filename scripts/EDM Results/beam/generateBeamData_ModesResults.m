clc; close all; clear all
% Create Geometry
%gm = multicuboid(0.1,0.005,0.00125);options = optimset('PlotFcns',@optimplotfval);

%Create Geometry
d=0.005/4;
xc_vec = 0.01:d:0.09;
L=0.1;
p = length(xc_vec);
Hs = ones(1, p+2)*d;
Hs(1) = (L-p*d)/2;
Hs(end) = (L-p*d)/2;
Z_offs = zeros(1, p+2);
for i=1:p
    Z_offs(i+1) = (L-p*d)/2 + (i-1)*d;
end
Z_offs(end) = (L-p*d)/2 + p*d;
gm = multicuboid(0.00125, 0.005, Hs ,ZOffset=Z_offs);
rotate(gm, 90 , [6.25E-4 -2.5E-3 0], [6.25E-4 2.5E-3 0])
% Vertex
l0 = 6.25e-4; % Distancia al origen de la geometria
%bottomVertex = addVertex(gm,'Coordinates',[0.03 + l0, 0.0025,0]); % Vertice Sensor

%x_sensor1 = l0 + L.*rand;
x_sensor1 = l0 + L/4;
x_vertex2 = l0 + 3*L/4;

%x_sensor2 = l0 + L.*rand;
%bottomVertex = addVertex(gm,'Coordinates',[x_vertex2, 0.0025, 0]); % Arista Cara Inferior
%Sensor2Vertex = addVertex(gm,'Coordinates',[x_sensor2, 0, 0.00125]);
%upperVertex = addVertex(gm,'Coordinates',[x_sensor1, -0.0025,0.00125]); % Arista Cara superior
% Properties
E = 210E9;
nu = 0.3;
rho = 7800;
sModel = createpde('structural','modal-solid');
sModel.Geometry = gm;
numCells = sModel.Geometry.NumCells;
generateMesh(sModel, 'Hmax', 0.002);
n = length(sModel.Mesh.Nodes);
% Plot Cells/(Faces
figure
%pdegplot(sModel,'VertexLabels','on', 'FaceAlpha', 0.5);
pdegplot(sModel);
axis off

% Plot Mesh
figure
pdeplot3D(sModel.Mesh)
%% Modes for sampled parameters
E2 = 10E9; % 10E9 Flaw Young Modulus

structuralProperties(sModel, 'Cell', [1:1 3:numCells], 'YoungsModulus',E,'PoissonsRatio',nu,'MassDensity',rho);
structuralProperties(sModel, 'Cell', 2, 'YoungsModulus', E2 ,'PoissonsRatio',nu,'MassDensity',rho); %5E8
structuralBC(sModel,'Face',1,'Constraint','fixed');
%sys_i = linearize(sModel);
%num_mechss = length(sys_i.M);

numExtractedModes = 5;
% arrays para modos Pde toolbox
modes_matrix = zeros(3*n,numExtractedModes,p);
freqs = zeros(p,numExtractedModes);

% arrays para modos eigs()
%modes_aug = zeros(num_mechss,numExtractedModes,p);
%omega_squared = zeros(p,numExtractedModes);
for i=1:p
    disp(['Computing parameter: x_c=', num2str(xc_vec(i)) ])
    delete(sModel.MaterialProperties.MaterialAssignments)
    delete(sModel.BoundaryConditions)
    structuralProperties(sModel, 'Cell', [1:i i+2:numCells], 'YoungsModulus',E,'PoissonsRatio',nu,'MassDensity',rho);
    structuralProperties(sModel, 'Cell', i+1, 'YoungsModulus', E2 ,'PoissonsRatio',nu,'MassDensity',rho); %5E8
    structuralBC(sModel,'Face',1,'Constraint','fixed');
    RModal_i = solve(sModel, 'FrequencyRange', [-0.001, 1e6]);
    modes_matrix(1:n, :, i) = RModal_i.ModeShapes.ux(:,1:numExtractedModes);
    modes_matrix(n+1:2*n, :, i) = RModal_i.ModeShapes.uy(:,1:numExtractedModes);
    modes_matrix(2*n+1:3*n, :, i) = RModal_i.ModeShapes.uz(:,1:numExtractedModes);
    freqs(i,:) = RModal_i.NaturalFrequencies(1:numExtractedModes);

    %sys_i = linearize(sModel);
    %[Vi,Di] = eigs(sys_i.K, sys_i.M, numExtractedModes, 'smallestabs');
    %modes_aug(:,:,i) = Vi;
    %omega_squared(i, :) = diag(Di);
    % Plot eigenvalues
    %figure
    %scatter(zeros(1,numExtractedModes), RModal_i.NaturalFrequencies(1:numExtractedModes), 'ro')
    %hold on;
    %scatter(zeros(1,numExtractedModes), sqrt(diag(Di)), 'bx')
    % Add legend
    %legend('Toolbox Modal Analysis', 'eigs(K,M)');
end
%% Sort based in dot product
for j=1:numExtractedModes
    for i=1:p-1
        %
        dot_prod_i = modes_matrix(:,j,i)'*modes_matrix(:,j,i+1);
        if dot_prod_i<0
            modes_matrix(:,j,i+1)=-modes_matrix(:,j,i+1);
        end
    end
end
%% Plot modes
close all

comp_idx = 1; % Ux,Uy,Uz -> 1,2,3
mode_idx = 3; %1:numExtractedModes

max_max = max(modes_matrix((comp_idx-1)*n+1:comp_idx*n, mode_idx, :), [], 'all');
min_min = min(modes_matrix((comp_idx-1)*n+1:comp_idx*n, mode_idx, :), [], 'all');

for i=1:4:p
    figure
    pdeplot3D(sModel, 'ColorMapData', modes_matrix((comp_idx-1)*n+1:comp_idx*n, mode_idx, i) )
    %clim([min_min max_max])
    title(['Position $x_c= $', num2str(xc_vec(i))] , 'Interpreter', 'latex')
end

save beamData_ModeResults.mat sModel modes_matrix  xc_vec n
%%
%FEM_M = assembleFEMatrices(sModel, 'M');
%M = FEM_M.M;
%[R,flag] = chol(M);
%database = modes_matrix(:,:,1:4:32);
%[def_modes, coeffs, singular_values] = themethod_noref(database, R, 2);
%%
%figure
%pdeplot3D(sModel, 'ColorMapData', def_modes((comp_idx-1)*n+1:comp_idx*n, 2, 1) )


%%
% Uxmode = permute(modes_matrix(1:n,1,:), [1 3 2]);
% Uymode = permute(modes_matrix(n+1:2*n,1,:), [1 3 2]);
% Uzmode = permute(modes_matrix(2*n+1:end,1,:), [1 3 2]);
% 
% stacked_modes = permute(modes_matrix(:,1,:), [1 3 2]);
% [U,S,V] = svd(stacked_modes, 'econ');
% figure
% semilogy(diag(S), 'o')
% figure
% plot(cumsum(diag(S))/sum(diag(S)), 'o')