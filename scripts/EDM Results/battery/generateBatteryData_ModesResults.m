%% Generate Data Base for battery
clear variables; close all; clc;
% Define the key geometric parameters of a prismatic battery cell.
cellWidth = 150/1000;
cellThickness = 15/1000;
tabThickness = 10/1000;
tabWidth = 15/1000;
cellHeight = 100/1000;
tabHeight = 5/1000;
connectorHeight = 3/1000;
%Define the number of cells in the battery module.
numCellsInModule = 20;
%Create the geometry of the battery module by using the supporting function createBatteryModuleGeometry, which is provided in a supporting file.
[geomModule,volumeIDs,boundaryIDs,volume,area,ReferencePoint] = ...
    createBatteryModuleGeometry(numCellsInModule, ...
                                cellWidth, ...
                                cellThickness, ...
                                tabThickness, ...
                                tabWidth, ...
                                cellHeight, ...
                                tabHeight, ...
                                connectorHeight);
%Plot the geometry.
pdegplot(geomModule)

% Create a thermal model for transient analysis, and assign the battery module geometry to the model.
model = createpde("thermal","transient");
model.SolverOptions.AbsoluteTolerance = 1e-8; % default es 1e-6
model.SolverOptions.RelativeTolerance = 1e-8; % default es 1e-3
model.Geometry = geomModule;
% Generate and plot the mesh
generateMesh(model);
pdemesh(model)

% Collect IDs for assigning material properties and boundary conditions.
cellIDs = [volumeIDs.Cell];
tabIDs = [volumeIDs.TabLeft,volumeIDs.TabRight];
connectorIDs = [volumeIDs.ConnectorLeft,volumeIDs.ConnectorRight];
bottomPlateFaces = [boundaryIDs.BottomFace];
% Specify the thermal conductivity of the battery, in W/(K*m).
cellThermalCond.inPlane = 80;
cellThermalCond.throughPlane = 2;
tabThermalCond = 386;
connectorThermalCond = 400;
% Specify the mass densities of the battery components, in kg/m³.
density.Cell = 780;
density.TabLeft = 2700;
density.TabRight = 2700;
density.ConnectorLeft = 540;
density.ConnectorRight = 540;
% Specify the specific heat values of the battery components, in J/(kg*K).
spHeat.Cell = 785;
spHeat.TabLeft = 890;
spHeat.TabRight = 890;
spHeat.ConnectorLeft = 840;
spHeat.ConnectorRight = 840;
% Assign the thermal material properties of the cell body, tabs, and connectors.
thermalProperties(model,Cell=cellIDs, ...
    ThermalConductivity=[cellThermalCond.throughPlane
                         cellThermalCond.inPlane
                         cellThermalCond.inPlane], ...
    MassDensity=density.Cell, ...
    SpecificHeat=spHeat.Cell);
thermalProperties(model,Cell=tabIDs, ...
    ThermalConductivity=tabThermalCond, ...
    MassDensity=density.TabLeft, ...
    SpecificHeat=spHeat.TabLeft);
thermalProperties(model,Cell=connectorIDs, ...
    ThermalConductivity=connectorThermalCond, ...
    MassDensity=density.ConnectorLeft, ...
    SpecificHeat=spHeat.ConnectorLeft);
% Define the ambient temperature, in K.
Tambient = 293;
% Apply natural convection cooling applied at the front and back faces of the module.
% thermalBC(model, ...
%           Face=[boundaryIDs(1).FrontFace,boundaryIDs(end).BackFace], ...
%           ConvectionCoefficient=15, ...
%           AmbientTemperature=Tambient);
% Apply nominal heat generation. Assume that the heat generation during normal operation is 15 W.
nominalHeatGen = 15/volume(1).Cell;
internalHeatSource(model,nominalHeatGen,Cell=cellIDs);
% Include the cooling effect on the bottom faces of the module
thermalBC(model,Face=bottomPlateFaces, ...
                ConvectionCoefficient=100, ...
                AmbientTemperature=Tambient);
% Switch the analysis type of the thermal model to "modal".
model.AnalysisType = "modal";
% Solve for modes of the thermal model in the specified decay range.
%% Parametros sampleados
[~,n] = size(model.Mesh.Nodes);
%h_vec = 0:40:120; %0:5:20
%h_vec = 0:1:5;
h_vec = 0:4:28;
%h_vec(1) = 0.001;
p=length(h_vec);
numExtractedModes=10;
modes_matrix = zeros(n,numExtractedModes,p);
lambdas = zeros(p,numExtractedModes);
fixed_points = zeros(n,p);
for i=1:p
    model.AnalysisType = "Modal";
    thermalBC(model, ...
          Face=[boundaryIDs(1).FrontFace,boundaryIDs(end).BackFace], ...
          ConvectionCoefficient=h_vec(i), ...
          AmbientTemperature=Tambient);
    RModal_i = solve(model, "DecayRange", [-Inf,5e-2]);
    size(RModal_i.DecayRates)
    modes_matrix(:,:,i) = RModal_i.ModeShapes(:,1:numExtractedModes);
    lambdas(i,:) = RModal_i.DecayRates(1:numExtractedModes);
    model.AnalysisType = "SteadyState";
    RSteady_i = solve(model);
    fixed_points(:,i) = RSteady_i.Temperature;
end

%Sort based in dot product
for j=1:numExtractedModes
    for i=1:length(h_vec)-1
        dot_prod_i = modes_matrix(:,j,i)'*modes_matrix(:,j,i+1);
        if dot_prod_i<0
            modes_matrix(:,j,i+1)=-modes_matrix(:,j,i+1);
        end
    end
end
%%
mode_plot=1;
max_max = max(modes_matrix(:,mode_plot,:), [], 'all');
min_min = min(modes_matrix(:,mode_plot,:), [], 'all');
for i=1:length(h_vec)
    figure
    pdeplot3D(model.Mesh, "ColorMapData",modes_matrix(:,mode_plot,i) )
    title(['h = ', num2str(round(h_vec(i)))], 'Interpreter', 'latex')
    colorbar
    %clim([min_min, max_max])
    %filename=strcat('third_mode_', num2str(round(h_vec(i))), '.png' );
    %saveas(gcf,filename)
    % ---- Uncomment to plot w/ black background ----
    %set(gca, 'Color','k', 'XColor','w', 'YColor','w')
    %set(gcf, 'Color','k')
    %colorbar('Color', [1 1 1])
    %set(gcf, 'InvertHardCopy', 'off');
end
%% Compute first mode shape ground truth
par_space = linspace(h_vec(1), h_vec(end)); 
first_modes_gt = zeros(n, length(par_space)); %ground truth first modes
third_modes_gt = zeros(n, length(par_space)); %ground truth third modes
model.AnalysisType = "Modal"; % Change model to Modal
for k=1:length(par_space)
    k
    % Change parameter to mu_k
    thermalBC(model, ...
          Face=[boundaryIDs(1).FrontFace,boundaryIDs(end).BackFace], ...
          ConvectionCoefficient=par_space(k), ...
          AmbientTemperature=Tambient);
    % Solve Modal
    RModal_i = solve(model, "DecayRange", [-Inf,5e-2]);
    %Store first mode shape for mu_k
    first_modes_gt(:,k) = RModal_i.ModeShapes(:,1);
    %Store third mode shape for mu_k
    third_modes_gt(:,k) = RModal_i.ModeShapes(:,3);
end
%%
%Sort grond truth based in dot product
for k=1:length(par_space)
    dot_prod_k = modes_matrix(:,1,1)'*first_modes_gt(:,k);
    if dot_prod_k<0
        first_modes_gt(:,k)=-first_modes_gt(:,k);
    end
    third_dot_prod_k = modes_matrix(:,3,1)'*third_modes_gt(:,k);
    if third_dot_prod_k<0
        third_modes_gt(:,k)=-third_modes_gt(:,k);
    end
end
%% Check if modes third modes are sorted
mode_plot=1;
max_max = max(modes_matrix(:,mode_plot,:), [], 'all');
min_min = min(modes_matrix(:,mode_plot,:), [], 'all');
for k=1:10:length(par_space)
    figure
    pdeplot3D(model.Mesh, "ColorMapData",third_modes_gt(:,k) )
    title(['h = ', num2str(round(par_space(k)))], 'Interpreter', 'latex')
    colorbar
end
%%
save BatteryData_EDMsResults.mat model h_vec modes_matrix first_modes_gt third_modes_gt par_space