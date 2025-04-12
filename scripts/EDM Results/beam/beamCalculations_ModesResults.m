clc; close all; clear all
load("beamData_ModeResults.mat")
%%
dataset = modes_matrix(:,:,1:8:end-1); %we use 8 parameters from the original data set
xc_data = xc_vec(1:8:end-1);

FEM_M = assembleFEMatrices(sModel);
M = FEM_M.M;
[R, flag] = chol(M);
[def_modes, coeffs, singular_values] = themethod_meanRef(dataset, R, 2);

save 'beamCalculationsE2_ModesResults.mat' sModel def_modes coeffs singular_values xc_data n

%%
for k=1:8:64
    figure
    pdeplot3D(sModel.Mesh, "ColorMapData", modes_matrix(1:n, 1, k))
    figfilename = strcat('Uxmode_xc', num2str(xc_vec(k)*100), '.png');
    exportgraphics(gcf,figfilename,'Resolution', 300)
end