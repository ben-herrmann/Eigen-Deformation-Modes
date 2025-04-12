clear all; close all; clc;
load("BatteryData.mat")
%% Mass matrix and method application

% We assemble mass matrix and compute its Cholesky decomposition
% Then we apply the method to obtain EDMs (computed using meaningful dot product) and EDM coefficients

FEM_M = assembleFEMatrices(model, 'M');
[R,flag] =  chol(FEM_M.M);
% Method
m=6; % Number of eigenmodes used for model reduction
modes_matrix_m = modes_matrix(:,1:m,:);
lambdas_m = lambdas(:,1:m);
[def_modes, coeffs] = themethod(modes_matrix_m, R, 2);
ref_basis = modes_matrix_m(:,:,1);

%% Comparison against ROM solution interpolation

% Comparison: direct interpolation of the local basis and interpolation of
% ROM soutions
% We apply initial condition to the model: Steady state equilibria for a parameter
% outside the training set

[n,m,p] = size(modes_matrix_m);

% Change Initial condition
h_ic = 100; %300 We set the initial condition as the fixed point for h_ic
model.AnalysisType = 'SteadyState';
thermalBC(model, ...
          Face=[boundaryIDs(1).FrontFace,boundaryIDs(end).BackFace], ...
          ConvectionCoefficient=h_ic, ...
          AmbientTemperature=Tambient);
RS_ic = solve(model);
thermalIC(model, RS_ic);
ic_fs = RS_ic.Temperature;

%% Compute some validation data

% Steady state equilibria and some validation trajectories
% The validation points will be the mid point between sampled parameters
% (p-1=3 validation points)

% Calculate validation data
t0 = 1/lambdas(1,1);
tlist=linspace(0,5*t0, 1000);

h_val_vec = zeros(p-1,1);
val_fixed_points = zeros(n,p-1); %array containing validation fixed points
val_ground_truth = zeros(n, length(tlist), p-1); %array containing validation trayectories
for j=1:p-1
    % Compute Validation Fixed Points
    model.AnalysisType = 'SteadyState';
    h_val_vec(j) = mean([h_vec(j) h_vec(j+1)]);
    thermalBC(model, ...
          Face=[boundaryIDs(1).FrontFace,boundaryIDs(end).BackFace], ...
          ConvectionCoefficient=h_val_vec(j), ...
          AmbientTemperature=Tambient);
    RS_vali = solve(model);
    val_fixed_points(:,j) = RS_vali.Temperature;
    % Compute Validation Trajectory
    model.AnalysisType='Transient';
    tic
    Rval = solve(model, tlist);
    toc
    val_ground_truth(:,:,j) = Rval.Temperature;
end

%% ROM estimation at training points and interpolation at validation points

% We calculate ROM estimation at training points and interpolate

% ROM estimation at training points and interpolation at validation points

roms_par = zeros(m,length(tlist),p); % Reduced State Snapshots for each p
perturbacion = zeros(n,length(tlist),p);
roms_recons_par = zeros(n,length(tlist),p); % Full State Snapshots for each p
for i=1:p
    Vi = modes_matrix_m(:,:,i);
    fp_i = fixed_points(:,i);
    %fp_i = val_fixed_points(:,1); %(borrar) Que los ROMs de dataset esten calculados con el eq de val 1 y asi se saca ese factor de la comparacion
    lambi = lambdas_m(i,:)';
    ic_i = Vi'*FEM_M.M*(ic_fs - fp_i);
    roms_par(:,:,i) = my_simulate_ROM(ic_i, lambi, tlist);
    perturbacion(:,:,i) = Vi*roms_par(:,:,i);% + fp_i;
end
interp_rom = zeros(m,length(tlist), p-1); % Reduced State Snapshots for each Validation Parameter
interp_pert = zeros(n,length(tlist), p-1); % Full State perturbation for each Validation Parameter
interp_full_state = zeros(n,length(tlist), p-1); %Full State Reconstruction for each val parameter
for j=1:p-1
    h_val = h_val_vec(j);
    fp_val_j = val_fixed_points(:,j);
    for k=1:length(tlist)
        coords_k = permute(roms_par(:,k,:), [1 3 2]); % Fixed k-th snapshot for each parameter
        interp_rom(:,k,j) = interp1(h_vec, coords_k', h_val, 'spline');
        pert_k = permute(perturbacion(:,k,:), [1 3 2]); % % Fixed k-th snapshot for each parameter
        interp_pert(:,k,j) = interp1(h_vec, pert_k', h_val, 'spline');
    end
    interp_full_state(:,:,j) = interp_pert(:,:,j) + fp_val_j;
end

%% Solution interpolation at by reduced state interpolation

% Calculate validation trajectories but by reduced state interpolation
% Here to "lift" the interpolated reduced state a local basis and
% equilibria is needed. Directly interpolated basis and validation
% equilibria is used. (This comparison is not included)

red_state_int = zeros(n,length(tlist), p-1); %
for j=1:p-1
    Vj = interpolate_whole(modes_matrix_m, h_vec, h_val_vec(j));
    fp_j = val_fixed_points(:,j);
    red_state_int(:,:,j) = Vj*interp_rom(:,:,j) + fp_j;
end

%% Method application: Estimated trajectories at validation points using EDM-based interpolation and direct interpolation

% Calculate the ROM estimation for validation parameters using EDM-based
% interpolated modes and directly interpolated modes

di_val = zeros(n,length(tlist), p-1);
ri_val = zeros(n,length(tlist), p-1);
for j=1:p-1
    Vdi = interpolate_whole(modes_matrix_m, h_vec, h_val_vec(j));
    Vri = interpolate_model(def_modes, coeffs, h_vec, ref_basis, h_val_vec(j));
    lamb_i = interp1(h_vec, lambdas_m, h_val_vec(j), 'spline')';
    fp_i = val_fixed_points(:,j);
    % ROM direct interpolation
    ic_di = Vdi'*FEM_M.M*(ic_fs - fp_i);
    rom_di = my_simulate_ROM(ic_di, lamb_i, tlist);
    di_val(:,:,j) = Vdi*rom_di + fp_i;
    % ROM reduced interpolation
    ic_ri = Vri'*FEM_M.M*(ic_fs - fp_i);
    rom_ri = my_simulate_ROM(ic_ri, lamb_i, tlist);
    ri_val(:,:,j) = Vri*rom_ri + fp_i;
end

%% Error comparison

% Error calculation in each case, against Ground Truth (Full Order Model
% (FOM) transient simulation)
% 1) ROM solutions interpolation by interpolating the full state
% 2) ROM solutions interpolation by interpolating the reduced state
% 3) Simulate a ROM using directly interpolated modes
% 4) Simulate a ROM using EDM-based interpolated modes

error_di = zeros(p-1, length(tlist)); %Base interpolada de forma directa
error_ri = zeros(p-1, length(tlist)); %Base interpolada de forma reducida
error_rom_int = zeros(p-1, length(tlist)); %Interpolacion del estado full
error_red_state_int = zeros(p-1, length(tlist)); %Interpolacion del estado reducido
for j=1:p-1
    fp_j = val_fixed_points(:,j); %Fixed point de validacion j
    gt_j = val_ground_truth(:,:,j); %Trayectoria FOM j
    di_j = di_val(:,:,j); %Trayectoria c/ Base DI j
    ri_j = ri_val(:,:,j); %Trayectoria c/ Base RI j
    rom_int_j = interp_full_state(:,:,j); %Trayectoria ROM full state interpolation
    rom_red_int_j = red_state_int(:,:,j); %Trayectoria ROM red. state interpolation
    % Error
    Rfp_j = fp_j; %agregar R Cholesky
    error_di(j,:) = vecnorm((gt_j - di_j))/norm(Rfp_j);
    error_ri(j,:) = vecnorm((gt_j - ri_j))/norm(Rfp_j);
    error_rom_int(j,:) = vecnorm((gt_j - rom_int_j))/norm(Rfp_j);
    error_red_state_int(j,:) = vecnorm((gt_j - rom_red_int_j))/norm(Rfp_j);
end
%% Time calculation and time-integrated error
% Time test
tlist=linspace(0,5*t0, 1000);
par_space = linspace(0.001, 120);

time_rom_int = zeros(1,length(par_space));
time_DI = zeros(1,length(par_space));
time_RI = zeros(1,length(par_space));
time_fom = zeros(1,length(par_space));

integral_rom_int = zeros(1,length(par_space)); % arrays for integrated error over time
integral_DI = zeros(1,length(par_space));
integral_RI = zeros(1,length(par_space));

integral_fom = zeros(1,length(par_space));

for j=1:length(par_space)
    j
    h_val = par_space(j);

    % Transient FOM for validation h
    model.AnalysisType = 'transient';
    thermalBC(model, ...
          Face=[boundaryIDs(1).FrontFace,boundaryIDs(end).BackFace], ...
          ConvectionCoefficient=h_val, ...
          AmbientTemperature=Tambient);
    tic
    RT_val = solve(model, tlist);
    time_fom(j) = toc; %Time fom
    val_gt = RT_val.Temperature;

    %Fixed point calculation
    model.AnalysisType = 'steadystate';
    RSteady = solve(model);
    fp = RSteady.Temperature;

    % ROM Full state interpolation

    perturbacion = zeros(n,length(tlist),p); % Perturbation for each p
    tic % Start stopwatch
    for i=1:p %We simulate for our data set
        Vi = modes_matrix_m(:,:,i);
        fp_i = fixed_points(:,i);
        lambi = lambdas_m(i,:)';
        ic_i = Vi'*FEM_M.M*(ic_fs - fp_i);
        rom_i = my_simulate_ROM(ic_i, lambi,tlist); %Simulate
        perturbacion(:,:,i) = Vi*rom_i; %Lift to full dim
    end
    
    interp_full_state = zeros(n,length(tlist)); % Full State Snapshots for 1 val Parameter
    for k=1:length(tlist)
        pert_k = permute(perturbacion(:,k,:), [1 3 2]); % % Fixed k-th snapshot for each parameter
        interp_full_state(:,k) = interp1(h_vec, pert_k', h_val, 'spline'); % Interpolate perturbation
    end
    interp_full_state = interp_full_state + fp;
    time_rom_int(j) = toc;

    % Direct Interpolation of basis
    tdi=0;
    for k=1:100
        tic %Start Stopwatch
        Vdi = interpolate_whole(modes_matrix_m, h_vec, h_val);
        lamb_aa = interp1(h_vec, lambdas_m, h_val, 'spline')';
        ic_di = Vdi'*FEM_M.M*(ic_fs - fp);
        rom_di = my_simulate_ROM(ic_di,lamb_aa,tlist);
        tdi=tdi+toc;
    end
    time_DI(j) = tdi/100; % Average time of 100 simulations
    recons_di = Vdi*rom_di + fp;

    % Reduced Interpolation of basis
    tri=0;
    for k=1:100
        tic %Start Stopwatch
        Vri = interpolate_model(def_modes, coeffs, h_vec,ref_basis,h_val);
        lamb_aa = interp1(h_vec, lambdas_m, h_val, 'spline')';
        ic_ri = Vri'*FEM_M.M*(ic_fs - fp);
        rom_ri = my_simulate_ROM(ic_ri,lamb_aa,tlist);
        tri=tri+toc;
    end
    time_RI(j) = tri/100; % Average time of 100 simulations
    recons_ri = Vri*rom_ri + fp;

    % Error calculation
    error_rom_int_j = vecnorm(val_gt - interp_full_state);
    error_di_j = vecnorm(val_gt - recons_di);
    error_ri_j = vecnorm(val_gt - recons_ri);
    norm_fom_j = vecnorm(val_gt-fp); %Norma en el tiempo del FOM

    % Integrate error over time
    integral_rom_int(j) = trapz(tlist, error_rom_int_j);
    integral_DI(j) = trapz(tlist, error_di_j);
    integral_RI(j) = trapz(tlist, error_ri_j);
    integral_fom(j) = trapz(tlist, norm_fom_j);
end

%% Calculate average FOM time over the parameter space
par_Fom = 0:12:120;
par_Fom(1) = 0.001;
model.AnalysisType='Transient';
tic
for i=1:length(par_Fom)
    i
    thermalBC(model, ...
          Face=[boundaryIDs(1).FrontFace,boundaryIDs(end).BackFace], ...
          ConvectionCoefficient=par_Fom(i), ...
          AmbientTemperature=Tambient);
    solve(model, tlist);
end
timeL_FOM = toc;
time_FOM=timeL_FOM/length(par_Fom);

%% Save calculations needed for plots
%save BatteryCalculations.mat  tlist error_di error_ri error_rom_int error_red_state_int h_val_vec par_space time_rom_int time_DI time_RI time_FOM time_fom integral_rom_int integral_DI integral_RI integral_fom