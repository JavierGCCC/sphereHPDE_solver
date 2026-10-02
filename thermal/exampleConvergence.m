clear;
clc;
close all;

% ============================================================
% PARAMETERS
% ============================================================

dom.a = 25e-9;
dom.R_sim = 40*dom.a;

mat.n_core   = 'Au';
mat.n_medium = 1.33;

mat.rho = [19300,1000];
mat.cp  = [129,4181];
mat.k   = [318,0.6];

mat.ITC = Inf;

mat.tau1 = 1.7e-12;

src.lda0 = linspace(500,600,100)*1e-9;

src.F = 1e8;

% Constant illumination
src.tau2 = Inf;


% ============================================================
% SPATIAL DISCRETIZATION
% ============================================================

mesh.N_p = 250;

mesh.N_e = round( ...
    (dom.R_sim-dom.a)/dom.a * mesh.N_p);

mesh.N_t = 10000;


% ============================================================
% BOUNDARY CONDITION
% ============================================================

bc.bound = 'Robin';
bc.T_ini = 0;


% ============================================================
% CW REFERENCE
%
% CW does NOT depend on simulation duration,
% therefore calculate it only ONCE.
% ============================================================

solver.mode = 'cw';

solution_cw = heatEqSolver( ...
    dom,mat,src,mesh,bc,solver);


% Reconstruct complete CW spectrum

Tcw = ...
    solution_cw.opticalBasis * ...
    solution_cw.spectralWeights + ...
    solution_cw.boundaryBasis;


% Maximum temperature for every wavelength

Tmax_cw = max(Tcw,[],1);


% ============================================================
% TRANSIENT SIMULATION DURATIONS
% ============================================================

duracion = linspace(100,10000,25)*1e-9;

residual = zeros(size(duracion));


% ============================================================
% TRANSIENT -> CW CONVERGENCE
% ============================================================

solver.mode = 'transient';
solver.showWaitbar = true;

for j = 1:length(duracion)

    % Simulation duration

    dom.t_max = duracion(j);


    % Solve transient problem

    solution_t = heatEqSolver( ...
        dom,mat,src,mesh,bc,solver);


    % Reconstruct temperature at FINAL time
    %
    % Dimensions:
    %
    % opticalBasis(:,end)    Nr x 1
    % spectralWeights         1 x Nlambda
    %
    % --> Nr x Nlambda

    Ttr_final = ...
        solution_t.opticalBasis(:,end) * ...
        solution_t.spectralWeights + ...
        solution_t.boundaryBasis(:,end);


    % Maximum spatial temperature for every wavelength

    Tmax_t = max(Ttr_final,[],1);


    % Maximum relative spectral difference

    relativeError = ...
        abs(Tmax_t-Tmax_cw) ./ ...
        abs(Tmax_cw);

    residual(j) = ...
        100*max(relativeError);


    % Progress

    fprintf( ...
        't_max = %8.2f ns | error = %.6f %%\n', ...
        dom.t_max*1e9, ...
        residual(j));

end


% ============================================================
% CONVERGENCE FIGURE
% ============================================================

figure;

plot( ...
    duracion*1e9, ...
    residual, ...
    '-o', ...
    'LineWidth',1.8);

xlabel('Simulation duration [ns]');
ylabel('Maximum CW relative error [%]');

grid on;
box on;

set(gca,'FontSize',22);