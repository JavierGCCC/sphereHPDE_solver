function result = check_transient_to_cw(showFigure)

% CHECK_TRANSIENT_TO_CW
%
% Verifies that the transient solution under continuous-wave
% illumination converges towards the stationary CW solution
% as the total simulation time increases.
%
% INPUT
%   showFigure  Logical value indicating whether the convergence
%               figure is displayed.
%
% OUTPUT
%   result       Structure containing validation information:
%
%               result.name
%               result.pass
%               result.error
%               result.threshold
%               result.simulationTime
%               result.residual


% ============================================================
% INPUT
% ============================================================

if nargin < 1
    showFigure = true;
end


% ============================================================
% REFERENCE CASE
% ============================================================

% Structure and simulation region
dom.a = 25e-9;
dom.R_sim = 40*dom.a;


% Material properties
mat.n_core   = 'Au';
mat.n_medium = 1.33;

mat.rho = [19300,1000];
mat.cp  = [129,4181];
mat.k   = [318,0.6];
mat.ITC = Inf;
mat.tau1 = 1.7e-12;

% Source parameters
src.lda0 = linspace(500,600,100)*1e-9;
src.F = 1e8;
src.tau2 = Inf; % Continuous-wave illumination

% Spatial and temporal discretization
mesh.N_p = 250;
mesh.N_e = round( ...
    (dom.R_sim-dom.a)/dom.a * mesh.N_p);

% Fixed number of temporal points for all transient calculations.
mesh.N_t = 10000;


% Boundary condition.
bc.bound = 'Robin';
bc.T_ini = 0;


% ============================================================
% SIMULATION TIMES
% ============================================================

% Total transient simulation times.
sim_time = linspace(100,10000,5)*1e-9;

% Allocate error array
residual = zeros(size(sim_time));


% ============================================================
% CW REFERENCE SOLUTION
% ============================================================

solver.mode = 'cw';
solver.preconditioner.enabled = false;

solution_cw = heatEqSolver( ...
    dom,mat,src,mesh,bc,solver);


% Reconstruct complete CW temperature distribution
%
% Dimensions:
%
% opticalBasis        Nr x 1
% spectralWeights      1 x Nlambda
%
% T_cw                Nr x Nlambda

T_cw = ...
    solution_cw.opticalBasis * ...
    solution_cw.spectralWeights + ...
    solution_cw.boundaryBasis;


% Maximum temperature for every wavelength
Tmax_cw = max(T_cw,[],1);


% ============================================================
% TRANSIENT -> CW CONVERGENCE
% ============================================================

solver.mode = 'transient';

% Disable transient waitbar during validation
solver.showWaitbar = false;

for j = 1:length(sim_time)
   
    dom.t_max = sim_time(j); % Set total transient simulation time
    solution_t = heatEqSolver( ...
        dom,mat,src,mesh,bc,solver); % Solve transient thermal problem

    % Reconstruct the temperature distribution at the
    % final simulation time for every wavelength
    %
    % opticalBasis(:,end)     Nr x 1
    % spectralWeights          1 x Nlambda
    %
    % T_t                      Nr x Nlambda

    T_t = ...
        solution_t.opticalBasis(:,end) * ...
        solution_t.spectralWeights + ...
        solution_t.boundaryBasis(:,end);

    % Maximum temperature for every wavelength
    Tmax_t = max(T_t,[],1);

    % Relative difference between transient and CW solutions
    relativeError = ...
        abs(Tmax_cw-Tmax_t) ./ ...
        abs(Tmax_cw);

    % Maximum error across the complete wavelength range
    residual(j) = ...
        100*max(relativeError);

end


% ============================================================
% VALIDATION CRITERION
% ============================================================

% Maximum accepted final relative error [%]
threshold = 1;

% Error at the longest simulation time
finalError = residual(end);

% PASS / FAIL
pass = finalError <= threshold;


% ============================================================
% OUTPUT STRUCTURE
% ============================================================

result.name = 'Transient -> CW';

result.pass = pass;

result.error = finalError;

result.threshold = threshold;

result.simulationTime = sim_time;

result.residual = residual;


% ============================================================
% COMMAND WINDOW OUTPUT
% ============================================================

if pass

    fprintf( ...
        'Transient -> CW: PASS | error = %.4f %% (threshold = %.2f %%)\n', ...
        finalError,threshold);

else

    fprintf( ...
        'Transient -> CW: FAIL | error = %.4f %% (threshold = %.2f %%)\n', ...
        finalError,threshold);

end


% ============================================================
% CONVERGENCE FIGURE
% ============================================================

if showFigure

    figure;

    plot( ...
        sim_time*1e9, ...
        residual, ...
        '-o', ...
        'LineWidth',1.5);

    hold on;

    yline( ...
        threshold, ...
        '--', ...
        'Threshold', ...
        'LineWidth',1.2);

    xlabel('Simulation time [ns]');
    ylabel('Maximum relative error [%]');

    title('Transient convergence towards CW solution');

    grid on;
    box on;

    set(gca,'FontSize',18);

end


end