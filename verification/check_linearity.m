function result = check_linearity(showFigure)

% CHECK_LINEARITY
%
% Verifies the linearity of the thermal solver with respect to the
% incident optical intensity.
%
% The checker tests both:
%
%   1) Continuous-wave heating
%   2) Transient heating
%
% For a linear heat equation:
%
%       T(alpha*I) = alpha*T(I)
%
% Here alpha = 2.
%
% OUTPUT
%   result.name
%   result.pass
%   result.error
%   result.threshold
%   result.errorCW
%   result.errorTransient


%% ============================================================
% INPUT
% ============================================================

if nargin < 1
    showFigure = true;
end


%% ============================================================
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


% Use one representative wavelength.
% Linearity does not depend on the spectral sampling.
src.lda0 = 550e-9;


% Reference intensity
I0 = 1e9;

% Scaling factor
alpha = 2;


% Spatial discretization
mesh.N_p = 250;

mesh.N_e = round( ...
    (dom.R_sim-dom.a)/dom.a * mesh.N_p);


% Boundary condition
bc.bound = 'Robin';
bc.T_ini = 0;


%% ============================================================
% CW LINEARITY
% ============================================================

solver.mode = 'cw';
solver.preconditioner.enabled = false;

src.tau2 = Inf;


%% Reference intensity

src.F = I0;

solution_cw_1 = heatEqSolver( ...
    dom,mat,src,mesh,bc,solver);


Tcw_1 = ...
    solution_cw_1.spectralWeights * ...
    solution_cw_1.opticalBasis + ...
    solution_cw_1.boundaryBasis;


%% Scaled intensity

src.F = alpha*I0;

solution_cw_2 = heatEqSolver( ...
    dom,mat,src,mesh,bc,solver);


Tcw_2 = ...
    solution_cw_2.spectralWeights * ...
    solution_cw_2.opticalBasis + ...
    solution_cw_2.boundaryBasis;


%% CW linearity error

errorCW = ...
    100 * max( ...
    abs(Tcw_2-alpha*Tcw_1) ./ ...
    max(abs(alpha*Tcw_1),eps));


%% ============================================================
% TRANSIENT LINEARITY
% ============================================================

solver.mode = 'transient';
solver.showWaitbar = false;


% Finite pulse to exercise the transient source implementation
src.tau2 = 100*mat.tau1;


% Simulate sufficiently beyond the pulse
dom.t_max = 10*src.tau2;


% Temporal discretization
mesh.N_t = 10000;


%% Reference intensity

src.F = I0;

solution_t_1 = heatEqSolver( ...
    dom,mat,src,mesh,bc,solver);


Tt_1 = ...
    solution_t_1.spectralWeights * ...
    solution_t_1.opticalBasis + ...
    solution_t_1.boundaryBasis;


%% Scaled intensity

src.F = alpha*I0;

solution_t_2 = heatEqSolver( ...
    dom,mat,src,mesh,bc,solver);


Tt_2 = ...
    solution_t_2.spectralWeights * ...
    solution_t_2.opticalBasis + ...
    solution_t_2.boundaryBasis;


%% Transient linearity error

errorTransient = ...
    100 * max( ...
    abs(Tt_2(:)-alpha*Tt_1(:)) ./ ...
    max(abs(alpha*Tt_1(:)),eps));


%% ============================================================
% VALIDATION CRITERION
% ============================================================

threshold = 1e-6;   % [%]


passCW = ...
    errorCW <= threshold;

passTransient = ...
    errorTransient <= threshold;


pass = ...
    passCW && passTransient;


%% ============================================================
% OUTPUT STRUCTURE
% ============================================================

result.name = 'Thermal linearity';

result.pass = pass;

result.error = ...
    max(errorCW,errorTransient);

result.threshold = threshold;

result.errorCW = errorCW;

result.errorTransient = errorTransient;

result.intensity = I0;

result.scalingFactor = alpha;


%% ============================================================
% COMMAND WINDOW OUTPUT
% ============================================================

if pass

    fprintf( ...
        ['Thermal linearity: PASS | ' ...
         'CW = %.3e %% | ' ...
         'Transient = %.3e %%\n'], ...
        errorCW,errorTransient);

else

    fprintf( ...
        ['Thermal linearity: FAIL | ' ...
         'CW = %.3e %% | ' ...
         'Transient = %.3e %%\n'], ...
        errorCW,errorTransient);

end


%% ============================================================
% FIGURES
% ============================================================

if showFigure

    %% CW

    figure;

    plot( ...
        solution_cw_1.nodes*1e9, ...
        alpha*Tcw_1, ...
        'LineWidth',1.8);

    hold on;

    plot( ...
        solution_cw_2.nodes*1e9, ...
        Tcw_2, ...
        '--', ...
        'LineWidth',1.8);

    xlabel('r [nm]');
    ylabel('\DeltaT [K]');

    title('CW linearity');

    legend( ...
        '\alpha T(I)', ...
        'T(\alpha I)', ...
        'Location','best');

    grid on;
    box on;

    set(gca,'FontSize',18);


    %% Transient maximum temperature

    Tmax_t1 = max(Tt_1,[],1);
    Tmax_t2 = max(Tt_2,[],1);


    figure;

    plot( ...
        solution_t_1.timeDomain*1e9, ...
        alpha*Tmax_t1, ...
        'LineWidth',1.8);

    hold on;

    plot( ...
        solution_t_2.timeDomain*1e9, ...
        Tmax_t2, ...
        '--', ...
        'LineWidth',1.8);

    xlabel('Time [ns]');
    ylabel('\DeltaT_{max} [K]');

    title('Transient linearity');

    legend( ...
        '\alpha T(I)', ...
        'T(\alpha I)', ...
        'Location','best');

    grid on;
    box on;

    set(gca,'FontSize',18);

end


end