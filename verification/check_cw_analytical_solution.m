function result = check_cw_analytical_solution(showFigure)

% CHECK_CW_ANALYTICAL_SOLUTION
%
% Verifies the stationary CW thermal solver against the analytical
% solution for a uniformly heated sphere embedded in an infinite medium.
%
% Two cases are tested:
%
%   1) Perfect thermal contact:
%          ITC = Inf
%
%   2) Finite interfacial thermal conductance:
%          ITC < Inf
%
% The optical absorption cross section is taken directly from the
% numerical solver. Therefore, this checker validates only the
% thermal formulation.
%
% OUTPUT
%   result.name
%   result.pass
%   result.error
%   result.threshold
%   result.errorPerfectITC
%   result.errorFiniteITC
%   result.lambda
%   result.TnumericalPerfect
%   result.TanalyticalPerfect
%   result.TnumericalFinite
%   result.TanalyticalFinite

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
% [particle, environment]
mat.k = [318,0.6];
mat.tau1 = 1.7e-12;

% Source parameters
src.lda0 = linspace(500,600,100)*1e-9;
src.F = 1e8;% CW intensity
src.tau2 = Inf;% Continuous-wave illumination

% Spatial/temporal discretization
mesh.N_p = 250;
mesh.N_e = round((dom.R_sim-dom.a)/dom.a * mesh.N_p);
mesh.N_t = 10000;

% Boundary condition
bc.bound = 'Robin';
bc.T_ini = 0;

% Solver
solver.mode = 'cw';
solver.preconditioner.enabled = false;

% ============================================================
% CONSTANTS
% ============================================================
a  = dom.a;
kp = mat.k(1);
km = mat.k(2);
I = src.F;
lambda = src.lda0;

% ============================================================
% CASE 1 — PERFECT THERMAL CONTACT
% ============================================================
mat.ITC = Inf;

% Numerical solution
solution_inf = heatEqSolver(dom,mat,src,mesh,bc,solver);

% Reconstruct full spectral temperature
T_inf = ...
    solution_inf.opticalBasis * ...
    solution_inf.spectralWeights + ...
    solution_inf.boundaryBasis;

% Numerical maximum temperature
Tmax_num_inf = max(T_inf,[],1);

% Absorbed optical power
Pabs_inf = ...
    I * solution_inf.spectralWeights;

% Analytical maximum temperature
R_particle = 1/(8*pi*kp*a);
R_medium = 1/(4*pi*km*a);
Tmax_ana_inf = Pabs_inf.*(R_particle+R_medium);

% ============================================================
% CASE 2 — FINITE INTERFACIAL THERMAL CONDUCTANCE
% ============================================================
G = 2e8;
mat.ITC = G;

% Numerical solution
solution_G = heatEqSolver(dom,mat,src,mesh,bc,solver);

% Reconstruct full spectral temperature
T_G = ...
    solution_G.opticalBasis * ...
    solution_G.spectralWeights + ...
    solution_G.boundaryBasis;

% Numerical maximum temperature
Tmax_num_G = max(T_G,[],1);

% Absorbed optical power
Pabs_G = I * solution_G.spectralWeights;

% Interface thermal resistance
R_interface = 1/(4*pi*a^2*G);

% Analytical maximum temperature
Tmax_ana_G = Pabs_G.*(R_particle+R_medium+R_interface);

% ============================================================
% ERROR CALCULATION
% ============================================================
relativeErrorInf = abs(Tmax_num_inf-Tmax_ana_inf)./max(abs(Tmax_ana_inf),eps);
relativeErrorG = abs(Tmax_num_G-Tmax_ana_G)./max(abs(Tmax_ana_G),eps);
errorInf = 100*max(relativeErrorInf);
errorG = 100*max(relativeErrorG);

% ============================================================
% VALIDATION CRITERION
% ============================================================
% Maximum accepted error [%]
threshold = 1;
passInf = errorInf <= threshold;
passG = errorG <= threshold;
pass = passInf && passG;

% ============================================================
% OUTPUT STRUCTURE
% ============================================================
result.name = 'CW analytical solution';
result.pass = pass;
result.error = max(errorInf,errorG);
result.threshold = threshold;
result.errorPerfectITC = errorInf;
result.errorFiniteITC = errorG;
result.lambda = lambda;
result.TnumericalPerfect = Tmax_num_inf;
result.TanalyticalPerfect = Tmax_ana_inf;
result.TnumericalFinite = Tmax_num_G;
result.TanalyticalFinite = Tmax_ana_G;
result.ITCfinite = G;

% ============================================================
% COMMAND WINDOW OUTPUT
% ============================================================
if pass
    fprintf( ...
        ['CW analytical solution: PASS | ' ...
         'G = Inf: %.4f %% | ' ...
         'G = %.2e: %.4f %%\n'], ...
        errorInf,G,errorG);
else
    fprintf( ...
        ['CW analytical solution: FAIL | ' ...
         'G = Inf: %.4f %% | ' ...
         'G = %.2e: %.4f %%\n'], ...
        errorInf,G,errorG);
end

% ============================================================
% FIGURE
% ============================================================
if showFigure
    figure; 
    subplot(1,2,1)% Perfect thermal contact
    plot(lambda*1e9, Tmax_num_inf, 'LineWidth',1.8);
    hold on;
    plot(lambda*1e9, Tmax_ana_inf, '--', 'LineWidth',1.8);
    xlabel('\lambda [nm]');
    ylabel('\DeltaT_{max} [K]');
    title('G = \infty');
    legend('Numerical', 'Analytical', 'Location','best');
    grid on;
    box on;
    subplot(1,2,2)% Finite ITC
    plot(lambda*1e9, Tmax_num_G, 'LineWidth',1.8);
    hold on;
    plot(lambda*1e9, Tmax_ana_G, '--', 'LineWidth',1.8);
    xlabel('\lambda [nm]');
    ylabel('\DeltaT_{max} [K]');
    title(sprintf('G = %.1e W m^{-2} K^{-1}',G));
    legend('Numerical', 'Analytical', 'Location','best');
    grid on;
    box on;
    set(findall(gcf,'Type','axes'), 'FontSize',18);
end
end