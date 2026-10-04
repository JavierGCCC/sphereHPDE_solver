%% example_cw_itc_gold_water.m
%
% CW heating of a gold nanosphere in water.
%
% Comparison between:
%
%   1) Perfect thermal contact: G = Inf
%   2) Finite interfacial thermal conductance:
%          G = 2e8 W m^-2 K^-1
%
% The radial temperature profiles are reconstructed using
% temperatureProfile().
% ============================================================

clear; clc; close all; 
initialize; %Folder routes.

% ============================================================
% DEFINITION OF PHYSICAL PARAMETERS
% ============================================================

% Geometry
dom.a = 25e-9;
dom.R_sim = 40*dom.a;

%Material properties.
mat.n_core = 'Au';
mat.n_medium = 1.33;
mat.rho = [19300 1000];
mat.cp = [129 4181];
mat.k = [318 0.6];
mat.tau1 = 1.7e-12;%electron-phonon relaxation time.

%Source properties.
src.lda0 = 532e-9; %Wavelength selection.
src.F = 1e9;% Irradiance.
src.tau2 = Inf; %cw source.

%Mesh properties.
mesh.N_p = 250;
mesh.N_e = round((dom.R_sim-dom.a)/dom.a*mesh.N_p);

%Boundary and initial conditions.
bc.bound = 'Robin';
bc.T_ini = 0;

%Solver options
solver.mode = 'cw';
solver.preconditioner.enabled = false;

% ============================================================
% CASE 1 - PERFECT THERMAL CONTACT
% ============================================================

mat.ITC =Inf;
solution_inf = heatEqSolver(dom, mat, src, mesh, bc, solver);
[r,T_inf,info_inf] = temperatureProfile(solution_inf);

% ============================================================
% CASE 2 - FINITE INTERFACIAL THERMAL CONDUCTANCE
% ============================================================

G = 2e8;
mat.ITC = G;
solution_G = heatEqSolver(dom, mat, src, mesh, bc, solver);
[~,T_G,info_G] = temperatureProfile(solution_G);

% ============================================================
% MAXIMUM RELATIVE DIFFERENCES. 
% ============================================================
Tcenter_inf = T_inf(1);
Tcenter_G = T_G(1);
epsilon_ITC = (Tcenter_G-Tcenter_inf)/Tcenter_inf;
epsilon_ITC_percent = 100*epsilon_ITC;

% ============================================================
% PRINT RESULTS
% ============================================================
fprintf('\n');
fprintf('=============================================\n');
fprintf(' CW GOLD SPHERE IN WATER\n');
fprintf('=============================================\n');
fprintf('Radius: %.1f nm\n', dom.a*1e9);
fprintf('Wavelength: %.1f nm\n', info_inf.lambda*1e9);
fprintf('Irradiance: %.3e W/m^2\n', src.F);
fprintf('Finite ITC: %.3e W/m^2/K\n', G);
fprintf('\n');
fprintf('T_center (G = Inf): %.6f K\n',Tcenter_inf);
fprintf('T_center (finite G): %.6f K\n', Tcenter_G);
fprintf( 'Relative ITC difference:    %.3f %%\n', epsilon_ITC_percent);
fprintf('=============================================\n\n');

% ============================================================
% PLOT
% ============================================================
figure('Color','w', 'Position',[250 150 850 600]);
hold on;

% Perfect thermal contact
plot(r/dom.a, T_inf, 'LineWidth',2.2);

% Finite ITC
plot(r/dom.a, T_G, '--', 'LineWidth',2.2);


% Particle-medium interface
xline(1, ':', 'LineWidth',1.4);
xlabel('r/a', 'FontSize',16);
ylabel('\DeltaT(r) [K]', 'FontSize',16);
title('CW heating of a gold nanosphere in water', 'FontSize',16);
legend('G = \infty', 'G = 2\times10^8 W m^{-2} K^{-1}', 'Particle surface', ...
    'Location','northeast');
box on;
grid off;
set(gca, 'FontSize',15, 'LineWidth',1.1);
xlim([0 dom.R_sim/dom.a]);
ylim([0 1.10*max(T_G)]);
annotationText = sprintf(['$\\frac{T_c^{G}-T_c^{\\infty}}' ...
     '{T_c^{\\infty}} = %.2f\\,\\%%$'], ...
    epsilon_ITC_percent);
text( ...
    0.58, ...
    0.80, ...
    annotationText, ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'FontSize',16, ...
    'BackgroundColor','w', ...
    'EdgeColor',[0.7 0.7 0.7], ...
    'Margin',8);
hold off;