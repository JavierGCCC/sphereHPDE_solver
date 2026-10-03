function result = check_quasistatic_limit(showFigure)

% CHECK_QUASISTATIC_LIMIT
%
% Verifies that the Mie solution and the numerical CW thermal
% solution converge towards the quasi-static approximation as
% the particle size decreases.
%
% OUTPUT
%   result.name
%   result.pass
%   result.error
%   result.threshold
%   result.radii
%   result.opticalError
%   result.temperatureError
%   result.thermalSolverError


%% ============================================================
% INPUT
% ============================================================

if nargin < 1
    showFigure = true;
end


%% ============================================================
% REFERENCE CASE
% ============================================================

% Particle radii used to approach the quasi-static limit
radii = [2.5 5 7.5 10 12.5 15 17.5 20]*1e-9;


% Material properties
mat.n_core   = 'Au';
mat.n_medium = 1.33;

mat.rho = [19300,1000];
mat.cp  = [129,4181];
mat.k   = [318,0.6];

mat.ITC = Inf;

% Only required by the general solver structure
mat.tau1 = 1.7e-12;


% Spectral range
src.lda0 = linspace(500,600,100)*1e-9;

% CW intensity
src.F = 1e9;

% Continuous-wave illumination
src.tau2 = Inf;


% Spatial discretization
mesh.N_p = 250;


% Boundary condition
bc.bound = 'Robin';
bc.T_ini = 0;


% Solver
solver.mode = 'cw';
solver.preconditioner.enabled = false;


%% ============================================================
% INITIALIZATION
% ============================================================

lambda = src.lda0;

Nradii = length(radii);
Nlambda = length(lambda);

sigma_mie = zeros(Nradii,Nlambda);
sigma_qa  = zeros(Nradii,Nlambda);

Tmax_num = zeros(Nradii,Nlambda);

Tmax_qa = zeros(Nradii,Nlambda);

% Analytical temperature calculated using the Mie absorption.
% This is useful for checking the thermal solver independently
% from the quasi-static optical approximation.
Tmax_mie_analytical = zeros(Nradii,Nlambda);


%% ============================================================
% MEDIUM OPTICAL PROPERTIES
% ============================================================

eps_m = mat.n_medium^2;

k_m = ...
    2*pi*mat.n_medium ./ lambda;


%% ============================================================
% RADIUS SWEEP
% ============================================================

for j = 1:Nradii

    %% --------------------------------------------------------
    % Geometry
    % ---------------------------------------------------------

    dom.a = radii(j);

    dom.R_sim = 40*dom.a;


    %% --------------------------------------------------------
    % Spatial mesh
    % ---------------------------------------------------------

    mesh.N_e = round( ...
        (dom.R_sim-dom.a)/dom.a * mesh.N_p);


    %% --------------------------------------------------------
    % Optical calculation
    % ---------------------------------------------------------

    for i = 1:Nlambda

        % Full Mie solution
        [n_particle,sigma_mie(j,i)] = ...
            mie_absorption( ...
                dom.a, ...
                mat.n_core, ...
                mat.n_medium, ...
                lambda(i));


        % Particle dielectric function
        eps_s = n_particle^2;


        % Quasi-static polarizability factor
        ratio = ...
            (eps_s-eps_m) / ...
            (eps_s+2*eps_m);


        % Quasi-static absorption cross section
        %
        % sigma_abs = sigma_ext - sigma_sca

        sigma_qa(j,i) = ...
            4*pi*k_m(i)*dom.a^3*imag(ratio) ...
            - ...
            (8*pi/3)*k_m(i)^4*dom.a^6*abs(ratio)^2;

    end


    %% --------------------------------------------------------
    % Numerical CW thermal solution
    % ---------------------------------------------------------

    solution = heatEqSolver( ...
        dom,mat,src,mesh,bc,solver);


    % Reconstruct complete spectral temperature
    T_num = ...
        solution.opticalBasis * ...
        solution.spectralWeights + ...
        solution.boundaryBasis;


    % Maximum temperature for each wavelength
    Tmax_num(j,:) = ...
        max(T_num,[],1);


    %% --------------------------------------------------------
    % Analytical stationary thermal resistance
    % ---------------------------------------------------------

    % Heat diffusion through surrounding medium
    Rth = ...
        1/(4*pi*mat.k(2)*dom.a);


    % Add interfacial thermal resistance when ITC is finite
    if ~isinf(mat.ITC)

        Rth = ...
            Rth + ...
            1/(4*pi*dom.a^2*mat.ITC);

    end


    %% --------------------------------------------------------
    % Analytical temperatures
    % ---------------------------------------------------------

    % Analytical thermal solution using full Mie absorption
    Tmax_mie_analytical(j,:) = ...
        src.F * sigma_mie(j,:) * Rth;


    % Full quasi-static prediction
    Tmax_qa(j,:) = ...
        src.F * sigma_qa(j,:) * Rth;

end


%% ============================================================
% ERROR CALCULATION
% ============================================================

opticalError = zeros(1,Nradii);

temperatureError = zeros(1,Nradii);

thermalSolverError = zeros(1,Nradii);


for j = 1:Nradii

    %% Mie vs quasi-static optics

    opticalError(j) = ...
        100 * max( ...
        abs(sigma_mie(j,:)-sigma_qa(j,:)) ./ ...
        max(abs(sigma_mie(j,:)),eps));


    %% Numerical temperature vs complete QA prediction

    temperatureError(j) = ...
        100 * max( ...
        abs(Tmax_num(j,:)-Tmax_qa(j,:)) ./ ...
        max(abs(Tmax_num(j,:)),eps));


    %% Thermal solver only:
    % numerical CW vs analytical CW using the SAME Mie sigma_abs

    thermalSolverError(j) = ...
        100 * max( ...
        abs(Tmax_num(j,:)-Tmax_mie_analytical(j,:)) ./ ...
        max(abs(Tmax_mie_analytical(j,:)),eps));

end


%% ============================================================
% VALIDATION CRITERION
% ============================================================

% The smallest radius represents the quasi-static limit
limitIndex = 1;

% Accepted error [%]
threshold = 1;


% Require both optical and temperature predictions to agree
% with the quasi-static approximation at the smallest radius.

limitOpticalError = opticalError(limitIndex);

limitTemperatureError = temperatureError(limitIndex);


% Also check that approaching smaller particles improves
% the agreement relative to the largest particle.

convergingOptics = ...
    opticalError(1) < opticalError(end);

convergingTemperature = ...
    temperatureError(1) < temperatureError(end);


pass = ...
    limitOpticalError <= threshold && ...
    limitTemperatureError <= threshold && ...
    convergingOptics && ...
    convergingTemperature;


%% ============================================================
% OUTPUT STRUCTURE
% ============================================================

result.name = 'Quasi-static limit';

result.pass = pass;

% Main reported error
result.error = max( ...
    limitOpticalError, ...
    limitTemperatureError);

result.threshold = threshold;

result.radii = radii;

result.opticalError = opticalError;

result.temperatureError = temperatureError;

result.thermalSolverError = thermalSolverError;

result.sigmaMie = sigma_mie;

result.sigmaQA = sigma_qa;

result.TmaxNumerical = Tmax_num;

result.TmaxQA = Tmax_qa;

result.lambda = lambda;


%% ============================================================
% COMMAND WINDOW OUTPUT
% ============================================================

if pass

    fprintf( ...
        ['Quasi-static limit: PASS | ' ...
         'optical error = %.4f %% | ' ...
         'temperature error = %.4f %%\n'], ...
        limitOpticalError, ...
        limitTemperatureError);

else

    fprintf( ...
        ['Quasi-static limit: FAIL | ' ...
         'optical error = %.4f %% | ' ...
         'temperature error = %.4f %%\n'], ...
        limitOpticalError, ...
        limitTemperatureError);

end


%% ============================================================
% FIGURES
% ============================================================

if showFigure

    %% --------------------------------------------------------
    % Error vs particle radius
    % ---------------------------------------------------------

    figure;

    plot( ...
        radii*1e9, ...
        opticalError, ...
        '-o', ...
        'LineWidth',1.5);

    hold on;

    plot( ...
        radii*1e9, ...
        temperatureError, ...
        '-s', ...
        'LineWidth',1.5);

    yline( ...
        threshold, ...
        '--', ...
        'Threshold');

    xlabel('Particle radius [nm]');
    ylabel('Maximum relative error [%]');

    legend( ...
        'Mie vs quasi-static', ...
        'Numerical T vs quasi-static T', ...
        'Location','best');

    title('Convergence towards the quasi-static limit');

    grid on;
    box on;

    set(gca,'FontSize',18);


    %% --------------------------------------------------------
    % Spectrum at the smallest particle radius
    % ---------------------------------------------------------

    figure;

    subplot(1,2,1)

    plot( ...
        lambda*1e9, ...
        sigma_mie(1,:), ...
        'LineWidth',1.8);

    hold on;

    plot( ...
        lambda*1e9, ...
        sigma_qa(1,:), ...
        '--', ...
        'LineWidth',1.8);

    xlabel('\lambda [nm]');
    ylabel('\sigma_{abs} [m^2]');

    legend( ...
        'Mie', ...
        'Quasi-static', ...
        'Location','best');

    grid on;
    box on;


    subplot(1,2,2)

    plot( ...
        lambda*1e9, ...
        Tmax_num(1,:), ...
        'LineWidth',1.8);

    hold on;

    plot( ...
        lambda*1e9, ...
        Tmax_qa(1,:), ...
        '--', ...
        'LineWidth',1.8);

    xlabel('\lambda [nm]');
    ylabel('\DeltaT_{max} [K]');

    legend( ...
        'Numerical', ...
        'Quasi-static', ...
        'Location','best');

    grid on;
    box on;

    set(findall(gcf,'Type','axes'), ...
        'FontSize',18);

end


end