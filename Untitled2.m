%% ============================================================
% RECALCULATE r_gamma MAPS + DOMAIN-SIZE CHECK
%
% New heatEqSolver implementation.
%
% Calculates:
%
%   1) r_gamma(gamma,tau2) for G = Inf
%   2) r_gamma(gamma,tau2) for G = 2e8 W m^-2 K^-1
%   3) epsilon_gamma between both ITC cases
%   4) Domain-size sensitivity:
%           R_sim = 40a  vs  R_sim = 80a
%
% ============================================================

clear;
clc;
close all;

initialize;

totalTimer = tic;


%% ============================================================
% USER SETTINGS
% ============================================================

% Original domain and enlarged domain
domainFactors = [40 60];

% The LAST domain is used for the main reproduced figure
mainDomainIndex = numel(domainFactors);

% Checkpoint file
checkpointFile = 'checkpoint_rgamma_domainStudy.mat';

% Save final results
outputFile = 'rgamma_domainStudy_newSolver.mat';


%% ============================================================
% PHYSICAL PARAMETERS
% ============================================================

% Particle
dom.a = 25e-9;


% Material
mat.n_core   = 'Au';
mat.n_medium = 1.33;

mat.rho = [19300 1000];
mat.cp  = [129 4181];
mat.k   = [318 0.6];

mat.tau1 = 1.7e-12;


% Optical source
src.lda0 = 532e-9;
src.F    = 38.3;


% Boundary condition
bc.bound = 'Robin';
bc.T_ini = 0;


%% ============================================================
% DISCRETIZATION
% ============================================================

% Original spatial resolution
mesh.N_p = 200;


%% ============================================================
% SOLVER
% ============================================================

solver.mode = 'transient';

solver.showWaitbar = false;


% For this sweep I would initially keep the preconditioner OFF.
%
% A different F_d is generated for each pulse duration and domain,
% so rebuilding the approximate inverse every time may cost more
% than the sparse direct solve.
solver.preconditioner.enabled = false;


%% ============================================================
% PARAMETER SWEEP
% ============================================================

% Pulse durations
pulses = ...
    logspace(-2,8,44) * mat.tau1;

Npulses = numel(pulses);


% Gamma
%
% Keep the same 0 -> 1 range as the original calculations.
gammas = ...
    linspace(0,1,42);

Ngamma = numel(gammas);


% ITC cases
Gvalues = [Inf 2e8];

NG = numel(Gvalues);


% Domains
Ndomains = numel(domainFactors);


%% ============================================================
% OUTPUT ARRAYS
%
% Dimensions:
%
%   pulse x gamma x G x domain
%
% ============================================================

Rgamma = ...
    nan(Npulses,Ngamma,NG,Ndomains);


%% ============================================================
% MAIN CALCULATION
% ============================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' r_gamma recalculation with the new solver\n');
fprintf('====================================================\n\n');


calculationCounter = 0;

Ncalculations = ...
    Npulses * NG * Ndomains;


for idom = 1:Ndomains

    %% --------------------------------------------------------
    % Domain
    % ---------------------------------------------------------

    dom.R_sim = ...
        domainFactors(idom) * dom.a;


    % Keep the SAME radial step when the domain increases.
    %
    % Therefore Ne increases proportionally with R_sim.
    mesh.N_e = round( ...
        (dom.R_sim-dom.a)/dom.a * ...
        mesh.N_p);


    fprintf('\n');
    fprintf('---------------------------------------------\n');
    fprintf(' Domain: R_sim = %g a\n',domainFactors(idom));
    fprintf(' Nr approximately = %d\n',mesh.N_p + mesh.N_e);
    fprintf('---------------------------------------------\n');


    for j = 1:Npulses

        pulseTimer = tic;


        %% ----------------------------------------------------
        % Pulse duration
        % -----------------------------------------------------

        src.tau2 = ...
            pulses(j);


        % Same criterion used in the original calculations
        dom.t_max = ...
            max(10*src.tau2,50e-9);


        % Same temporal discretization strategy
        if dom.t_max == 50e-9

            mesh.N_t = 60000;

        else

            mesh.N_t = 10000;

        end


        %% ----------------------------------------------------
        % Both ITC cases
        % -----------------------------------------------------

        for iG = 1:NG

            mat.ITC = ...
                Gvalues(iG);


            calculationCounter = ...
                calculationCounter + 1;


            fprintf( ...
                ['[%3d/%3d] R = %2ga | pulse %2d/%2d | ' ...
                 'tau2 = %.3e s | G = %s | Nt = %d\n'], ...
                calculationCounter, ...
                Ncalculations, ...
                domainFactors(idom), ...
                j, ...
                Npulses, ...
                src.tau2, ...
                Glabel(mat.ITC), ...
                mesh.N_t);


            %% ================================================
            % THERMAL SOLUTION
            % ================================================

            solution = heatEqSolver( ...
                dom, ...
                mat, ...
                src, ...
                mesh, ...
                bc, ...
                solver);


            %% ================================================
            % MAXIMUM-TEMPERATURE ENVELOPE
            %
            % Do this in temporal chunks so that another huge
            % Nr x Nt temperature matrix is NOT created.
            % ================================================

            Tenv = temperatureEnvelope( ...
                solution, ...
                1000);


            %% ================================================
            % r_gamma
            % ================================================

            Rgamma(j,:,iG,idom) = ...
                calculateRgamma( ...
                    solution.nodes, ...
                    Tenv, ...
                    dom.a, ...
                    mesh.N_p, ...
                    gammas);


            %% Free the large transient solution ASAP

            clear solution Tenv

        end


        %% ----------------------------------------------------
        % Progress information
        % ---------------------------------------------------------

        elapsedPulse = toc(pulseTimer);


        fprintf( ...
            '    pulse completed in %.2f min\n', ...
            elapsedPulse/60);


        %% ----------------------------------------------------
        % CHECKPOINT
        % ---------------------------------------------------------

        save( ...
            checkpointFile, ...
            'Rgamma', ...
            'pulses', ...
            'gammas', ...
            'domainFactors', ...
            'Gvalues', ...
            'idom', ...
            'j', ...
            '-v7.3');

    end

end


%% ============================================================
% MAIN DOMAIN RESULTS
% ============================================================

R_inf = ...
    squeeze(Rgamma(:,:,1,mainDomainIndex));

R_fin = ...
    squeeze(Rgamma(:,:,2,mainDomainIndex));


%% ============================================================
% ITC RELATIVE DIFFERENCE
% ============================================================

epsilon_gamma = ...
    100 * ...
    abs(R_fin-R_inf) ./ ...
    R_inf;


%% ============================================================
% DOMAIN-SIZE EFFECT
% ============================================================

% Compare the first and last domains.
%
% By default:
%
%       40a vs 80a

R_inf_small = ...
    squeeze(Rgamma(:,:,1,1));

R_inf_large = ...
    squeeze(Rgamma(:,:,1,end));


R_fin_small = ...
    squeeze(Rgamma(:,:,2,1));

R_fin_large = ...
    squeeze(Rgamma(:,:,2,end));


domainError_inf = ...
    100 * ...
    abs(R_inf_small-R_inf_large) ./ ...
    R_inf_large;


domainError_fin = ...
    100 * ...
    abs(R_fin_small-R_fin_large) ./ ...
    R_fin_large;


%% ============================================================
% PRINT DOMAIN-CONVERGENCE INFORMATION
% ============================================================

% Only use the region shown in the paper:
%
% gamma >= 10 %

gammaMask = ...
    gammas >= 0.1;


validInf = ...
    domainError_inf(:,gammaMask);

validFin = ...
    domainError_fin(:,gammaMask);


validInf = ...
    validInf(isfinite(validInf));

validFin = ...
    validFin(isfinite(validFin));


fprintf('\n');
fprintf('====================================================\n');
fprintf(' DOMAIN-SIZE CHECK\n');
fprintf('====================================================\n');

fprintf( ...
    '%ga vs %ga\n', ...
    domainFactors(1), ...
    domainFactors(end));


fprintf( ...
    'G = Inf   | max error = %.3f %% | median = %.3f %%\n', ...
    max(validInf), ...
    median(validInf));


fprintf( ...
    'G = 2e8   | max error = %.3f %% | median = %.3f %%\n', ...
    max(validFin), ...
    median(validFin));


%% ============================================================
% SAVE FINAL DATA
% ============================================================

data = struct();

data.Rgamma = ...
    Rgamma;

data.pulses = ...
    pulses;

data.gammas = ...
    gammas;

data.domainFactors = ...
    domainFactors;

data.Gvalues = ...
    Gvalues;

data.epsilonGamma = ...
    epsilon_gamma;

data.domainErrorInf = ...
    domainError_inf;

data.domainErrorFinite = ...
    domainError_fin;

data.radius = ...
    dom.a;


save( ...
    outputFile, ...
    'data', ...
    '-v7.3');


%% ============================================================
% PLOTTING VARIABLES
% ============================================================

gammaPct = ...
    100*gammas;

logTau = ...
    log10(pulses);


contourLevels = ...
    [1.1 1.2 1.3 1.5 2 3 5 8];


%% ============================================================
% FIGURE 1
%
% Reproduction of the paper figure using the LARGE domain.
% ============================================================

figure( ...
    'Color','w', ...
    'Position',[100 100 1350 430]);


tiledlayout( ...
    1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ------------------------------------------------------------
% PANEL A — G = Inf
% ------------------------------------------------------------

ax1 = nexttile;


imagesc( ...
    gammaPct, ...
    logTau, ...
    R_inf);


set(ax1,'YDir','normal');

xlim([10 100]);
ylim([-14 -4]);

clim([1 10]);


xlabel('\gamma (%)');
ylabel('log_{10}(\tau_2 [s])');


if exist('inferno','file')

    colormap(ax1,inferno(512));

else

    colormap(ax1,turbo(512));

end


hold(ax1,'on');


[C1,h1] = contour( ...
    ax1, ...
    gammaPct, ...
    logTau, ...
    R_inf, ...
    contourLevels, ...
    '--', ...
    'LineColor','w', ...
    'LineWidth',1.4);


clabel( ...
    C1,h1, ...
    'Color','w', ...
    'FontSize',10);


cb1 = colorbar(ax1,'northoutside');

title(cb1,'\bar{r}_{\gamma} (a.u.)');


title( ...
    ax1, ...
    sprintf('G = \\infty, R_{sim} = %ga', ...
    domainFactors(mainDomainIndex)));


set(ax1,'FontSize',16,'Box','on');


%% ------------------------------------------------------------
% PANEL B — FINITE G
% ------------------------------------------------------------

ax2 = nexttile;


imagesc( ...
    gammaPct, ...
    logTau, ...
    R_fin);


set(ax2,'YDir','normal');

xlim([10 100]);
ylim([-14 -4]);

clim([1 10]);


xlabel('\gamma (%)');


if exist('inferno','file')

    colormap(ax2,inferno(512));

else

    colormap(ax2,turbo(512));

end


hold(ax2,'on');


[C2,h2] = contour( ...
    ax2, ...
    gammaPct, ...
    logTau, ...
    R_fin, ...
    contourLevels, ...
    '--', ...
    'LineColor','w', ...
    'LineWidth',1.4);


clabel( ...
    C2,h2, ...
    'Color','w', ...
    'FontSize',10);


cb2 = colorbar(ax2,'northoutside');

title(cb2,'\bar{r}_{\gamma} (a.u.)');


title( ...
    ax2, ...
    sprintf('G = 2\\times10^8, R_{sim} = %ga', ...
    domainFactors(mainDomainIndex)));


set(ax2,'FontSize',16,'Box','on');


%% ------------------------------------------------------------
% PANEL C — ITC EFFECT
% ------------------------------------------------------------

ax3 = nexttile;


imagesc( ...
    gammaPct, ...
    logTau, ...
    epsilon_gamma);


set(ax3,'YDir','normal');

xlim([10 100]);
ylim([-14 -4]);

clim([0 20]);


xlabel('\gamma (%)');


if exist('cividis','file')

    colormap(ax3,cividis(512));

else

    colormap(ax3,parula(512));

end


cb3 = colorbar(ax3,'northoutside');

title(cb3,'\epsilon_{\gamma} (%)');


title(ax3,'Finite-ITC effect');


set(ax3,'FontSize',16,'Box','on');


%% ============================================================
% FIGURE 2 — DOMAIN-SIZE SENSITIVITY
% ============================================================

figure( ...
    'Color','w', ...
    'Position',[150 150 950 430]);


tiledlayout( ...
    1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ------------------------------------------------------------
% G = Inf
% ------------------------------------------------------------

ax4 = nexttile;


imagesc( ...
    gammaPct, ...
    logTau, ...
    domainError_inf);


set(ax4,'YDir','normal');

xlim([10 100]);
ylim([-14 -4]);


xlabel('\gamma (%)');
ylabel('log_{10}(\tau_2 [s])');


title( ...
    sprintf( ...
    'Domain effect: G = \\infty (%ga vs %ga)', ...
    domainFactors(1), ...
    domainFactors(end)));


if exist('cividis','file')

    colormap(ax4,cividis(512));

else

    colormap(ax4,parula(512));

end


colorbar(ax4);

set(ax4,'FontSize',16,'Box','on');


%% ------------------------------------------------------------
% Finite G
% ------------------------------------------------------------

ax5 = nexttile;


imagesc( ...
    gammaPct, ...
    logTau, ...
    domainError_fin);


set(ax5,'YDir','normal');

xlim([10 100]);
ylim([-14 -4]);


xlabel('\gamma (%)');


title( ...
    sprintf( ...
    'Domain effect: G = 2\\times10^8 (%ga vs %ga)', ...
    domainFactors(1), ...
    domainFactors(end)));


if exist('cividis','file')

    colormap(ax5,cividis(512));

else

    colormap(ax5,parula(512));

end


colorbar(ax5);

set(ax5,'FontSize',16,'Box','on');


%% ============================================================
% TOTAL TIME
% ============================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' TOTAL TIME: %.2f min\n',toc(totalTimer)/60);
fprintf('====================================================\n');


%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function Tenv = temperatureEnvelope(solution,chunkSize)

% Calculates:
%
%       Tenv(r) = max_t T(r,t)
%
% without constructing another full Nr x Nt temperature matrix.

sigma = ...
    solution.spectralWeights(1);


Nr = ...
    size(solution.opticalBasis,1);

Nt = ...
    size(solution.opticalBasis,2);


Tenv = ...
    -inf(Nr,1);


for i1 = 1:chunkSize:Nt

    i2 = ...
        min(i1+chunkSize-1,Nt);


    ind = ...
        i1:i2;


    Tchunk = ...
        sigma * solution.opticalBasis(:,ind) + ...
        solution.boundaryBasis(:,ind);


    Tenv = ...
        max( ...
            Tenv, ...
            max(Tchunk,[],2));

end

end


function rgamma = calculateRgamma( ...
    nodes,Tenv,a,Np,gammas)

% Calculates the dimensionless confinement radius:
%
%       rbar_gamma = r_gamma/a
%
% in the surrounding medium.
%
% The reference temperature is the maximum-temperature envelope
% evaluated at the first external-medium node.


% First node immediately outside the particle.
%
% This follows the indexing used in the original calculations:
%
%       T(Np+1,:) = immediate external medium.

firstMedium = ...
    Np + 1;


r = ...
    nodes(firstMedium:end);


T = ...
    Tenv(firstMedium:end);


% Reference temperature at a+
Tref = ...
    T(1);


% Normalize envelope
Tnorm = ...
    T/Tref;


rgamma = ...
    nan(size(gammas));


for k = 1:numel(gammas)

    gamma = ...
        gammas(k);


    % gamma = 0 generally lies outside a finite simulation domain
    if gamma <= 0
        continue;
    end


    % Locate first crossing
    idx = ...
        find(Tnorm <= gamma,1,'first');


    if isempty(idx)

        % The requested confinement radius lies outside R_sim.
        rgamma(k) = NaN;

        continue;

    end


    if idx == 1

        rgamma(k) = ...
            r(1)/a;

        continue;

    end


    % Linear interpolation between neighboring radial nodes
    r1 = r(idx-1);
    r2 = r(idx);

    T1 = Tnorm(idx-1);
    T2 = Tnorm(idx);


    if T2 == T1

        rCross = r2;

    else

        rCross = ...
            r1 + ...
            (gamma-T1) * ...
            (r2-r1)/(T2-T1);

    end


    rgamma(k) = ...
        rCross/a;

end

end


function txt = Glabel(G)

if isinf(G)

    txt = 'Inf';

else

    txt = sprintf('%.2e',G);

end

end