function solution = heatEqSolver(dom,mat,src,mesh,bc,solver)
% ============================================================
% HEAT EQUATION SOLVER
%
% solver.mode = 'cw'
% solver.mode = 'transient'
%
% Optional inverse-based preconditioning:
%
% solver.preconditioner.enabled    = true / false
% solver.preconditioner.targetCond = 1e3
% solver.preconditioner.iterMax    = 10
% solver.preconditioner.timeMax    = 360
%
% Common spectral reconstruction:
%   T(lambda) =
%       spectralWeights(lambda) * opticalBasis
%       + boundaryBasis
% CW:
%   opticalBasis  -> N_r x 1
%   boundaryBasis -> N_r x 1
%
% TRANSIENT:
%   opticalBasis  -> N_r x N_t
%   boundaryBasis -> N_r x N_t
% ===========================================================

total_ref = tic;

% ============================================================
% CHECK SOLVER MODE
% ============================================================
if ~isfield(solver,'mode')
    error('solver.mode must be defined as ''cw'' or ''transient''.');
end
mode = lower(solver.mode);
if ~ismember(mode,{'cw','transient'})
    error( ...
        'Unknown solver mode "%s". Use ''cw'' or ''transient''.', ...
        solver.mode);
end

% ============================================================
% PRECONDITIONER OPTIONS
% ============================================================

% Default: no preconditioning
usePreconditioner = false;

% Default preconditioner parameters
target_cond = 1e3;
iter_max    = 10;
time_max    = 360;

% Read user-defined options
if isfield(solver,'preconditioner')
    if isfield(solver.preconditioner,'enabled')
        usePreconditioner = solver.preconditioner.enabled;
    end
    if isfield(solver.preconditioner,'targetCond')
        target_cond = solver.preconditioner.targetCond;
    end
    if isfield(solver.preconditioner,'iterMax')
        iter_max = solver.preconditioner.iterMax;
    end
    if isfield(solver.preconditioner,'timeMax')
        time_max = solver.preconditioner.timeMax;
    end
end

% Initialize preconditioner information
info_pre = [];

% ============================================================
% COMMON PARAMETERS
% ============================================================

% Simulation domain
a     = dom.a;
R_sim = dom.R_sim;

% Optical properties
n_core   = mat.n_core;
n_medium = mat.n_medium;

% Thermal properties
rho = mat.rho;
cp  = mat.cp;
k   = mat.k;
ITC = mat.ITC;

% Source
lambda     = src.lda0;
irradiance = src.F;

% Boundary conditions
boundary = bc.bound;
T_ini    = bc.T_ini;

% ============================================================
% OPTICAL CALCULATION
% ============================================================
optical_ref = tic;
nLambda = numel(lambda);
sigma = zeros(1,nLambda);
for i = 1:nLambda
    [~,sigma(i)] = mie_absorption(a, n_core, n_medium, lambda(i));

end
optical_time = toc(optical_ref);

% ============================================================
% THERMAL SOLVER
% ============================================================
thermal_ref = tic;
switch mode
    case 'cw' % CONTINUOUS WAVE  
        N = [mesh.N_p, mesh.N_e ]; % Spatial discretization 
      
        % Build stationary thermal system
        [A,Q_unit,nodes] = kernelAssemblerCW(a, R_sim, N, k, T_ini, ...
            boundary, ITC);
        
        % Separate optical and boundary contributions
        Q_unit = Q_unit(:);

        % Unit volumetric heating contribution
        b_optical = Q_unit;
        b_optical(end) = 0;

        % Boundary contribution
        b_boundary = zeros(size(Q_unit));
        b_boundary(end) = Q_unit(end);

        
        % Include illumination amplitude
        b_optical = irradiance * b_optical;

        % Combined RHS
        RHS = [b_optical, b_boundary ];

        if usePreconditioner
            [M,info_pre] = preconditioner(A, target_cond, iter_max, time_max);
            if info_pre.time_limit_reached
                warning( ...
                    ['Preconditioner construction was interrupted. ' ...
                     'Falling back to the unpreconditioned CW system.']);
                dA = decomposition(A);
                X = dA\RHS;
            else                
                A_pre = M*A;
                % Factorize preconditioned matrix once
                dA = decomposition(A_pre);
                
                % Solve both RHS simultaneously
                X = dA\(M*RHS);
            end
        else          
            % STANDARD DIRECT SOLUTION
            dA = decomposition(A);
            X = dA\RHS;
        end
        
      
        % Thermal bases
        T_optical = X(:,1);
        T_boundary = X(:,2);

        % ----------------------------------------------------
        % Output
        % -----------------------------------------------------
        solution.opticalBasis = T_optical;
        solution.boundaryBasis = T_boundary;

        % CW has no temporal domain
        solution.timeDomain = [];

    % ========================================================
    % TRANSIENT
    % =========================================================
    case 'transient'
        if ~isfield(dom,'t_max')
            error('dom.t_max is required for transient simulations.');
        end
        if ~isfield(mesh,'N_t')
            error('mesh.N_t is required for transient simulations.');
        end
        if ~isfield(src,'tau2')
            error('src.tau2 is required for transient simulations.');
        end
        if ~isfield(mat,'tau1')
            error('mat.tau1 is required for transient simulations.');
        end
        t_max = dom.t_max;
        tau1 = mat.tau1;
        tau2 = src.tau2;
        N = [mesh.N_p, mesh.N_e, mesh.N_t ];

        % Temporal illumination profile
        Qfunc = illumination(irradiance, tau1, tau2);

        % ----------------------------------------------------
        % Build Crank-Nicolson system
        % -----------------------------------------------------
        [F_d,F_a, ...
         Q_vals,Q_spatial,Q_boundary, ...
         nodes,t] = kernelAssembler( ...
            a, ...
            R_sim, ...
            t_max, ...
            N, ...
            rho, ...
            cp, ...
            k, ...
            Qfunc, ...
            T_ini, ...
            boundary, ...
            ITC);  
        
        % PREPARE LINEAR SOLVER
        if usePreconditioner
            [M,info_pre] = preconditioner(F_d, target_cond, iter_max, ...
                time_max);          
            if info_pre.time_limit_reached
                warning( ...
                    ['Preconditioner construction was interrupted. ' ...
                     'Falling back to the unpreconditioned transient system.']);
                M = [];
                dFd = decomposition(F_d);
                preconditionerActive = false;
            else
                Fd_pre = M*F_d; % Factorize only once
                dFd = decomposition(Fd_pre);
                preconditionerActive = true;
            end
        else
            % STANDARD FACTORIZATION
            M = [];
            dFd = decomposition(F_d);
            preconditionerActive = false;
        end
     
        % Allocate thermal bases
        Nr = length(nodes);
        Nt = length(t);
        T_optical = zeros(Nr,Nt);
        T_boundary = zeros(Nr,Nt);

        
        % Initial condition
        T_boundary(:,1) = T_ini;
                
        % Progress bar
        showWaitbar = true;
        if isfield(solver,'showWaitbar')
            showWaitbar = solver.showWaitbar;
        end
        if showWaitbar
            h = waitbar(0, 'Solving transient heat transfer...');
            updateEvery = max(1, floor((Nt-1)/100));
        end

        % ----------------------------------------------------
        % Crank-Nicolson time evolution
        %
        % Optical basis:
        %
        % F_d*T_o^(n+1) =
        %
        %       F_a*T_o^n
        %       + Q_spatial*Q_vals(n)
        %
        %
        % Boundary basis:
        %
        % F_d*T_b^(n+1) =
        %
        %       F_a*T_b^n
        %       + Q_boundary
        %
        % -----------------------------------------------------
        for n = 1:Nt-1
            % Optical contribution
            rhs_optical = F_a*T_optical(:,n) + Q_spatial*Q_vals(n);

            % Boundary contribution
            rhs_boundary = F_a*T_boundary(:,n) + Q_boundary;

            % Combined RHS
            RHS = [rhs_optical, rhs_boundary ];

            % Solve both RHS simultaneously
            if preconditionerActive
                X = dFd\(M*RHS);
            else
                X = dFd\RHS;
            end

            % Store next timestep
            T_optical(:,n+1) = X(:,1);
            T_boundary(:,n+1) = X(:,2);

            % Progress
            if showWaitbar && (mod(n,updateEvery)==0 || n==Nt-1)
                waitbar(n/(Nt-1), h);
            end
        end

        % ----------------------------------------------------
        % Close progress bar
        % -----------------------------------------------------
        if showWaitbar
            close(h);
        end

        % ----------------------------------------------------
        % Output
        % -----------------------------------------------------
        solution.opticalBasis = T_optical;
        solution.boundaryBasis = T_boundary;
        solution.timeDomain = t;
end

% ============================================================
% COMMON OUTPUT
% ============================================================
thermal_time = toc(thermal_ref);
solution.spectralWeights = sigma;
solution.lambda = lambda;
solution.nodes = nodes;
solution.mode = mode;

% ============================================================
% SOURCE INFORMATION
% ============================================================
solution.source.irradiance = irradiance;
if strcmp(mode,'transient')
    solution.source.tau1 = tau1;
    solution.source.tau2 = tau2;
end

% ============================================================
% COMPUTATIONAL INFORMATION
% ============================================================
solution.info.opticalTime = optical_time;
solution.info.thermalTime = thermal_time;
solution.info.totalTime = toc(total_ref);

% ============================================================
% PRECONDITIONER INFORMATION
% ============================================================
solution.info.preconditionerRequested = usePreconditioner;
if usePreconditioner && ~isempty(info_pre)
    solution.info.preconditioner = info_pre;
else
    solution.info.preconditioner = [];
end

% ============================================================
% RECONSTRUCTION
%
% For wavelength index iLambda:
%
% T_lambda =
%
%   solution.spectralWeights(iLambda)
%       * solution.opticalBasis
%
%   + solution.boundaryBasis;
%
% This works identically for CW and transient.
%
% ============================================================

end