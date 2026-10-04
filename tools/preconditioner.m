function [M,info] = preconditioner(A,target_cond,iter_max,time_max)

% PRECONDITIONER
% Builds a sparse truncated approximation of the inverse matrix: M ~=
% A^(-1),such that the preconditioned matrix: M*A reaches a target 
% condition number whenever possible.
%
% INPUTS
%   A             System matrix.
%   target_cond   Target condition number for M*A.
%   iter_max      Maximum number of truncation iterations.
%   time_max      Maximum allowed construction time [s].
%
% OUTPUTS
%   M             Sparse approximate inverse preconditioner.
%   info          Structure containing preconditioner diagnostics.

% ============================================================
% INPUT CHECK
% ============================================================

if nargin < 1
    error('System matrix A is required.'); %Error statement.
end

if nargin < 2 || isempty(target_cond)
    target_cond = 1e3; %Default target_cond.
end

if nargin < 3 || isempty(iter_max)
    iter_max = 10; %Default maximum number of iterations.
end

if nargin < 4 || isempty(time_max)
    time_max = Inf; %Defoult: there is no time restriction.
end

%Check matrix is square.
[n,m] = size(A);
if n ~= m 
    error('System matrix A must be square.');
end

%Check condition number of A.
condition_original = condest(A);
if condition_original > 1e13
    warning( ...
        'Matrix A appears to be near-singular (condest = %.2e).', ...
        condition_original); %Warning: A is ill-conditioned.
end

% ============================================================
% CONDITIONING INITIALIZATION
% ============================================================

% Initial number of coefficients retained per inverse column
k_max = 5;

% Check sampling.
if n > 1000
    check_every = max(1,floor(n/100));
else
    check_every = 10;
end

% Print information.
fprintf('\n');
fprintf('=========================================\n');
fprintf('   INVERSE PRECONDITIONER STARTED\n');
fprintf('=========================================\n\n');
fprintf('Original condition number: %.3e\n', ...
    condition_original);
fprintf('Target condition number:   %.3e\n\n', ...
    target_cond);

% ============================================================
% UNIQUE FACTORIZATION OF A
% ============================================================

dA = decomposition(A); %Decomposition of A.

% ============================================================
% BUILD APPROXIMATE INVERSE
% ============================================================

t_ref = tic; % Reference simulation time. 
h = waitbar(0,'Building truncated inverse preconditioner...');

% Default information in case the process stops early.
f_cond = condition_original;

% Simulation time limit nor reached by default. 
time_limit_reached = false;

for t = 1:iter_max %Loop.

    fprintf( ...
        'Iteration %d | retained coefficients per column: %d\n', ...
        t, k_max); %Print information.
    M = spalloc(n, n, k_max*n); % Allocate sparse approximate inverse

    % --------------------------------------------------------
    % Construct inverse column by column
    % ---------------------------------------------------------
    for j = 1:n
        % Check maximum execution time
        if mod(j,check_every)==0 && ...
                toc(t_ref) > time_max
            %Print information.
            fprintf('\n');
            fprintf( ...
                'Maximum preconditioning time exceeded (%.1f s).\n', ...
                time_max);
            time_limit_reached = true; %Change default.
            break; % Stops loop.

        end

        % Solve inverse column.
        e_j = sparse(j, 1, 1, n, 1); % Unit vector e_j.
        x = dA\e_j; % Column j of the exact inverse.

        % Sort coefficients by magnitude.
        [~,idx_sorted] = sort(abs(x), 'descend');

        % Keep only the largest coefficients.
        keep_idx = idx_sorted(1:min(k_max, nnz(x)));

        % Store sparse truncated inverse column.
        M(:,j) = sparse(keep_idx, 1, x(keep_idx), n, 1);

        % Update waitbar.
        if mod(j,check_every)==0 || j==n
            progress = ((t-1)*n+j)/(iter_max*n);
            waitbar(progress, h);
        end
    end

    % Stop if time limit was reached.
    if time_limit_reached
        break;
    end

    % Evaluate preconditioned matrix.
    A_pre = M*A;
    f_cond = condest(A_pre);
    fprintf('cond(M*A) = %.3e\n\n', f_cond);
    
    % Target reached.
    if f_cond <= target_cond
        break;
    end
    
    % Increase inverse density.
    k_max = k_max + 5*t;
end

% Close progress window.
if isvalid(h)
    close(h);
end

% ============================================================
% FINAL DIAGNOSTICS
% ============================================================

elapsed_time = toc(t_ref);

% If execution was interrupted during construction, evaluate the
% available approximate inverse when possible.
if time_limit_reached
    try
        f_cond = condest(M*A);
    catch
        f_cond = Inf;
    end
end

% ============================================================
% OUTPUT INFORMATION
% ============================================================

info.iter = t;
info.conditionOriginal = condition_original;
info.conditioning = f_cond;
info.targetCondition = target_cond;
info.k_max = k_max;
info.check_every = check_every;
info.elapsed_time = elapsed_time;
info.time_limit_reached = time_limit_reached;
info.success = (~time_limit_reached) && (f_cond <= target_cond);
info.method = "truncated-inverse";

% Sparsity diagnostics.
info.nnzOriginal = nnz(A);
info.nnzPreconditioner = nnz(M);
info.nnzPreconditioned = nnz(M*A);
info.fillRatio = info.nnzPreconditioned/max(info.nnzOriginal,1);

% ============================================================
% FINAL MESSAGE
% ============================================================
fprintf('\n');
fprintf('=========================================\n');
fprintf('Preconditioner finished in %.2f s\n', elapsed_time);
fprintf('cond(A)   = %.3e\n', condition_original);
fprintf('cond(M*A) = %.3e\n', f_cond);
fprintf( 'Fill ratio = %.2f\n', info.fillRatio);

if info.success
    fprintf('Target conditioning reached successfully.\n');
elseif time_limit_reached
    fprintf('Stopped before reaching target: time limit reached.\n');
else
    fprintf('Target conditioning was not reached.\n');
end
fprintf('=========================================\n\n');
end