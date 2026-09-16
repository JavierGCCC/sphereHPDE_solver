function [M,info]=preconditioner(U,target_cond,iter_max,time_max)
% Validación mínima de U
if nargin < 1
    error('Matrix U is required.');
end

% Chequear triangularidad
is_upper_tri = istriu(U);

% En caso de que no se pase target_cond etc., se pasa al modo diagonal
use_diag = (~is_upper_tri || nargin < 4);

singularity = condest(U);
if singularity > 1e13
    warning('⚠ Matrix U appears to be near-singular (condest = %.2e)', singularity);
end

n = size(U,1);
k_max = 5;
estimated_nnz = n * k_max;

if n > 2e4
    fprintf('\n  Estimated preconditioner size: %d elements (very large).\n', estimated_nnz);
    fprintf('    This may take significant time for %d columns.\n', n);
    fprintf('\nHow do you want to proceed?\n');
    fprintf('  1) Continue with full truncated inverse preconditioner\n');
    fprintf('  2) Use simple diagonal preconditioner based on row norms\n');
    
    choice = input('Select an option (1 or 2): ');
elseif use_diag
    choice = 2;
else
    choice = 1;
end

fprintf('\n==============================\n');
fprintf('   PRECONDITIONER STARTED\n');
fprintf('==============================\n\n');

switch choice
    case 1
        fprintf('🚀 Proceeding with full preconditioner...\n');
        if n > 1000
            check_every = floor(n / 100);
        else
            check_every = 10;
        end

        M = spalloc(n, n, k_max * n);
        t_ref = tic;
        h = waitbar(0, '⏳ Building truncated inverse preconditioner...');

        for t = 1:iter_max
            for j = 1:n
                if mod(j, check_every) == 0 && toc(t_ref) > time_max
                    fprintf('\n⏱ Maximum time exceeded (%.1f s).\n', time_max);
                    fprintf('⛔ Execution stopped early due to time limit.\n');
                    close(h);
                    
                    info.iter = t;
                    info.conditioning = condest(M*U);
                    info.k_max = k_max;
                    info.check_every = check_every;
                    info.elapsed_time = toc(t_ref);
                    info.time_limit_reached = true;
                    info.success = false;
                    info.method = "truncated-inverse";
                    info.rcond_U = singularity;
                    return;
                end

                % Vector columna j de la identidad
                e_j = sparse(j, 1, 1, n, 1);
                x = U \ e_j;

                [~, idx_sorted] = sort(abs(x), 'descend');
                keep_idx = idx_sorted(1:min(k_max, nnz(x)));
                M(:,j) = sparse(keep_idx, 1, x(keep_idx), n, 1);

                if mod(j, check_every) == 0 || j == n
                    waitbar(j/n, h);
                end
            end

            f_cond = condest(M * U);
            if f_cond < target_cond
                break;
            end
            k_max = k_max + 5 * t;
            M = spalloc(n, n, k_max * n);
        end

        close(h);
        info.iter = t;
        info.conditioning = f_cond;
        info.k_max = k_max;
        info.check_every = check_every;
        info.elapsed_time = toc(t_ref);
        info.time_limit_reached = false;
        info.success = f_cond < target_cond;
        info.method = "truncated-inverse";
        info.rcond_U = singularity;

    case 2
        fprintf('✅ Using simple diagonal preconditioner based on row norms.\n');
        t_ref = tic;
        row_norms = sum(abs(U), 2);
        M = spdiags(1 ./ row_norms, 0, n, n);
        
        info.iter = 0;
        info.conditioning = condest(M * U);
        info.k_max = k_max;
        info.check_every = NaN;
        info.elapsed_time = toc(t_ref);
        info.time_limit_reached = false;
        info.success = true;
        info.method = "diagonal";
        info.rcond_U = singularity;

    otherwise
        error('⛔ Invalid option selected. Execution aborted.');
end

fprintf('\n✅ Preconditioner finished in %.2f s with cond(M·U) ≈ %.2e\n', ...
    info.elapsed_time, info.conditioning);
end
