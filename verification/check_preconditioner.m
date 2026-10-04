function result = check_preconditioner(showFigure)

% CHECK_PRECONDITIONER
%
% Verifies that the truncated-inverse preconditioner:
%
%   1) Reduces the condition number of the system matrix.
%   2) Preserves the numerical solution.
%   3) Reaches the requested conditioning target for a
%      controlled sparse test problem.
%
% OUTPUT
%   result.name
%   result.pass
%   result.error
%   result.threshold
%   result.conditionOriginal
%   result.conditionPreconditioned
%   result.improvementFactor
%   result.residualDirect
%   result.residualPreconditioned
%   result.preconditionerInfo
% ============================================================
if nargin < 1
    showFigure = true;
end

% ============================================================
% REFERENCE SYSTEM
% ============================================================

% Matrix size.
n = 200;

% ------------------------------------------------------------
% Build a diffusion-like sparse matrix. T is well behaved and 
% trongly diagonally dominant.
% ------------------------------------------------------------
e = ones(n,1);
T = spdiags([-0.20*e, 1.00*e, -0.20*e], -1:1, n,n);

% ------------------------------------------------------------
% Introduce a strong row scaling to deliberately generate an
% ill-conditioned system.
%
% The resulting matrix remains sparse and representative of the
% type of scaling problems that can appear after discretization.
% ------------------------------------------------------------

scale = logspace(0,8,n)';

D = spdiags(scale,0,n,n);

A = D*T;


% ============================================================
% RIGHT-HAND SIDE
% ============================================================

% Deterministic reference solution.
x_exact = sin(linspace(0,2*pi,n))' + 1;


% Generate RHS from the known solution.
b = A*x_exact;

% ============================================================
% DIRECT SOLUTION
% ============================================================
x_direct = A\b;

% ============================================================
% PRECONDITIONER SETTINGS
% ============================================================
target_cond = 1e3;
iter_max = 10;
time_max = 60;

% ============================================================
% BUILD PRECONDITIONER
% ============================================================
[M,info] = preconditioner(A, target_cond, iter_max, time_max);

% ============================================================
% PRECONDITIONED SYSTEM
% ============================================================
A_pre = M*A;
b_pre = M*b;
x_pre = A_pre\b_pre;

% ============================================================
% CONDITIONING
% ============================================================
conditionOriginal = condest(A);
conditionPreconditioned = condest(A_pre);
improvementFactor = conditionOriginal/conditionPreconditioned;

% ============================================================
% SOLUTION ERROR
% ============================================================
relativeSolutionError = norm(x_pre-x_direct)/max(norm(x_direct),eps);

% ============================================================
% RESIDUALS
% ============================================================
residualDirect = norm(A*x_direct-b)/max(norm(b),eps);
residualPreconditioned = norm(A*x_pre-b)/max(norm(b),eps);

% ============================================================
% PASS CRITERIA
% ============================================================
% Solution preservation tolerance.
threshold = 1e-8;
conditioningImproved = conditionPreconditioned<conditionOriginal;
solutionPreserved = relativeSolutionError<=threshold;
residualAcceptable = residualPreconditioned<=threshold;
targetReached = conditionPreconditioned<=target_cond;
timeOK = ~info.time_limit_reached;
pass = ...
    conditioningImproved && ...
    solutionPreserved && ...
    residualAcceptable && ...
    targetReached && ...
    timeOK;

% ============================================================
% RESULT STRUCTURE
% ============================================================
result.name = 'Inverse preconditioner';
result.pass = pass;
result.error = relativeSolutionError;
result.threshold = threshold;
result.conditionOriginal = conditionOriginal;
result.conditionPreconditioned = conditionPreconditioned;
result.targetCondition = target_cond;
result.improvementFactor = improvementFactor;
result.residualDirect = residualDirect;
result.residualPreconditioned = residualPreconditioned;
result.conditioningImproved = conditioningImproved;
result.solutionPreserved = solutionPreserved;
result.targetReached = targetReached;
result.preconditionerInfo = info;

% ============================================================
% PRINT RESULTS
% ============================================================

fprintf('\n');
fprintf('=========================================\n');
fprintf(' PRECONDITIONER VERIFICATION\n');
fprintf('=========================================\n');
fprintf('cond(A) = %.3e\n', conditionOriginal);
fprintf( 'cond(M*A) = %.3e\n', conditionPreconditioned);
fprintf( 'Improvement = %.3e x\n', improvementFactor);
fprintf( 'Solution error = %.3e\n', relativeSolutionError);
fprintf( 'Direct residual = %.3e\n', residualDirect);
fprintf( 'Precond. residual = %.3e\n', residualPreconditioned);
if pass
    fprintf('\nRESULT: PASS\n');
else
    fprintf('\nRESULT: FAIL\n');
end
fprintf('=========================================\n\n');

% ============================================================
% OPTIONAL FIGURE
% ============================================================
if showFigure
    figure('Color','w','Name','Preconditioner verification');
    values = [conditionOriginal, conditionPreconditioned];
    bar(values);
    set(gca, 'YScale','log', 'XTick',[1 2], 'XTickLabel', {'A','M A'}, ...
        'FontSize',14);
    ylabel('Condition number');
    title(sprintf(['Preconditioner verification\n' ...
        'solution error = %.2e'], relativeSolutionError));
    grid on;
end
end