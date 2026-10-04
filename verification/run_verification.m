function results = run_verification()

% RUN_VERIFICATION
% Runs the complete sphereHPDE_solver verification suite.

% ============================================================
% HEADER
% ============================================================
fprintf('\n');
fprintf('=============================================\n');
fprintf(' sphereHPDE_solver verification suite\n');
fprintf('=============================================\n\n');

% ============================================================
% RUN CHECKS
% ============================================================
results = cell(1,4);
results{1} = check_transient_to_cw(false);
results{2} = check_quasistatic_limit(false);
results{3} = check_cw_analytical_solution(false);
results{4} = check_linearity(false);

% ============================================================
% SUMMARY
% ============================================================
fprintf('\n');
fprintf('---------------------------------------------\n');
fprintf(' Verification summary\n');
fprintf('---------------------------------------------\n');
for i = 1:length(results)
    if results{i}.pass
        status = 'PASS';
    else
        status = 'FAIL';
    end
    fprintf('%-28s %s\n', results{i}.name, status);
end

% ============================================================
% GLOBAL RESULT
% ============================================================

nPassed = sum(cellfun(@(x) x.pass, results));
nTotal = length(results);
fprintf('---------------------------------------------\n');
fprintf('Result: %d/%d checks passed.\n', nPassed,nTotal);
if nPassed == nTotal
    fprintf('\n');
    fprintf('All verification checks passed successfully.\n');
else
    fprintf('\n');
    fprintf('WARNING: One or more verification checks failed.\n');
end
fprintf('=============================================\n\n');
end