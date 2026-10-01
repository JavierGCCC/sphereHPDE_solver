function initialize()
% INITIALIZE  Configure sphereHPDE_solver for the current MATLAB session.
%
% Adds the required sphereHPDE_solver and MatScat directories
% to the MATLAB search path.

    % Repository root
    rootDir = fileparts(mfilename('fullpath'));

    % sphereHPDE_solver
    addpath(fullfile(rootDir,'thermal'));
    addpath(fullfile(rootDir,'optical'));
    addpath(fullfile(rootDir,'tools'));

    % MatScat
    matscatDir = fullfile(rootDir,'external','matscat-1.4.0');

    addpath(matscatDir);
    addpath(fullfile(matscatDir,'util'));
    addpath(fullfile(matscatDir,'expcoeff'));
    addpath(fullfile(matscatDir,'bessel'));

    fprintf('[sphereHPDE_solver] Initialization completed.\n');

end