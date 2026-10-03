function initialize()

% INITIALIZE
% Configure sphereHPDE_solver for the current MATLAB session.
%
% Adds the required sphereHPDE_solver directories and external
% dependencies to the MATLAB search path.


% ============================================================
% REPOSITORY ROOT
% ============================================================

rootDir = fileparts(mfilename('fullpath'));


% ============================================================
% sphereHPDE_solver
% ============================================================

addpath(fullfile(rootDir,'thermal'));
addpath(fullfile(rootDir,'optical'));
addpath(fullfile(rootDir,'tools'));
addpath(fullfile(rootDir,'verification'));


% ============================================================
% EXTERNAL DEPENDENCIES
% ============================================================

% MatScat
matscatDir = fullfile( ...
    rootDir, ...
    'external', ...
    'matscat-1.4.0');

addpath(matscatDir);
addpath(fullfile(matscatDir,'util'));
addpath(fullfile(matscatDir,'expcoeff'));
addpath(fullfile(matscatDir,'bessel'));


% ============================================================
% INITIALIZATION MESSAGE
% ============================================================

fprintf('[sphereHPDE_solver] Initialization completed.\n');

end