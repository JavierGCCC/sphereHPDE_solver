function [ns, sigma_abs] = mie_absorption(a, nCore, nMedium, lda0)
% MIE_ABSORPTION
% Computes the absorption cross-section of a spherical particle
% using Mie theory.
%
% INPUTS:
%   a        - particle radius [m]
%   nCore    - complex refractive index (numeric) OR string with material
%              name 'Au'. See /optical/optical_properties/.
%   nMedium  - refractive index of surrounding medium
%   lda0     - wavelength [m]
%
% OUTPUTS:
%   ns        - complex refractive index at lda0
%   sigma_abs - absorption cross-section [m^2]
%
% REQUIREMENTS:
%   calcmie.m and its dependencies must be available in the MATLAB path.
%   Run initialize.m before using this function.
%   Optical data files must be stored inside
%   /optical/optical_properties/
%   with format: [lambda, Re(epsilon), Im(epsilon)]
%
% OBSERVATIONS:
%   - Optical properties are interpolated at the input wavelength lda0.
%
% EXAMPLES OF CALL:
%   [ns, sigma] = mie_absorption(20e-9, 0.2+3.5i, 1.33, 532e-9);
%   [ns, sigma] = mie_absorption(20e-9, 'Au',    1.33, 532e-9);


% ========== Optical data path ==========
baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir, 'optical_properties');


% ========== SECTION I: Interpret nCore ==========
if isnumeric(nCore)

    % Case 1: user gave directly nCore (complex number)
    ns = nCore;

elseif ischar(nCore) || isstring(nCore)

    % Case 2: user gave a material name -> load optical data
    nCore = char(nCore);

    dataFileMat = fullfile(dataDir, ['data_' nCore '.mat']);
    dataFileTxt = fullfile(dataDir, ['data_' nCore '.txt']);

    if exist(dataFileMat, 'file') == 2

        data = importdata(dataFileMat);

    elseif exist(dataFileTxt, 'file') == 2

        data = importdata(dataFileTxt);

    else

        error(['Optical data for %s not found in ' ...
            'optical/optical_properties/'], nCore);

    end

    % Expect columns: [lambda, Re(epsilon), Im(epsilon)]
    r_eps = interp1(data(:,1), data(:,2), lda0, 'spline');
    i_eps = interp1(data(:,1), data(:,3), lda0, 'spline');

    epsMat = r_eps + 1i*i_eps;
    ns     = sqrt(epsMat);

    fprintf('[Mie] lambda = %.1f nm (interpolated)\n', lda0*1e9);

else

    error('nCore must be complex or string (material name).');

end


% ========== SECTION II: Absorption cross-section ==========
nang = 2000; % angular resolution
conv = 1;    % convergence factor

[~, C, ~] = calcmie(a, ns, nMedium, lda0, nang, ...
    'ConvergenceFactor', conv);

sigma_abs = C.abs;

fprintf('[Mie] sigma_abs = %.3e m^2\n', sigma_abs);

end