function [ns, sigma_abs] = mie_absorption(a, nCore, nMedium, lda0)
% MIE_ABSORPTION
% Computes the absorption cross-section of a spherical particle
% using Mie theory.
%
% INPUTS:
%   a        - particle radius [m]
%   nCore    - complex refractive index (numeric) OR string with material
%              name ('Au', 'Ag', etc.). See /mie/permittivity folder.
%   nMedium  - refractive index of surrounding medium
%   lda0     - wavelength [m]
%
% OUTPUT:
%   sigma_abs - absorption cross-section [m^2]
%
% REQUIREMENTS:
%   calcmie.m and its dependencies inside /mie/
%   optical data files inside /mie/dat_mat/ (format: [λ, Re(ε), Im(ε)])
%
% OBSERVATIONS:
%   -Interpolation is applied for input lda0. Check the number of
%   interpolated elements in the definition of dom_I that suits your
%   specifics.
%
% EXAMPLES OF CALL:
%   sigma = mie_absorption(20e-9, 0.2+3.5i, 1.33, 532e-9);
%   sigma = mie_absorption(20e-9, 'Au',    1.33, 532e-9);

    % ========== SECTION I: Ensure Mie paths ==========
    baseDir = fileparts(mfilename('fullpath'));
    mieDir  = fullfile(baseDir, 'mie');
    addpath(mieDir);
    addpath(fullfile(mieDir, 'util'));
    addpath(fullfile(mieDir, 'expcoeff'));
    addpath(fullfile(mieDir, 'permittivity'));
    addpath(fullfile(mieDir, 'bessel'));

    % ========== SECTION II: Interpret nCore ==========
    if isnumeric(nCore)
        % Case 1: user gave directly nCore (complex number)
        ns = nCore;

    elseif ischar(nCore) || isstring(nCore)
        % Case 2: user gave a material name → load optical data
        nCore = char(nCore);
        dataFileMat = fullfile(mieDir, 'permittivity', ['data_' nCore '.mat']);
        dataFileTxt = fullfile(mieDir, 'permittivity', ['data_' nCore '.txt']);

        if exist(dataFileMat, 'file') == 2
            data = importdata(dataFileMat);
        elseif exist(dataFileTxt, 'file') == 2
            data = importdata(dataFileTxt);
        else
            error('Optical data for %s not found in mie/dat_mat/', nCore);
        end

        % Expect columns: [λ, Re(ε), Im(ε)]
        dom_I = linspace(min(data(:,1)), max(data(:,1)), 4000);
        r_eps = interp1(data(:,1), data(:,2), dom_I, 'spline');
        i_eps = interp1(data(:,1), data(:,3), dom_I, 'spline');
        [~, idx_lda0] = min(abs(dom_I - lda0));
        lambda = dom_I(idx_lda0);
        epsMat = r_eps(idx_lda0) + 1i*i_eps(idx_lda0);
        ns     = sqrt(epsMat);

        fprintf('[Mie] λ = %.1f nm (interpolated)\n', lambda*1e9);

    else
        error('nCore must be complex or string (material name).');
    end

    % ========== SECTION III: Absorption cross-section ==========
    nang = 2000; % angular resolution
    conv = 1;    % convergence factor

    [~, C, ~] = calcmie(a, ns, nMedium, lda0, nang, ...
                        'ConvergenceFactor', conv);

    sigma_abs = C.abs;

    fprintf('[Mie] σ_abs = %.3e m²\n', sigma_abs);

end
