function [lambda, sigma_abs, sigma_abs_qe, deltaT, deltaT_qe]=check_QA(dom,mat,range_lam,I,mesh,bc, cond, doPlot)

%ABSORPTION
% ============================================================
% INPUT CHECK
% ============================================================

% nCore can be:
%   1) Material name: 'Au', 'Ag', ...
%   2) Scalar complex refractive index
%   3) Complex refractive-index spectrum n(lambda)

if ischar(mat.n_core) || (isstring(mat.n_core) && isscalar(mat.n_core))

    % Material database mode
    core_mode = 'material';

elseif isnumeric(mat.n_core)

    % Numerical refractive-index mode
    core_mode = 'numeric';

    if isscalar(mat.n_core)

        % Constant refractive index over the whole spectrum
        mat.n_core = repmat(mat.n_core, size(range_lam));

    elseif isvector(mat.n_core) && numel(mat.n_core) == numel(range_lam)

        % Ensure same orientation as lda0
        mat.n_core = reshape(mat.n_core, size(range_lam));

    else
        error(['If mat.n_core is numeric, it must be either a scalar ' ...
            'or a vector with the same number of elements as range_lam.']);
    end

else
    error('mat.n_core must be a material name or a numeric refractive index.');
end


% ============================================================
% INITIALIZATION
% ============================================================

lambda       = range_lam;
sigma_abs_qe = zeros(size(range_lam));
sigma_abs    = zeros(size(range_lam));
deltaT       = zeros(size(range_lam));
deltaT_qe    = zeros(size(range_lam));

eps_m        = mat.n_medium^2;
k_m          = 2*pi*mat.n_medium ./ range_lam;
mat.tau1     = 1.7e-12; %Characteristic time for gold.
src.tau2     = Inf;     %CW illumination.
src.F        = I;       %Intensity of CW.
ind=1;                  %Defoult: finite Interfacial Thermal Conductance.

% ============================================================
% SPECTRAL OPTICAL CALCULATION
% ============================================================

for i = 1:numel(range_lam)
    % Select optical input at this wavelength
    switch core_mode
        case 'material'
            core_i = mat.n_core;
        case 'numeric'
            core_i = mat.n_core(i);
    end

    % MIE CALCULATION
    [ns, sigma_abs(i)] = mie_absorption(dom.a, core_i, mat.n_medium,...
        range_lam(i));

    % QUASI-STATIC CALCULATION
    eps_s = ns^2;
    ratio = (eps_s - eps_m) / ...
        (eps_s + 2*eps_m);
    sigma_abs_qe(i) = ...
        4*pi*k_m(i)*dom.a^3*imag(ratio) ...
        - (8*pi/3)*k_m(i)^4*dom.a^6*abs(ratio)^2;
end

if isinf(mat.ITC) %Check for Interfacial Thermal Conductance.
    ind=0;
end

% ============================================================
% SPECTRAL THERMAL CALCULATION
% ============================================================

for i=1:numel(range_lam)
    src.lda0=range_lam(i);

    % NUMERICAL CALCULATION
    [~, ~, T, ~, ~, ~]=heatEqCNFD_solver(dom,mat,src,mesh,bc, cond);
    deltaT(i)=max(max(T));

    % QUASI-STATIC CALCULATION
    deltaT_qe(i)=src.F*sigma_abs_qe(i)/(4*pi*mat.k(2)*dom.a)+...
        ind*src.F*sigma_abs_qe(i)/(4*pi*dom.a^2*mat.ITC)  ;
end

% ========================================================
% PLOT RESULTS
% ========================================================

if doPlot
    figure;

    % ABSORPTION CROSS SECTION
    subplot(1,2,1)
    plot(lambda*1e9, sigma_abs, ...
        'k-', 'LineWidth', 1.8);
    hold on
    plot(lambda*1e9, sigma_abs_qe, ...
        'r--', 'LineWidth', 1.8);
    xlabel('\lambda (nm)')
    ylabel('\sigma_{abs} (m^2)')
    legend('Mie', 'Quasi-static', ...
        'Location','best')
    grid on
    box on

    % TEMPERATURE INCREASE
    subplot(1,2,2)
    plot(lambda*1e9, deltaT, ...
        'k-', 'LineWidth', 1.8);
    hold on
    plot(lambda*1e9, deltaT_qe, ...
        'r--', 'LineWidth', 1.8);
    xlabel('\lambda (nm)')
    ylabel('\DeltaT_{max} (K)')
    legend('Numerical', 'Quasi-static', ...
        'Location','best')
    grid on
    box on

end
end
