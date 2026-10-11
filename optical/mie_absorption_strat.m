function [ns,sigma_abs,sigma_layers,info] = ...
    mie_absorption_strat(a,nSphere,nMedium,lda0)

% MIE_ABSORPTION
%
% Computes total and layer-resolved absorption cross-sections
% of a homogeneous or concentric multilayer sphere using MatScat.
%
% INPUTS
%   a
%       Outer radius of each particle layer [m].
%       Homogeneous sphere:
%           a = 25e-9;
%       Core-shell:
%           a = [10e-9 25e-9];
%       Three-layer particle:
%           a = [10e-9 20e-9 25e-9];
%   nSphere
%       Optical properties of each particle layer.
%       Numerical refractive indices:
%           nSphere = [n1 n2 n3];
%       Material database:
%           nSphere = {'Fe','SiO2','Au'};
%       Mixed:
%           nSphere = {2.4+0.1i,'SiO2','Au'};
%   nMedium
%       Refractive index of the surrounding medium.
%       The surrounding medium is assumed non-absorbing.
%   lda0
%       Vacuum wavelength [m].
%
% OUTPUTS
%   ns
%       Complex refractive index of each particle layer.
%   sigma_abs
%       Total absorption cross-section [m^2].
%   sigma_layers
%       Layer-resolved absorption cross-sections [m^2].
%       sigma_layers(j) corresponds to particle layer j.
%   info
%       Diagnostic information:
%           info.sigmaExt
%           info.sigmaSca
%           info.sigmaCumulativeFlux
%           info.sigmaOuterFlux
%           info.closureError
%           info.absorptionFraction
%
% PHYSICAL INTERPRETATION
%
%   CW:
%       P_abs,j = I * sigma_layers(j)
%   Transient:
%       E_abs,j = F * sigma_layers(j)
%
% ASSUMPTION FOR THE THERMAL MODEL
%
%   mie_absorption determines the TOTAL absorbed power associated
%   with each material layer.
%
%   The thermal solver can subsequently assume that this power is
%   uniformly distributed inside thecorresponding layer.
% ============================================================


% ============================================================
% INPUT ORGANIZATION
% ============================================================

a = a(:).';
nParticleLayers = numel(a);

if any(a <= 0)
    error('All particle radii must be positive.');
end

if any(diff(a) <= 0)
    error('Particle radii must be strictly increasing.');
end

if ~isscalar(lda0)
    error('lda0 must currently contain one wavelength.');
end

if ~isscalar(nMedium)
    error('nMedium must be scalar.');
end

if abs(imag(nMedium)) > 1e-12
    error([ ...
        'Layer-resolved absorption currently assumes a ' ...
        'non-absorbing surrounding medium.']);
end

% ============================================================
% OPTICAL DATA PATH
% ============================================================
baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir, 'optical_properties');

% ============================================================
% REFRACTIVE INDICES OF PARTICLE LAYERS
% ============================================================
ns = resolveLayerIndices(nSphere, nParticleLayers, lda0, dataDir);

% ============================================================
% MIE PARAMETERS
% ============================================================
conv = 1;
k = 2*pi/lda0 * nMedium; % Wavenumber in surrounding medium
x = k*a; % Size parameter at every interface
m = ns/nMedium; % Relative refractive index of every layer

% ============================================================
% EXTERNAL MIE COEFFICIENTS
%
% an, bn completely describe the field outside the particle.
%
% For a multilayer sphere MatScat calculates these directly
% without requiring any spatial/angular mesh.
% ============================================================
if nParticleLayers == 1
    [an,bn] = expcoeff_mie(x, m, conv); % Homogeneous sphere.
else
    [an,bn] = expcoeff_mie_strat(x, m, conv); % Stratified sphere.
end

% ============================================================
% TOTAL ABSORPTION
%
% MatScat's asmmie uses:
%   C_ext =
%       2*pi/k^2 sum_n (2n+1) Re(an+bn)
%   C_sca =
%       2*pi/k^2 sum_n (2n+1)(|an|^2+|bn|^2)
%   C_abs = C_ext - C_sca
%
% We calculate these quantities directly.
% Therefore we do NOT need:
%       calcmie(...)
% and we do NOT calculate the angular scattering matrix.
% ============================================================

orders = 1:numel(an);
weights = 2*orders + 1;
sigma_ext = 2*pi/k^2*sum(weights.*real(an + bn));
sigma_sca = 2*pi/k^2*sum(weights.*(abs(an).^2 + abs(bn).^2));
sigma_abs = real(sigma_ext-sigma_sca);

% ============================================================
% HOMOGENEOUS SPHERE
%
% There is only one absorbing domain.
% No internal decomposition is necessary.
% ============================================================
if nParticleLayers == 1
    sigma_layers = sigma_abs;
    info.sigmaExt = sigma_ext;
    info.sigmaSca = sigma_sca;
    info.sigmaCumulativeFlux = sigma_abs;
    info.sigmaOuterFlux = sigma_abs;
    info.closureError = 0;
    
    if abs(sigma_abs) > 0
        info.absorptionFraction = 1;
    else
        info.absorptionFraction = 0;
    end
    return
end

% ============================================================
% INTERNAL MIE COEFFICIENTS
% MatScat uses these four matrices to reconstruct the
% electromagnetic field inside every particle layer.
% Dimensions:
%       number of Mie orders x number of layers
% ============================================================
[fn,gn,vn,wn] = expcoeff_mie_strat_int(an, bn, x, m);

% ============================================================
% MATSCAT CONSTANTS
% ============================================================
MatScat_const;

% ============================================================
% INCIDENT POYNTING FLUX
%
% MatScat uses E0 = 1 and defines:
%       P = real(E x conj(H))
% The incident magnetic field is:
%       H0 = nMedium/(c0*mu0)
% Therefore, using exactly MatScat's normalization:
%       Iinc = nMedium/(c0*mu0)
% ============================================================
Iinc = real(nMedium)/(c0*mue0);

% ============================================================
% ANALYTICALLY INTEGRATED POYNTING FLUX
% At every interface a(j), calculate:
%       Phi(a_j) =
%           integral_S P . n dA
% Phi > 0:
%       net electromagnetic power leaving the enclosed volume.
% Therefore:
%       P_abs(<a_j) = -Phi(a_j)
% and:
%       sigma_abs(<a_j) =
%           -Phi(a_j)/Iinc
% The angular integral is performed analytically using the
% orthogonality of the vector spherical harmonics.
% No theta/phi mesh is required.
% ============================================================
M = numel(an);
n = (1:M)';
fluxInterfaces =zeros(1,nParticleLayers);
sigmaCumulativeFlux = zeros(1,nParticleLayers);

for j = 1:nParticleLayers
    r = a(j);
    % --------------------------------------------------------
    % Optical size parameter inside layer j
    %
    % rho = k_j r
    %     = m_j k_medium r
    % ---------------------------------------------------------
    rho = m(j)*x(j);

    % --------------------------------------------------------
    % SPHERICAL BESSEL FUNCTIONS
    % ---------------------------------------------------------
    jn = sbesselj(n,rho);
    yn = sbessely(n,rho);

    % --------------------------------------------------------
    % DERIVATIVES OF THE RICCATI-BESSEL FUNCTIONS
    %
    % MatScat definitions:
    %   dricbesj =
    %       d[rho*j_n(rho)]/d rho
    % ricbesy is defined as:
    %       -rho*y_n(rho)
    % therefore:
    %       d[rho*y_n]/d rho =
    %           -dricbesy
    % --------------------------------------------------------
    dPsiJ = dricbesj(n,rho);
    dPsiY = -dricbesy(n,rho);

    % --------------------------------------------------------
    % RADIAL COMBINATIONS USED BY MATSCAT
    %
    % These are precisely the j_n + y_n combinations appearing
    % in nfmie.m for a stratified sphere.
    % ---------------------------------------------------------
    F = fn(:,j).*jn+vn(:,j).*yn;
    G = gn(:,j).*jn+wn(:,j).*yn;

    %--------------------------------------------------------
    % Corresponding Riccati-Bessel derivative terms
    % ---------------------------------------------------------

    dF = (fn(:,j).*dPsiJ+vn(:,j).*dPsiY)./rho;
    dG = (gn(:,j).*dPsiJ+wn(:,j).*dPsiY)./rho;

    % --------------------------------------------------------
    % MAGNETIC-FIELD PREFACTOR
    % MatScat uses internally:
    %       k_j/(omega*mu0)
    % and:
    %       k_j/omega = n_j/c0
    % ---------------------------------------------------------
    Hprefactor = ns(j)/(c0*mue0);

    %--------------------------------------------------------
    % INTEGRATED RADIAL POYNTING FLUX
    %
    % Orthogonality removes all cross terms between different
    % Mie orders after integration over the full sphere.
    % What remains is a simple sum over n.
    % ---------------------------------------------------------
    modalFlux = (2*n + 1).*1i.*(F.*conj(dF)-dG.*conj(G));
    fluxInterfaces(j) = 2*pi*r^2*real(conj(Hprefactor)*sum(modalFlux));

    % --------------------------------------------------------
    % CUMULATIVE ABSORPTION
    % ---------------------------------------------------------
    sigmaCumulativeFlux(j) = -fluxInterfaces(j)/Iinc;
end

% ============================================================
% CLOSURE CHECK
%
% There are now two independent ways of obtaining the total
% absorption:
%
%   1) external Mie coefficients:
%
%           sigma_abs
%
%   2) Poynting flux through outer particle surface:
%
%           sigmaCumulativeFlux(end)
%
% They must agree.
% ============================================================
sigmaOuterFlux = sigmaCumulativeFlux(end);
closureError =abs(sigmaOuterFlux-sigma_abs)/max(abs(sigma_abs),realmin);

% ============================================================
% ABSORPTION OF EACH LAYER
%
% For example:
%   sigmaCum(1) = sigma_core
%   sigmaCum(2) = sigma_core + sigma_shell1
%   sigmaCum(3) = sigma_core + sigma_shell1 + sigma_shell2
%
% Therefore:
%       sigma_j =
%           sigmaCum(j) - sigmaCum(j-1)
% For the outermost cumulative absorption we use the total
% farfield cross-section.
%
% This guarantees exact global energy conservation, while the
% independently calculated outer Poynting flux remains available
% through closureError as a validation check.
% ============================================================
sigmaCumulative =sigmaCumulativeFlux;
sigmaCumulative(end) = sigma_abs;
sigma_layers = diff([0 sigmaCumulative]);

% ============================================================
% ABSORPTION FRACTIONS
% ============================================================
if abs(sigma_abs) > realmin
    absorptionFraction =sigma_layers/sigma_abs;
else
    absorptionFraction = zeros(size(sigma_layers));
end

% ============================================================
% DIAGNOSTIC OUTPUT
% ============================================================
info.sigmaExt = sigma_ext;
info.sigmaSca = sigma_sca;
info.sigmaCumulativeFlux = sigmaCumulativeFlux;
info.sigmaOuterFlux = sigmaOuterFlux;
info.closureError = closureError;
info.absorptionFraction = absorptionFraction;
end



