function ns = resolveLayerIndices(nSphere, nParticleLayers, lda0, dataDir)
%============================================================
% NUMERIC INPUT
% Example:
%   nSphere = [2.3+0.1i 1.45 0.5+2i];
% ============================================================
if isnumeric(nSphere)
    if numel(nSphere) ~= nParticleLayers
        error([ ...
            'nSphere must contain one refractive index ' ...
            'for each particle layer.']);
    end
    ns = reshape(nSphere,1,[]);
    return
end
% ============================================================
% SINGLE MATERIAL NAME
% ============================================================
if (ischar(nSphere) || isstring(nSphere)) && nParticleLayers == 1
    nSphere = {char(nSphere)};
end

% ============================================================
% STRING ARRAY
% ============================================================
if isstring(nSphere)
    nSphere = cellstr(nSphere);
end
% ============================================================
% CELL ARRAY
%
% Example:
%
%   nSphere = {'Fe','SiO2','Au'};
%
% or:
%
%   nSphere = {2.4+0.1i,'SiO2','Au'};
% ============================================================
if ~iscell(nSphere)
    error([ ...
        'nSphere must be either a numeric vector or a cell ' ...
        'array containing one optical definition per layer.']);
end

if numel(nSphere) ~= nParticleLayers
    error([ ...
        'nSphere must contain one optical definition ' ...
        'for every particle layer.']);
end

ns = zeros(1,nParticleLayers);

for j = 1:nParticleLayers

    material = nSphere{j};
    % --------------------------------------------------------
    % Direct numerical refractive index
    % ---------------------------------------------------------
    if isnumeric(material)
        if ~isscalar(material)
            error( ...
                'Numerical refractive indices must be scalar.');
        end
        ns(j) = material;
        
    % --------------------------------------------------------
    % Material database
    % ---------------------------------------------------------
    elseif ischar(material) || isstring(material)
        material = char(material);
        dataFileMat = fullfile(dataDir,['data_' material '.mat']);
        dataFileTxt = fullfile(dataDir,['data_' material '.txt']);
        if exist(dataFileMat,'file') == 2
            data = importdata(dataFileMat);
        elseif exist(dataFileTxt,'file') == 2
            data = importdata(dataFileTxt);
        else
            error( ...
                ['Optical data for "%s" not found in ' ...
                 'optical/optical_properties/.'], material);
        end
        % Expected format:
        %
        %   [lambda, Re(epsilon), Im(epsilon)]
        r_eps = interp1(data(:,1), data(:,2), lda0, 'spline');
        i_eps = interp1(data(:,1), data(:,3), lda0, 'spline');
        if isnan(r_eps) || isnan(i_eps)
            error( ...
                ['Wavelength %.2f nm lies outside the ' ...
                 'optical data range for "%s".'], lda0*1e9, material);
        end
        epsMat = r_eps + 1i*i_eps;
        ns(j) = sqrt(epsMat);
    else
        error( ...
            'Invalid optical definition for particle layer %d.', j);
    end
end
end