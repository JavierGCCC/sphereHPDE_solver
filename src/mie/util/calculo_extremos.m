%% Define test parameters
N_pol=7;
N_au=1;
N_core=1;
nCore=1.51;
nPol=1.41;
nMedium=1.33;
A=importdata('datos_Au.mat');
dom=A(A(:,1)>448e-9 & A(:,1)<1202e-9,1);
B=A(A(:,1)>448e-9 & A(:,1)<1202e-9,2:3);
epsAu = B(:,1) + 1i*B(:,2); 
nAu = sqrt(epsAu);

rCore = linspace(36e-9, 36e-9, N_core);
deltaAu = linspace(5e-9, 5e-9, N_au);
deltaPol = linspace(0e-9, 60e-9, N_pol);
%% 
sigmaA = cell(length(rCore), length(deltaAu));
sigmaApart = zeros(length(dom), length(deltaPol)+1);
sigmaApart(:,1) = dom*1e9;

QA = cell(length(rCore), length(deltaAu));
QApart = zeros(length(dom), length(deltaPol)+1);
QApart(:,1) = dom*1e9;

totalSteps = N_core * N_au * N_pol;
stepCounter = 0;

hWait = waitbar(0, 'Calculando espectros...');
%% 
for i = 1:length(rCore)
    for j = 1:length(deltaAu)
        for k = 1:length(deltaPol)
            stepCounter = stepCounter + 1;

            for l = 1:length(dom)
                % stratified sphere
                rad = [rCore(i), rCore(i)+deltaAu(j), rCore(i)+deltaAu(j)+deltaPol(k)];
                ns = [nCore, nAu(l), nPol];
                lambda = dom(l);
                nang = 3000;
                conv = 1;

                % Calcular matriz de dispersión
                [S, C, ang] = calcmie(rad, ns, nMedium, lambda, nang, ...
                    'ConvergenceFactor', conv);
                A = pi * (rCore(i)+deltaAu(j)+deltaPol(k))^2;
                sigmaApart(l, k+1) = C.abs;
                QApart(l, k+1) = C.abs / A;

                % Descomentar para debug: [i j k l]
            end

            sigmaA{i,j} = sigmaApart;
            QA{i,j} = QApart;

            % Actualización eficiente del waitbar
            if mod(stepCounter, round(0.05*totalSteps)) == 0 || stepCounter == totalSteps
                waitbar(stepCounter / totalSteps, hWait, ...
                    sprintf('Progreso: %.1f %%', 100 * stepCounter / totalSteps));
            end
        end
    end
end

close(hWait); % Cerrar la barra de progreso


%% Análisis espectral: búsqueda del pico más al infrarrojo
lambdaMaxMap = cell(length(rCore), 1);
absMaxMap = cell(length(rCore), 1);
QMaxMap = cell(length(rCore), 1);
umbral_resonancia = 1e-16;

for i = 1:length(rCore)
    lambdaMap_i = NaN(length(deltaAu), length(deltaPol));
    absMap_i = NaN(length(deltaAu), length(deltaPol));
    QMap_i = NaN(length(deltaAu), length(deltaPol));

    for j = 1:length(deltaAu)
        specAbs = sigmaA{i,j}; % (N_lambda x N_pol+1)
        specQ = QA{i,j};
        lambda = specAbs(:,1); % vector lambda

        for k = 1:length(deltaPol)
            s = specAbs(:,k+1); % sigma_abs
            q = specQ(:,k+1);   % eficiencia Q_abs

            if max(s) < umbral_resonancia
                lambdaMax = NaN;
                absMax = NaN;
                qMax = NaN;
            else
                [pks, locs] = findpeaks(s, lambda);

                if isempty(pks)
                    lambdaMax = NaN;
                    absMax = NaN;
                    qMax = NaN;
                else
                    [lambdaMax, idxMax] = max(locs);
                    absMax = pks(idxMax);
                    qMax = q(lambda == lambdaMax); % Interp si es necesario
                end
            end

            lambdaMap_i(j,k) = lambdaMax;
            absMap_i(j,k) = absMax;
            QMap_i(j,k) = qMax;
        end
    end

    lambdaMaxMap{i} = lambdaMap_i;
    absMaxMap{i} = absMap_i;
    QMaxMap{i} = QMap_i;
end