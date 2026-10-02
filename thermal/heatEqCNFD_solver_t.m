function solution=heatEqCNFD_solver_t(dom,mat,src,mesh,bc)

% ==SIMULATION DOMAIN==
a=dom.a;
R_sim=dom.R_sim;
t_max=dom.t_max;

% ==MATERIAL THERMOPHYSICAL PROPERTIES==
n_core=mat.n_core;
n_medium=mat.n_medium;
rho=mat.rho;
cp=mat.cp;
k=mat.k;
ITC=mat.ITC;
tau1=mat.tau1;

% ==SOURCE CHARACTERISTICS==
lda0=src.lda0;
F=src.F;
tau2=src.tau2;

% ==MESH CONSTRUCTION==
N=[mesh.N_p, mesh.N_e, mesh.N_t];

% ==BOUNDARY CONDITIONS==
bound=bc.bound;
T_ini=bc.T_ini;

% ==MATRIX PRECONDITIONING==
% target_cond=cond.target_cond;
% iter_max=cond.iter_max;
% time_max=cond.time_max;

% ==SOLVER==
sigma=zeros(1,numel(lda0));
for i=1:numel(lda0)
[~,sigma(i)]=mie_absorption(a, n_core, n_medium, lda0(i)); %Calculate absorption.
end

Qfunc = illumination(F, tau1, tau2); %Set illumination function.

[F_d,F_a,Q_vals, Q_spatial, Q_boundary, nodes, t]=kernelAssembler(a, R_sim, t_max, N, rho, cp,...
    k, Qfunc, T_ini,bound, ITC); %Build core matrix.

t_ref = tic;
dFd = decomposition(F_d);

T_optical  = zeros(length(nodes),length(t));
T_boundary = zeros(length(nodes),length(t));
T_boundary(:,1) = T_ini;

Nt = length(t);


h = waitbar(0,'Solving heat transfer...');

updateEvery = max(1,floor((Nt-1)/100));

for n = 1:Nt-1

    % Optical contribution
    rhs_optical = ...
        F_a*T_optical(:,n) + ...
        Q_spatial*Q_vals(n);

    % Boundary contribution
    rhs_boundary = ...
        F_a*T_boundary(:,n) + ...
        Q_boundary;

    % Solve both RHS simultaneously
    RHS = [rhs_optical, rhs_boundary];
    X = dFd\RHS;
    T_optical(:,n+1)  = X(:,1);
    T_boundary(:,n+1) = X(:,2);

    if mod(n,updateEvery)==0 || n==Nt-1
        waitbar(n/(Nt-1),h);
    end

end
elapsed_time = toc(t_ref);
close(h);

solution.opticalBasis   = T_optical;
solution.boundaryBasis  = T_boundary;
solution.spectralWeights = sigma;

solution.timeDomain = t;
solution.nodes      = nodes;
solution.lambda     = lda0;

solution.fluence = F;
solution.tau2    = tau2;

solution.info.elapsedTime = elapsed_time;

%Reconstruction.
% T_lambda = ...
%     solution.spectralWeights(iLambda) * ...
%     solution.opticalBasis + ...
%     solution.boundaryBasis;


% [L,U]=ilu(F_d); %LU factorization.
% 
% [M,info_cond]=preconditioner(U,target_cond,iter_max,time_max); %Precond.
% 
% [T, info_solver]=CN_linearSolver(F_d, F_a, Q, T_ini, M,L,U); %Solver.

end