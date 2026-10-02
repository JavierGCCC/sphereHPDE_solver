function solution_cw=heatEqCNFD_solverCW(dom,mat,src,mesh,bc)

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

% % ==MATRIX PRECONDITIONING==
% target_cond=cond.target_cond;
% iter_max=cond.iter_max;
% time_max=cond.time_max;

% ==SOLVER==
Q=zeros(1,numel(lda0));
for i=1:numel(lda0)
[~,sigma]=mie_absorption(a, n_core, n_medium, lda0(i)); %Calculate absorption.

Qfunc = illuminationCW(sigma, F, tau1, Inf); %Set illumination function.
Q(i)=Qfunc(1);
end

[K,Q_unit, nodes]=kernelAssemblerCW(a, R_sim,...
    N, k, T_ini, bound, ITC); %Build core matrix.


Q_unit_a = Q_unit; Q_unit_a(end) = 0;
Q_unit_b = zeros(size(Q_unit)); Q_unit_b(end) = Q_unit(end);

X = K\[Q_unit_a, Q_unit_b];

T_unit = X(:,1);
T_0    = X(:,2);

T=T_unit*Q+T_0;


solution_cw.T   = T;
solution_cw.Q=Q;
solution_cw.nodes      = nodes;
solution_cw.lambda     = lda0;

solution_cw.fluence = F;
solution_cw.tau2    = tau2;


end