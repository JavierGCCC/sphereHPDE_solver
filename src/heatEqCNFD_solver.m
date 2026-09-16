function [t, nodes, T, Qfunc, info_solver, info_cond]=heatEqCNFD_solver(dom,mat,src,mesh,bc,cond)
%% READING
%% ==SIMULATION DOMAIN==
a=dom.a;
R_sim=dom.R_sim;
t_max=dom.t_max;

%% ==MATERIAL THERMOPHYSICAL PROPERTIES==
n_core=mat.n_core;
n_medium=mat.n_medium;
rho=mat.rho;
cp=mat.cp;
k=mat.k;
ITC=mat.ITC;
tau1=mat.tau1;

%% ==SOURCE CHARACTERISTICS==
lda0=src.lda0;
F=src.F;
tau2=src.tau2;

%% ==MESH CONSTRUCTION==
N=[mesh.N_p, mesh.N_e, mesh.N_t];

%% ==BOUNDARY CONDITIONS==
bound=bc.bound;
T_ini=bc.T_ini;

%% ==MATRIX PRECONDITIONING==
target_cond=cond.target_cond;
iter_max=cond.iter_max;
time_max=cond.time_max;

%% ==SOLVER==

sigma=mie_absorption(a, n_core, n_medium, lda0);

[Qfunc] = illumination(sigma, F, tau1, tau2);

[F_d,F_a,Q, nodes, t]=kernelAssembler(a, R_sim, t_max, N, rho, cp,...
    k, Qfunc, T_ini,bound, ITC);

[L,U]=ilu(F_d);

[M,info_cond]=preconditioner(U,target_cond,iter_max,time_max);

[T, info_solver]=CN_linearSolver(F_d, F_a, Q, T_ini, M,L,U);

end