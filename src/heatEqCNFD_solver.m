function [t, nodes, T, info_solver]=heatEqCNFD_solver(dom,mat,src,mesh,bc)
%% READING
%% ==SIMULATION DOMAIN==
a=dom.a;
R_sim=dom.R_sim;
t_max=dom.t_max;

%% ==MATERIAL THERMOPHYSICAL PROPERTIES==
rho=mat.rho;
cp=mat.cp;
k=mat.k;

%% ==SOURCE CHARACTERISTICS==
sigma=src.sigma;
F=src.F;

%% ==MESH CONSTRUCTION==
N=[mesh.N_p, mesh.N_e, mesh.N_t];

%% ==BOUNDARY CONDITIONS==
bound=bc.bound;
T_ini=bc.T_ini;

%% ==SOLVER==

[Qfunc] = illumination(sigma, F);

[F_d,F_a,Q, nodes, t]=kernelAssembler(a, R_sim, t_max, N, rho, cp,...
    k, Qfunc, T_ini,bound);

[T, info_solver]=CN_linearSolver(F_d, F_a, Q, T_ini);
end