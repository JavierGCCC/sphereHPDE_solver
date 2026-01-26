dom.a=50e-9;
dom.R_sim=40*dom.a;
mat.rho=[19300, 1000];
mat.cp=[129, 4181];
mat.k=[318,0.6];
src.sigma=1e-15;
src.F=1e9;
dom.t_max=100e-9;
mesh.N_p=200;
mesh.N_e=round((dom.R_sim-dom.a)/dom.a*mesh.N_p);
mesh.N_t=10000;
bc.bound='Dirichlet';
bc.T_ini=0;

[t, nodes, T, info_solver]=heatEqCNFD_solver(dom,mat,src,mesh,bc);

