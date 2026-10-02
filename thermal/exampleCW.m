dom.a=25e-9;
dom.R_sim=40*dom.a;


mat.n_core='Au';
mat.n_medium=1.33;
mat.rho=[19300, 1000];
mat.cp=[129, 4181];
mat.k=[318,0.6];
mat.ITC=2e8; 
%mat.ITC=Inf;
mat.tau1=1.7e-12;

src.lda0=linspace(500,600,110)*1e-9;
src.F=1e8;
src.tau2=10000000*mat.tau1;
src.tau2=1e5*mat.tau1;
dom.t_max=max([10*src.tau2,100e-9]);


% %Discretization.
% N=zeros(1,dim_rad);
% N(1)=200;
% for i=2:dim_rad
%     N(i) = round((a(i)-a(i-1))/a(1)*N(1));
% end
% mesh.N=N;
mesh.N_p=250;
mesh.N_e=round((dom.R_sim-dom.a)/dom.a*mesh.N_p);
%mesh.N_t=max([10000,10*dom.t_max/src.tau2])
mesh.N_t=10000;

bc.bound='Robin';
bc.T_ini=0;

cond.target_cond=1000;
cond.iter_max=10;
cond.time_max=360;

[nodes, T, Q, F]=heatEqCNFD_solverCW(dom,mat,src,mesh,bc,cond);

plot(src.lda0, max(T));