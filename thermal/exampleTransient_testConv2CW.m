dom.a=25e-9;
dom.R_sim=40*dom.a;


mat.n_core='Au';
mat.n_medium=1.33;
mat.rho=[19300, 1000];
mat.cp=[129, 4181];
mat.k=[318,0.6];
mat.ITC=2e8; 
mat.ITC=Inf;
mat.tau1=1.7e-12;

src.lda0=linspace(500,600,100)*1e-9;
src.F=1e8;
src.tau2=1000*mat.tau1;
%src.tau2=1e5*mat.tau1;
src.tau2=Inf;
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

% cond.target_cond=1000;
% cond.iter_max=10;
% cond.time_max=360;
duracion=linspace(100,10000,25)*1e-9;
residual=zeros(1,length(duracion));
for j=1:length(duracion)
dom.t_max=duracion(j);

solution_t=heatEqCNFD_solver_t(dom,mat,src,mesh,bc)

solution_cw=heatEqCNFD_solverCW(dom,mat,src,mesh,bc);

%%

nodes_t=solution_t.nodes;N_nodes=length(nodes);
lam_t=solution_t.lambda; N_lam=length(lam);
timeDomain=solution_t.timeDomain; N_t=length(timeDomain);
T_t=zeros(1,N_lam);
for i=1:N_lam
 T_t(i) = max (...
    solution_t.spectralWeights(i) * ...
    solution_t.opticalBasis(:,end) + ...
    solution_t.boundaryBasis(:,end));
end

nodes_cw=solution_cw.nodes;N_nodes=length(nodes);
lam_cw=solution_cw.lambda; N_lam=length(lam);
T_cw=solution_cw.T;

residual(j)=max(abs(max(T_cw)-T_t)./max(T_cw)*100);
% figure(1)
% plot(src.lda0, max(T_cw));hold on;
% plot(lam_t, T_t);


end

figure(2);
plot(duracion, residual);
xlabel('pulse duration (ns)');
ylabel('Error (%)');
set(gca,'FontSize',22);