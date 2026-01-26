
function [F_d,F_a,Q, nodes, t]=kernelAssembler(a, L, t_max, N, rho, cp, k, Qfunc, Tamb, boundary)
%rho=[rho_p, rho_e];
%cp=[cp_p, cp_e];
%k=[k_p, k_e];
%G=[G];
%Tamb=Tamb;
%a=a;
%L=L;
%t_max=t_max;
%N=[N_p,N_e,N_t];
%Q_func=Q_func;
Tinf=Tamb;
%Check external parameters.

switch boundary
    case 'Dirichlet'
        fprintf('>> Using boundary condition: Dirichlet\n');
        fprintf('   Equation: T(R,t) = T_amb\n');

    case 'Neumann'
        fprintf('>> Using boundary condition: Neumann (insulated/reflecting)\n');
        fprintf('   Equation: dT/dr |_(R) = 0\n');

    otherwise
        error('Invalid boundary condition "%s". Please choose Dirichlet, Neumann, or Robin.', boundary);
end

%Thermophysical properties.
rho_p=rho(1); cp_p=cp(1); k_p=k(1); alpha_p=k_p/(rho_p*cp_p);
rho_e=rho(2); cp_e=cp(2); k_e=k(2); alpha_e=k_e/(rho_e*cp_e);

%Space & Time Discretization.
N_p=N(1); %Number of nodes in particle.
N_e=N(2); %Number of nodes in environment.
N_t=N(3); %Number of temporal nodes.
nodes_p=linspace(0,a,N_p); dx_p=nodes_p(2)-nodes_p(1);
nodes_e=linspace(a,L,N_e);  dx_e=nodes_e(2)-nodes_e(1);
nodes=[nodes_p(1:end-1) nodes_e]; N=length(nodes);
t=linspace(0,t_max,N_t); dt=t(2)-t(1);

%DEFINITION OF MATRIX ELEMENTS:
%Simmetry BC.
k0=1.5*dt/(dx_p*dx_p)*alpha_p;

%Crank-Nicolson matrix elements inside particle and environment.
k1_p=0.5*dt/(dx_p*dx_p)*alpha_p;
k2_p=0.5*dt/dx_p*alpha_p./nodes_p(2:end-1);
k1_e=0.5*dt/(dx_e*dx_e)*alpha_e;
k2_e=0.5*dt/dx_e*alpha_e./nodes_e(2:end-1);

%Flux conservation BC (no ITC).
s1=k_p*dx_e;
s2=k_e*dx_p;
s=s1/s2;

%Source.
q=.5*dt/(rho_p*cp_p); V=4/3*pi*a*a*a; q=q/V;
Q_vals = q*(Qfunc(t(1:end-1))+Qfunc(t(2:end)));        % (1 × N_t-1)
Q_spatial = [ones(N_p-1, 1); zeros(N_e, 1)];            % (N × 1)
Q = Q_spatial * Q_vals;


%Matrix assembly:
%Lower diagonals.
F1_d=[-k1_p+k2_p,-s,-k1_e+k2_e,pi,0]; F1_a=-F1_d;

%Upper diagonals.
F3_d=[0, -2*k0, -k1_p-k2_p, -1, -k1_e-k2_e]; F3_a=-F3_d;

%Principal diagonal.
Aux=[2*k0, 2*k1_p*ones(1,N_p-2), pi, 2*k1_e*ones(1,N_e-2),pi];
F2_d=1+Aux; F2_d(N_p)=s+1;
F2_a=1-Aux; F2_a(N_p)=-s-1;


switch boundary
    case 'Dirichlet'
        F2_d(N)=1; F1_d(N-1)=0;
        F2_a(N)=0; F1_a(N-1)=0;
        Q(N,:)=Tinf*ones(1,N_t-1); %Dirichlet condition.
    case 'Neumann'
        F2_d(N)=1+2*k1_e; F1_d(N-1)=-2*k1_e;
        F2_a(N)=1-2*k1_e; F1_a(N-1)=2*k1_e;
        Q(N,:)=zeros(1,N_t-1); %Dirichlet condition.
end

%Kernel matrix.
F_a = spdiags([F1_a(:), F2_a(:), F3_a(:)], [-1, 0, 1], N, N);
F_d = spdiags([F1_d(:), F2_d(:), F3_d(:)], [-1, 0, 1], N, N);

end

