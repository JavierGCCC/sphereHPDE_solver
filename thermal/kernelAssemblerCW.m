
function [F,Q_unit, nodes]=kernelAssemblerCW(a, L, N, k, Tinf, boundary, ITC)


if isinf(ITC)
    use_ITC = false;
    disp('Infinite ITC-->Ill-conditioned core matrix...');
    disp('Switching to continuous solver');
else
    use_ITC = true;
end

switch boundary
    case 'Dirichlet'
        fprintf('>> Using boundary condition: Dirichlet\n');
        fprintf('   Equation: T(R,t) = T_amb\n');

    case 'Neumann'
        fprintf('>> Using boundary condition: Neumann (insulated/reflecting)\n');
        fprintf('   Equation: dT/dr |_(R) = 0\n');

    case 'Robin'
        fprintf('>> Using boundary condition: Robin (convective/absorbing)\n');
        fprintf('   Equation: -k * dT/dr |_(R) = h * (T(R,t) - T_amb)\n');

    otherwise
        error('Invalid boundary condition "%s". Please choose Dirichlet, Neumann, or Robin.', boundary);
end




%Thermophysical properties.
k_p=k(1);
k_e=k(2);
if use_ITC, G=ITC; end

%Space & Time Discretization.
N_p=N(1); %Number of nodes in particle.
N_e=N(2); %Number of nodes in environment.
nodes_p=linspace(0,a,N_p); dx_p=nodes_p(2)-nodes_p(1);
nodes_e=linspace(a,L,N_e);  dx_e=nodes_e(2)-nodes_e(1);
nodes=[nodes_p(1:end-1) nodes_e]; N=length(nodes);

%DEFINITION OF MATRIX ELEMENTS:

%Crank-Nicolson matrix elements inside particle and environment.
aux_p=dx_p./nodes(2:N_p-1);
aux_e=dx_e./nodes(N_p+1:N-1);


%Flux conservation BC (no ITC).
s1=k_p*dx_e;
s2=k_e*dx_p;



%Source.
V=4/3*pi*a*a*a;
Q_unit=[-dx_p^2/(6*k_p), -dx_p^2/k_p*ones(size(aux_p)),0, ...
    zeros(size(aux_e)),pi]/V; Q_unit=Q_unit';



%Matrix assembly:
%Lower diagonals.
F1=[1-aux_p,-s1/(s1+s2),1-aux_e,pi,0];

%Upper diagonals.
F3=[0, 1, 1+aux_p, -s2/(s1+s2), 1+aux_e];

%Principal diagonal.
F2=[-1, -2*ones(1,N_p-2), 1, -2*ones(1,N_e-2),pi];

if use_ITC
    alpha_p = k_p/(k_p + G*dx_p);
    beta_p  = G*dx_p/(k_p + G*dx_p);

    alpha_e = G*dx_e/(G*dx_e + k_e);
    beta_e  = k_e/(G*dx_e + k_e);

    % Particle side
    F1(N_p-1) = -alpha_p;
    F2(N_p)   = 1;
    F3(N_p+1) = -beta_p;

    % Environment side
    F1(N_p)   = alpha_e;
    F2(N_p+1) = -1;
    F3(N_p+2) = beta_e;
end

switch boundary
    case 'Dirichlet'
        F2(N)=-2; F1(N-1)=1-dx_e/nodes(N);

        Q_unit(N)=-(1+dx_e/nodes(N))*Tinf; %Dirichlet condition.
    case 'Neumann'
        F2(N)=-1; F1(N-1)=1;

        Q_unit(N)=0; %Neumann condition.
    case 'Robin'
        aux_robin=dx_e/nodes(N)*(1+dx_e/nodes(N));
        F2(N)=-(1+aux_robin); F1(N-1)=1;
        Q_unit(N)=-aux_robin*Tinf; %Robin condition.
end

%Kernel matrix.
F = spdiags([F1(:), F2(:), F3(:)], [-1, 0, 1], N, N);


end

