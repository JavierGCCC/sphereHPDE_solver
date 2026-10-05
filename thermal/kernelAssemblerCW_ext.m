
function [F,Q_unit, nodes, Q_boundary]=kernelAssemblerCW_ext(a, L, N, k, Tinf, boundary, ITC)
%Define the input and usage. For later...
%a=[a1, a2, a3,...].
%L\propto max(a);
%N=[N1,N2,N3, N4, ...]--> length(N)=length(a)+1;
%k=[k1,k2,k3,k4,...]----> length(k)=length(a)+1;
%Tinf from outside. size(Itinf)=1x1. does not change. 
%Boundary: does not change.
%ITC=[ITC1, ITC2, ITC3, ...]-->length(ITC)=length(a).

nInterfaces=numel(a);
nRegions=nInterfaces+1;

%Check input self-consistency.
if numel(N)~=nRegions
    error('length(N) must be length(a)+1.');
end

if numel(k)~=nRegions
    error('length(k) must be length(a)+1.');
end

if numel(ITC)~=nInterfaces
    error('length(ITC) must be length(a).');
end

if any(diff(a)<=0)
    error('The list of radii must be strictly increasing.');
end

if L<=a(end)
    error('L must be larger than the particle size.');
end

if any(N<3)
    error('At least three nodes are required in every region.');
end

if any(N~=round(N))
    error('The number of nodes N must contain integer values.');
end

if any(ITC<=0 & ~isinf(ITC))
    error('ITC values must be positive or Inf.');

end


%Define type of boundary.
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

bounds = [0 a L]; %Unify bounds from input data.

%&Preallocate variables. 
Ntot=sum(N)-nInterfaces;
nodes=zeros(1,Ntot);
dr=zeros(1,nRegions);
regionStart=zeros(1,nRegions);
regionEnd=zeros(1,nRegions);
cursor=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%PREVIOUS VERSION %%%%%%%%%%%%%%%%%
% %Thermophysical properties.
% k_p=k(1);
% k_e=k(2);
% if use_ITC, G=ITC; end

% 
% %Space & Time Discretization.
% N_p=N(1); %Number of nodes in particle.
% N_e=N(2); %Number of nodes in environment.
% nodes_p=linspace(0,a,N_p); dx_p=nodes_p(2)-nodes_p(1);
% nodes_e=linspace(a,L,N_e);  dx_e=nodes_e(2)-nodes_e(1);
% nodes=[nodes_p(1:end-1) nodes_e]; N=length(nodes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Definition of radial grid.
for j=1:nRegions
    node_aux=linspace(bounds(j),bounds(j+1),N(j)); %auxiliary radii.
    dr(j)=node_aux(2)-node_aux(1); %Radial step.
    if j<nRegions
        node_aux=node_aux(1:end-1);
    end
    idx=cursor:(cursor+numel(node_aux)-1);
    nodes(idx)=node_aux;
    regionStart(j)=idx(1);
    regionEnd(j)=idx(end);
    cursor=idx(end)+1;
end
interfaceIdx=regionStart(2:end);

%Preallocate matrix diagonals.
lower=zeros(Ntot,1);
main=zeros(Ntot,1);
upper=zeros(Ntot,1);

%Preallocate source matrices. 
Q_unit=zeros(Ntot, nInterfaces);
Q_boundary=zeros(Ntot,1);

%Layer volumes.
V=zeros(1,nInterfaces);
for j=1:nInterfaces
V(j)=4/3*pi*(bounds(j+1)^3-bounds(j)^3);
end

%Ghost-node treatment.
main(1)=-1;
upper(1)=1;
Q_unit(1,1)=-dr(1)^2/(6*k(1)*V(1));

for j = 1:nRegions
    if j == 1
        idx = (regionStart(j)+1):regionEnd(j);
    elseif j < nRegions
        idx = (regionStart(j)+1):regionEnd(j);
    else
        idx = (regionStart(j)+1):(regionEnd(j)-1);
    end

    if isempty(idx)
        continue
    end

    % Radial finite-difference coefficients
    aux = dr(j)./nodes(idx);
    lower(idx) = 1-aux;
    main(idx) = -2;
    upper(idx) = 1+aux;

    % Unit absorbed-power source.
    if j <= nInterfaces
        Q_unit(idx,j) = -dr(j)^2/(k(j)*V(j));
    end
end

% Thermal interfaces.
for j = 1:nInterfaces

    % Global interface index.
    i = interfaceIdx(j);

    % Properties on both sides.
    kL = k(j);
    kR = k(j+1);
    drL = dr(j);
    drR = dr(j+1);
    G = ITC(j);

    %Perfect thermal contact.
    if isinf(G)
        s1 = kL*drR;
        s2 = kR*drL;
        lower(i) = -s1/(s1+s2);
        main(i) = 1;
        upper(i) = -s2/(s1+s2);
        Q_unit(i,:) = 0;

        %Finite thermal interfacial conductance.
    else
        %Left.
        alphaL = kL/(kL+G*drL);
        betaL = G*drL/(kL+G*drL);
        lower(i) = -alphaL;
        main(i) = 1;
        upper(i) = -betaL;      
        
        %Right.
        alphaR = G*drR/(G*drR+kR);
        betaR = kR/(G*drR+kR);
        lower(i+1) = alphaR;
        main(i+1) = -1;
        upper(i+1) = betaR;

        %No heat source at the interfaces.
        Q_unit(i,:) = 0;
        Q_unit(i+1,:) = 0;
    end
end

% Ecternal boundary conditions. 
i = Ntot;
dx_e = dr(end);
R = nodes(end);
switch boundary
    case 'Dirichlet'
        main(i) = -2;
        lower(i) = 1-dx_e/R;
        Q_boundary(i) = -(1+dx_e/R)*Tinf;
    case 'Neumann'
        main(i) = -1;
        lower(i) = 1;
        Q_boundary(i) = 0;
    case 'Robin'
        aux_robin = dx_e/R*(1+dx_e/R);
        main(i) = -(1+aux_robin);
        lower(i) = 1;
        Q_boundary(i) = -aux_robin*Tinf;
end

%MATRIX ASSEMBLY.

F1 = [lower(2:end); 0];
F2 = main;
F3 = [0; upper(1:end-1)];
F = spdiags([F1 F2 F3], [-1 0 1], Ntot, Ntot);
end

