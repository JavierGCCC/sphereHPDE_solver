
%DEFINITION OF VARIABLES.
doSave = true;

%Domain
dom.a            = 15e-9;%NP radius.
dom.R_sim        = 40*dom.a;%Simulation domain radius.

%Materials.
mat.n_core       = 'Au';%Particle RI.
mat.n_medium     = 1.33;%Environment RI.
mat.rho          = [19300, 1000];%Density.
mat.cp           = [129, 4181];%Specific heat.
mat.k            = [318,0.6];%Thermal cond.
mat.ITC          = 2e8;%Interfacial thermal conductance.

%illumination
range_lam        = linspace(500,600,10)*1e-9;%Wavelenght range.
I                = 1e9; %Intensity of CW illumination.
dom.t_max        = 100e-9;%Simulated time.

%Space/time discretization.
mesh.N_p         = 250;%Number of spatial elements inside NP.
mesh.N_e         = round(...
    (dom.R_sim-dom.a)/dom.a*mesh.N_p);%Number of spatial elements outside.
mesh.N_t         = 10000;%Number of temporal nodes.

%Thermal boundary condition.
bc.bound         ='Robin';%BC type.
bc.T_ini         = 0;%Initial temperature.

%Core matrix conditioning
cond.target_cond = 1000;%Maximum condition number allowed.
cond.iter_max    = 10;%Maximum number of iterations allowed.
cond.time_max    = 360;%Maximum time for conditioning in seconds.

%Visualization.
doPlot = true;

%CALCULATION.
[lambda, sigma_abs, sigma_abs_qe, deltaT, deltaT_qe] ...
    = check_QA(dom, mat, range_lam, I, mesh, bc, cond, doPlot);

%SAVEDATA.
if doSave 
    out.lambda       = lambda;
    out.sigma_abs    = sigma_abs;
    out.sigma_abs_qe = sigma_abs_qe;
    out.deltaT       = deltaT;
    out.deltaT_qe    = deltaT_qe;
    save('check_QA_results.mat','out');
end