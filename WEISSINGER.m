% this program computes the velocity around a finite 3-D wing
% the wing is invested by an external flow of V_inf = 1
% the wing can have multiple displacement wrt the airstream
% the wing position is described by aerodynamic angles
%
% INPUT: 
%   WING properties:
%       sweep angle            -- lambda                    [deg]
%       dihedral angle         -- delta                     [deg]
%       root chord             -- root
%       semi-wing length       -- L
%       taper ratio            -- taper
%       AOA                    -- wing/tail angle of attack [deg]
%       # of spanwise panels   -- M
%       # of chordwise panels  -- N
%
%   AIRSTREAM properties:
%       airflow incidence      -- alpha [deg]
%       sideslip angle         -- beta  [deg]
%       U                         = 1   [m/s]
%

%% computing coefficients from initial conditions
clc
clear 
close all

tic

% AERODYNAMIC properties
alpha  = 0;
beta   = 10;
% GEOMETRIC properties
delta  = 5;
lambda = 0;
root   = 8;
L      = 30;
taper  = 1;
AOA    = 0;

M = 15;
N = 7;

flag = "plot";

% panel creation function 
[PANELwing] = PANELING(delta,lambda,AOA,root,taper,L,M,N,flag,[0,0,0]);

% system matrix generation
% setting tollerance --> useful to avoid singular MATRIX 
toll        = 1e-4;
[MATRIX]    = BS(PANELwing,AOA,M,N,L,toll);

% system known vector 
[b]         = compute_vector(PANELwing,alpha,beta,M,N);

% solve system
GAMMA       = MATRIX\b;

% plotting spanwise GAMMA distribution 
plotGAMMA(GAMMA,M,N);

% allocating GAMMA values in PANELwing class array
for i=1:N*2*M
    PANELwing(i).GAMMA = GAMMA(i);
end 

% computing LIFT
U            = 1;
rho          = 1;
S            = (root + root/taper) * L*cos(lambda/180*pi);
[~,L_vec,Cl] = compute_LIFT(GAMMA,PANELwing,lambda,delta,M,N,rho,U,S,"yes");

% computing induced velocity 
[v_ind,alpha_ind] = compute_INDUCEDvel(GAMMA,PANELwing,alpha,M,N,U,"yes"); 

% computing DRAG
[D,D_vec,Cd]      = compute_DRAG(L_vec,-alpha_ind,rho,U,S,M);

% computing Cl and Cd wrt alpha
%coeff_PLOT(MATRIX,PANELwing,beta,lambda,delta,AOA,M,N,S,linspace(-10,10,100),"yes");

toc

%% computing coefficients varying wing geometry 
clc
clear 
close all

tic

% varying taper ratio
beta   = 0;
delta  = 0;
lambda = 0;
root   = 8;
L      = 30;
AOA    = 0;

M = 25;
N = 10;

flag = "noplot";

TAPERvec  = 1:1:5;
alpha_vec = linspace(-10,10,30);

figure 

subplot(3,1,1)
hold on

for taper = TAPERvec
    
    % panel creation function 
    [PANELwing] = PANELING(delta,lambda,AOA,root,taper,L,M,N,flag,[0,0,0]);

    % system matrix generation
    % setting tollerance --> useful to avoid singular MATRIX 
    toll        = 1e-4;
    [MATRIX]    = BS(PANELwing,AOA,M,N,L,toll);
    
    % computing surface
    S               = (root + root/taper) * L*cos(lambda/180*pi);
    
    % compute Cl and Cd wrt alpha
    [Cl_vec,Cd_vec] = coeff_PLOT(MATRIX,PANELwing,beta,lambda,delta,AOA,M,N,S,alpha_vec,flag);
    
    % plotting results
    plot(Cd_vec,Cl_vec,'LineWidth',3);
    drawnow
    
end

grid on
grid minor

xlabel('$C_{D}$','Interpreter','latex');
ylabel('$C_{L}$','Interpreter','latex');
TEXT = "$\Lambda \ , \ \beta = " + string(beta) + "$"; 
title(TEXT,'Interpreter','latex');
TEXT = "$\Lambda = " + string(TAPERvec) + "$";
legend(TEXT,'Interpreter','latex');

% varying delta
lambda = 0;
root   = 8;
L      = 30;
taper  = 1;
AOA    = 0;

M = 25;
N = 10;

flag = "noplot";

subplot(3,1,2)
hold on

DELTAvec  = 0:1:5;
alpha_vec = linspace(-10,10,30);

for delta = DELTAvec
    
    % panel creation function 
    [PANELwing] = PANELING(delta,lambda,AOA,root,taper,L,M,N,flag,[0,0,0]);

    % system matrix generation
    % setting tollerance --> useful to avoid singular MATRIX 
    toll        = 1e-4;
    [MATRIX]    = BS(PANELwing,AOA,M,N,L,toll);
    
    % computing surface
    S               = (root + root/taper) * L*cos(lambda/180*pi);
    
    % compute Cl and Cd wrt alpha
    [Cl_vec,Cd_vec] = coeff_PLOT(MATRIX,PANELwing,beta,lambda,delta,AOA,M,N,S,alpha_vec,flag);
    
    % plotting results
    plot(Cd_vec,Cl_vec,'LineWidth',3);
    drawnow
    
end

grid on
grid minor

xlabel('$C_{D}$','Interpreter','latex');
ylabel('$C_{L}$','Interpreter','latex');
TEXT = "$\Delta \ , \ \beta = " + string(beta) + "$"; 
title(TEXT,'Interpreter','latex');
TEXT = "$\Delta = " + string(DELTAvec) + "$";
legend(TEXT,'Interpreter','latex');

% varying lambda
delta  = 0;
root   = 8;
L      = 30;
taper  = 1;
AOA    = 0;

M = 25;
N = 10;

flag = "noplot";

subplot(3,1,3)
hold on

LAMBDAvec = 0:5:30;
alpha_vec = linspace(-10,10,30);

for lambda = LAMBDAvec
    
    % panel creation function 
    [PANELwing] = PANELING(delta,lambda,AOA,root,taper,L,M,N,flag,[0,0,0]);

    % system matrix generation
    % setting tollerance --> useful to avoid singular MATRIX 
    toll        = 1e-4;
    [MATRIX]    = BS(PANELwing,AOA,M,N,L,toll);
    
    % computing surface
    S               = (root + root/taper) * L*cos(lambda/180*pi);
    
    % compute Cl and Cd wrt alpha
    [Cl_vec,Cd_vec] = coeff_PLOT(MATRIX,PANELwing,beta,lambda,delta,AOA,M,N,S,alpha_vec,flag);
    
    % plotting results
    plot(Cd_vec,Cl_vec,'LineWidth',3);
    drawnow
    
end

grid on
grid minor

xlabel('$C_{D}$','Interpreter','latex');
ylabel('$C_{L}$','Interpreter','latex');
TEXT = "$\lambda \ , \ \beta = " + string(beta) + "$"; 
title(TEXT,'Interpreter','latex');
TEXT = "$\lambda = " + string(LAMBDAvec) + "$";
legend(TEXT,'Interpreter','latex');

toc

%% wing interactions
close all
clear
clc

tic

% AERODYNAMIC ANGLES
    alpha   = 0;
    beta    = 0;
    
% 1ST WING GEOMETRY
    delta1  = 0;
    lambda1 = 0;
    root1   = 8;
    L1      = 30;
    taper1  = 1;
    AOA1    = 10;
    transl1 = [0,0,0]; 
    M1      = 15;
    N1      = 7;
    
% 2ND WING GEOMETRY
    delta2  = 0;
    lambda2 = 0;
    root2   = 4;
    L2      = 15;
    taper2  = 1;
    AOA2    = 10;
    transl2 = [25,0,0]; 
    M2      = 10;
    N2      = 5;

% toggling plotting
flag = "plot";

% panel creation function 
[PANELwing1] = PANELING(delta1,lambda1,AOA1,root1,taper1,L1,M1,N1,flag,transl1);
[PANELwing2] = PANELING(delta2,lambda2,AOA2,root2,taper2,L2,M2,N2,flag,transl2);

% assemblying MATRIX
tol          = 1e-4;
PANELwing    = [PANELwing1,PANELwing2];
[MATRIX]     = BS_multi(PANELwing,[AOA1,AOA2],[M1,M2],[N1,N2],(L1+L2)/2,tol);    

% assembling vector 
[b]          = compute_vector_multi(PANELwing,alpha,beta,[M1,M2],[N1,N2]);

% solving system
GAMMA        = MATRIX\b; 

% computing surfaces
S1           = (root1 + root1/taper1) * L1*cos(lambda1/180*pi);
S2           = (root2 + root2/taper2) * L2*cos(lambda2/180*pi);

% computing Cl vs alpha 
[Cl_vec1,Cd_vec1,Cl_vec2,Cd_vec2] = coeff_PLOT_multi(MATRIX,PANELwing,beta,[AOA1,AOA2],[lambda1,lambda2],[delta1,delta2],[M1,M2],[N1,N2],[S1,S2],linspace(-10,0,100),"yes");

%% computing coefficients varying tail AOA
close all
clear
clc

% AERODYNAMIC ANGLES
    alpha     = 0;
    beta      = 0;
    alpha_vec = linspace(-10,10,50);
    
% 1ST WING GEOMETRY
    delta1  = 0;
    lambda1 = 0;
    root1   = 8;
    L1      = 30;
    taper1  = 1;
    AOA1    = 10;
    transl1 = [0,0,0]; 
    M1      = 15;
    N1      = 7;
    
% 2ND WING GEOMETRY
    delta2  = 0;
    lambda2 = 0;
    root2   = 4;
    L2      = 15;
    taper2  = 1;
    AOA2vec = -10:5:10;
    transl2 = [25,0,0]; 
    M2      = 10;
    N2      = 5;

% toggling plotting
flag = "noplot";

figure 

for AOA2 = AOA2vec
   
    % panel creation function 
    [PANELwing1] = PANELING(delta1,lambda1,AOA1,root1,taper1,L1,M1,N1,flag,transl1);
    [PANELwing2] = PANELING(delta2,lambda2,AOA2,root2,taper2,L2,M2,N2,flag,transl2);

    % assemblying MATRIX
    tol          = 1e-4;
    PANELwing    = [PANELwing1,PANELwing2];
    [MATRIX]     = BS_multi(PANELwing,[AOA1,AOA2],[M1,M2],[N1,N2],(L1+L2)/2,tol);    

    % assembling vector 
    [b]          = compute_vector_multi(PANELwing,alpha,beta,[M1,M2],[N1,N2]);

    % solving system
    GAMMA        = MATRIX\b; 

    % computing surfaces
    S1           = (root1 + root1/taper1) * L1*cos(lambda1/180*pi);
    S2           = (root2 + root2/taper2) * L2*cos(lambda2/180*pi);

    % computing Cl vs alpha 
    [Cl_vec1,Cd_vec1,Cl_vec2,Cd_vec2] = coeff_PLOT_multi(MATRIX,PANELwing,beta,[AOA1,AOA2],[lambda1,lambda2],[delta1,delta2],[M1,M2],[N1,N2],[S1,S2],linspace(-10,0,100),"noplot");

    
    subplot(2,2,1)
    plot(AOA1 + alpha_vec,Cl_vec1,'LineWidth',3);
    hold on
    drawnow

    subplot(2,2,2)
    plot(Cd_vec1,Cl_vec1,'LineWidth',3);
    hold on
    drawnow

    subplot(2,2,3)
    plot(AOA2 + alpha_vec,Cl_vec2,'LineWidth',3);
    hold on
    drawnow

    subplot(2,2,4)
    plot(Cd_vec2,Cl_vec2,'LineWidth',3);
    hold on
    drawnow
    
end

TEXT = "$\alpha_{TAIL} = " + string(AOA2vec) + "$";

subplot(2,2,1)
legend(TEXT,'Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('WING','Interpreter','latex');
grid on 
grid minor

subplot(2,2,2)
legend(TEXT,'Interpreter','latex');
xlabel('$C_D$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('WING','Interpreter','latex');
grid on 
grid minor

subplot(2,2,3)
legend(TEXT,'Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('TAIL','Interpreter','latex');
grid on 
grid minor

subplot(2,2,4)
legend(TEXT,'Interpreter','latex');
xlabel('$C_D$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('TAIL','Interpreter','latex');
grid on 
grid minor

toc
