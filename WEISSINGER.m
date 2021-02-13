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
%       rho                       = 1   [kg/m**3]
%
%
% !!! the program gives different Cd results wrt XFLR5. The main problem
% could hide behind the computation of the induced velocity. The whole
% procedure follows the aerodynamics' paper on the WEISSINGER method !!!
%
%% computing coefficients from initial conditions
close all
clear 
clc

% setting path
flpath = pwd;
addpath(append(flpath,'/src/'));

tic

% AERODYNAMIC properties
    alpha     = 0;
    beta      = 0;
% GEOMETRIC properties
    delta     = 0;
    lambda    = 0;
    root      = 8;
    L         = 30;
    taper     = 1;
    AOA       = 0;
% DISCRETIZATION properties
    M         = 10;
    N         = 5;
    alpha_vec = linspace(-10,10,30);

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

% computing LIFT
U                = 1;
rho              = 1;
S                = (root + root/taper) * L*cos(lambda/180*pi);
[~,L_vec,Cl,~,~] = compute_LIFT(GAMMA,PANELwing,lambda,M,N,rho,U,S,"yes");

% computing induced velocity 
[v_ind,alpha_ind] = compute_INDUCEDvel(GAMMA,PANELwing,AOA,alpha,M,N,U,"yes"); 

% computing DRAG
[D,D_vec,Cd]      = compute_DRAG(L_vec,-alpha_ind,alpha,rho,U,S,M);

% computing Cl and Cd wrt alpha
[Cl_vec,Cd_vec]   = coeff_PLOT(MATRIX,PANELwing,beta,lambda,AOA,M,N,S,alpha_vec,"yes");

toc

% removing path
flpath = pwd;
rmpath(append(flpath,'/src/'));

%% computing coefficients varying wing geometry 
close all
clear 
clc

% setting path
flpath = pwd;
addpath(append(flpath,'/src/'));

tic

% TAPER RATIO study
% GEOMETRY properties
    delta     = 0;
    lambda    = 0;
    root      = 8;
    L         = 30;
    AOA       = 0;
% AERODYNAMIC properties
    beta      = 0;
% DISCRETIZATION properties
    M         = 10;
    N         = 5;
    alpha_vec = linspace(0,5,30);
% TAPER RATIO discretization properties
    TAPERvec  = 1:1:5;

% plotting options
flag = "noplot";

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
    [Cl_vec,Cd_vec] = coeff_PLOT(MATRIX,PANELwing,beta,lambda,AOA,M,N,S,alpha_vec,flag);
    
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
TEXT = "$\lambda = " + string(TAPERvec) + "$";
legend(TEXT,'Interpreter','latex','location','best');

% DIHEDRAL study
% GEOMETRY properties
    lambda    = 0;
    root      = 8;
    L         = 30;
    taper     = 1;
    AOA       = 0;
% DIHEDRAL discretization properties
    DELTAvec  = 0:1:5;

% plotting options
flag = "noplot";

subplot(3,1,2)
hold on

for delta = DELTAvec
    
    % panel creation function 
    [PANELwing] = PANELING(delta,lambda,AOA,root,taper,L,M,N,flag,[0,0,0]);

    % system matrix generation
    % setting tollerance --> useful to avoid singular MATRIX 
    toll        = 1e-4;
    [MATRIX]    = BS(PANELwing,AOA,M,N,L,toll);
    
    % computing surface
    S = (root + root/taper) * L*cos(lambda/180*pi);
    
    % compute Cl and Cd wrt alpha
    [Cl_vec,Cd_vec] = coeff_PLOT(MATRIX,PANELwing,beta,lambda,AOA,M,N,S,alpha_vec,flag);
    
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
legend(TEXT,'Interpreter','latex','location','best');

% SWEEP ANGLE study
    delta     = 0;
    root      = 8;
    L         = 30;
    taper     = 1;
    AOA       = 0;
% DISCRETIZATION properties
    LAMBDAvec = 0:10:30;

% plotting options
flag = "noplot";

subplot(3,1,3)
hold on

for lambda = LAMBDAvec
    
    % panel creation function 
    [PANELwing] = PANELING(delta,lambda,AOA,root,taper,L,M,N,flag,[0,0,0]);

    % system matrix generation
    % setting tollerance --> useful to avoid singular MATRIX 
    toll        = 1e-4;
    [MATRIX]    = BS(PANELwing,AOA,M,N,L,toll);
    
    % computing surface
    S = (root + root/taper) * L*cos(lambda/180*pi);
    
    % compute Cl and Cd wrt alpha
    [Cl_vec,Cd_vec] = coeff_PLOT(MATRIX,PANELwing,beta,lambda,AOA,M,N,S,alpha_vec,flag);
    
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
TEXT = "$\Lambda = " + string(LAMBDAvec) + "$";
legend(TEXT,'Interpreter','latex','location','best');

toc

% removing path
flpath = pwd;
rmpath(append(flpath,'/src/'));

%% wing interactions
close all
clear
clc

% setting path
flpath = pwd;
addpath(append(flpath,'/src/'));

tic

% AERODYNAMIC ANGLES
    alpha     = 0;
    beta      = 0;
    alpha_vec = linspace(-10,10,30);
    
% 1ST WING GEOMETRY
    delta1    = 10;
    lambda1   = 5;
    root1     = 8;
    L1        = 30;
    taper1    = 1;
    AOA1      = 5;
    transl1   = [0,0,0]; 
    M1        = 7;
    N1        = 3;
    
% 2ND WING GEOMETRY
    delta2    = 15;
    lambda2   = 5;
    root2     = 4;
    L2        = 10;
    taper2    = 1;
    AOA2      = 5;
    transl2   = [35,0,0]; 
    M2        = 7;
    N2        = 3;

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

% plotting gamma distribution for the 2 wings
plotGAMMA(GAMMA(1:N1*2*M1)    ,M1,N1);
plotGAMMA(GAMMA(N1*2*M1+1:end),M2,N2);

% computing surfaces
S1           = (root1 + root1/taper1) * L1*cos(lambda1/180*pi);
S2           = (root2 + root2/taper2) * L2*cos(lambda2/180*pi);

% computing Cl vs alpha and Cl vd Cd
coeff_PLOT_multi(MATRIX,PANELwing,beta,[AOA1,AOA2],[lambda1,lambda2],[M1,M2],[N1,N2],[S1,S2],alpha_vec,"yes",70);

% AERODYNAMIC PROPERTIES FOR UNDISTURBED WINGS
% WING STUDY %
    flag = "noplot";

    % panel creation function 
    [PANELwing]       = PANELING(delta1,lambda1,AOA1,root1,taper1,L1,M1,N1,flag,[0,0,0]);

    % system matrix generation
    % setting tollerance --> useful to avoid singular MATRIX 
    toll              = 1e-4;
    [MATRIX]          = BS(PANELwing,AOA1,M1,N1,L1,toll);

    % computing Cl and Cd wrt alpha --> 1st wing
    [Cl_vec1,Cd_vec1] = coeff_PLOT(MATRIX,PANELwing,beta,lambda1,AOA1,M1,N1,S1,alpha_vec,"no");

% TAIL STUDY %
    flag = "noplot";

    % panel creation function 
    [PANELwing]       = PANELING(delta2,lambda2,AOA2,root2,taper2,L2,M2,N2,flag,[0,0,0]);

    % system matrix generation
    % setting tollerance --> useful to avoid singular MATRIX 
    toll              = 1e-4;
    [MATRIX]          = BS(PANELwing,AOA2,M2,N2,L2,toll);

    % computing Cl and Cd wrt alpha --> 2nd wing
    [Cl_vec2,Cd_vec2] = coeff_PLOT(MATRIX,PANELwing,beta,lambda2,AOA2,M2,N2,S2,alpha_vec,"no");

% PLOTTING DATA --> this is for comparing the different plots -- [WING+TAIL; WING; TAIL]
figure(70)
hold on

subplot(2,2,1)
hold on
plot(AOA1+alpha_vec,Cl_vec1,'or','LineWidth',5);
legend("WING + TAIL","WING ALONE",'location','best')

subplot(2,2,2)
hold on
plot(Cd_vec1,Cl_vec1,'or','LineWidth',5);
legend("WING + TAIL","WING ALONE",'location','best')

subplot(2,2,3)
hold on
plot(AOA2+alpha_vec,Cl_vec2,'or','LineWidth',5);
legend("WING + TAIL","TAIL ALONE",'location','best')

subplot(2,2,4)
hold on
plot(Cd_vec2,Cl_vec2,'or','LineWidth',5);
legend("WING + TAIL","TAIN ALONE",'location','best')

toc

% removing path
flpath = pwd;
rmpath(append(flpath,'/src/'));

%% computing coefficients varying tail AOA
close all
clear
clc

% setting path
flpath = pwd;
addpath(append(flpath,'/src/'));

tic

% AERODYNAMIC ANGLES
    alpha     = 0;
    beta      = 0;
    alpha_vec = 0;
    
% 1ST WING GEOMETRY
    delta1  = 0;
    lambda1 = 0;
    root1   = 8;
    L1      = 30;
    taper1  = 1;
    AOA1    = 0;
    transl1 = [0,0,0]; 
    M1      = 7;
    N1      = 3;
    
% 2ND WING GEOMETRY
    delta2  = 0;
    lambda2 = 0;
    root2   = 4;
    L2      = 15;
    taper2  = 1;
    AOA2vec = -10:0.2:10;
    transl2 = [25,0,0]; 
    M2      = 5;
    N2      = 2;

% toggling plotting
flag = "noplot";

figure(70) 

% the 1st WING isn't inside the for loop because its geometry doesn't vary
% with the angle of attack of the TAIL
[PANELwing1] = PANELING(delta1,lambda1,AOA1,root1,taper1,L1,M1,N1,flag,transl1);
    
for AOA2 = AOA2vec
   
    % panel creation function 
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
    [Cl_vec1,Cd_vec1,Cl_vec2,Cd_vec2] = coeff_PLOT_multi(MATRIX,PANELwing,beta,[AOA1,AOA2],[lambda1,lambda2],[M1,M2],[N1,N2],[S1,S2],alpha_vec,"noplot");

    subplot(2,2,1)
    plot(AOA1 + alpha_vec,Cl_vec1,'ok','LineWidth',3);
    drawnow
    hold on

    subplot(2,2,2)
    plot(Cd_vec1,Cl_vec1,'ok','LineWidth',3);
    drawnow
    hold on

    subplot(2,2,3)
    plot(AOA2 + alpha_vec,Cl_vec2,'ok','LineWidth',3);
    drawnow
    hold on

    subplot(2,2,4)
    plot(Cd_vec2,Cl_vec2,'ok','LineWidth',3);
    drawnow
    hold on
    
end

TEXT = "$\alpha_{TAIL} = " + string(AOA2vec) + "$";

subplot(2,2,1)
%legend(TEXT,'Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('$WING \Rightarrow \ \alpha_{TAIL}$','Interpreter','latex');
grid on 
grid minor

subplot(2,2,2)
%legend(TEXT,'Interpreter','latex');
xlabel('$C_D$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('$WING \Rightarrow \ \alpha_{TAIL}$','Interpreter','latex');
grid on 
grid minor

subplot(2,2,3)
%legend(TEXT,'Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('TAIL','Interpreter','latex');
grid on 
grid minor

subplot(2,2,4)
%legend(TEXT,'Interpreter','latex');
xlabel('$C_D$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('TAIL','Interpreter','latex');
grid on 
grid minor

toc

% removing path
flpath = pwd;
rmpath(append(flpath,'/src/'));

%% computing ground effect
close all
clear
clc

% setting path
flpath = pwd;
addpath(append(flpath,'/src/'));

tic

% AERODYNAMIC ANGLES
    alpha   = 0;
    beta    = 0;
    AOA_vec = -1:0.05:1;
    
% 1ST WING GEOMETRY
    delta1  = 0;
    lambda1 = 0;
    root1   = 8;
    L1      = 30;
    taper1  = 1;
    AOA1    = 10;
    transl1 = [0,0,5]; 
    M1      = 7;
    N1      = 3;
    
% 2ND WING GEOMETRY
    delta2  = -delta1;
    lambda2 = lambda1;
    root2   = root1;
    L2      = L1;
    taper2  = taper1;
    transl2 = - transl1; 
    M2      = M1;
    N2      = N1;

% toggling plotting
flag = "noplot";

for AOA = AOA_vec

    [PANELwing1] = PANELING(delta1,lambda1, AOA,root1,taper1,L1,M1,N1,flag,transl1);
    [PANELwing2] = PANELING(delta2,lambda2,-AOA,root2,taper2,L2,M2,N2,flag,transl2);
    
    % assemblying MATRIX
    tol          = 1e-4;
    PANELwing    = [PANELwing1,PANELwing2];
    [MATRIX]     = BS_multi(PANELwing,[AOA,-AOA],[M1,M2],[N1,N2],transl1(3)*3/4,tol);  
    
    % assembling vector 
    [b]          = compute_vector_multi(PANELwing,alpha,beta,[M1,M2],[N1,N2]);
    
    % computing surfaces
    S1           = (root1 + root1/taper1) * L1*cos(lambda1/180*pi);
    S2           = (root2 + root2/taper2) * L2*cos(lambda2/180*pi);

    % computing Cl vs alpha 
    [Cl_vec1,Cd_vec1,Cl_vec2,Cd_vec2] = coeff_PLOT_multi(MATRIX,PANELwing,beta,[AOA,-AOA],[lambda1,lambda2],[M1,M2],[N1,N2],[S1,S2],0,"noplot");
    
    figure(70)
    
    subplot(2,1,1)
    GE1 = plot(AOA,Cl_vec1,'ok','LineWidth',3);
    drawnow
    hold on

    subplot(2,1,2)
    GE2 = plot(Cd_vec1,Cl_vec1,'ok','LineWidth',3);
    drawnow
    hold on

end

% AERODYNAMIC PROPERTIES FOR UNDISTURBED WINGS
% WING STUDY %
    flag = "noplot";

    % panel creation function 
    [PANELwing]       = PANELING(delta1,lambda1,0,root1,taper1,L1,M1,N1,flag,[0,0,0]);

    % system matrix generation
    % setting tollerance --> useful to avoid singular MATRIX 
    toll              = 1e-4;
    [MATRIX]          = BS(PANELwing,0,M1,N1,L1,toll);

    % computing Cl and Cd wrt alpha --> 1st wing
    [Cl_vec1,Cd_vec1] = coeff_PLOT(MATRIX,PANELwing,beta,lambda1,0,M1,N1,S1,AOA_vec,"no");

% adding UNDISTURBED WING data    
subplot(2,1,1)
hold on
P1 = plot(AOA_vec,Cl_vec1,'-r','LineWidth',3);
xlabel('$\alpha$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('$WING_{real}$','Interpreter','latex');
legend([GE1,P1],'$GROUND \ EFFECT$','$UNDISTURBED$','Interpreter','latex','location','best');
grid on 
grid minor

subplot(2,1,2)
hold on
P2 = plot(Cd_vec1,Cl_vec1,'-r','LineWidth',3);
xlabel('$C_D$','Interpreter','latex');
ylabel('$C_L$','Interpreter','latex');
title('$WING_{real}$','Interpreter','latex');
legend([GE2,P2],'$GROUND \ EFFECT$','$UNDISTURBED$','Interpreter','latex','location','best');
grid on 
grid minor  

toc 

% removing path
flpath = pwd;
rmpath(append(flpath,'/src/'));
    