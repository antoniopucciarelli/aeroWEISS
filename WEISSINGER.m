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
%       airflow incidence -- alpha [deg]
%       sideslip angle    -- beta  [deg]
%       U                    = 1   [m/s]
%

clc
clear 
close all

alpha  = 0;
beta   = 0;
delta  = 0;
lambda = 0;
root   = 8;
L      = 30;
taper  = 1;
AOA    = 10;

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

for i=1:N*2*M
    PANELwing(i).GAMMA = GAMMA(i);
end 

% computing LIFT
U            = 1;
rho          = 1;
S            = (root + root/taper) * L*cos(lambda/180*pi);
[L,L_vec,Cl] = compute_LIFT(GAMMA,PANELwing,lambda,M,N,rho,U,S,"yes");

% computing Cl vs alpha
coeff_PLOT(MATRIX,PANELwing,beta,lambda,M,N,S,linspace(-10,10,100));

% computing induced velocity 
[v_ind,alpha_ind] = compute_INDUCEDvel(GAMMA,PANELwing,M,N,U,"yes"); 

% computing DRAG
[D,D_vec,Cd]      = compute_DRAG(L_vec,-alpha_ind,rho,U,S,M);

% computing aerodynamic coeff


% doing some tests on different configurations

%% wing interactions
close all
clear
clc

% AERODYNAMIC COEFFICIENTS
    alpha  = 0;
    beta   = 0;
    
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

flag = "plot";

% panel creation function 
[PANELwing1] = PANELING(delta1,lambda1,AOA1,root1,taper1,L1,M1,N1,flag,transl1);
[PANELwing2] = PANELING(delta2,lambda2,AOA2,root2,taper2,L2,M2,N2,flag,transl2);

% assemblying MATRIX
tol       = 1e-4;
PANELwing = [PANELwing1,PANELwing2];
[MATRIX]  = BS_multi(PANELwing,[AOA1,AOA2],[M1,M2],[N1,N2],[L1+L2]/2,tol);    

% assembling vector 
[b]       = compute_vector_multi(PANELwing,alpha,beta,[M1,M2],[N1,N2]);

% solving system
GAMMA     = MATRIX\b; 

% computing surfaces
S1            = (root1 + root1/taper1) * L1*cos(lambda1/180*pi);
S2            = (root2 + root2/taper2) * L2*cos(lambda2/180*pi);

% computing Cl vs alpha 
[Cl_vec1,Cd_vec1,Cl_vec2,Cd_vec2] = coeff_PLOT_multi(MATRIX,PANELwing,beta,[lambda1,lambda2],[M1,M2],[N1,N2],[S1,S2],linspace(-10,10,100));




