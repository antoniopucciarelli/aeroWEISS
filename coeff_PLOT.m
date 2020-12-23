function [Cl_vec,Cd_vec] = coeff_PLOT(MATRIX,PANELwing,beta,lambda,M,N,S,alpha_vec,flag)
% this function computes the Cl alpha plot of a 3D wing given geometry and
% flow conditions
% 
% INPUT: 
%   MATRIX    : system matrix -- describes the non penetration condition
%   PANELwing : PANEL class array
%   beta      : sideslip angle
%   M         : spanwise # of discretization points
%   N         : chordwise # of discretization points
%   S         : wing surface
%   alpha_vec : set AOA
%

% initializing variables
Cl_vec = zeros(size(alpha_vec));
Cd_vec = zeros(size(alpha_vec));

for i=1:length(alpha_vec)
        
    alpha = alpha_vec(i);
    
    % system known vector 
    [b]         = compute_vector(PANELwing,alpha,beta,M,N);

    % solve system
    GAMMA       = MATRIX\b;

    for j=1:N*2*M
        PANELwing(j).GAMMA = GAMMA(j);
    end 

    % computing LIFT
    U               = 1;
    rho             = 1;
    [~,L_vec,Cl_vec(i)] = compute_LIFT(GAMMA,PANELwing,lambda,M,N,rho,U,S,"no");
    
    % computing induced velocity 
    [~,alpha_ind]   = compute_INDUCEDvel(GAMMA,PANELwing,M,N,U,"no"); 

    % computing DRAG
    [~,~,Cd_vec(i)] = compute_DRAG(L_vec,-alpha_ind,rho,U,S,M);
    
end 

if(flag == "yes")
    % plotting procedure
    figure

    subplot(2,1,1)
    plot(alpha_vec,Cl_vec,'k-','LineWidth',2);
    xlabel('$\alpha$','Interpreter','latex');
    ylabel('$C_L$','Interpreter','latex');
    grid on 
    grid minor

    subplot(2,1,2)
    plot(Cd_vec,Cl_vec,'k-','LineWidth',2);
    xlabel('$C_D$','Interpreter','latex');
    ylabel('$C_L$','Interpreter','latex');
    grid on 
    grid minor
end
end