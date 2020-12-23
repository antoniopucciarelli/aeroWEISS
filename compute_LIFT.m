function [L,L_vec,Cl] = compute_LIFT(GAMMA,PANELwing,lambda,M,N,rho,U,S,flag)
% this function compute the total circulation of the wing 
%
% INPUT:
%   GAMMA     : filament circulation array
%   PANELwing : PANEL class array 
%   M         : # of discretization points in the spanwise direction
%   N         : # of discretization points in the chordwise direction
%   lambda    : sweep angle 
%   alpha     : AOA
%   beta      : sideslip angle
%
% OUTPUT:
%   L     : 3D wing total lift
%   L_vec : 3D wing lift distribution spanwise
%

tic

% initializing values
L_vec  = zeros(2*M,1);
lambda = lambda/180*pi;

for i=1:2*M
    % computing lift distribution spanwise
    for j=1:N
        L_vec(i) = L_vec(i) + rho * U * cos(lambda) * GAMMA(i+(j-1)*M) * norm(PANELwing(i+(j-1)*M).C4(1,:) - PANELwing(i+(j-1)*M).C4(2,:));
    end 
end 

% summing lift distribution spanwise
L  = sum(L_vec);

% computing Cl 
Cl = L/(0.5 * rho * U^2 * S); 

if(flag == "yes")
    
    K = zeros(length(L_vec),1);
    
    for i=1:M
       K(i) = L_vec(M+1-i); 
    end
    
    for i=1:M
       K(M+i) = L_vec(M+i);
    end
    
    figure
    x = linspace(-1,1,length(L_vec));
    plot(x,K,'ok-','LineWidth',2.5);
    grid on
    grid minor
    xlabel('$\% SPAN$','Interpreter','latex');
    ylabel('$L_{(i)}$','Interpreter','latex');
    
end

toc

end