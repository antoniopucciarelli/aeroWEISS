function [L,L_vec,Cl] = compute_LIFT(GAMMA,PANELwing,M,N,rho,U,S)
% this function compute the total circulation of the wing 
%
% INPUT:
%   GAMMA     : filament circulation array
%   PANELwing : PANEL class array 
%   M         : # of discretization points in the spanwise direction
%   N         : # of discretization points in the chordwise direction
%   alpha     : AOA
%   beta      : sideslip angle
%
% OUTPUT:
%   L     : 3D wing total lift
%   L_vec : 3D wing lift distribution spanwise
%

tic

% initializing values
L_vec = zeros(2*M,1);

for i=1:2*M
    % computing lift distribution spanwise
    for j=1:N
        L_vec(i) = L_vec(i) + rho * U * GAMMA(i+(j-1)*2*M) * norm(PANELwing(i+(j-1)*2*M).C4(1,:) - PANELwing(i+(j-1)*2*M).C4(2,:));
    end 
end 

% summing lift distribution spanwise
L  = sum(L_vec);

% computing Cl 
Cl = L/(0.5 * rho * U^2 * S); 

toc

end