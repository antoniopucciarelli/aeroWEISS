function [v_ind,alpha_ind] = compute_INDUCEDvel(GAMMA,PANEL,M,N,U)
% this function computes the induced velocity of a circulation
% discribution GAMMA over a finite 3D wing described by PANEL
%
% INPUT: 
%   GAMMA : vortex filamente circulation
%   PANEL : PANEL class array 
%   M     : # of discretization panels spanwise
%   N     : # of discretization panels chordwise
%
%

tic

% induced velocity vector
v_ind = zeros(N*2*M,1);

for i=1:2*M*N
    % !!!PAY ATTENTION!!! the contribute of induced velocity is computed 
    % only through the lateral filament of each panel that induced a 
    % vertical velocity at the midpoint of the studied panel 
    for j=1:N*2*M
        
        v_ind1 = + 1/(4*pi) * GAMMA(j) / (PANEL(i).MIDPOINT(2) - PANEL(j).C4(2,2));
        
        v_ind2 = - 1/(4*pi) * GAMMA(j) / (PANEL(i).MIDPOINT(2) - PANEL(j).C4(1,2));
        
        v_ind(i) = v_ind(i) + (v_ind1 + v_ind2);
     
    end 
end

% initializing alpha_ind vec
alpha_ind = zeros(2*M,1);

for i=1:2*M
    for j=1:N
        alpha_ind(i) = alpha_ind(i) + v_ind(i+(j-1)*2*M);
    end
    
    alpha_ind(i) = atan(alpha_ind(i)/U);
    
end 

toc

end 