function [v_ind] = compute_INDUCEDvel(GAMMA,PANEL,M,N)
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

% induced velocity vector
v_ind = zeros(N*2*M,1);

for i=1:2*M
    % !!!PAY ATTENTION!!! the contribute of induced velocity is computed 
    % only through the lateral filament of each panel that induced a 
    % vertical velocity at the midpoint of the studied panel 
    for j=1:N
        v_ind(i) = v_ind(i) + 1/(4*pi)*;
    end 
end 

end 