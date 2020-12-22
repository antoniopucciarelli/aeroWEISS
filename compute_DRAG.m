function [D,D_vec] = compute_DRAG(L_vec,alpha_ind,M)
% This function computes the induced drag of a 3D wing 
% After have computed the distribution of lift over the wing and the
% induced velocity/angle, one can compute the variation of the LIFT vector
% that produce the induced drag of the wing.
% The induced DRAG is a kind of resistence that is related to the fact that
% the vortex' filaments are horse-shoe vortices. Their curlig induced a
% downward velocity that alters the flow direction from infinity. So this,
% the KUTTA-JOUKOWSKY theorem predicts a component of LIFT that is
% essentially a DRAG force. 
% It is interesting to understand that this kind of DRAG is produced
% primarly by a difference of pressure over the 3D wing but it is not a
% pressure DRAG because the flow is not deatached.
%
% INPUT: 
%   L_vec     : spanwise distribution of lift
%   alpha_ind : induced angle of attack deriving from the decomposition  
%               of the spanwise distribution of induced velocity
%   M         : # of spanwise discretization points
%
% OUTPUT:
%   D         : total induced DRAG   
%   D_vec     : spanwise ditribution of induced DRAG 
%

tic

% initializing D_vec
D_vec = zeros(2*M,1);

for i=1:2*M
    D_vec(i) = L_vec(i) * sin(alpha_ind(i));
end

D = sum(D_vec);

toc

end 