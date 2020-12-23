function [b] = compute_vector_multi(PANELwing,alpha,beta,M,N)
% this function compute the known vector of the system
%
% INPUT: 
%   PANELwing : PANEL class array            [PANEL class]
%   alpha     : AOA of the wing              [deg]
%   beta      : side slip angle of the wing  [deg]
%   M         : longitudinal discretization  [M1,M2]
%   N         : horizontal discretization    [N1,N2]
%
% OUTPUT:
%   b         : known vector --> describes the non penetration condition for
%               the know velocity at infinity
%
% PANEL DESCRIPTION BY PANELING FUNCTION 
%
%     X1  C42        X2
%     o---o----------o  
%     |   |          |
%     |   |          |
%     |   |   MID    |
%     |   |    o     | 
%     |   |          |    
%     |   |          |
%     o---o----------o
%     X4  C41       X3
%
%  | y -- spanwise
%  |
%  |_____ x -- chordwise
%

alpha = alpha/180*pi;
beta  = beta/180*pi;

vel = [cos(beta)*cos(alpha); - sin(beta); cos(beta)*sin(alpha)];

VECsize = N(1)*2*M(1) + N(2)*2*M(2);

b = ones(VECsize,1);

for i=1:VECsize
    b(i) = - PANELwing(i).normal' * vel;
end 

end