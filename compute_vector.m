function [b] = compute_vector(PANELwing,alpha,beta,M,N)
% this function compute the known vector of the system
%
% INPUT: 
%   PANELwing : PANEL class array            [PANEL class]
%   alpha     : AOA of the wing              [deg]
%   beta      : side slip angle of the wing  [deg]
%   M         : longitudinal discretization
%   N         : horizontal discretization
%
% OUTPUT:
%   b         : known vector --> describes the non penetration condition for
%               the know velocity at infinity
%
%

tic

alpha = alpha/180*pi;
beta  = beta/180*pi;

PANEL = PANELwing;

vel = [cos(alpha),sin(beta),sin(alpha)];

b = ones(N*2*M,1);

for i=1:N*2*M
    b(i) = PANEL(i).normal' * vel';
end 

toc

end 