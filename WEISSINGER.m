% this program computes the velocity around a finite 3-D wing
% the wing is invested by an external flow of V_inf = 1
% the wing can have multiple displacement wrt the airstream
% the wing position is described by aerodynamic angles
%
% INPUT: 
%   WING properties:
%       sweep angle            -- gamma [deg]
%       dihedral angle         -- delta [deg]
%       root chord             -- root
%       semi-wing length       -- L
%       taper ratio            -- taper
%       # of spanwise panels   -- M
%       # of chordwise panels  -- N
%
%   AIRSTREAM properties:
%       AOA            -- alpha [deg]
%       sideslip angle -- beta  [deg]
%       U                 = 1   [m/s]
%

clc
clear 
close all

alpha = 5;
beta  = 0;
delta = 0;
gamma = 0;
root  = 8;
L     = 30;
taper = 0.7;

M = 20;
N = 10;

flag = "plot";

% panel creation function 
[PANELwing] = PANELING(delta,gamma,root,taper,L,M,N,flag);

% system matrix generation
[MATRIX]    = BS(PANELwing,M,N,L);

% system known vector 
[b]         = compute_vector(PANELwing,alpha,beta,M,N);

% solve system
GAMMA       = MATRIX\b;

GAMMA