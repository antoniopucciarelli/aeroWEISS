%function [] = BS()
% this function computes the induced velocity in a point by the panel of
% unit circulation value with the BIOT SAVART law.
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
%
% BIOT-SAVART law
%   
%  R2 o beta2            
%     | *   ^       
%     |   * | Vc     
%   R |     o 
%     |   * C
%     | *         
%  R1 o beta1
%

clc
close all
clear
 
[PANELwing,MIDwing,C4wing,VERTEXwing] = PANELING(2,7,7.5,0.6,25,25,5,"noplot");

size(PANELwing)
% 1
% PANELwing(1,1,:)
% 2
% PANELwing(1,2,:)
% 3
% PANELwing(1,3,:)
% 4
% PANELwing(1,4,:)

r_c = [1,1,1]
L = 25;

NORM = @(x) sqrt(x(1)^2 + x(2)^2 + x(3)^2);
DOT  = @(x11,x12,x13,x21,x22,x23) x11*x21 + x12*x22 + x13*x23;

% computing C/4 line induction 
% setting up variables
C41(:) = C4wing(1,1,:);
C42(:) = C4wing(1,2,:);

r0_vec(:) = C42 - C41;
r1_vec(:) = r_c - C41;
r2_vec(:) = r_c - C42;
r1        = NORM(r1_vec);
r2        = NORM(r2_vec);
pr_vec    = cross(r1_vec,r2_vec);
pr        = NORM(pr_vec);

r0_vec
r1_vec
r2_vec
pr_vec
r1
r2
pr

A = r1_vec /r1 - r2_vec/r2;
B = pr_vec/pr^2;

A
B

C = DOT(A(1),A(2),A(3),B(1),B(2),B(3))

Vc_c4  = 1/(4*pi) * r0_vec * C;

% % computing lateral line induction -- UPPER line from C42 to infinity 
% X1(:) = PANELwing(1,1,:);
% X2(:) = PANELwing(1,2,:);
% X3(:) = PANELwing(1,3,:);
% X4(:) = PANELwing(1,4,:);
% 
% r0_vec(:) = (X2 - X1) * L * 100;
% r1_vec(:) = r_c - C42;
% r2_vec(:) = r_c - (C42 + r0_vec);
% r1        = NORM(r1_vec);
% r2        = NORM(r2_vec);
% pr_vec    = cross(r1_vec,r2_vec);
% pr        = NORM(pr_vec);
% 
% Vc_inf2 = 1/(4*pi) * r0_vec * dot((r1_vec/r1 - r2_vec/r2),(pr_vec / pr^2));
% 
% % computing lateral line induction -- LOWER line from C41 to infinity
% r0_vec = (X3 - X4) * L * 100;
% r1_vec = r_c - C41;
% r2_vec = r_c - (C41 + r0_vec);
% r1     = NORM(r1_vec);
% r2     = NORM(r2_vec);
% pr_vec = cross(r1_vec,r2_vec);
% pr     = NORM(pr_vec);
% 
% Vc_inf1 = 1/(4*pi) * r0_vec * dot((r1_vec/r1 - r2_vec/r2),(pr_vec / pr^2));
% 
% % computing velocity at control point 
% VC = Vc_c4 + Vc_inf1 + Vc_inf2;




% end