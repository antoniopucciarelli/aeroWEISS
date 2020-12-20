%function [VERTEX,MID,C4] = PANELING(delta,gamma,root,taper,L,M,N)
% this function computes the panlization of the wing
% it allows delta wings and sweept angled wings
%
% INPUT: 
%   delta  = diedre angle 
%   gamma  = sweep angle 
%   root   = root chord
%   taper  = wing taper ratio
%   L      = wing length --> it's the actual wing length not the projection
%            on the x axes
%   M      = # of points (-1) used to discretize the longitudinal lines 
%   N      = # of points (-1) used to discretize the horizontal lines 
%
% OUTPUT:  
%   VERTEX = panel vertices matrix --> stores [X1,X2,X3,X4] vertices of the
%            panel 
%   MID    = matrix that stores the control of each panel   
%   C4     = matrix that stores the C/4 points that end at the edge of the
%            panel
%   
%     X1            X2
%     o-------------o  
%     |             |
% C41 o-------------o C42 -- c/4
%     |             |
%     |     MID     | 
%     |      o      |     -- 3c/4
%     |             |
%     o-------------o
%     X4            X3
%
%  | y -- spanwise
%  |
%  |_____ x -- chordwise
%
%        

clear 
close all
clc

% allocating dimensions
root  = 5;
gamma = 5/180*pi;
delta = 2/180*pi;
L     = 25;
taper = 0.5;
M     = 10;
N     = 5;

% allocating extremes 
ext1 = [0   ,0];
ext2 = [root,0];

ext3 = ext2 + L * [ - sin(gamma), + cos(gamma)];
ext4 = ext3 - taper * root * [1,0];

EXT = [ext1;ext2;ext3;ext4;ext1];

% allocating points 
% --- split root and tip with N+1 points       --> N chordwise panels
% --- split longitudinal lines with M+1 points --> M tipwise panels

VERTEX = zeros(N+1,M+1,3);

deltax_root = root / N;
deltax_tip  = taper * root / N; 

for j = 0:N
        
    X1 =        [j*deltax_root, 0];
    X2 = ext4 + [j*deltax_tip,  0];
    
    % compute the vector that links the root with tape
    vec = X2 - X1;
    % dividing vec into N points
    vec_l     = norm(vec);
    vec       = vec / vec_l;
    delta_vec = vec_l / M;

    for i = 0:M

        VERTEX(j+1,i+1,1:2) = X1 + delta_vec * i * vec;

    end 
    
end

% incline points procedure
R = ROT(0,0,delta,'rad');
% rotating points wrt the x axes
for j=1:N+1
    for i=1:M+1
        X = [VERTEX(j,i,1),VERTEX(j,i,2),VERTEX(j,i,3)];
        X = R*X';
        VERTEX(j,i,:) = X;
    end
end 

figure(1)
plot3(EXT(:,1),EXT(:,2),zeros(5),'-ok','LineWidth',3)
hold on
plot3(VERTEX(:,:,1),VERTEX(:,:,2),VERTEX(:,:,3),'or','LineWidth',4)
axis('equal')

% generate panels 
PANEL = zeros(N*M,4,3);

for j=0:N-1
    for i=1:M
        PANEL(i+j*M,1,:) = VERTEX(j+1,i+1,:);
        PANEL(i+j*M,2,:) = VERTEX(j+2,i+1,:);
        PANEL(i+j*M,3,:) = VERTEX(j+2,i,:);
        PANEL(i+j*M,4,:) = VERTEX(j+1,i,:);
    end 
end

% computing control points
% --- 3/4 of the panel
%     --- extreme [1]> [X1,X4]
%     --- extreme [2]> [X2,X3]
MID = zeros(N*M,3);

for j=1:N*M
    X1       = PANEL(j,1,:);
    X2       = PANEL(j,2,:);
    X3       = PANEL(j,3,:);
    X4       = PANEL(j,4,:);
    MID(j,:) = ( X1 - (X1 - X2)*3/4 + X4 - (X4 - X3)*3/4 ) / 2;
end 

% computing C4 points
C4 = zeros(N*M,2,3);

for j=1:N*M
    X1        = PANEL(j,1,:);
    X2        = PANEL(j,2,:);
    X3        = PANEL(j,3,:);
    X4        = PANEL(j,4,:);
    C4(j,1,:) = X1 - (X1 - X2)/4;
    C4(j,2,:) = X4 - (X4 - X3)/4;
end 



figure(2)
hold on
for j=1:N*M
    P = [ PANEL(j,1,:); 
          PANEL(j,2,:);
          PANEL(j,3,:);
          PANEL(j,4,:);
          PANEL(j,1,:)];
      
    plot3(P(:,1),P(:,2),P(:,3),'ok-','LineWidth',4);
    plot3(MID(j,1),MID(j,2),MID(j,3),'or','LineWidth',3);
    plot3(C4(j,1,1),C4(j,1,2),C4(j,1,3),'*b','LineWidth',3);
    plot3(C4(j,2,1),C4(j,2,2),C4(j,2,3),'oc','LineWidth',3);
end 
axis('equal')



%end /