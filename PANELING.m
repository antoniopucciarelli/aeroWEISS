function [PANELwing,MIDwing,C4wing,VERTEXwing] = PANELING(delta,gamma,root,taper,L,M,N,flag)
% this function computes the panlization of the wing
% it allows delta wings and sweept angled wings
%
% INPUT: 
%   delta  = dihedral angle 
%   gamma  = sweep angle 
%   root   = root chord
%   taper  = wing taper ratio
%   L      = wing length --> it's the actual wing length not the projection
%            on the x axes
%   M      = # of points (-1) used to discretize the longitudinal lines 
%   N      = # of points (-1) used to discretize the horizontal lines 
%   flag   = describes if the fuction should print the goemetry
%
%
% OUTPUT:  
%   VERTEX = panel vertices matrix --> stores [X1,X2,X3,X4] vertices of the
%            panel 
%   MID    = matrix that stores the control of each panel   
%   C4     = matrix that stores the C/4 points that end at the edge of the
%            panel
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

% converting input in rad
delta = delta/180*pi;
gamma = gamma/180*pi;

% allocating extremes 
ext4 = [0   ,0];
ext3 = [root,0];

ext1 = ext4 + L * [ + sin(gamma), + cos(gamma)];
ext2 = ext1 + taper * root * [1,0];

% allocating points 
% --- split root and tip with N+1 points       --> N chordwise panels
% --- split longitudinal lines with M+1 points --> M tipwise panels

VERTEX = zeros(N+1,M+1,3);

deltax_root = root / N;
deltax_tip  = taper * root / N; 

for j = 0:N
        
    X1 =        [j*deltax_root, 0];
    X2 = ext1 + [j*deltax_tip,  0];
    
    % compute the vector that links the root with tape
    vec       = X2 - X1;
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

% rotating extremes
ext1 = R*[ext1,0]';
ext2 = R*[ext2,0]';
ext3 = R*[ext3,0]';
ext4 = R*[ext4,0]';
EXT = [ext4';ext3';ext2';ext1';ext4'];

% rotating points wrt the x axes
for j=1:N+1
    for i=1:M+1
        X = [VERTEX(j,i,1),VERTEX(j,i,2),VERTEX(j,i,3)];
        X = R*X';
        VERTEX(j,i,:) = X;
    end
end 

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

% computing control points and C4 points
% --- control points: 3/4 of the panel
%     --- extreme [1]> [X1,X4]
%     --- extreme [2]> [X2,X3]
% --- C4 points at the c/4 of the chord
%
MID = zeros(N*M,3);
C4 = zeros(N*M,2,3);

for j=1:N*M
    X1        = PANEL(j,1,:);
    X2        = PANEL(j,2,:);
    X3        = PANEL(j,3,:);
    X4        = PANEL(j,4,:);
    MID(j,:)  = ( (X1 + (X2 - X1)*3/4) + (X4 + (X3 - X4)*3/4) ) / 2;
    C4(j,1,:) = X4 + (X3 - X4)/4;
    C4(j,2,:) = X1 + (X2 - X1)/4;
end 

% generating the rest of the wing
LEFT_WING        = PANEL;
LEFT_WING(:,:,2) = - LEFT_WING(:,:,2);
RIGTH_WING       = PANEL;
PANELwing        = [LEFT_WING; RIGTH_WING];

% generating C4 points
LEFT_C4        = C4;
LEFT_C4(:,:,2) = - LEFT_C4(:,:,2);
C4wing         = [LEFT_C4; C4];

% generating MID points
LEFT_MID      = MID;
LEFT_MID(:,2) = - LEFT_MID(:,2);
MIDwing       = [LEFT_MID; MID];

% generating vertices 
LEFT_VERTEX        = VERTEX;
LEFT_VERTEX(:,:,2) = - LEFT_VERTEX(:,:,2); 
VERTEXwing         = [LEFT_VERTEX; VERTEX];

if(flag == "plot")
    
    figure(1)
    % extremes
    plot3(EXT(:,1),EXT(:,2),EXT(:,3),'-ok','LineWidth',3)
    hold on
    % internal points
    plot3(VERTEX(:,:,1),VERTEX(:,:,2),VERTEX(:,:,3),'or','LineWidth',4)
    axis('equal')

    figure(2)
    hold on
    for j=1:2*N*M
        P = [ PANELwing(j,1,:); 
              PANELwing(j,2,:);
              PANELwing(j,3,:);
              PANELwing(j,4,:);
              PANELwing(j,1,:) ];

        plot3(P(:,1),P(:,2),P(:,3),         'ok-','LineWidth',4);
        plot3(MIDwing(j,1),MIDwing(j,2),MIDwing(j,3),   'or', 'LineWidth',3);
        plot3(C4wing(j,1,1),C4wing(j,1,2),C4wing(j,1,3),'*b', 'LineWidth',3);
        plot3(C4wing(j,2,1),C4wing(j,2,2),C4wing(j,2,3),'oc', 'LineWidth',3);
    end 
    axis('equal')

    figure(3)
    surf(VERTEX(:,:,1),VERTEX(:,:,2),VERTEX(:,:,3));
    hold on 
    surf(LEFT_VERTEX(:,:,1),LEFT_VERTEX(:,:,2),LEFT_VERTEX(:,:,3));
    axis('equal')
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$z$','Interpreter','latex')
    
end 

end