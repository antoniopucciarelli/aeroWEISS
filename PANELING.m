function [PANELwing] = PANELING(delta,gamma,root,taper,L,M,N,flag)
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
%   PANELwing = PANEL class array --> it stores all the necessary info 
%               to compute the solution  
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

tic

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
EXT  = [ext4';ext3';ext2';ext1';ext4'];

% rotating points wrt the x axes
for j=1:N+1
    for i=1:M+1
        X = [VERTEX(j,i,1),VERTEX(j,i,2),VERTEX(j,i,3)];
        X = R*X';
        VERTEX(j,i,:) = X;
    end
end 

% generate PANEL array
PANELarray = PANEL.empty(N*M,0);
% allocating PANEL array class variables
for j=0:N-1
    for i=1:M
        PANELarray(i+j*M).VERTEX(1,:) = VERTEX(j+1,i+1,:);    
        PANELarray(i+j*M).VERTEX(2,:) = VERTEX(j+2,i+1,:);
        PANELarray(i+j*M).VERTEX(3,:) = VERTEX(j+2,i,:);        
        PANELarray(i+j*M).VERTEX(4,:) = VERTEX(j+1,i,:);
    end 
end

% computing control points and C4 points
% --- control points: 3/4 of the panel
%     --- extreme [1]> [VERTEX_1,VERTEX_4]
%     --- extreme [2]> [VERTEX_2,VERTEX_3]
% --- C4 points at the c/4 of the chord
% --- normal vector
normal_vec = R*[0;0;1];

for j=1:N*M
    X1                     = PANELarray(j).VERTEX(1,:);
    X2                     = PANELarray(j).VERTEX(2,:);
    X3                     = PANELarray(j).VERTEX(3,:);
    X4                     = PANELarray(j).VERTEX(4,:);
    PANELarray(j).MIDPOINT = ( (X1 + (X2 - X1)*3/4) + (X4 + (X3 - X4)*3/4) ) / 2; 
    PANELarray(j).C4(1,:)  = X4 + (X3 - X4)/4;
    PANELarray(j).C4(2,:)  = X1 + (X2 - X1)/4;
    PANELarray(j).normal   = normal_vec;
end 

% generating the rest of the wing
LEFT_WING = PANELarray;
for i=1:N*M
    LEFT_WING(i).VERTEX(:,2) = - LEFT_WING(i).VERTEX(:,2);
    LEFT_WING(i).MIDPOINT(2) = - LEFT_WING(i).MIDPOINT(2);
    LEFT_WING(i).C4(:,2)     = - LEFT_WING(i).C4(:,2);
    LEFT_WING(i).normal(2)   = - LEFT_WING(i).normal(2);
    
    temp                     = LEFT_WING(i).VERTEX(1,:);
    LEFT_WING(i).VERTEX(1,:) = LEFT_WING(i).VERTEX(4,:);
    LEFT_WING(i).VERTEX(4,:) = temp;
    temp                     = LEFT_WING(i).VERTEX(3,:);
    LEFT_WING(i).VERTEX(3,:) = LEFT_WING(i).VERTEX(2,:);
    LEFT_WING(i).VERTEX(2,:) = temp;
    
    temp                 = LEFT_WING(i).C4(1,:);
    LEFT_WING(i).C4(1,:) = LEFT_WING(i).C4(2,:);
    LEFT_WING(i).C4(2,:) = temp;
end 

% allocating PANELwing
PANELwing = [LEFT_WING, PANELarray];

% plotting process
if(flag == "plot")
    
    figure(1)
    % extremes
    plot3(EXT(:,1),EXT(:,2),EXT(:,3),'-ok','LineWidth',3)
    hold on
    % internal points
    for i=1:N*M
        plot3(PANELarray(i).VERTEX(:,1), ...
              PANELarray(i).VERTEX(:,2), ...
              PANELarray(i).VERTEX(:,3),'or','LineWidth',4)
    end 
    axis('equal')
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$z$','Interpreter','latex')

    figure(2)
    hold on
    for j=1:2*N*M
        plot3(PANELwing(j).MIDPOINT(1), ... 
              PANELwing(j).MIDPOINT(2), ... 
              PANELwing(j).MIDPOINT(3), 'or', 'LineWidth',3);
        plot3(PANELwing(j).C4(1,1), ...
              PANELwing(j).C4(1,2), ... 
              PANELwing(j).C4(1,3),'*b', 'LineWidth',3);
        plot3(PANELwing(j).C4(2,1), ...
              PANELwing(j).C4(2,2), ... 
              PANELwing(j).C4(2,3),'*b', 'LineWidth',3);
    end 
    axis('equal')
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$z$','Interpreter','latex')

    figure(3)
    hold on
    for j=1:N*2*M
        PANELwing(j).PANELplot("c","yes");
    end  
    axis('equal')
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$z$','Interpreter','latex')
    
end

toc

end