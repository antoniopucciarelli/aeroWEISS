function [v_ind,alpha_ind] = compute_INDUCEDvel(GAMMA,PANEL,alpha,M,N,U,flag)
% this function computes the induced velocity of a circulation
% discribution GAMMA over a finite 3D wing described by PANEL
%
% INPUT: 
%   GAMMA                 : vortex filamente circulation
%   PANEL                 : PANEL class array 
%   flow stream incidence : alpha
%   dihedral angle        : delta
%   M                     : # of discretization panels spanwise
%   N                     : # of discretization panels chordwise
% 
% WING -- schematic (* = panels)
%             0 >-------------------> M          0 >-------------------> M
%  0  k = 0   ************************* semiwing ************************* k = 0
%  |  k = 1   *************************          ************************* k = 1
%  |  k = ... *************************          ************************* k = ...
%  |          *************************          *************************
%  N  k = N   *************************          ************************* k = N
% 

% induced velocity vector
v_ind = zeros(N*2*M,3);
V     = zeros(N*2*M);

% determining vector normal to alpha --> inclination of the streamflow
alpha        = alpha/180*pi;
ROT_alpha    = ROT(0,-(alpha+pi/2),0);
alpha_normal = ROT_alpha * [1;0;0];

for i=1:2*M*N
    % !!!PAY ATTENTION!!! the contribute of induced velocity is computed 
    % only through the lateral filament of each panel that induced a 
    % vertical velocity at the midpoint of the studied panel 
    
    % this check is made for studying the 2 semiwing parts
    if(i == 1)
        k = 0;
    elseif(i == (M*N+1))
        k = 0;
    end
    
    % induced velocity computation
    for j=1:M
        
        v_ind1     = + 1/(4*pi) * GAMMA(j + k*M) / (PANEL(i).MIDPOINT(2) - PANEL(j + k*M).C4(2,2));
        
        v_ind2     = - 1/(4*pi) * GAMMA(j + k*M) / (PANEL(i).MIDPOINT(2) - PANEL(j + k*M).C4(1,2));
        
        v_ind(i,:) = v_ind(i,:) + (v_ind1 + v_ind2) * PANEL(j + k*M).normal';
        
        v_ind1     = + 1/(4*pi) * GAMMA(j + (N + k)*M) / (PANEL(i).MIDPOINT(2) - PANEL(j + (N + k)*M).C4(2,2));
        
        v_ind2     = - 1/(4*pi) * GAMMA(j + (N + k)*M) / (PANEL(i).MIDPOINT(2) - PANEL(j + (N + k)*M).C4(1,2));
        
        v_ind(i,:) = v_ind(i,:) + (v_ind1 + v_ind2) * PANEL(j + (N + k)*M).normal';
        
    end 
    
    V(i) = dot(v_ind(i,:),alpha_normal); 
    
    % if at k is adding 1 ==> means the program is studying the next wing line
    if(mod(i,M) == 0)
        k = k + 1;
    end
        
end

% initializing alpha_ind vec
alpha_ind = zeros(2*M,1);

for i=1:2*M 
    
    for j=1:N
        alpha_ind(i) = alpha_ind(i) + V(i+(j-1)*M);
    end
  
    % computing induced AOA
    alpha_ind(i) = atan(alpha_ind(i)/U);
    
end 

% plotting induced velocity 
if(flag == "yes")
    figure
    for i=1:2*M*N

        PANEL(i).PANELplot("c","no");
        quiver3(PANEL(i).MIDPOINT(1),PANEL(i).MIDPOINT(2),PANEL(i).MIDPOINT(3), ...
                0,0,100*v_ind(i),'b','LineWidth',3);
        hold on

    end 
    axis('equal');
    xlabel('x','Interpreter','latex');
    ylabel('y','Interpreter','latex');
    zlabel('x','Interpreter','latex');
    
    figure
    x = linspace(-1,1,length(alpha_ind));
    K = zeros(length(alpha_ind),1);
    for i=1:M
        K(i) = alpha_ind(M+1-i);
    end
    
    K(M+1:end) = alpha_ind(M+1:end);
    
    plot(x,K,'or-','LineWidth',2.5);
    grid on
    grid minor
    xlabel('$\% SPAN$','Interpreter','latex');
    ylabel('$\alpha_{ind}$','Interpreter','latex');

end 

end 