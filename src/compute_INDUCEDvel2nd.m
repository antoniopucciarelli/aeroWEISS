function [v_ind,alpha_ind] = compute_INDUCEDvel2nd(GAMMA,PANEL,AOA,alpha,M,N,U,flag)
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

% induced velocity vector
v_ind = zeros(N*2*M,1);

% rotation matrix
AOA       = AOA/180*pi;
ROTnormal = ROT(0,-(AOA+pi/2),0);
normal    = ROTnormal * [1;0;0];

% determining vector normal to alpha --> inclination of the streamflow
alpha        = alpha/180*pi;
ROT_alpha    = ROT(0,-(alpha-pi/2),0);
alpha_normal = ROT_alpha * [1;0;0];

% determining wing span lenght
L = abs(PANEL(N*M).VERTEX(3,1));
%L = norm(PANEL(1).VERTEX(1,:)-PANEL(1).VERTEX(2,:));

for i=1:2*M*N
    % !!!PAY ATTENTION!!! the contribute of induced velocity is computed 
    % only through the lateral filament of each panel that induced a 
    % vertical velocity at the midpoint of the studied panel 
    
    k = 0;
    
    % initializing V
    V = zeros(1,3);
    
    for j=1:M    
        
        % PANEL kth -- right filament -- 1st part
        R0_vec = PANEL(j + k*M).VERTEX(2,:) - PANEL(j + k*M).C4(2,:);
        R1_vec = PANEL(i).MIDPOINT - PANEL(j + k*M).C4(2,:);
        R2_vec = PANEL(i).MIDPOINT - (PANEL(j + k*M).C4(2,:) + R0_vec * L);
        PR_vec = cross(R1_vec,R2_vec);
        
        % computing length
        R1     = norm(R1_vec);
        R2     = norm(R2_vec);
        PR     = norm(PR_vec);
        
        % computing induced velocity from the PANEL kth right filament -- 1st part
        V_1_K = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2 * GAMMA(j + k*M);
        
        % PANEL kth -- right filament -- 2nd part (to infinity)
        R1_vec = PANEL(i).MIDPOINT - (PANEL(j + k*M).C4(2,:) + R0_vec * L);
        R2_vec = PANEL(i).MIDPOINT - (PANEL(j + k*M).C4(2,:) + R0_vec * L + [1,0,0]*1e+3);
        R0_vec = [1,0,0]*1e+3;
        PR_vec = cross(R1_vec,R2_vec);
        
        % computing length
        R1     = norm(R1_vec);
        R2     = norm(R2_vec);
        PR     = norm(PR_vec);
        
        % computing induced velocity from the PANEL kth right filament -- 2nd part
        V_2_K = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2 * GAMMA(j + k*M);
        
        % PANEL kth -- left filament -- 1st part
        R0_vec = PANEL(j + k*M).C4(1,:) - PANEL(j + k*M).VERTEX(3,:);
        R2_vec = PANEL(i).MIDPOINT - PANEL(j + k*M).C4(1,:);
        R1_vec = PANEL(i).MIDPOINT - (PANEL(j + k*M).C4(1,:) - R0_vec * L);
        PR_vec = cross(R1_vec,R2_vec);
        
        % computing length
        R1     = norm(R1_vec);
        R2     = norm(R2_vec);
        PR     = norm(PR_vec);
        
        % computing induced velocity from the PANEL kth left filament -- 1st part
        V_3_K = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2 * GAMMA(j + k*M);
        
        % PANEL kth -- left filament -- 2nd part (to infinity)
        R2_vec = PANEL(i).MIDPOINT - (PANEL(j + k*M).C4(1,:) - R0_vec * L);
        R1_vec = PANEL(i).MIDPOINT - (PANEL(j + k*M).C4(1,:) - R0_vec * L + [1,0,0]*1e+3);
        R0_vec = - [1,0,0]*1e+3;
        PR_vec = cross(R1_vec,R2_vec);
        
        % computing length
        R1     = norm(R1_vec);
        R2     = norm(R2_vec);
        PR     = norm(PR_vec);
        
        % computing induced velocity from the PANEL kth left filament -- 2nd part
        V_4_K = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2 * GAMMA(j + k*M);
        
        % PANEL (k + N)th -- right filament -- 1st part
        R0_vec = PANEL(j + k*M).VERTEX(2,:) - PANEL(j + k*M).C4(2,:);
        R1_vec = PANEL(i).MIDPOINT - PANEL(j + k*M).C4(2,:);
        R2_vec = PANEL(i).MIDPOINT - (PANEL(j + k*M).C4(2,:) + R0_vec * L);
        PR_vec = cross(R1_vec,R2_vec);
        
        % computing length
        R1     = norm(R1_vec);
        R2     = norm(R2_vec);
        PR     = norm(PR_vec);
        
        % computing induced velocity from the PANEL (N + k)th right filament -- 1st part
        V_1_NK = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2 * GAMMA(j + (N + k)*M);
        
        %  PANEL (k + N)th -- right filament -- 2nd part (to infinity)
        R1_vec = PANEL(i).MIDPOINT - (PANEL(j + (N + k)*M).C4(2,:) + R0_vec * L);
        R2_vec = PANEL(i).MIDPOINT - (PANEL(j + (N + k)*M).C4(2,:) + R0_vec * L + [1,0,0]*1e+3);
        R0_vec = [1,0,0]*1e+3;
        PR_vec = cross(R1_vec,R2_vec);
        
        % computing length
        R1     = norm(R1_vec);
        R2     = norm(R2_vec);
        PR     = norm(PR_vec);
        
        % computing induced velocity from the PANEL (N + k)th right filament -- 2nd part
        V_2_NK = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2 * GAMMA(j + (N + k)*M);
       
        %  PANEL (k + N)th -- left filament -- 1st part
        R0_vec = PANEL(j + (N + k)*M).C4(1,:) - PANEL(j + k*M).VERTEX(3,:);
        R2_vec = PANEL(i).MIDPOINT - PANEL(j + (N + k)*M).C4(1,:);
        R1_vec = PANEL(i).MIDPOINT - (PANEL(j + (N + k)*M).C4(1,:) - R0_vec * L);
        PR_vec = cross(R1_vec,R2_vec);
        
        % computing length
        R1     = norm(R1_vec);
        R2     = norm(R2_vec);
        PR     = norm(PR_vec);
        
        % computing induced velocity from the PANEL (N + k)th left filament -- 1st part
        V_3_NK = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2 * GAMMA(j + (N + k)*M);
        
        % PANEL (k + N)th -- left filament -- 2nd part (to infinity)
        R2_vec = PANEL(i).MIDPOINT - (PANEL(j + (N + k)*M).C4(1,:) - R0_vec * L);
        R1_vec = PANEL(i).MIDPOINT - (PANEL(j + (N + k)*M).C4(1,:) - R0_vec * L + [1,0,0]*1e+3);
        R0_vec = - [1,0,0]*1e+3;
        PR_vec = cross(R1_vec,R2_vec);
        
        % computing lenght
        R1     = norm(R1_vec);
        R2     = norm(R2_vec);
        PR     = norm(PR_vec);
        
        % computing induced velocity from the PANEL (N + k)th left filament -- 2nd part 
        V_4_NK = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2 * GAMMA(j + (N + k)*M);
        
        % updating V
        V = V + V_1_K + V_1_NK + V_3_K + V_3_NK;
        V = V + V_2_K + V_2_NK + V_4_K + V_4_NK;
        
    end 
    
    V = - V;
    
    v_ind(i) = dot(V,alpha_normal); 
    
    k = k + 1;
        
end

% initializing alpha_ind vec
alpha_ind = zeros(2*M,1);

for i=1:2*M 
    
    for j=1:N
        alpha_ind(i) = alpha_ind(i) + v_ind(i+(j-1)*M);
    end
    
    % induced AoA for the ith airfoil of the wing
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