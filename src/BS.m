function [MATRIX] = BS(PANEL,AOA,M,N,L,toll)
% this function computes the induced velocity in a point by the panel of
% unit circulation value with the BIOT SAVART law.
%
% INPUT:
%   PANEL : PANEL array class 
%   M     : # of spanwise discretization points
%   N     : # of chordwise discetization points
%   toll  : tollerance used to avoid singular MATRIX caused by short
%           distance between inducing filament and induced point
%
% OUTPUT: 
%   MATRIX : system matrix -- describes the non penetration condition for
%            the vortex filament inducing flow stream
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
%  R0 |     o 
%     |   * C
%     | *         
%  R1 o beta1
%

MATRIX = ones(N*2*M,N*2*M);

% rotating the first part of the vortex filament
AOA     = AOA/180*pi;
R       = ROT(0, AOA, 0, 'noprint');
R0_vec1 = [1;0;0];
R0_vec1 = R * R0_vec1;
R0_vec1 = R0_vec1';

for i=1:N*2*M
    for j=1:N*2*M
        
        % studying PANEL-jth induction on PANEL-ith
        
        % computing C/4 line induction 
        R0_vec = PANEL(j).C4(2,:)  - PANEL(j).C4(1,:);
        R1_vec = PANEL(i).MIDPOINT - PANEL(j).C4(1,:);
        R2_vec = PANEL(i).MIDPOINT - PANEL(j).C4(2,:);
        R0     = norm(R0_vec);
        R1     = norm(R1_vec);
        R2     = norm(R2_vec);
        PR_vec = cross(R1_vec,R2_vec);
        PR     = norm(PR_vec);
        
        % C/4 line induction : INDUCED VELOCITY CONDITION
        if(PR/R0 < toll)
            Vc_C4 = 0.0;
        else
            Vc_C4  = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2;
        end
        
        % computing LATERAL2 induced velocity
        R0_vec       = R0_vec1 * 0.5 * L; % the 1st part of the vortex filament extends for X times the wing span       
        R1_vec       = PANEL(i).MIDPOINT - PANEL(j).C4(2,:);
        R2_vec       = PANEL(i).MIDPOINT - (PANEL(j).C4(2,:) + R0_vec);
        START2ndPART = PANEL(j).C4(2,:)  + R0_vec; 
        R0           = norm(R0_vec);
        R1           = norm(R1_vec);
        R2           = norm(R2_vec);
        PR_vec       = cross(R1_vec,R2_vec);
        PR           = norm(PR_vec);
        
        % second filament : INDUCED VELOCITY CONDITION --> 1st part
        if(PR/R0 < toll)
            Vc_inf2_start = 0.0;
        else
            Vc_inf2_start = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2;
        end
        
        % computing LATERAL2 [second filament] induced velocity
        % !!! PAY ATTENTION !!! the second filament starts at the point
        % indicated by START2ndPART 
        R0_vec  = [1,0,0] * 1000 * L; % this time the R0_vec points to infinity parallel to x-axes
        R1_vec  = PANEL(i).MIDPOINT - START2ndPART;
        R2_vec  = PANEL(i).MIDPOINT - (START2ndPART + R0_vec);
        R0      = norm(R0_vec);
        R1      = norm(R1_vec);
        R2      = norm(R2_vec);
        PR_vec  = cross(R1_vec,R2_vec);
        PR      = norm(PR_vec);
        
        % second filament : INDUCED VELOCITY CONDITION --> 2nd part
        if(PR/R0 < toll)
            Vc_inf2_end = 0.0;
        else
            Vc_inf2_end = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2;
        end        
        
        % computing LATERAL1 induced velocity
        R0_vec       = - R0_vec1 * 0.5 * L;
        R1_vec       = PANEL(i).MIDPOINT - (PANEL(j).C4(1,:) - R0_vec);
        R2_vec       = PANEL(i).MIDPOINT - PANEL(j).C4(1,:);
        START2ndPART = PANEL(j).C4(1,:) - R0_vec;
        R0           = norm(R0_vec);
        R1           = norm(R1_vec);
        R2           = norm(R2_vec);
        PR_vec       = cross(R1_vec,R2_vec);
        PR           = norm(PR_vec);
        
        % first filament : INDUCED VELOCITY CONDITION --> 1st part
        if(PR/R0 < toll)
            Vc_inf1_start = 0.0;
        else
            Vc_inf1_start = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2;   
        end
        
        % computing LATERAL1 induced velocity
        % !!! PAY ATTENTION !!! the second filament starts at the point
        % indicated by START2ndPART 
        R0_vec  = - [1,0,0] * 1000 * L;
        R1_vec  = PANEL(i).MIDPOINT - (START2ndPART - R0_vec);
        R2_vec  = PANEL(i).MIDPOINT - START2ndPART;
        R0      = norm(R0_vec);
        R1      = norm(R1_vec);
        R2      = norm(R2_vec);
        PR_vec  = cross(R1_vec,R2_vec);
        PR      = norm(PR_vec);
        
        % first filament : INDUCED VELOCITY CONDITION --> 2nd part
        if(PR/R0 < toll)
            Vc_inf1_end = 0.0;
        else
            Vc_inf1_end = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2;   
        end
        
        % allocating normal induced (UNITY GAMMA) velocity
        MATRIX(i,j) = dot(Vc_C4 + Vc_inf1_start + Vc_inf1_end + Vc_inf2_start + Vc_inf2_end ,PANEL(i).normal);
    
    end 
end 

end