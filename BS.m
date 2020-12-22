function [MATRIX] = BS(PANEL,M,N,L,toll)
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

tic

MATRIX = ones(N*2*M,N*2*M);

for i=1:N*2*M
    for j=1:N*2*M
        
        % computing C/4 line induction 
        R0_vec = PANEL(j).C4(2,:)  - PANEL(j).C4(1,:);
        R1_vec = PANEL(i).MIDPOINT - PANEL(j).C4(1,:);
        R2_vec = PANEL(i).MIDPOINT - PANEL(j).C4(2,:);
        R0     = norm(R0_vec);
        R1     = norm(R1_vec);
        R2     = norm(R2_vec);
        PR_vec = cross(R1_vec,R2_vec);
        PR     = norm(PR_vec);
        
        if(PR/R0 < toll)
            Vc_C4 = 0.0;
        else
            Vc_C4  = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2;
        end
        
        % computing LATERAL2 induced velocity
        R0_vec  = (PANEL(j).VERTEX(2,:) - PANEL(j).C4(2,:)) * 1000 * L;
        R1_vec  = PANEL(i).MIDPOINT     - PANEL(j).C4(2,:);
        R2_vec  = PANEL(i).MIDPOINT     - (PANEL(j).C4(2,:) + R0_vec);
        R0      = norm(R0_vec);
        R1      = norm(R1_vec);
        R2      = norm(R2_vec);
        PR_vec  = cross(R1_vec,R2_vec);
        PR      = norm(PR_vec);
        
        if(PR/R0 < toll)
            Vc_inf2 = 0.0;
        else
            Vc_inf2 = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2;
        end
        
        % computing LATERAL1 induced velocity
        R0_vec  = - (PANEL(j).VERTEX(3,:) - PANEL(j).C4(1,:)) * 1000 * L;
        R2_vec  = PANEL(i).MIDPOINT       - PANEL(j).C4(1,:);
        R1_vec  = PANEL(i).MIDPOINT       - (PANEL(j).C4(1,:) - R0_vec);
        R0      = norm(R0_vec);
        R1      = norm(R1_vec);
        R2      = norm(R2_vec);
        PR_vec  = cross(R1_vec,R2_vec);
        PR      = norm(PR_vec);
        
        if(PR/R0 < toll)
            Vc_inf1 = 0.0;
        else
            Vc_inf1 = 1/(4*pi) * dot(R0_vec,(R1_vec/R1 - R2_vec/R2)) * PR_vec/PR^2;   
        end
        
        MATRIX(i,j) = dot(Vc_C4 + Vc_inf1 + Vc_inf2,PANEL(i).normal);
    
    end 
end 

toc

end