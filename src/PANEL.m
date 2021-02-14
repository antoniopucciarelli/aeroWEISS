classdef PANEL
% this class describes the surface of the wing
% it describes the wing as a flat surface 
% the PANEL class hasn't any thickness
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
    
    properties
        VERTEX   = zeros(4,3)
        MIDPOINT = zeros(1,3)
        C4       = zeros(2,3)
        normal   = zeros(1,3)
        Vc       = 0.0
    end
    
    methods        
        function PANELplot(obj,flag,flag1)
            X = [obj.VERTEX; obj.VERTEX(1,:)];
            fill3(X(:,1),X(:,2),X(:,3),flag);

            hold on
            
            if(flag1 == "yes")
                plot3(obj.MIDPOINT(1),obj.MIDPOINT(2),obj.MIDPOINT(3),'*r','LineWidth',4)
                plot3(obj.C4(:,1),obj.C4(:,2),obj.C4(:,3),'*b','LineWidth',4)
                quiver3(obj.MIDPOINT(1),obj.MIDPOINT(2),obj.MIDPOINT(3), ...
                        obj.normal(1),obj.normal(2),obj.normal(3),'k','LineWidth',2);
            end 
        end
        
    end
end


