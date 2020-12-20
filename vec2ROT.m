function [ROT_BASE, ROT_BASE2] = vec2ROT(BASE)
% this function computes the rotation matrix from a rotation vector

if(or(size(BASE) == size(zeros(1,3)), size(BASE) == size(zeros(3,1))))

    ROT_BASE = [ 0      , -BASE(3),  BASE(2);
                 BASE(3), 0       , -BASE(1);
                -BASE(2),  BASE(1),        0 ];
    if(nargout > 1)
        ROT_BASE2 = ROT_BASE * ROT_BASE;
    end 
         
else
    
    error('INVALID VECTOR INPUT');
    
end 
end 