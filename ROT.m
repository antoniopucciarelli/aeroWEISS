function [R,RA,RB,RC,RBI,RCI,RCJ] = ROT(phi, theta, psi,flag1)
% this function computes the rotation matrix for EULER angles
% !!! PAY ATTENTION => IN ORDER TO COMPUTE THE ROTATION MATRIX 
%     IT IS NECESSARY YOU PUT THE PROIECTION OF THE TRANSFORMED 
%     VECTORS IN COLUMNS                                   !!!
% INPUT: ANGLES MUST BE IN RADIANTS


RA = [  cos(phi), -sin(phi), 0; ...
        sin(phi),  cos(phi), 0; ...
        0       ,  0       , 1 ];

RB = [ cos(theta), 0, sin(theta);
       0         , 1, 0         ;
      -sin(theta), 0, cos(theta) ];
    
RC = [ 1,         0,         0;
       0,  cos(psi), -sin(psi);
       0,  sin(psi),  cos(psi) ];
 
RBI = RA * RB * RA';
RCJ = RB * RC * RB';
RCI = RA * RCJ * RA'; 
R   = RA * RB * RC;

if(nargin == 4 && flag1 == "print")

        I = eye(3);

        figure

        J = RA  * I;
        K = RBI * I;
        L = RCI * I;

        for i=1:3
            O1 = quiver3(0,0,0,I(i,1),I(i,2),I(i,3),0,'k');
            hold on
        end
        for i=1:3
            O2 = quiver3(0,0,0,J(i,1),J(i,2),J(i,3),0,'r');
        end 
        for i=1:3
            O3 = quiver3(0,0,0,K(i,1),K(i,2),K(i,3),0,'b');
        end
        for i=1:3
            O4 = quiver3(0,0,0,L(i,1),L(i,2),L(i,3),0,'g');
        end

        title('$\mathcal{F_{I}} \Rightarrow \ \phi \ \theta \ \psi$','Interpreter','latex');
        xlabel('x');
        ylabel('y');
        zlabel('z');
        legend([O1,O2,O3,O4],'$ \mathcal{F_{I}} $','$ \phi $','$ \theta $','$ \psi $','Interpreter','latex');
        axis('equal')
        fprintf(' φ = %f \n Θ = %f \n Φ = %f \n',phi,theta,psi);
        view([1,1,1]);

    end 
end 
