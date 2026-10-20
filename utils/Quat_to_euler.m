function [phi, theta, psi] = Quat_to_euler(qq)
% QUAT_TO_EULER  Convert quaternion to Z-Y-X Euler angles.
%
%   Utility for comparing A01 states to A00 states in test scripts.
%   Returns [phi; theta; psi] in radians (Z-Y-X, psi in [-pi,pi]).
 
      q0 = qq(1);  q1 = qq(2);  q2 = qq(3);  q3 = qq(4);
 
      phi   = atan2(2*(q0*q1+q2*q3), q0^2-q1^2-q2^2+q3^2);
      theta = asin (max(-1, min(1, 2*(q0*q2-q3*q1))));
      psi   = atan2(2*(q0*q3+q1*q2), q0^2+q1^2-q2^2-q3^2);
end
 
   

