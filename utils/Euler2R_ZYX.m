function R = Euler2R_ZYX(phi, theta, psi)
% Rotation Matrix (FRD) -> NED  Z-Y-X
%
%   phi   : roll  (rotation around x_B)
%   theta : pitch (rotation around y_B)
%   psi   : yaw   (rotation around z_B)

    % trigonometric functions
    cphi = cos(phi);   sphi = sin(phi);
    cthe = cos(theta); sthe = sin(theta);
    cpsi = cos(psi);   spsi = sin(psi);

    % Rotation Matrix
    R = [cthe*cpsi,  sphi*sthe*cpsi - cphi*spsi,  cphi*sthe*cpsi + sphi*spsi;
         cthe*spsi,  sphi*sthe*spsi + cphi*cpsi,  cphi*sthe*spsi - sphi*cpsi;
         -sthe,      sphi*cthe,                   cphi*cthe               ];
end