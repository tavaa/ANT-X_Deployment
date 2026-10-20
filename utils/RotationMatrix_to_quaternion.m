function qq = RotationMatrix_to_quaternion(R)
        % ROTM_TO_QUAT  Convert rotation matrix to unit quaternion.
        %
        %   Uses Shepperd's method for numerical stability.
        %   Returns qq = [q0; q1; q2; q3] with q0 >= 0 (scalar first).
        %   Satisfies: v_I = R * v_b  where R = R(qq).
 
            tr = R(1,1) + R(2,2) + R(3,3);
 
            if tr > 0
                s  = 0.5 / sqrt(tr + 1.0);
                q0 = 0.25 / s;
                q1 = (R(3,2) - R(2,3)) * s;
                q2 = (R(1,3) - R(3,1)) * s;
                q3 = (R(2,1) - R(1,2)) * s;
            elseif (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
                s  = 2.0 * sqrt(1.0 + R(1,1) - R(2,2) - R(3,3));
                q0 = (R(3,2) - R(2,3)) / s;
                q1 = 0.25 * s;
                q2 = (R(1,2) + R(2,1)) / s;
                q3 = (R(1,3) + R(3,1)) / s;
            elseif R(2,2) > R(3,3)
                s  = 2.0 * sqrt(1.0 + R(2,2) - R(1,1) - R(3,3));
                q0 = (R(1,3) - R(3,1)) / s;
                q1 = (R(1,2) + R(2,1)) / s;
                q2 = 0.25 * s;
                q3 = (R(2,3) + R(3,2)) / s;
            else
                s  = 2.0 * sqrt(1.0 + R(3,3) - R(1,1) - R(2,2));
                q0 = (R(2,1) - R(1,2)) / s;
                q1 = (R(1,3) + R(3,1)) / s;
                q2 = (R(2,3) + R(3,2)) / s;
                q3 = 0.25 * s;
            end
 
            % Canonical form: q0 >= 0 
            if q0 < 0
                q0 = -q0;  q1 = -q1;  q2 = -q2;  q3 = -q3;
            end
 
            % Normalize to correct floating-point drift
            nrm = sqrt(q0^2 + q1^2 + q2^2 + q3^2);
            qq = [q0; q1; q2; q3] / nrm;
        end