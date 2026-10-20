function xref_out = wrap_yaw_reference(xref_seq, psi_current)
% Sequential yaw unwrapping along the MPC horizon.
% Prevents 2*pi jumps in the QP cost function.
    xref_out   = xref_seq;
    psi_anchor = psi_current;
    for j = 1:size(xref_seq, 2)
        delta           = atan2(sin(xref_seq(9,j) - psi_anchor), cos(xref_seq(9,j) - psi_anchor));
        psi_adj         = psi_anchor + delta;
        xref_out(9,j)   = psi_adj;
        psi_anchor      = psi_adj;
    end
end