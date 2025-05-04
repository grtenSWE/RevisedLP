function [z, x, pivalues, indices, exitflag] = revisedsimplex(A, b, c, ~, n, B, indices, phase)
%REVISEDSIMPLEX   Revised simplex method (Phase I or Phase II)

%   Inputs:
%     A         — m×n_total constraint matrix (augmented if phase II)
%     b         — m×1 right-hand side vector
%     c         — n_total×1 cost vector
%     m         — number of constraints
%     n         — total number of variables (original + artificials)
%     B         — m×m current basis matrix
%     indices   — m×1 current basic-variable index list
%     phase     — 1 for Phase I (feasibility), 2 for Phase II (optimality)
%
%   Outputs:
%     z         — objective value (Phase I sum of artificials or Phase II original cᵀx)
%     x         — solution vector (length n_total)
%     pivalues  — m×1 simplex multipliers π
%     indices   — m×1 vector of final basic-variable indices
%     exitflag  — 0 if optimal; 1 if infeasible (Phase I); -1 if unbounded

    cb = c(indices);
    isbasic = false(n,1);
    isbasic(indices) = true;
    exitflag = 0;

    while true
        % compute multipliers π (m×1):  B' * π = c_B
        pivalues = B' \ cb;
        % find entering variable
        [as, cs, s] = revisedfindenter(A, pivalues, c, isbasic, phase);
        if s == 0
            % optimal (no negative reduced cost)
            exitflag = 0;
            break;
        end
        % compute current basic solution x_B = B^{-1} b
        xb = B \ b;
        % find leaving row
        leave = revisedfindleave(B, as, xb, phase, n, indices);
        if leave == 0
            % unbounded
            exitflag = -1;
            break;
        end
        % pivot update
        [B, indices, cb] = revisedbupdate(B, indices, cb, cs, s, as, leave);
        % update isbasic map
        isbasic(:) = false;
        isbasic(indices) = true;
    end

    % assemble full solution vector
    x = zeros(n,1);
    xb = B \ b;
    x(indices) = xb;
    % objective value
    z = c' * x;
end
