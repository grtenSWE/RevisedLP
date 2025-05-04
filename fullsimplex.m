function [z, x, pivalues, indices, exitflag] = fullsimplex(A, b, c, m, n)
% FULLSIMPLEX   Two-phase revised simplex algorithm
%
%   [z, x, pivalues, indices, exitflag] = fullsimplex(A, b, c, m, n)
%
%   Inputs:
%     A         — m×n constraint matrix
%     b         — m×1 right-hand side
%     c         — n×1 cost vector
%     m         — number of constraints
%     n         — number of original variables
%
%   Outputs:
%     z         — objective value
%     x         — solution vector (length n, originals only)
%     pivalues  — m×1 final simplex multipliers π
%     indices   — m×1 basic-variable indices (one per row of B, in order)
%     exitflag  — 0 = optimal, 1 = infeasible, –1 = unbounded

    % Phase I: setup
    A1       = [A, eye(m)];               
    c1       = [zeros(n,1); ones(m,1)];   
    B1       = eye(m);                    
    indices1 = (n+1 : n+m)';              

    [z1, x1_aug, pivalues1, indices1, flag1] = ...
        revisedsimplex(A1, b, c1, m, n+m, B1, indices1, 1);

    if flag1 == -1
        % unbounded in Phase I (degenerate LP), shouldn’t happen
        z = []; x = []; pivalues = []; indices = []; exitflag = -1;
        return;
    end

    if z1 > 1e-8
        % infeasible
        z        = z1;
        x        = x1_aug;     % only original variables
        pivalues = pivalues1;
        indices  = indices1;        % keep the row-order from Phase I
        exitflag = 1;
        return;
    end

    % Phase II: true objective on same augmented system
    A_aug  = A1;
    c2     = [c; zeros(m,1)];
    B2     = A_aug(:, indices1);

    [z2, x2_aug, pivalues, indices, exitflag] = ...
        revisedsimplex(A_aug, b, c2, m, n+m, B2, indices1, 2);

    % Always return only the first n (non-artificial) variables
    x = x2_aug(1:n);
    z = z2;
end