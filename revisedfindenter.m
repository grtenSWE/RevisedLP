function [as, cs, s] = revisedfindenter(A, pivalues, c, isbasic, phase)
% REVISED_FINDENTER   Choose entering variable with Bland's anti-cycling rule
%
%   [as, cs, s] = revisedfindenter(A, pivalues, c, isbasic, phase)
%
%   Inputs:
%     A         — m×n_total (possibly augmented) constraint matrix
%     pivalues  — m×1 simplex multipliers π
%     c         — n_total×1 cost vector
%     isbasic   — n_total×1 logical mask of basic variables
%     phase     — 1 or 2
%
%   Outputs:
%     as        — m×1 entering-column of A (empty if none)
%     cs        — cost c(s) (empty if none)
%     s         — index of entering var (0 if optimal / none enter)

    tol = 1e-8;
    [m, ~] = size(A);

    % Compute reduced costs for all variables
    r = c.' - (pivalues.' * A);   % 1×n_total

    % Mask out all basic variables
    r(isbasic.') = +Inf;

    % In Phase 2, also mask out artificial columns
    orig_n = length(c) - m;
    if phase == 2
        r(orig_n+1 : end) = +Inf;
    end

    % Bland's Rule: select the smallest index among those with negative reduced cost
    eligible = find(r < -tol);
    if isempty(eligible)
        % no negative reduced cost ⇒ optimal
        s  = 0;
        as = [];
        cs = [];
    else
        % pick smallest index to avoid cycling
        s = min(eligible);
        as = A(:, s);
        cs = c(s);
    end
end