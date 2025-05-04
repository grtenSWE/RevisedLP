function leave = revisedfindleave(Bmatrix, as, xb, phase, ~, indices)
% REVISEDFINDLEAVE   Minimum‐ratio test with Phase II tie‐breaking
%
%   leave = revisedfindleave(Bmatrix, as, xb, phase, n, indices)
%
%   Inputs:
%     Bmatrix — m×m current basis matrix
%     as      — m×1 column of A for the entering variable
%     xb      — m×1 current basic solution (B^{-1} b)
%     phase   — 1 or 2 (Phase 2 applies Bland‐style tie‐breaking)
%     n       — total number of original vars (not used here)
%     indices — m×1 vector of indices of the current basic variables
%
%   Output:
%     leave — row (1…m) of the leaving variable in the basis
%             0 if the problem is unbounded

    tol = 1e-8;

    % 1) compute direction
    d = Bmatrix \ as;      

    % 2) compute full ratio vector, invalid entries as Inf
    ratio = xb ./ d;
    invalid = (d <= tol);
    ratio(invalid) = Inf;

    % 3) check for unbounded
    rmin = min(ratio);
    if isinf(rmin)
        leave = 0;
        return;
    end

    % 4) gather all rows within tolerance of the min
    mask = abs(ratio - rmin) <= tol * max(1, abs(rmin));
    candidates = find(mask);   % guaranteed non‐empty

    % 5) tie‐break in Phase II
    if phase == 2 && numel(candidates) > 1
        [~, k] = min(indices(candidates));
        leave = candidates(k);
    else
        leave = candidates(1);
    end
end