function [newBmatrix, newindices, newcb] = revisedbupdate(Bmatrix, indices, cb, cs, s, as, leave)
%   Inputs:
%     Bmatrix  — m×m matrix whose columns are the current basic columns of A
%     indices  — m×1 vector of the variable indices corresponding to each column of Bmatrix
%     cb       — m×1 vector of costs for the current basic variables
%     cs       — scalar cost of the entering (nonbasic) variable
%     s        — scalar index (identifier) of the entering variable
%     as       — m×1 column vector of constraint coefficients for the entering variable
%     leave    — integer position (1 ≤ leave ≤ m) of the column in Bmatrix to be replaced
%
%   Outputs:
%     newBmatrix — updated m×m basis matrix after replacing column ‘leave’
%     newindices — updated m×1 basic‐variable index vector
%     newcb      — updated m×1 basic‐cost vector
    
    newBmatrix = Bmatrix;
    newindices = indices;
    newcb = cb;

    newBmatrix(:,leave) = as;
    newindices(leave) = s;
    newcb(leave) = cs;
    
end