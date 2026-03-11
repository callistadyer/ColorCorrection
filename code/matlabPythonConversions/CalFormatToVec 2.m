function [xVec] = CalFormatToVec(calFormat)
% Converts a CalFormat matrix into one long vector using MATLAB ordering
%
% Syntax:
%    [xVec] = CalFormatToVec(calFormat)
%
% Inputs:
%   calFormat:   Input CalFormat matrix of size 3 x N
%
% Outputs:
%   xVec:        Output vector of size (3*N) x 1
%
% Description:
%   Converts a CalFormat matrix into one long vector using
%   column-major ordering.
%
%   For a CalFormat matrix of size 3 x N, the output is exactly what is
%   returned by calFormat(:).
%
%   Thus if the columns of calFormat correspond to pixels, the output
%   vector is ordered as
%       [R1; G1; B1; R2; G2; B2; ... ].
%
% History:
%   03/10/2026  cmd    Wrote it.
%
% Examples:
%{
calFormat = [1 2 3; 10 20 30; 100 200 300];
xVec = CalFormatToVec(calFormat);
%}

xVec = calFormat(:);

end