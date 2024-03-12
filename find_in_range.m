function [in,out, inPos, outPos] = find_in_range(f,r)
% find_in_range - Find indices of values within a specified range.
%
%   [in, out, inPos, outPos] = find_in_range(f, r)
%
%   This function finds the indices of values within a specified range 'r' in the input vector 'f'.
%
%   Parameters:
%   - f (vector): Input vector of values.
%   - r (vector): Range specified by a two-element vector [min_value, max_value].
%
%   Returns:
%   - in (vector): Indices of values within the specified range.
%   - out (vector): Indices of values outside the specified range.
%   - inPos (vector): Indices of positive values within the specified range.
%   - outPos (vector): Indices of positive values outside the specified range.
%
% Additional comments for clarity
% in: indices of all values within the specified range (including negative values)
% out: indices of all values outside the specified range
% inPos: indices of positive values within the specified range
% outPos: indices of positive values outside the specified range

% Sort the range in ascending order to avoid errors or NaN outputs
r = sort(r, 'ascend');

inPosRange = find(f < r(2) & f > r(1));
inNegRange = find(f < -r(1) & f > -r(2));
in = sort([inNegRange' inPosRange'], 'ascend');
out = setdiff(1:length(f), in);

% Only positive values
inPos = find(f >=0 & f < r(2) & f > r(1))';
outPos = find(f > 0 & (f <= r(1) | f >= r(2)))';
end