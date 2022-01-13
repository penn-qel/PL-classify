function SSE = WeightedSSE(fdata,ydata,weights)
% Function to compute weighted sum of squared error.
% Useage:
%   SSE = WeightedSSE(fdata,ydata,weights)
% computes the sum over the weighted squared errors 
%   wi*(yi-fi)^2
% where wi are weights, and yi are observations.
%
% All arguments may be matrices, but they must be the same size.
%
% If no weights are given the unweighted SSE is computed.

% First check inputs

weighted = nargin==3;
datasize = size(fdata);
if any(size(ydata)~=datasize)
    error('All inputs must be the same size!');
end
if weighted && any(size(weights)~=datasize)
    error('Weights array must be the same size as fdata,ydata.');
end

if ~weighted
    weights = ones(datasize);
end

weights = 1./ydata; 
weights(isinf(weights)) = 0; 
SSE = sum((weights(:).*(ydata(:)-fdata(:)).^2));


