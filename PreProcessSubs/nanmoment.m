function sigma = nanmoment(x,order)
% MOMENT Central moments of all orders.
%   SIGMA = MOMENT(X,ORDER) returns the ORDER-th central sample moment of
%   the values in X.  For vector input, SIGMA is MEAN((X-MEAN(X)).^ORDER).
%   For a matrix input, MOMENT(X,ORDER) returns a row vector containing the
%   central moment of each column of X.  For N-D arrays, MOMENT operates
%   along the first non-singleton dimension.
%
%   MOMENT(X,ORDER,DIM) takes the moment along dimension DIM of X.
%
%   The first central moment is exactly zero. The second central moment is
%   the variance, using a divisor of N instead of N-1, where N is the
%   sample size.
%
%   See also MEAN, STD, VAR, SKEWNESS, KURTOSIS.

%   Copyright 1993-2009 The MathWorks, Inc.
%   $Revision: 1.8.2.4 $  $Date: 2009/03/09 19:49:54 $

if nargin < 2
    error('stats:moment:TooFewInputs', 'Moment requires two inputs.');
elseif ~isscalar(order)
    error('stats:moment:BadOrder', 'ORDER must be a scalar.');
end

% Return the first moment exactly.
if order == 1
    sigma = nanmean(zeros(size(x)));

else
    
    % Center X and compute the specified moment.
    sigma = nanmean((x - nanmean(x)).^order);
end
