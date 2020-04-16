function YLag = hankel_matrix(Y,lags)
% HANKEL_MATRIX  Create matrix of lagged time series.
%   YLAG = HANKEL_MATRIX(Y,LAGS) given the time-series Y and the lags, 
%   creates the matrix of lagged time series.


% Check for a vector:
if numel(Y) == length(Y)
    Y = Y(:); % Ensure a column vector
end


lags = lags(:); % Ensure a column vector

default_value = NaN;  % Assign default missing value

numLags = length(lags); % Number of lags to apply to each time series

[numObs,numSeries] = size(Y);

YLag = default_value(ones(numObs,numSeries*numLags)); % Preallocation

for k = 1:numLags
    
    L       = lags(k);
    columns = (numSeries*(k-1)+1):k*numSeries; % Columns to fill, this lag
    
    if L > 0 % Time delays
        YLag((L + 1):end,columns) = Y(1:(end - L), :);
    else
        if L < 0 % Time leads
            YLag(1:(end + L),columns) = Y((1 - L):end, :);
        else % No shifts
            YLag(:,columns) = Y;
        end
        
    end
    
end

end