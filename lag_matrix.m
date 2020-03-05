%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create matrix of lagged time series
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function YLag = lag_matrix(Y,lags)

% Check for a vector:
if numel(Y) == length(Y)
   Y = Y(:); % Ensure a column vector
end


lags = lags(:); % Ensure a column vector

missingValue = NaN;  % Assign default missing value

numLags = length(lags); % Number of lags to apply to each time series

[numObs,numSeries] = size(Y);

YLag = missingValue(ones(numObs,numSeries*numLags)); % Preallocate

for c = 1:numLags

    L       = lags(c);
    columns = (numSeries*(c-1)+1):c*numSeries; % Columns to fill, this lag

    if L > 0 % Time delays

       YLag((L + 1):end,columns) = Y(1:(end - L), :);

    elseif L < 0 % Time leads

       YLag(1:(end + L),columns) = Y((1 - L):end, :);

    else % No shifts

       YLag(:,columns) = Y;

    end

end