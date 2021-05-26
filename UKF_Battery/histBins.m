function [yvals] = histBins(xvals,bin,xmax)
% FUNCTION divides x interval into 100 bins and computes the frequency of
% occurrence for each bin

% INPUTS: "raw" data, number of bins, max val (upper limit)

% OUPUTS: y-range (freq of occurrence)

rng = 2*xmax / bin;  % width of bins
low = -xmax;  % initialize lower bound
tot = length(xvals);   % number of sample points
for i = 1:1:bin
    % for each bin, count frequency of occurrence
    count = 0;   % counts frequency
    for j = 1:1:length(xvals)
        % find and count frequency from specified range
        if xvals(j) >= low && xvals(j) < low + rng
            count = count + 1;
        end
    end
    yvals(i) = (count / tot) / rng;
    low = low + rng;   % increment range to next bin
end

end

