function [CI, avg] = confidenceinterval(M, alpha)
%confidenceinterval calculates the average and confidence interval of each
%column in matrix M using CI = t*s/sqrt(n).
%   INPUT:
%   > M     - Matrix containing the data, where samples are in the rows and
%             variables in the columns
%   > alpha - Limit for 2-sided t statistic (if alpha is 0.05 the upper and
%             lower limit are 0.025 and 0.975)
%   OUTPUT:
%   > CI    - vector containing confidenceinterval of each column in M
%   > avg   - vector containing average value of each column in M

% calculate the standard error:
SEM = std(M)/sqrt(size(M,1));
% calculate two sided t-statistic
t_value = tinv(1-alpha./2, size(M,1)-1);

CI = t_value.*SEM;
avg = mean(M);

end

