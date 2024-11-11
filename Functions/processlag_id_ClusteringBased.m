function shift_pred = processlag_id_ClusteringBased(X, measurement_times, maxshift, estimation_method, varargin)

%DESCRIPTION:
%This function will determine the lag between sensors of a dataset by first
%applying a clustering method and then applying a time delay estimation
%method for each link in the cluster moving from the bottom to the top of
%the tree. It wil thus always shift two variables at a time (a 
%representative is selected when shifting two groups to each other).

%INPUT:
%-X: The data for the sensors for which the lag should be identified 
%  (sensor as columns, time points as rows).
%-measurement_times: The time axis for the sensor data (as column vector).
%-maxshift: A scalar indicating the maximum lag allowed between two sensors 
%  in either direction (symetric). THIS SHOULD BE A NUMBER OF DATAPOINTS!
%-estimation_method: Method used to estimate the shifts
%   - 'biv-seq': Sequentially adds sensors to a reference sensor, while
%                optimizing its (bivariate) correlation to that reference 
%                through shifting.
%   - 'pca-con': Adds shifting copies of all sensors (except the reference),
%                makes a PCA model and selects per sensor the most 
%                contributing shift.
%   - 'pls-con': Adds shifting copies of all sensors (except the reference),
%                regresses the reference sensor on all others, and selects 
%                per (other) sensor the most contributing shift. (Default)
%   - 'pls-seq': Sequentially adds sensors to a growing block of reference
%                sensors, while optimizing its shift in terms of how well it
%                can by regressed on the reference sensors with PLS.

% Additional name-value pairs:
%- 'contribution': When using the 'pls-con' method for identifying sensor
%  lags, this is the sensor contribution measure used. This can be:
%  - 'loading':     The sum of squared PC loadings, weighted according to 
%                   the fraction of variance explained by the PC.
%  - 'vip':         The variable importance in projetion.
%  - 'coefficient': The absolute regression coefficients (B).

%- 'optimisation': Only applicanle for 'pca-con' and 'pls-con'
%  - 'all':         optimise number of PCs in created PCA model with all of
%                   the possible PCs, not recomended.
%  - '#var':        optimise number of PCs in created PCA model with a 
%                   maximum equal to the number of variables, better for
%                   larger datasets (>30 variables) (default)
%  - 'non':         skip optimisation and only use one PC, better for
%                   smaller datasets (<30 variables)

% OUTPUT
% - shift_pred: Contains the actual lags determined, in the same units as
%   the time axis was inputted.

% Check the funtion input
if nargin < 4
    estimation_method = 'pls-con';
end
if nargin < 3
    error('Some funtion input missing, please check description')
end

contribution = 'loading';
optimisation = '#var';
for i=1:2:length(varargin)
    eval([varargin{i} ' = varargin{i+1};']);
end
if strcmpi(estimation_method, 'pca-con')
    contribution = 'loading';
end

% calculate sampling frequency
sf = hours( measurement_times(2) - measurement_times(1) );


%% Perform clustering with the data
shifts = -maxshift:1:maxshift;
% Create the correlation matrix with which the clustering is performed.
maxcorr = NaN(size(X, 2), size(X, 2));
maxcorr_shifts = NaN(size(X, 2), size(X, 2));
for r = 1:size(X,2)
    for c = setdiff(1:size(X,2), r)
        [cor,lags] = xcorr(X(:,r), X(:,c), maxshift, 'coeff');
        cor = cor.^2;
        [M,I] = max(cor);
        maxcorr(r,c) = M;
        maxcorr_shifts(r,c) = shifts(I);
    end
end
% Transform to distance-like matrix
maxcorr_inv = 1-maxcorr;
% Set the diagonal (NaN values) to 0
maxcorr_inv(isnan(maxcorr_inv)) = 0;

% Now the important information is transformed using squareform. Because
% only this transformed for work with the linkage funtion
maxcorr_inv_sq = squareform(maxcorr_inv);
Z = linkage(maxcorr_inv_sq, 'average');


%% Create a cell array that contains the variables present for each number in the three (Z)
var_numb = size(X, 2);

% initialise
VarsInLink = cell(max(Z(:)), 1);
% loop through all of the links
for i = 1:max(Z(:))
    % If the number is less than 95 (var_numb) it will be the same as the 
    % original variable, otherwise it is a combination of variables
    if i <= var_numb
        temp_output = i;    
    else % otherwise it is the 2 variables put together in position-95
        temp_output = Z(i-var_numb, 1:2); 
        % while there is still numbers higher than 95 they will be reduced
        % to their originally put together variables. 
        while any(temp_output > var_numb)
            temp_lessthen = temp_output(temp_output <= var_numb);
            temp_morethen = temp_output(temp_output > var_numb);
            temp_output = temp_lessthen;
            for n = temp_morethen
                temp_insidemorethen = Z(n-var_numb, 1:2);
                temp_output = [temp_output temp_insidemorethen];
            end
        end    
    end  
    VarsInLink{i} = temp_output;
end


%% loop through the rows in Z to determine shifts cluster by cluster

% we want to save the shifts in an array and the actually shifted variables
% in a matrix. NaN values are added in front and after each variable to be
% able to use circshift() and such that no excess information is removed. 

% initialise both the shifted x matrix and the array containing the shift
% values
X_shiftinfo = [nan(maxshift, var_numb); X; nan(maxshift, var_numb)];
shifts_est = zeros(1,var_numb);

% loop through all of the links in the tree
for i = 1:size(Z,1)
    % find the variables that are shifted to each other and put them into a
    % matrix
    group1 = VarsInLink{Z(i,1)};
    group2 = VarsInLink{Z(i,2)};
    
    all_current_vars = [X_shiftinfo(:,group1) X_shiftinfo(:,group2)];
    columns1 = 1:length(group1);
    columns2 = length(group1)+1:length(group1)+length(group2);
    
    % Remove all the rows containing NaN values and determine a
    % representative variable (for both groups) if multiple variables are 
    % present in the groups.
    rows_removed = any(isnan(all_current_vars), 2);
    all_current_vars_nonan = all_current_vars;
    all_current_vars_nonan(rows_removed, :) = [];
    current_times = [hours(-maxshift*sf:sf:0-sf)'; measurement_times; hours(length(measurement_times)*sf+sf:sf:maxshift*sf+length(measurement_times)*sf)'];
    current_times(rows_removed, :) = [];
    
    % Determine the representative variables for each group (this is the
    % variable with the highest sum of correlations to the variables in the
    % other group
    corr_matrix = abs(corr(all_current_vars_nonan(:, columns1), all_current_vars_nonan(:, columns2)));
    [~, max_corr_I] = max(corr_matrix(:));
    [represent1_I, represent2_I] = ind2sub(size(corr_matrix), max_corr_I);
    represent1 = all_current_vars_nonan(:, columns1(represent1_I));
    represent2 = all_current_vars_nonan(:, columns2(represent2_I));


    % determine the shifts that are allowed (cannot shift more than 24 or
    % less then -24 hours). Since the first group is always used as the
    % reference, only the allowed shifts of the second group have to be
    % determined
    nan_rows = all(isnan(all_current_vars(:, columns2)), 2);
    neg_shifts = sum(nan_rows(1:ceil(length(nan_rows)/2)));
    pos_shifts = sum(nan_rows(ceil(length(nan_rows)/2):end));
    shifts_allowed = hours(-neg_shifts*sf:4:pos_shifts*sf);
    
    % Add the representatives together and determine the time delay between
    % the two groups
    X = [represent1 represent2];
    [lags_pred, ~] = processlag_id(X, current_times, shifts_allowed, estimation_method, 'reference', 'pearson-shift', 'contribution', contribution, 'optimisation', optimisation);
    lags_pred_h = hours(lags_pred);  
     
    % Change the X_shiftinfo matrix to contain the newly found shift
    X_shiftinfo(:,group1) = circshift(X_shiftinfo(:,group1), lags_pred_h(1)/sf, 1);
    X_shiftinfo(:,group2) = circshift(X_shiftinfo(:,group2), lags_pred_h(2)/sf, 1);
    
    % Save the found shifts in the shifts_est array
    shifts_est(group1) = shifts_est(group1) + lags_pred_h(1);
    shifts_est(group2) = shifts_est(group2) + lags_pred_h(2);
    
end    

% give the output
shift_pred = hours(shifts_est);

end

