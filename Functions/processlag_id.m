function [lags, meta] = processlag_id(data, time, maxlag, method, varargin)

%DESCRIPTION:
%This function will determine the lag between sensors based on the data
%measured for them. Multiple methods for doing this are supported, as given
%below. The data needs to be sampled at a regular time frequency.
%
%INPUT:
%- data: The data for the sensors for which the lag should be identified 
%  (sensor as columns, time points as rows).
%- time: The time axis for the sensor data (as column vector).
%- maxlag: This can either be a scalar (size=1x1) indicating the maximum 
%  lag allowed between two sensors in either direction (symetric). Or it 
%  can be an array (1xn) containing every possible shift, where each variable 
%  has the same maxlag (asymetric). Lastly, the maxlag can be a cell array
%  where each cell contains the allowed shift of that variable. Any of 
%  these should be in the same units as the time axis.
%- method: The method for determining the lag, which can be:
%  - 'BIV-PEAR-SEQ': Sequentially adds sensors to a reference sensor, 
%                    while optimizing its (bivariate) Pearson correlation 
%                    to that reference through shifting.
%  - 'BIV-DCOR-SEQ': Sequentially adds ssensors to a reference sensor,
%                    while optimizing its (bivariate) distance correlation 
%                    to that reference through shifting.
%  - 'PCA-CON':      Adds shifting copies of all sensors (except the 
%                    reference), makes a PCA model and selects per sensor 
%                    the most contributing shift.
%  - 'PLS-CON':      Adds shifting copies of all sensors (except the 
%                    reference), regresses the reference sensor on all 
%                    others, and selects per (other) sensor the most 
%                    contributing shift.
%  - 'PLS-SEQ':      Sequentially adds sensors to a growing block of 
%                    reference sensors, while optimizing its shift in terms 
%                    of how well it can by regressed on the reference 
%                    sensors with PLS.
%  - 'RF-CON':       Adds shifting copies of all sensors (except the 
%                    reference), make a Random Forest-model and selects
%                    per sensor the most contributing shift. Contribution
%                    is estimated through out-of-bag predictor importance
%                    by permutation. That is, each predictor is iteratively
%                    permuted to see how much it being randomized hurts or
%                    improves the model.
%
%Some additional settings that can be set as name-value pairs:
%- 'reference': The index of sensor to which all other sensors are lagged 
%  with respect to. The following methods can also be used to
%  auto-determine the reference:
% - 'forward':       Iterates from first to last sensor.
% - 'backward':      Iterates from last to first sensor.
% - 'pearson':       Orders sensors based on sum of absolute correlations 
%                    to all other sensors (without shifting them).
% - 'pearson-shift': Similar to 'pearson', but it will select per sensor
%                    pair the maximum correlation obtained when shifting 
%                    them. This can take quite long.
%  For the sequential methods ('biv-seq' and 'pls-seq'), this function
%  input can also be a vector of indices, which is then the order in which
%  the sensors are added.
%
%- 'contribution': When using the 'pls-con' method for identifying sensor
%  lags, this is the sensor contribution measure used. This can be:
%  - 'loading':     The sum of squared PC loadings, weighted according to 
%                   the fraction of variance explained by the PC.
%  - 'vip':         The variable importance in projetion.
%  - 'coefficient': The absolute regression coefficients (B).
%  For the method 'pca-con', the 'loading' option is always used.
%
%- 'optimisation': Only applicanle for 'pca-con' and 'pls-con'
%  - 'all':         optimise number of PCs in created PCA model with all of
%                   the possible PCs, not recomended.
%  - '#var':        optimise number of PCs in created PCA model with a 
%                   maximum equal to the number of variables, better for
%                   larger datasets (>30 variables) (default)
%  - 'non':         skip optimisation and only use one PC, better for
%                   smaller datasets (<30 variables)
%
%OUTPUT:
%- lags: Will contain the actual lags determined, in the same units as the
%  time axis was inputted.
%- meta: Contains additional information about the lag-identification that
%  is method-specific. Note that the calculation-time only reflects the
%  time of the lag identification, not the reference identification.
%
%AUTHOR INFORMATION:
%- Tim Offermans, Radboud University (tim.offermans@ru.nl), March 2023,
%  tested on Matlab R2021a.

%Check function input:
if nargin<4
    method = 'BIV-PEAR-SEQ';
end
if nargin<3
    error('Some function input is missing, please check description.');
end
reference = 'pearson-shift';
contribution = 'loading';
optimisation = '#var';
for i=1:2:length(varargin)
    eval([varargin{i} ' = varargin{i+1};']);
end
if strcmpi(method, 'PCA-CON')
    contribution = 'loading';
end

%Initialize meta-output:
meta = [];
meta.lags = [];
meta.data = data;
meta.time = time;
meta.maxlag = maxlag;
meta.method = method;
meta.calculationtime = [];
if strcmpi(method, 'PCA-CON') || strcmpi(method, 'PLS-CON')
    meta.contribution = contribution;
end

%Check if the data is sampled at a regular time series:
if length(unique(time(2:end)-time(1:end-1)))~=1
    error('Data is not measured on a regular time frequency, which is a requirement for this function.');
end

tic

%Find the to-be considered lags, and the sample indices that should be
%used for the reference and shifting sensor for each of these lags:
numvars = size(data,2);
if isscalar(maxlag)
    chklags = -maxlag:(time(2)-time(1)):maxlag;
    [~, idref] = intersect(time, intersect(time-maxlag, time+maxlag), 'sorted');
    idsft = repmat(idref, 1, length(chklags)) + repmat(chklags./(time(2)-time(1)), length(idref), 1);
    % Create a cell array that contains copies of chklags, idref and idsft
    % for each variable
    chklags_cell = mat2cell(repmat(chklags,numvars,1), ones(1,numvars));
    idref_cell = mat2cell(repmat(idref,1,numvars), size(idref,1), ones(1,numvars))';
    idsft_cell = mat2cell(repmat(idsft,1,1,numvars), size(idsft,1), size(idsft,2), ones(1,numvars));

elseif iscell(maxlag)
    % Here maxlag is a cell array when each variable has its own allowed 
    % shift. Here we thus need to loop through all variables since each is 
    % different.
    % first we initiate the cell arrays
    chklags_cell = cell(numvars,1);
    idref_cell = cell(numvars,1);
    idsft_cell = cell(1,1,numvars);
    
    % For later calculations it is important that for each variable the
    % time points are of the same length (same number of rows). For
    % this we need to determine the max shift and the min shift of all shifts.
    max_shift = max(cellfun(@max,maxlag));
    min_shift = min(cellfun(@min,maxlag));

    for var = 1:numvars
        chklags = sort(maxlag{var});
        [~, idref] = intersect(time, intersect(time-min_shift, time-max_shift), 'sorted');
        idsft = repmat(idref, 1, length(chklags)) + repmat(chklags./(time(2)-time(1)), length(idref), 1);
        % save the chklags, idref and idsft in the correct spot of the cell
        % arrays
        chklags_cell{var,1} = chklags;
        idref_cell{var,1} = idref;
        idsft_cell{1,1,var} = idsft;  
    end    

    
elseif ismatrix(maxlag)
    chklags = sort(maxlag);
    [~, idref] = intersect(time, intersect(time-min(chklags), time-max(chklags)), 'sorted');
    idsft = repmat(idref, 1, length(chklags)) + repmat(chklags./(time(2)-time(1)), length(idref), 1);
    % Create a cell array that contains copies of chklags, idref and idsft
    % for each variable
    chklags_cell = mat2cell(repmat(chklags,numvars,1), ones(1,numvars));
    idref_cell = mat2cell(repmat(idref,1,numvars), size(idref,1), ones(1,numvars))';
    idsft_cell = mat2cell(repmat(idsft,1,1,numvars), size(idsft,1), size(idsft,2), ones(1,numvars));
end


%Identify reference / order of adding sensors:
if ~isnumeric(reference)
    if strcmpi(reference, 'forward')
        reference = 1:size(data, 2);
    elseif strcmpi(reference, 'backward')
        reference = size(data, 2):-1:1;
    elseif strcmpi(reference, 'pearson')
        reference = abs(corr(data));
        reference = reference.^2;
        reference = reference.*~eye(size(reference, 1));
        reference = sum(reference, 1);
        [~, reference] = sort(reference, 'descend');
    elseif strcmpi(reference, 'pearson-shift')
        reference = zeros(size(data, 2), size(data, 2), length(chklags));
        for isenref = 1:size(data, 2)
            for isensft = setdiff(1:size(data, 2), isenref)
                for ichklags=1:length(chklags)
                    reference(isenref, isensft, ichklags) = abs(corr(data(idref, isenref), data(idsft(:, ichklags), isensft), 'rows', 'complete'));
                end
            end
        end
        reference = reference.^2;
        reference = max(reference, [], 3);
        reference = sum(reference, 1);
        [~, reference] = sort(reference, 'descend');
        clearvars isenref isensft ichklags
    else
        error(['Method for selecting reference sensor (''' reference ''') not recognized.']);
    end
elseif length(reference)==1
    reference = [reference setdiff(1:size(data, 2), reference)];
end
if strcmpi(method, 'PCA-CON') || strcmpi(method, 'PLS-CON')
    reference = reference(1);
end
meta.reference = reference;

%Identify the process lags with either of the methods:
% first assume that the lag is 0 for every variable
lagsid = nan(1,numvars);
for i = 1:numvars
    lagsid(i) = find(~hours(chklags_cell{i}));
end

%Bivariate sequential method:
if strcmpi(method, 'BIV-PEAR-SEQ')
    %Loop through non-reference sensors:
    for isensft = reference(2:end)
        c = NaN(1, length(chklags_cell{isensft}));
        %Find maximum correlation to reference while considering all lags:
        for ichklags = 1:length(chklags_cell{isensft})
            c(ichklags) = abs(corr(data(idref_cell{isensft}, reference(1)), data(idsft_cell{:,:,isensft}(:, ichklags), isensft), 'rows', 'complete'));
        end
        [~, lagsid(isensft)] = max(c);
    end
    clearvars isensft c ichklags

elseif strcmpi(method, 'BIV-DCOR-SEQ')
    %Loop through non-reference sensors:
    for isensft = reference(2:end)
        c = NaN(1, length(chklags_cell{isensft}));
        %Find maximum correlation to reference while considering all lags:
        for ichklags = 1:length(chklags_cell{isensft})
            c(ichklags) = distcorr(data(idref_cell{isensft}, reference(1)), data(idsft_cell{:,:,isensft}(:, ichklags), isensft));
        end
        [~, lagsid(isensft)] = max(c);
    end
    clearvars isensft c ichklags
    
%PCA and PLS contribution-based method (they have similar steps):
elseif strcmpi(method, 'PCA-CON') || strcmpi(method, 'PLS-CON') || strcmpi(method, 'RF-CON')
    %Copy and add each non-reference sensor for each lag:
    SC1 = [];
    SC2 = [];
    for i = 1:numvars
        SC1 = [SC1 ones(size(chklags_cell{i}))*i];
        SC2 = [SC2 1:length(chklags_cell{i})];
    end
    SC = [SC1; SC2];
    SC(:, SC(1, :)==reference(1)) = [];
    SC = [[reference(1); lagsid(reference(1))], SC];
    
    clearvars SC1 SC2
    
    XC = NaN(size(idref_cell{1}, 1), size(SC, 2));
    for i=1:size(SC, 2)
        XC(:, i) = data(idsft_cell{:,:,SC(1, i)}(:, SC(2, i)), SC(1, i));
    end

    %Create calibration and testing set:
    if strcmpi(method, 'PCA-CON') || strcmpi(method, 'PLS-CON')
        test = mod((1:size(XC, 1)), 5)'==0;
    elseif strcmpi(method, 'RF-CON')
        test = zeros(size(XC, 1),1);
    end
    %Mean-center data (on calibration set):
    XC = XC - repmat(mean(XC(~test, :)), size(XC, 1), 1);
    %Autoscale data (on calibration set):
    XC = XC ./ repmat(std(XC(~test, :)), size(XC, 1), 1);
    clearvars i
    %Find per sensor and lag the contribution to a PCA model:
    if strcmpi(method, 'PCA-CON')
        %Find optimal number of PCs:
        if strcmpi(optimisation, '#var')
            E = NaN(1, size(data, 2));
        elseif strcmpi(optimisation, 'all')
            E = NaN(1, size(XC, 2));
        end
        
        [P, ~, V] = pca(XC(~test, :));
        if strcmpi(optimisation, 'non')
            pc = 1;  
        else
            for pc=1:length(E)-1
                T = XC(test, :) * P(:, 1:pc);
                Xp = T * P(:, 1:pc)'; 
                E(pc) = sqrt(mean((XC(test, :) - Xp).^2, [1 2]));
                if pc>5 && ~any(((E(pc-5:pc-1) - E(pc-4:pc)) ./ E(pc-5:pc-1)) > 0.01)
                    break
                end
            end
            [~, pc] = min(E); 
        end
        %Calculate and add contributions per sensor and lag:
        C = sum((P(:, 1:pc).^2) .* (repmat(V(1:pc)'./100, size(P(:, 1:pc), 1), 1)), 2)';
        %C = sum((P(:, 1:pc).^2) .* (repmat(V(1:pc), size(P(:, 1:pc), 1), 1)), 2)';
        SC = [SC; C];
        clearvars E P V T Xp pc C
        
	%Find per sensor and lag the contribution to a PLS model:
    elseif strcmpi(method, 'PLS-CON')
        %Find optimal number of LVs:
        if strcmpi(optimisation, '#var')
            E = NaN(1, size(data, 2)-1);
        elseif strcmpi(optimisation, 'all')
            E = NaN(1, size(XC, 2)-1);
        end
        
        if strcmpi(optimisation, 'non')
            lv = 1;  
        else
            for lv=1:length(E)
                [~, ~, ~, ~, B] = plsregress(XC(~test, 2:end), XC(~test, 1), lv);
                E(lv) = sqrt(mean(((XC(test, 2:end) * B(2:end) + B(1)) - XC(test, 1)).^2));
                if lv>5 && ~any(((E(lv-5:lv-1) - E(lv-4:lv)) ./ E(lv-5:lv-1)) > 0.01)
                    break
                end
            end
            [~, lv] = min(E);
        end
        
        [P, ~, T, U, B, V, ~, stats] = plsregress(XC(~test, 2:end), XC(~test, 1), lv);
        %Calculate and add contributions per sensor and lag:
        if strcmpi(contribution, 'loading')
            C = sum((P.^2) .* repmat(V(1, :), size(P, 1), 1), 2)';
        elseif strcmpi(contribution, 'vip')
            W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
            SS = sum(T.^2, 1) .* sum(U.^2, 1);
            vip = sqrt(size(P, 1) * sum(SS.*(W0.^2), 2) ./ sum(SS, 2));
            C = vip';
        elseif strcmpi(contribution, 'coefficient')
            C = abs(B(2:end))';
        end
        SC = [SC; [1, C]];
        clearvars E lv B P T U V stats C W0 SS vip

    elseif strcmpi(method, 'RF-CON')
        %Train the Random Forest model.
        model = fitrensemble(XC(:, 2:end), XC(:, 1), 'Method', 'Bag');
        %Apply the variable importance estimated (by permutation).
        C = oobPermutedPredictorImportance(model);
        SC = [SC; [1, C]];
        clearvars model C
    end 
    
    %Find per sensor the lag with maximum contribution to either PCA or
    %PLS:
    for i=1:size(data, 2)
        Sp = SC(:, SC(1, :) == i);
        [~, j] = max(Sp(3, :));
        lagsid(i) = Sp(2, j);
    end
    clearvars SC XC i test Sp j
    
%PLS sequential method:
elseif strcmpi(method, 'pls-seq')
    X = [];
    %Loop through non-reference sensors:
    for isensft=2:length(reference)
        X = [X data(idsft_cell{:,:,reference(isensft-1)}(:, lagsid(reference(isensft-1))), reference(isensft-1))];
        %Create mean-centered and autoscaled calibration and testing sets:
        test = mod((1:size(X, 1)), 5)'==0;
        Xpp = X;
        Xpp = Xpp - repmat(mean(Xpp(~test, :)), size(Xpp, 1), 1);
        Xpp = Xpp ./ repmat(std(Xpp(~test, :)), size(Xpp, 1), 1);
        %Minimize optimal PLS prediction error while considering all possible lags:
        E = NaN(length(chklags_cell{reference(isensft)}), size(Xpp, 2));
        for ichklags = 1:size(E, 1)
            Y = data(idsft_cell{:,:,reference(isensft)}(:, ichklags), reference(isensft));
            for lv=1:size(E, 2)
                %PLS can crash is too many LVs are asked for extremely covarying X and Y. This 'try' deals with that.
                try
                    [~, ~, ~, ~, B] = plsregress(Xpp(~test, :), Y(~test, 1), lv);
                    E(ichklags, lv) = sqrt(mean((((Xpp(test, :) * B(2:end))+B(1)) - Y(test, 1)).^2));
                end
                if lv>5 && ~any(((E(lv-5:lv-1) - E(lv-4:lv)) ./ E(lv-5:lv-1)) > 0.01)
                    break
                end
            end
        end
        [lagsid(reference(isensft)), ~] = find(E == min(E(:)), 1);
    end
end

%Get the determine lags in original unit:
lags = nan(1,numvars);
for i = 1:numvars
    lags(i) = hours(chklags_cell{i}(lagsid(i)));
end
lags = hours(lags);
%Copy the identified lags to the meta-output:
meta.lags = lags;
meta.calculationtime = seconds(toc);