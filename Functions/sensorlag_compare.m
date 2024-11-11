function difference = sensorlag_compare(lags1, lags2, normalize)

%DESCRIPTION:
%This function will compare two sets of sensor lags, and can be used to
%validate a set of identified lags with the reference lags, if they are
%known.
%
%The difference between the two sets of lags is calculated as follows.
%First, per set all cross-distances between all sensors are calculated.
%Then, per cross-distance the difference between the two sets is
%calculated, and summed. This sum is then divided by the number of
%sensors-1 (because otherwise the displacement of one sensor is counted
%multiple times by its cross-distances to all other sensor), and is then
%(again) divided by the number of sensors (to make the measure relative to
%the number of sensors).
%
%In this set of functions, a lower (negative) lag means that the sensor is
%earlier (upstream) in the process. A higher (positive) lag means that the
%sensor is measured later in the process, and that its measurements are
%thus 'lagging behind' to what is already measured.
%
%INPUT:
%- lags1: The first set of sensor lags to compare.
%- lags2: The second set of sensor lags to compare.
%- order: Whether to normalize the lags of both sets to have a minimum of 0
%  and maximum of 1, before calculating their cross-distances.
%
%AUTHOR INFORMATION:
%- Tim Offermans, Radboud University (tim.offermans@ru.nl), April 2023,
%  tested on Matlab R2021a.

%Check function input:
if nargin<3
    normalize = false;
end

%Ensure both lag sets are row vectors:
lags1 = reshape(lags1, 1, length(lags1));
lags2 = reshape(lags2, 1, length(lags2));

%Normalize the lags from 0 (start) to 1 (end):
if normalize
    lags1 = lags1 - min(lags1);
    lags1 = lags1 ./ max(lags1);
    lags2 = lags2 - min(lags2);
    lags2 = lags2 ./ max(lags2);
end

%Calculate the difference between the two sensor lag sets, as specified
%above:
d1 = repmat(lags1, length(lags1), 1) - repmat(lags1', 1, length(lags1));
d2 = repmat(lags2, length(lags2), 1) - repmat(lags2', 1, length(lags2));
difference = (d1-d2).^2;
difference = sum(difference(:));
difference = difference ./ (length(lags1)-1);
difference = difference ./ length(lags1);
difference = sqrt(difference);

end