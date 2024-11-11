function [data, time] = processlag_apply(data, time, lags)

%DESCRIPTION:
%This function will shift data according to the sensor lags as identified
%by 'processlagid.m'.
%
%INPUT:
%- data: The data for the sensors to which the lags should be applied.
%- time: The time axis for the sensor data.
%- lags: The lags for each sensor as identified.
%
%OUTPUT:
%- data: The lag-corrected data.
%- time: The lag-corrected time-axis.
%
%AUTHOR INFORMATION:
%- Tim Offermans, Radboud University (tim.offermans@ru.nl), March 2023,
%  tested on Matlab R2021a.

%Check function input:
if nargin<3
	error('Some function input is missing, please check description.');
end

%Get sample indices that should be used for the time-axis, and for each 
%sensor:
[~, idtime] = intersect(time, intersect(time-min(lags), time-max(lags)), 'sorted');
idsft = repmat(idtime, 1, length(lags)) + repmat(lags./(time(2)-time(1)), length(idtime), 1);

%Apply the lags to the data:
newdata = NaN(size(idsft, 1), size(data, 2));
for i=1:size(newdata, 2)
    newdata(:, i) = data(idsft(:, i), i);
end
data = newdata;

%Apply the lags to the time-axis:
time = time(idtime);

end