clear 
close all

addpath('.\Functions');

%% Generate random dataset
% rng(37)
% 
% n_datapoints = 2659;
% sf = 4; %sampling frequency in hours
% n_sensors = 20; %in the paper this was 95, however this can take a while to run. 
% maxshift = 96; %maximum possible shift in datapoints
% 
% data = randn(n_datapoints, n_sensors);
% measurement_times = hours(0:sf:(n_datapoints*sf)-1)';
% 
% possibleshifts = -maxshift:1:maxshift; %in datapoints
% references_shifts = hours(randsample(possibleshifts,n_sensors).*sf);


%% Generate data with idinput()
rng(37)

n_datapoints = 2659;
sf = 4; %sampling frequency in hours
n_sensors = 10; %in the paper this was 95, however this can take a while to run. 
n_periods = 3; %number of periods to use in the sin function
maxshift = 96; %maximum possible shift in datapoints

% Generate data
data = idinput([floor(n_datapoints./n_periods), n_sensors, n_periods], 'rgs', [0 0.5]);

% give data artifical delays
possibleshifts = -maxshift:1:maxshift; %in datapoints
references_shifts = randsample(possibleshifts, n_sensors, true);
data_withshift = circshift(data, references_shifts);

measurement_times = hours(0:sf:(n_datapoints*sf)-2*sf)';
references_shifts = hours((references_shifts.*sf));

% Plot the randomly generated data
figure;
subplot(4,1,1)
plot(measurement_times, data)
xlabel('Time in hours')
ylabel('Y')
title('Generated data with random time shifts applied')

% we now have random data (data) with time points (measurement_times) and a
% vector conating the actual shifts (references_shifts) which we can use for
% validation. 

%% Applying one of our methods
% To apply one of our methods to the full dataset we can use the
% processlag_id() function. Here you can change the method to the desired
% one. 
Method = 'PLS-CON'; % BIV-PEAR-SEQ, BIV-DCOR-SEQ, PCA-CON, PLS-CON, RF-CON, PLS-SEQ

% possible_shifts = -hours(maxshift*sf):hours(sf):hours(maxshift*sf);
% predictedshifts = processlag_id(data, measurement_times, possible_shifts, Method);

predicted_shifts = processlag_id(data, measurement_times, hours(maxshift*sf), Method);
% If the method is PLS-CON you can set the contribution metric with the
% 'contribution' input (default='loading'). If the method is PCA-CON or PLS-CON 
% you can set the 'optimisation' input (default='#var')
% NOTE: this can take a while to run, especially for RF-CON, BIV-DCOR-SEQ 
% PLS-SEQ and if there are many sensors

%% Validation
% With the predicted shifts and the reference shifts we can calculate the 
% Time Delay Estimation error using the sensorlag_compare() function
TDEerror = sensorlag_compare(hours(predicted_shifts), hours(references_shifts), false);
disp(['The TDE error obtained with ', Method ' is: ', num2str(TDEerror)])

% We can also correct the estimated shifts in the original data by using
% the processlag_apply() function. 
[data_corrected, time_corrected] = processlag_apply(data, measurement_times, predicted_shifts);

% we can compare the data with shifts, the data where the shifts have been
% perfectly corrected and the data corrected by the shifts estimated by one
% of our methods:
subplot(4,1,2)
[data_corrected_ref, time_corrected_ref] = processlag_apply(data, measurement_times, references_shifts);
plot(time_corrected_ref, data_corrected_ref)
xlabel('Time in hours')
ylabel('Y')
title('Shifts perfectly corrected (using reference shifts)')

subplot(4,1,3)
plot(time_corrected, data_corrected)
xlabel('Time in hours')
ylabel('Y')
title(['Shifts corrected with ', Method])

%% Clustering-based approach
% in the paper a clustering based approach was described for estimating the
% shifts. We can apply this by using the processlag_id_ClusteringBased()
% function.
% IMPORTANT NOTE: WHILE PROCESS_ID TAKES THE MAXLAG IN THE UNIT OF THE MEASURMENT
% TIMES, PROCESS_ID_CLUSTERINGBASED TAKES THE MAXSHIFT IN DATAPOINTS!!!!!
predictedshifts_cluster = processlag_id_ClusteringBased(data, measurement_times, maxshift, Method);

TDEerror_clust = sensorlag_compare(hours(predictedshifts_cluster), hours(references_shifts), false);
disp(['The TDE error obtained with clustering-based ', Method ' is: ', num2str(TDEerror_clust)])

subplot(4,1,4)
[data_corrected_clust, time_corrected_clust] = processlag_apply(data, measurement_times, predictedshifts_cluster);
plot(time_corrected_clust, data_corrected_clust)
xlabel('Time in hours')
ylabel('Y')
title(['Shifts corrected with clustering-based ', Method])