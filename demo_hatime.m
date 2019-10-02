clear all
close all

addpath('utility/')
addpath('hatime_func/')
addpath('hb_func/')

%% config simulation parameter
sim_parm.fs_raw = 20000;        % raw sampling rate [Hz]
sim_parm.fs = 50;               % acquired sampling rate [Hz]

sim_parm.p = 1;                 % nonlinearity parameter: 1-'linear', 
                                % <1: 'saturation', >1: 'superlinear'
                                
sim_parm.SNR = 5;               % signal-to-noise ratio

sim_parm.amp = 1;               % spike amplitude
sim_parm.tauOn = 0.01;          % rise time constant [s]
sim_parm.tauOff = 0.2;          % decay time constant [s]

sim_parm.Nspike = 200;          % number of generated spikes
sim_parm.fs_spike = 1;          % firing rate the Poisson spike train
sim_parm.MinISI = 0.1;          % minimum ISI [s]

%% instructions for the data structure
% data.t [1 x T]: sampling time vector
% data.y [1 x T]: fluorescence signal
% data.fs [scalar]: sampling rate
% data.dt [scalar]: sampling interval (=1/fs)
% data.spike_time [Nspike x 1]: ground-truth spike (required for training)
% data.Nspike [scalar]: number of spike time (=numel(spike_time));

% generate training data
train_data = generate_ca_response(sim_parm);

%% Train HA_time
% config the algorithm parameter (cf. the function for detailed explainations)
parm = set_default_param;

% The most influent (and difficult to set) parameters are the ones related
% to sampling the spike candidates, which may vary across conditions (i.e,
% firing rate, sampling rate, nonlinearity). It is thus best to optimize
% those parameters by the training data. For that, just set the 'optimize'
% flag to 1 (default, 0). The HA_time algorithm automatically determines
% the values by maximizing the F1-score for the training data. 
% Warning: it's time consuming and thus may be optimized for each data set
parm.optimize = 0;

% train the HA_time
[classifier, model, parm] = ha_time_train(train_data, parm);

%% Test HA_time
% generate test data
test_data = generate_ca_response(sim_parm);

% apply the HA_time
% the estimated spike train is in the 'est_spike_time' field
test_data = ha_time_test(test_data, classifier, model, parm);

%% evaluate performance
accept_window = 0.05;           % window for accepting 'hit'
sd_cost = 1;                    % cost for computing spike distance (= firing rate)

perf = eval_performance(test_data.spike_time, test_data.est_spike_time, ...
    accept_window, sd_cost);
perf.f1_score
perf.spk_distance
