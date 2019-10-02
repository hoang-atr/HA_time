function	[post_info, spike_info] = evaluate_posterior( ...
			data, Model, parm)
%  evaluate posterior
%
%	[post_info, spike_info] = evaluate_posterior(data, Model, parm)
%
% --- Estimated spike info for overlap window
% post_info.Pspike(ns,n) = P(ns | yk) of n-th overlap window
% post_info.twin{n} = time index of n-th overlap window
% post_info.tsec = sample time of observed data [sec]
% post_info.fs   = sampling frequency [Hz]
%
% post_info.spike_num(n)  = estimated spike number of n-th overlap window
% ----- spike_info
%	spike_info.spike_num 
%	spike_info.spike_time
%
% data.y   = observed spike data with noise 
% data.t   = sample time [sec]
% data.fs  = sampling rate for observed data [Hz]
% data.spike_time  = time of spike onset [sec] with (fs_raw) resolution 
%
%	parm.fs_est  = 100;   % estimation freq [Hz]
%	parm.Twin    = 1.5;   % (15) sample
%	parm.Tpre    = 0.4;   
%	parm.decay    = 1; % decay time [sec]
%	                   % using estimation spike configulation 
%	
%	parm.threshold = 1;
%	parm.threshold_max = 1.5;
%	
%	parm.max_freq = 10;  % max firing frequency [Hz]
%	parm.max_spike = 2;
%	parm.Ntau = 3;
%
%	Model.tau 
%	Model.a   
%	Model.b   
%	Model.sx  
%	Model.b0
			
ypred = [];
spike_info = [];

if nargin < 4, post_sw = 1; end;

% ----- Extract_spike_window;
Data = extract_spike_window(parm, data.y);

% ----- make multiple_spike_state for estimation
spike_state = multiple_spike_state(parm);

% ----- get spike_information in extracted window 
overlap = find_ovelap_window(data, Data, spike_state);

if isfield(data,'spike_time')
	spike_info = spike_overlap_info(data, overlap);
end

% ----- Spike Estimation by EM-algorithm
[Model, state] = ...
		estimate_multispike_em(Model,Data,spike_state);

% ----- Evaluate posterior

% --- estimate posterior of overlap window based on each window estimation
post_info = estimate_posterior_state(state,spike_state,overlap,Data);
