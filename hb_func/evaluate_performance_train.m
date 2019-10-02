function	[post_info, spike_info, ypred] = ...
			evaluate_performance_train( ...
			data, Model, parm, post_sw)
%  evaluate_performance
%
%	[post_info, spike_info, ypred] = ...
%		evaluate_performance(data, Model, parm)
%
% post_info.spike_num(n)  = estimated spike number
%                           of n-th overlap window
% post_info.spike_time{n} = estimated onset time
%                           of n-th overlap window
%                           absolute time [sec] 
% post_info.Pspike(ns,n) = P(ns | yk) 
%                          of n-th overlap window
% post_info.twin{n} = time index of n-th overlap window
% post_info.tsec = sample time of observed data [sec]
% post_info.fs   = sampling frequency [Hz]
%
%	spike_info.spike_num 
%	spike_info.spike_time
%
% data.y   = observed spike data with noise 
% data.t   = sample time [sec]
% data.fs  = sampling rate for observed data [Hz]
% data.spike_time  = time of spike onset [sec] 
%                    with (fs_raw) resolution 
%
%	parm.fs_est     % estimation freq [Hz]
%	parm.Twin       % 
%	parm.Tpre       
%	parm.decay   % decay time [sec]
%	             % using estimation spike configulation 
%	
%	parm.threshold 
%	parm.threshold_max
%	
%	parm.max_freq  % max firing frequency [Hz]
%	parm.max_spike
%	parm.Ntau 
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

Model_ini = Model;

% ----- Spike Estimation by EM-algorithm
[Model, state] = ...
		estimate_multispike_em(Model_ini,Data,spike_state);


if isfield(Model_ini,'W')
	Model.W = Model_ini.W;
end

% ----- Evaluate estimation performance

% get onset information from estimated model
switch post_sw
case	1
	% --- estimate posterior of overlap window
	% based on each window estimation
	post_info = estimate_posterior_prob(...
		state,spike_state,overlap,Data);
	
	post_info = estimate_spike_num_weight(post_info,Model);
	
	% number of samples in one window
	post_info.Twin = ceil(parm.Twin*parm.fs);
	
	[pred_info, ypred] = estimate_posterior_spike(...
			data, post_info, Model, parm.fs_est );
	
	[post_spike, ypred] = post_spike_search(...
		data, ypred, Model, parm);
	
	ns1 = length(pred_info.spike_num);
	ns2 = length(post_spike.spike_num);
	
%	fprintf('estimated spike num = %d\n',...
%		sum(pred_info.spike_num))
%	fprintf('estimated spike num = %d (post)\n',...
%		sum(post_spike.spike_num))
	
	post_info.spike_num  = zeros(ns1+ns2,1);
	post_info.spike_time = cell(1,ns1+ns2);
	
	for n = 1:ns1
		post_info.spike_num(n)  = pred_info.spike_num(n);
		post_info.spike_time{n} = pred_info.spike_time{n};
	end
	
	for n = 1:ns2
		post_info.spike_num(n+ns1)  = ...
			post_spike.spike_num(n);
		post_info.spike_time{n+ns1} = ...
			post_spike.spike_time{n};
	end
	
%	fprintf('estimated spike num = %d (total)\n',...
%		sum(post_info.spike_num))

end
