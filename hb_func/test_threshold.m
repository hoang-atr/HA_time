function	threshold_max = test_threshold(...
				data, Model, spike_mode, plot_mode)
%
if nargin < 3, spike_mode  = 2; end;
if nargin < 4, plot_mode  = 0; end;

print_mode = 0;

% test_threshold
% ---- Extract spike window changing threshold

% Model
tau  = Model.tau;
a  = Model.a; % spike waveform amplitude
b  = Model.b; % bias
if isfield(Model,'b0')
	b0 = Model.b0; % bias
else
	b0 = b;
end
sx = Model.sx; % noise variance

[tpeak, gpeak] = spike_func_peak(tau);
% Model peak amplitude
peak  = gpeak * a(1);
noise = max(b,b0) + sqrt(sx);

% mean & STD
ym = mean(data.y);
sd = std(data.y);

threshold  = 0.3:0.2:2.5;
threshold = threshold(:);
NT = length(threshold);

TP = zeros(NT,1);
FP = zeros(NT,1);
% amplitude ratio
AP = zeros(NT,1);

Nspike = length(data.spike_time);

% ------------------------------ 
% Basic spike parameter
% ------------------------------ 
parm.fs = data.fs;
parm.plot_mode = 0;

parm.max_freq = 10;  % max firing frequency [Hz]
parm.decay    = 1; % decay time [sec]
	               % using estimation spike configulation 

% ----- set test parameter setting
parm = set_spike_parm_cs(spike_mode, parm);

% ---- Extract spike window changing threshold
for n=1:NT
	% threshold value
	th = ym + sd*threshold(n);
	% amplitude ratio
	% th = noise + peak * AP(n);
	AP(n) = (th - noise) / peak;
	
	parm.threshold = threshold(n);
	parm.threshold_max = threshold(n);
	
	% ----- Extract_spike_window;
	Data = extract_spike_window(parm, data.y, print_mode);
	
	% ----- make multiple_spike_state for estimation
	spike_state = multiple_spike_state(parm, print_mode);
	
	% ----- get spike_information in extracted window 
	overlap = find_ovelap_window(...
				data, Data, spike_state, print_mode);
	
	spike_info = spike_overlap_info(...
				data, overlap, print_mode);
	
	% spike number in windows
	TP(n) = sum(spike_info.spike_num);
	
	% number of windows without spike
	iz =  find( spike_info.spike_num == 0);
	FP(n) = length(iz);
end

% noise
noise_id = find( AP < 0);
if ~isempty(noise_id), noise_id = noise_id(end);end;
amp_id = find( AP > 0.5);
if ~isempty(amp_id), amp_id = amp_id(1);end;

id = find(TP >= Nspike);

if isempty(id)
	max_id = [];
else
	max_id = id(end);
end

if isempty(max_id)
	if ~isempty(amp_id),
		threshold_max = threshold(amp_id);
	else
		threshold_max = threshold(NT - 5);
	end
end

if ~isempty(max_id)
	if ~isempty(amp_id),
		threshold_max = threshold( max( amp_id, max_id) );
	else
		threshold_max = threshold(max_id);
	end
end

if threshold_max > 1
	threshold_max = threshold_max - 0.1;
end

if plot_mode==0, return; end

if isempty(max_id)
	fprintf('empty max_id \n')
else
	fprintf('max_id = %d\n',max_id)
end

if isempty(amp_id),
	fprintf('empty amp_id \n')
else
	fprintf('amp_id = %d\n',amp_id)
end

figure;
plot(threshold,TP,'-')
hold on
plot([threshold(1); threshold(NT)], ...
	[Nspike;Nspike],'--r')

plot(threshold,FP,'--c')

plot([threshold_max; threshold_max], ...
	[0;Nspike],'--r')

if ~isempty(noise_id)
	plot([threshold(noise_id); threshold(noise_id)], ...
	[0;2*Nspike],'-c')
end

if ~isempty(amp_id)
	plot([threshold(amp_id); threshold(amp_id)], ...
	[0;2*Nspike],'-r')
end

	
%ylim([0 3*Nspike])

