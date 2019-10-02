function	data = generate_spike_signal(parm, spike_time)
% ----- Poisson spike signal generation 
%    data = generate_spike_signal(parm, spike_time)
% ----- Time of Poisson spike 
% spike_time = time of spike onset [sec]  (Nspike x 1)
% ----- spike generation setting
% parm.tau      = [onset time; decay time] [sec]
% parm.fs       = sampling rate for observed data [Hz]
% parm.SN       = SN ratio of noise
% ----- optional
% parm.fs_raw   = sampling rate for raw signal [Hz]
% parm.normalize_on = on/off of 
%                     normalization to [0 1] range [ = 1]
% --- Output
% data.y   = observed spike data with noise 
% data.t   = sample time [sec]
% data.fs  = sampling rate for observed data [Hz]
% data.spike_time  = time of spike onset [sec] 
%                    with (fs_raw) resolution 
%
% data.y_org      = clean spike data without noise
% spike amplitude & bias
% data.a  = amplitude;
% data.b  = bias;
% data.sx = noise variance;

fprintf('----- Spike signal generation\n')

plot_mode=0;
if isfield(parm,'plot_mode'),
	plot_mode = parm.plot_mode;
end;

% fs_raw   = sampling rate for raw signal [Hz]
fs_raw = 20000;
if isfield(parm,'fs_raw'), fs_raw = parm.fs_raw ;end;

% normalization on/off
normalize_on = 0;
if isfield(parm,'normalize_on'), 
	normalize_on = parm.normalize_on; 
end

% peak time of spike response
decayRate = 1000;
[tpeak, gpeak, tdecay] = ...
	spike_func_peak(parm.tau , decayRate);

% parm.tau = [0.01; 0.5; 1.0]; 
%            [onset ; decay_small; decay_long] [sec]
% amplitude of two decay component
% parm.amp = [1.0; 0.2];   
% parm.fs  = sampling rate for observed data [Hz]

Ntau = length(parm.tau);

fprintf('sampling freq = %8.3f [Hz]\n', parm.fs)
fprintf('spike response time constants ');

switch	Ntau
case	2
	fprintf('= [%8.3f, %8.3f] [sec]\n',parm.tau)
case	3
	fprintf('= [%8.3f, %8.3f, %8.3f] [sec]\n',parm.tau)
end

% 1. time index of spike onset (int) at fs_raw
Nspike = length(spike_time);
% add decay time for 1st spike
spike_time = spike_time + tdecay; 

spike_id = ceil(spike_time * fs_raw);
peaks_id = ceil((spike_time + tpeak) * fs_raw);
Tmax   = max(peaks_id) + fix(tdecay*fs_raw);

% 2. spike time with sample rate = fs_raw

% time step for observed data [sec]
dt = 1/parm.fs;
% time step for raw signal [sec]
dt_raw = 1/fs_raw;

% time sample for raw signal
% signal at t_raw(m) = (m * dt_raw) 
t_raw = (1:(Tmax))*dt_raw;

% time of spike onset with (dt_raw) resolution [sec] 
spike_time  = t_raw(spike_id);      % = m * dt_raw
% time of spike onset with dt resolution [sec] 
spike_flame = ceil(spike_time/dt);  
% n = m * (dt_raw/dt)

% 3. spike response curve with fs_raw [Hz]

% time sample with fs_raw [Hz] 
tt = 0:dt_raw:tdecay;
T  = length(tt);

data0.y  = zeros(1,T);
data0.Dt = 0;
data0.dt = dt_raw;

%  tau = [É—1; É—2]    or  [É—1; É—2; É—3]   (T2 < T3)
%  g   : (Nconf x Nt)  or  (Nconf x Nt x 2)
g  = spike_func_evaluate2(parm.tau, data0);

switch	Ntau
case	2
	g0 = parm.amp * g;
case	3
	g0 = parm.amp(1) * g(:,:,1) + parm.amp(2) * g(:,:,2);
end

gpeak = max(g0(:));

% 4. generation of spike response (# of spike = Nspike)
% time length of spike response at sampling freq fs [Hz]
tmax = ceil(tdecay/dt) + 1;

% spike onset time from flame start time
Dts = spike_time - (spike_flame*dt - dt);

data_n.Dt = Dts(:);
data_n.dt = dt;
data_n.y = zeros(1,tmax);

%  tau = [É—1; É—2]    or  [É—1; É—2; É—3]   (T2 < T3)
%  g   : (Nconf x Nt)  or  (Nconf x Nt x 2)
g = spike_func_evaluate2(parm.tau, data_n);

% 5. generation of spike signal

% sample length at fs_raw
Tmax  = Tmax + T; 
% sample length at fs
Nsample = fix(Tmax*dt_raw/dt);

y = zeros(1,Nsample);

% add each spike response
for n=1:Nspike
	tn = (1:tmax) + spike_flame(n) - 1;
	
	switch	Ntau
	case	2
		y(tn) = y(tn) + parm.amp * g(n,:);
	case	3
		y(tn) = y(tn) + parm.amp(1) * g(n,:,1)...
	                  + parm.amp(2) * g(n,:,2);
	end
	
end

% 6. noise generation
%    SNR = y_max / SD_noise

noise = (gpeak/parm.SN) * randn(1, Nsample);

% 7. make observed data
t = (1:Nsample)*dt;

ydata = y + noise;

% 8. normalization to [0 1]
%	 ydata = (ydata - ymin)/(max(ydata) - ymin) ;
if normalize_on == 1
	ymin  = min(ydata);
	scale = 1/(max(ydata) - ymin);
	ydata = (ydata - ymin) * scale ;
	noise = (noise - ymin) * scale;
	bias  = mean( noise );
	sx    = mean((noise-bias).^2);
else
	% variance of observation noise
	scale = 1;
	bias  = mean(noise);
	sx    = mean((noise-bias).^2);
end

% spike amplitude & bias
data.a  = scale * parm.amp;
data.b  = bias;
data.sx = sx;

fprintf('number of spikes = %d\n',Nspike)
fprintf('number of samples = %d\n',length(t))

% ---- Output variable
data.y = ydata;
data.t = t;
data.y_org = y;
data.fs = parm.fs;

% ----- spike info
% data.spike_time = time of spike onset [sec] 
%                   with (fs_raw) resolution 
data.spike_time = spike_time;

if plot_mode==0, return; end

fs_txt = sprintf('[fs=%5.1f]',parm.fs);

if isfield(parm,'fs_spike')
	if length(parm.fs_spike) > 1
		fmax = max(parm.fs_spike);
		fmin = min(parm.fs_spike(parm.fs_spike>0));
		fire_txt = sprintf('[firing=%5.1f -%5.1f]',...
			fmax,fmin);
	else
		fire_txt = sprintf('[firing=%5.1f]',parm.fs_spike);
	end
else
	fire_txt = [];
end

if plot_mode == -1
	figure
	plot(Dt,g0,'lineWidth',3)
	hold on
	xlabel('Time [sec]')
	title('Spike waveform')
end
	
NX = 2;
NY = 2;
nspike = min(15,Nspike);
T1  = spike_time( nspike );

ym  = mean(ydata);
ysd = std(ydata);
th  = ym + ysd;

if plot_mode==1
	figure;
	plot(tt,g0)
	hold on
	xlim([0 tt(end)])
	title(sprintf('Spike waveform [fs(raw data)=%6.1f]',...
		fs_raw))
	xlabel('Time [sec]')
	
	figure;
	plot(t,ydata)
	xlim([0 T1])
	title_txt = sprintf('Spike signal [SN=%5.1f]',parm.SN);
	title([title_txt, fs_txt, fire_txt])
end

NX = 1;
NY = 2;
T  = spike_time( end);

if plot_mode == 2
	figure;
	subplot(NY,NX,1)
	plot(t,y)
	xlim([0 T])
	title_txt = ['Clean Spike signal '];
	title([title_txt, fs_txt, fire_txt])
	
	subplot(NY,NX,2)
	plot(t,ydata)
	hold on 
	plot([0 T], [ym ym], 'r-')
	plot([0 T], [th th], 'r--')
	xlim([0 T])
	title(sprintf('Spike signal with noise [SN=%5.1f]',...
		parm.SN))

	figure
	plot(t,y)
	hold on
	
	for n=1:nspike
		%t0 = spike_time(n);
		tn = (1:tmax) + spike_flame(n) - 1;
		
		switch	Ntau
		case	2
			plot(t(tn), g(n,:),'r--')
		case	3
			plot(t(tn), g(n,:,1),'r--')
			hold on
			plot(t(tn), g(n,:,2),'r--')
		end
		
	end
	xlim([0 T1])
	title(['Spike signal ', fs_txt, fire_txt ])
end
