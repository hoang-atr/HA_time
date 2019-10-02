function	[post_spike, ypred] = post_spike_search(...
		data, ypred, Model, parm, plot_mode);

if nargin < 5, plot_mode = 0; end

tau  = Model.tau;
a  = Model.a; % spike waveform amplitude
b  = Model.b; % bias
if isfield(Model,'b0')
	b0 = Model.b0; % bias
else
	b0 = b;
end
sx = Model.sx; % noise variance
p = Model.p;

% number of samples in one window
Twin = ceil(parm.Twin*parm.fs);
Tpre = ceil(parm.Tpre*parm.fs);

% time of data [sec]
tsec = data.t; 
% sampling time step
dt   = 1/data.fs;
% estimation time step
dt_est = 1/parm.fs_est; 

% observed data
ydata = data.y;
T = length(ydata);
t = 1:T;

% peak time of spike response
decayRate = 1000;
[tpeak, gpeak, tdecay] = spike_func_peak(tau , decayRate);
tdecay = fix(tdecay*data.fs);

% Model peak amplitude
peak  = gpeak * a(1);
noise = max(b,b0) + sqrt(sx);
% threshold 
threshold = 0.5;
th    = b + peak * threshold;

if noise > th,
	th = noise + peak * threshold;
end

% residual signal
y = (ydata(:) - ypred(:))';

if plot_mode > 0,
	figure;
	
	plot(t,y)
	hold on
	plot([1;T],[th; th], 'c-')
	plot([1;T],[noise; noise], 'y--')
end

% 1. find start and end point over threshold
ton  = find( ([0, y(1:T-1)] < th) & (y >= th));
toff = find( ([y(2:T), 0]   < th) & (y >= th));
Nwin = length(ton);

tlen = toff - ton;

if plot_mode > 0,
	for n= 1:Nwin
		fprintf('length of %d-th win = %d\n', n, tlen(n))
	end
end

% window onset time variable
onset  = zeros(1,Nwin); % onset time of n-th window
offset = zeros(1,Nwin); % end   time of n-th window

% window index
tt    = 0:(Twin-1);

post_spike.twin = cell(Nwin,1);

% 2. set window onset time
for n=1:Nwin
	% start time index of n-th region over threshold
	t1 = ton(n);
	
	onset(n)  = t1 - Tpre;
	offset(n) = onset(n) + Twin - 1;
	
	if onset(n)  < 1,  onset(n) = 1; end;
	if offset(n) > T, offset(n) = T; end;
	
	post_spike.twin{n} = onset(n):offset(n);
end

post_spike.spike_num = zeros(Nwin,1);
post_spike.spike_time = cell(1,Nwin);

data_tk.dt  = dt;% sampling time step [sec]
data_opt.dt = dt;% sampling time step [sec]

for n=1:Nwin
	% time index of this period (n-th overlap window)
	t  = post_spike.twin{n} ; % time index
	ts = tsec(t);   % absolute time [sec]
	% evaluation time
	Tlen = length(t);
	teval = (1:Tlen) + t(1) - 1;
	
	% absolute start/end time of this period [sec]
	t0 = ts(1) - dt;
	te = ts(end) - dt;
	
	tt = (t0:dt_est:te) - t0;
	
	% spike onset pattarn for nspike = 1
	Tspike = tt(:);
	
	% n-th window data
	data_tk.Dt = Tspike;
	% evaluation time
	data_tk.y  = zeros(1,Tlen);

	% --- Evaluate spike response
	% (# of pattarn) x (# of spike)
	g0 = spike_func_evaluate2(tau, data_tk);
	
	switch	 length(a)
	case	1
		g  = a * g0 + b;
	case	2
		g  = a(1) * g0(:,:,1) + a(2) * g0(:,:,2) + b;
	end
	
	% --- Error for spike response and observed data
	err = repadd( g , - y(teval) ) ;
	
	% sum over time
	E  = sum(err.^2, 2);
	% minimum error for spike pattern
	[Emin ,imin]= min(E);
	
	% optimal spike onset pattern
	Ts  = Tspike(imin);
	
	post_spike.spike_num(n)  = 1; 
	post_spike.spike_time{n} = Ts + t0; 
	% absolute timw [sec]
	
	% evaluate prediction response ypred
	% optimal spike onset
	data_opt.Dt = Ts;
	data_opt.y  = zeros(1,Twin + tdecay);
	
	% evaluate optimal spike response
	gopt0 = spike_func_evaluate2(tau, data_opt);
	
	switch	 length(a)
	case	1
		%gopt  = a * gopt0 + b;
        gopt  = a * gopt0 ;
        gopt = gopt + p(1)*(gopt.^2-gopt)+p(2)*(gopt.^3-gopt) + b;    % Huu - Mar 6th 2019
	case	2
		gopt  = a(1) * gopt0(:,:,1) ...
		      + a(2) * gopt0(:,:,2) + b;
	end
	
	tend  = min(T, t(end) + tdecay);
	tpred = t(1):tend;
	Tpred = length(tpred);
	
	% spike response of estimated onsets in this period
	ypred(tpred) = ypred(tpred) + gopt(1:Tpred);
	
end

