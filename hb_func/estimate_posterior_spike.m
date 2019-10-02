function	[pred_info, ypred] = ...
			estimate_posterior_spike(...
			data, post_info, Model, fs_est)
% --- spike estimation based on posterior for overlap window 
%   pred_info = estimate_posterior_spike(...
%				data, post_info, Model, fs_est)
% ----- Output
% pred_info.spike_num(n)  = estimated spike number 
%                           of n-th overlap window
% pred_info.spike_time{n} = estimated onset time 
%                           of n-th overlap window
%                           absolute time [sec] 
% pred_info.twin{n} = time index of n-th overlap window
% ----- Input
% fs_est : estimation frequency [Hz]
% --- observed data
% data.y   = observed spike data with noise 
% data.t   = sample time [sec]
% data.fs  = sampling rate for observed data [Hz]
% --- Estimated spike info for overlap window
% post_info.spike_num(n)  = estimated spike number 
%                           of n-th overlap window
% post_info.twin{n} = time index of n-th overlap window
% --- Model parameters
% Model.tau
% Model.a
% Model.b
% Model.sx

tau  = Model.tau;
a = Model.a; % spike waveform amplitude
b = Model.b; % bias
p = Model.p;

% # of overlap window
Nwin = length(post_info.twin);

% estimation window length
Twin = post_info.Twin;
pred_info.Twin = Twin;

% time of data [sec]
tsec = data.t; 
% sampling time step
dt   = 1/data.fs;
% estimation time step
dt_est = 1/fs_est; 

% observed data
ydata = data.y;
T = length(ydata);
ypred = zeros(1,T);

% peak time of spike response
decayRate = 1000;
[tpeak, gpeak, tdecay] = spike_func_peak(tau , decayRate);
tdecay = fix(tdecay*data.fs);

data_tk.dt  = dt;% sampling time step [sec]
data_opt.dt = dt;% sampling time step [sec]

pred_info.twin = post_info.twin;
pred_info.spike_num = post_info.spike_num;
pred_info.spike_time = cell(1,Nwin);

for n=1:Nwin
	% time index of this period (n-th overlap window)
	t  = post_info.twin{n} ; % time index
	ts = tsec(t);   % absolute time [sec]
	tw = length(t); % length of period
	if tw < Twin, tw = Twin; end
	% evaluation time
	teval = (1:tw) + t(1) - 1;
	
	% absolute start/end time of this period [sec]
	t0 = ts(1) - dt;
	te = ts(end) - dt;
	%	% find spike within this period
	%	ix = find( (Tk(nc,:) >= (ts(1) - dt)) ...
	%			&  (Tk(nc,:) < ts(end)));
	
	% estimated spike number
	nspike = post_info.spike_num(n);
	
	% if no spike state, skip to next
	if nspike==0, continue; end;
	
	tt = (t0:dt_est:te) - t0;
	
	if nspike > 1
		[t1,t2] = meshgrid(tt);
		t1 = t1(:);
		t2 = t2(:);
		ix = find(t1 < t2);
		t1 = t1(ix);
		t2 = t2(ix);
	end
	
	if nspike == 1,
		Tspike = tt(:);
	elseif nspike == 2
		Tspike = [t1, t2];
	else
		Nconf  = length(t1);
		Tspike = zeros(Nconf,nspike);
		for m = 1:Nconf
			tlist = 0:(nspike-1);
			tstep = (t2(m) - t1(m))/(nspike-1);
			Tspike(m,:) = tlist * tstep + t1(m);
		end
	end
	
	Nconf = size(Tspike,1);
	
	% n-th window data
	data_tk.Dt = Tspike;
	% evaluation time
	data_tk.y  = zeros(1,tw);

	% --- Evaluate spike response
	% (# of pattarn) x (# of spike)
	g0 = spike_func_evaluate2(tau, data_tk);
	
	switch	 length(a)
	case	1
		%g  = a * g0 + b;
        g  = a * g0;
        g = g + p(1)*(g.^2-g)+p(2)*(g.^3-g) + b;    % Huu - Mar 6th 2019
	case	2
		g  = a(1) * g0(:,:,1) + a(2) * g0(:,:,2) + b;
	end
	% subtract previous spike contribution
	dy  = ydata(teval) - ypred(teval) ; 
	
	% temporal difference of waveform
	%	dy  = diff([0 dy]);
	%	dg  = diff([zeros(size(g,1),1) g],1, 2);
	
	% --- Error for spike response and observed data
	err = repadd( g , - dy ) ;
	
	% sum over time
	E  = sum(err.^2, 2);
	% minimum error for multi-spike pattern
	[Emin ,imin]= min(E);
	
	% error for no spike state
	E0 = sum(dy.^2);
	
	if Emin > E0
		% --- No spike state is less error
		% modify spike_num 
		pred_info.spike_num(n) = 0;
	else
		% --- Spike state is less error
		% optimal spike onset pattern
		Ts  = Tspike(imin,:);
		pred_info.spike_time{n} = Ts + t0; 
		% absolute timw [sec]
		
		% optimal spike onset
		data_opt.Dt = Ts;
		data_opt.y  = zeros(1,tw + tdecay);
	
		% evaluate optimal spike response
		gopt0 = spike_func_evaluate2(tau, data_opt);
		
		switch	 length(a)
		case	1
			%gopt  = a * gopt0 + b;
            gopt  = a * gopt0;
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
end

return

post_spike = post_spike_search(...
		data, pred_info, ypred, Model, fs_est);
