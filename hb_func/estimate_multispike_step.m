function  [Model,data] = estimate_multispike_step(...
						Model,Data,Spike_state)
% Estimate multispike state
%    [Model,data] = estimate_multispike_step(...
%					Model,Data,Spike_state)
% --- Model parameters (input: fixed parameter)
% Model.tau
% --- Model parameters (input/output: updated parameter)
% Model.a
% Model.b
% Model.b0
% Model.sx
% Model.LogPostY 
%
% --- extracted window data
% Data.twin  : (# of windows) x (time sample in one window)
% Data.y     : data im windows
% Data.yrest = data outside of windows (rest period)
% Data.trest = time in rest period
% --- possible range of spike state
% Spike_state.Nstate = number of states for (spike number) = nspike
% Spike_state.Tspike{nspike} = time list of (spike number) = nspike
% Spike_state.Dt{nspike} = t - Tspike{nspike}
%   Dt{nspike} : (# of pattarn) x (# of time sample) x (# of spike)
%
% ----- Output info
% --- estimated spike state for window
% Model.spike_state = number of spike
% Model.ix_spike = spike on index
% Model.ix_rest  = spike off index
% Model.Tspike = spike timing list
% --- histry
% Model.sx_hist = 
% Model.LogPostY_hist
% Model.sx_ini

print_mode=0;
% max value for log
logMax = log(1.0e+10);

% spike waveform time constants (this value is fixed in this step)
tau = Model.tau;
Model.sx_ini = Model.sx;

% Dim of window data : (# of windows) x (# of time sample in one windows)
[Nwin,Nt] = size(Data.y);
Nrest = length(Data.yrest);

if print_mode==1
	fprintf('Nwindow = %d, Ntime = %d, Nrest = %d\n',Nwin,Nt,Nrest)
end

% Nstate(n) = number of states for (spike number) = n
Nstate = Spike_state.Nstate;
% Max number of spikes in one window
Nspike = length(Nstate);

% Iteration number of amplitude & variance update
Nupdate = Model.Nupdate_shape;

% non-informative prior for multiple spike state
P0  = 1/(Nspike+1); % P0(nspike)
PS0 = 1./Nstate;  % P0(Tspike{nspike}|nspike)

%  Dt{n} : spike onset in n-th window
Dt = cell(Nwin,1);

% working variable
spike_state = zeros(Nwin,1);
spike_id = zeros(Nwin,Nspike);
logPost  = zeros(Nspike,1);
logPostY = zeros(Nwin,1);
LogPostY = zeros(Nupdate,1);

% update history
sx_hist = zeros(Nupdate,1);
sx_old  = Model.sx;

model_hist = cell(Nupdate,1);
state_hist = cell(Nupdate,1);

% Iteration of amplitude & variance update
for n_itr = 1:Nupdate
	% -------------------------------------------
	% ---- Optimal spike state estimate step
	% -------------------------------------------
	
	% Model parameters fixed in one iteration
	tau = Model.tau;
	b0  = Model.b0; % bias of rest state
	sx  = Model.sx; % noise variance
	
	% data structure for k-th window error
	data_tk.a = Model.a; % spike waveform amplitude
	data_tk.b = Model.b; % bias of spike state
	data_tk.dt = Data.dt;% sampling time step [sec]
	
	% loop for windows
	for k=1:Nwin
		% k-th window data
		yk = Data.y(k,:);
		
		% log-likelihood of spike state
		for nspike = 1:Nspike
			% # of configurations of 'nspike' state 
			Nconf = Nstate(nspike); 
			
			% spike_state.Tspike{nspike} = 
			%       spike time onset list
			%     : (# of pattarn) x (# of spike)
			data_tk.Dt = Spike_state.Tspike{nspike};
			data_tk.y  = yk;
			
			% error for multiple spike state 
			Evec = sq_error_spike_vec(tau, data_tk); 
			% (Nconf x 1)
			
			% find best multi-spike pattern
			% that minimize Evec
			[Emin ,Imin]= min(Evec);
			spike_id(k,nspike) = Imin;
			
			%     - (log-likelihood) 
			logP = (Evec - Emin)/(2*sx); % >= 0  (Nconf x 1)
			logP = min( logP , logMax );
			
			% Posterior for nspike state: P(nspike|yk)
			% P0  = 1/(Nspike+1); % P0(nspike)
			% PS0 = 1./Nstate;  % P0(Tspike{nspike}|nspike)
			P     = exp( - logP ) * PS0(nspike) * P0;
			Post  = sum(P);
			logPost(nspike) = log(Post) - Emin/(2*sx);
		end
		
		% Find optimal spike number
		[logPmax ,nspike_opt]= max(logPost);
		
		% PostOff = Post(off|Y) = exp(logPostOff)
		% PostOn  = Post(on|Y)  = sum_{nspike} Post(nspike|Y)
		logPost = logPost - logPmax; % <= 0
		logPost = max( logPost , - logMax );
		PostOn  = exp( logPost );
		
		logPostOn = log(sum(PostOn)) + logPmax;
		
		% error for no spike state
		E0 = sum((yk - b0).^2);
		logPostOff = - E0/(2*sx) + log(P0);
		
		% P(yk|Model) = P(On|yk) + P(Off|yk)
		logPostY(k) = log( (exp(logPostOn) ...
					+ exp(logPostOff)) );
		
%		if logPostOn >= logPostOff,
		if logPmax >= logPostOff,
			% Spike on state is selected
			
			% estimated number of spikes
			spike_state(k) = nspike_opt;     
			
			% best spike configulation index
			n_best = spike_id(k,nspike_opt); 
			
			% Spike timing list
			% spike_state.Tspike{nspike} = 
			%   spike time onset list
			%  : (# of pattarn) x (# of spike)
			Dt{k} = Spike_state.Tspike{nspike_opt}(n_best,:);
		else
			% No spike
			spike_state(k) = 0;
			Dt{k}   = [];
		end
	end
	
	% Rest data outside window
	if Nrest==0
		Erest = 0;
	else
		Erest = sum((Data.yrest(:) - b0).^2);
	end
	
	% Marginal log-likelihood  log( P(Y|Model) )
	LogPostY(n_itr) = sum(logPostY) - Erest/(2*sx) ...
	                - (Nwin*Nt+Nrest) * log(sx)/2;
	
	% -------------------------------------------
	% ---- Set optimal spike state
	% -------------------------------------------
	
	% set spike on/off index
	ix_rest  = find( spike_state == 0);
	ix_spike = find( spike_state >  0);
	
	% spike state data
	data.dt = Data.dt;% sampling time step [sec]
	data.y  = Data.y(ix_spike,:);
	data.Dt = cell(length(ix_spike),1);
	for m = 1:length(ix_spike)
		data.Dt{m} = Dt{ix_spike(m)};
	end

	% -------------------------------------------
	% ---- Parameter estimate step using optimal spike state
	% -------------------------------------------
	% update amplitude
	[Model.a, Model.b] = ...
		update_spike_amplitude(tau, data);
	
	% no spike state data
	data.yrest = Data.y(ix_rest,:);
	if Nrest > 0
		data.yrest = [data.yrest(:); Data.yrest(:)];
	end
	
	% update noise variance
	data.a   = Model.a;
	data.b   = Model.b;
	[Model.sx, Model.b0] = ...
		update_model_variance(Model.tau, data);
	
	% noise variance history
	sx_hist(n_itr) = Model.sx;
	
	sx_dif = (sx_old - Model.sx);
	if (sx_dif > 0) && ( sx_dif < Model.conv),
		break;
	end
	sx_old = Model.sx;
end

% estimated number of spikes
% spike_state(k) = nspike_opt;     

Model.spike_num = spike_state;

% 95 percentile of spike number
[spnum, indx] = sort(spike_state);
pos = fix(length(spike_state) * 0.95);
Model.max_spike = spnum(pos);

sx_hist  = sx_hist(1:n_itr);
LogPostY = LogPostY(1:n_itr);

Model.sx_hist  = sx_hist;
Model.LogPostY = LogPostY;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
