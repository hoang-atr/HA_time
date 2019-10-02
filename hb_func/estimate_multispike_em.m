function  [Model, state] = estimate_multispike_em(Model,Data,Spike_state)
% Estimate multispike state
%    [Model, state] = estimate_multispike_state(Model,Data,Spike_state)
% --- Model parameters (input: fixed parameter)
% Model.tau
% --- Model parameters (input/output: update parameter)
% Model.a
% Model.b
% Model.sx
% Model.LogPY 
%
% --- extracted window data
% Data.twin  : (# of windows) x (time sample in one window)
% Data.y     : data im windows
% Data.yrest = data outside of windows (rest period)
% Data.trest = time in rest period
% --- possible range of spike state
% Spike_state.Nstate = number of states for (spike number) = nspike
% Spike_state.Tspike{nspike} = spike peak time list of (spike number) = nspike
%
% ----- Output info
% --- estimated spike state for window
% --- histry
% Model.sx_hist = 
% Model.LogPY

print_mode=0;
% max value for log
logMax = log(1.0e+15);

% spike waveform time constants (this value is fixed in this step)
tau  = Model.tau;
Ntau = length(tau);

p = Model.p;    % Huu - Mar 6th 2019

% Iteration number of amplitude & variance update
Nupdate = Model.Nupdate;
Conv    = Model.conv;

% Dim of window data :
% (# of windows) x (# of time sample in one windows)
[Nwin,Nt] = size(Data.y);
Nrest = length(Data.yrest);
Ndata = (Nwin*Nt+Nrest);

if print_mode==1
	fprintf('Nwindow = %d, Ntime = %d, Nrest = %d\n',Nwin,Nt,Nrest)
end

% Nstate(n) = number of states (configuration) for (spike number) = n
Nstate = Spike_state.Nstate;
% Max number of spikes in one window
Nspike = length(Nstate);
% Max number of configuration in spike states
NconfMax = max(Nstate);

if isfield(Model, 'prior')
	% prior histogram for multiple spike state
	% histogram of spike states
	nhist = Model.prior(2:end); % hist(Sk=nspike)
	n_off = Model.prior(1);     % hist(Sk=0)
	nall  = sum(nhist) + n_off;
	
	P00 = n_off/nall; % P0(Sk = 0)
	PS0 = nhist/nall; % P0(Sk = nspike) 
else
	% non-informative prior for multiple spike state
	P00 = 1/(Nspike+1);                     % P0(Sk = 0)
	PS0 = repmat(1/(Nspike+1), [1 Nspike]); % P0(Sk = nspike) 
end

PT0 = 1./Nstate;  % P0(Tk|Sk)


% update history
LogPY = zeros(Nupdate,1);
sx_hist = zeros(Nupdate,1);
a_hist  = zeros(Nupdate,1);
b_hist  = zeros(Nupdate,1);

% statics variable for spike state
stat.P    = zeros(NconfMax,Nspike);
stat.logP = zeros(Nspike,1);
stat.Emin = zeros(Nspike,1);

% Posterior: P(Tk,Sk,Ak=1|yk)
state.P  = zeros(NconfMax,Nspike,Nwin);
state.E  = zeros(NconfMax,Nspike,Nwin);

% Posterior: P(Ak=0|yk)
state.P0 = zeros(Nwin,1);
state.E0 = zeros(Nwin,1);
state.logP0 = zeros(Nwin,1);

% statics value of spike state in k-th window
state.y  = zeros(NconfMax,Nspike,Nwin);

switch	Ntau
case	2
	state.g  = zeros(NconfMax,Nspike,Nwin,1);
	state.gg = zeros(NconfMax,Nspike,Nwin,1);
	state.yg = zeros(NconfMax,Nspike,Nwin,1);
case	3
	state.g  = zeros(NconfMax,Nspike,Nwin,2);
	state.gg = zeros(NconfMax,Nspike,Nwin,3);
	state.yg = zeros(NconfMax,Nspike,Nwin,2);
	a_hist  = zeros(Nupdate,2);
end
%state.dgg = zeros(NconfMax,Nspike,Nwin);
%state.dyg = zeros(NconfMax,Nspike,Nwin);

% mean of window data
%state.y0 = mean([Data.y(:); Data.yrest(:)]);
state.y0 = mean([Data.y(:)]);

% Rest data outside window
if Nrest > 0
	state.yrest = mean([Data.yrest(:)]);
	state.Erest = sum((Data.yrest(:) - state.yrest).^2);
else
	state.yrest = 0;
	state.Erest = 0;
end

% # of data samples
state.NTwin = Nwin*Nt;
state.Nrest = Nrest;
state.Ndata = Ndata;

% log of marginal prob.
state.logPY = zeros(Nwin,1);

% Iteration of amplitude & variance update
for n_itr = 1:Nupdate
	% update spike amplitude & variance
	if n_itr > 1
		Model = update_spike_model(state);
	end
	
	% Model parameters fixed in one iteration
	b   = Model.b;  % bias
	sx  = Model.sx; % noise variance
	% data structure for k-th window 
	data_tk.a = Model.a; % spike waveform amplitude
	data_tk.b = Model.b; % bias
	data_tk.dt = Data.dt;% sampling time step [sec]        
	
	% loop for k-th windows
	for k=1:Nwin
		% k-th window data
		yk = Data.y(k,:);
		
		% Spike state loop
		for nspike = 1:Nspike
			% # of configurations of 'nspike' state 
			Nconf = Nstate(nspike); 
			
			% spike_state.Tspike{nspike} = spike time onset list
			%                            : (# of pattarn) x (# of spike)
			data_tk.Dt = Spike_state.Tspike{nspike};
			data_tk.y  = yk;
			
			% spike state statics
			spike_stat = spike_state_statics(tau, data_tk, p); 
			
			% minimum error for multi-spike pattern
			[Emin ,Imin]= min(spike_stat.E);
			
			stat.Emin(nspike) = Emin;
			
			% Joint prob. for nspike state: 
			% PS0 = P0(Sk = nspike) 
			% PT0 = P0(Tk|Sk) 
			% P(yk,Tk,Sk=nspike) = P(yk|Tk,Sk) P0(Tk|Sk) P0(Sk=nspike) 
			%                    = P(yk|Tk,Sk) * PT0 * PS0
			
			% To avoid numerical overflow, subtract Emin
			logP =  (spike_stat.E - Emin)/(2*sx); % >= 0  (Nconf x 1)
			ix = find( logP < logMax );
			
			P     = zeros(NconfMax,1);
			P(ix) = exp( - logP(ix) ) * PT0(nspike) * PS0(nspike);
			
			% stat.P = P(yk,Tk,Sk=nspike) * exp( Emin/(2*sx) )
			stat.P(:,nspike) = P;
			
			% stat.logP = log(P(yk,Sk)) = log( sum_Tk P(yk,Tk,Sk) )
			stat.logP(nspike) = log( sum(P) ) - Emin/(2*sx);
			
			% statics
			state.E(1:Nconf,nspike,k)  = spike_stat.E  ;
			state.y(1:Nconf,nspike,k)  = spike_stat.y  ;
			
			state.g(1:Nconf,nspike,k,:)  = spike_stat.g  ;
			state.gg(1:Nconf,nspike,k,:) = spike_stat.gg ;
			state.yg(1:Nconf,nspike,k,:) = spike_stat.yg ;
			
%			state.dgg(1:Nconf,nspike,k) = spike_stat.dgg ;
%			state.dyg(1:Nconf,nspike,k) = spike_stat.dyg ;
		end
		% END of Spike state loop (nspike)
		
		% Maximum value of logP(yk,Sk)
		[logPmax ,nspike_opt]= max( stat.logP );
		
		% PSk = exp( logP - logPmax ) = P(yk,Sk) * exp(-logPmax)
		PSk = zeros(Nspike,1);
		
		logPSk = stat.logP - logPmax; % <= 0
		ix = find( - logMax < logPSk );
		PSk(ix) = exp( logPSk(ix) );  % <= 1
		
		% --- Joint prob. for spike on state
		% Pon    = P(yk,Sk>0) = sum_{Sk>0} P(yk,Sk) 
		% logPOn = log( sum_{Sk>0} P(yk,Sk) ) = log( sum(PSk) ) + logPmax
		logPon  = log( sum(PSk) ) + logPmax;
		
		% --- Error for no spike state
		E0 = sum((yk - b).^2);
		% statics for no spike state
		state.E0(k) = E0 ;
		
		% P(yk,Sk=0) = exp( - E0/(2*sx) ) * P0(Sk=0)
		logP0 = - E0/(2*sx) + log(P00);
		
		% --- Joint prob. for spike on/off states
		% Pon   = P(yk,Sk>0) = exp( logPon )
		% P0    = P(yk,Sk=0) = exp( logP0 )
		% --- Marginal prob.
		% P(yk) = Pon + P0 = exp( logPon ) + exp( logP0 )
		logPMAX = max( logPon, logP0 );
		
		% --- Scaled Marginal prob.
		% Psum = P(yk) * exp( -logPMAX )
		%      = exp( logP0 - logPMAX ) + exp( logPon - logPMAX )
		% --- Posterior for No spike state
		% Post0 = P(Sk=0|yk) = P(yk,Sk=0) / P(yk)
		%       = exp(logP0 - logPMAX) / Psum
		Poff  = exp(logP0  - logPMAX);
		Psum  = exp(logPon - logPMAX) + Poff;
		Post0 = Poff/Psum;
		
		% state.P0  = P(Sk=0|yk) = P(yk,Sk=0)/P(yk)
		state.P0(k) = Post0;
		
		% log(P(yk))   = log(Psum) + logPMAX
		state.logPY(k) = log(Psum) + logPMAX;
		
		% --- Posterior for spike state
		% stat.P = P(yk,Tk,Sk=ns) * exp( Emin/(2*sx) )
		% P(Tk,Sk=ns|yk) = P(yk,Tk,Sk=ns)/P(yk)
		%     = stat.P * exp( - Emin(k))/(2*sx) ) / P(yk)
		%     = stat.P * exp( - Emin(k))/(2*sx) - logPMAX)/Psum
		Pcoef = zeros(Nspike,1);
		logc  = stat.Emin/(2*sx) + logPMAX; 
		ix    = find( logc < logMax );
		Pcoef(ix) = exp( - logc(ix) )/Psum;
		
		% state.P = P(Tk,Sk=ns|yk) = P(yk,Tk,Sk=ns)/P(yk)
		state.P(:,:,k)  = repmultiply(stat.P , Pcoef');
	end
	% END of window loop (k-loop)
	
%	state.Erest = sum((Data.yrest(:) - b).^2);
	
	% Marginal log-likelihood  log( P(Y|Model) )
	LogPY(n_itr) = (sum(state.logPY) - state.Erest/(2*sx))/Ndata ...
	             - log(sx)/2;
	
	% history
	sx_hist(n_itr) = Model.sx;
	a_hist(n_itr,:)  = Model.a;
	b_hist(n_itr)  = Model.b;
	
	if n_itr > 1
		LP_dif = (LogPY(n_itr) - LogPY(n_itr-1));
		if (LP_dif > 0) && ( LP_dif < Conv),
			break;
		end
	end
end

if print_mode==1
	fprintf('# of spike_state = %d\n',length(spike_state))
end

Model.LogPY   = LogPY(1:n_itr);
Model.sx_hist = sx_hist(1:n_itr);
Model.a_hist  =  a_hist(1:n_itr,:);
Model.b_hist  =  b_hist(1:n_itr);

Model.tau = tau;
Model.Nupdate = Nupdate;
Model.conv    = Conv;
Model.p = p;
