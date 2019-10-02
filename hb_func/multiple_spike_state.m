function	spike_state = multiple_spike_state(...
			parm, print_mode)
% multiple_spike_state
%   spike_state = multiple_spike_state(parm)
%
% spike_state.Nstate(nspike) = number of states
%                              for (# of spike) = nspike
% spike_state.Tspike{nspike} = spike time onset list 
%                              of (# of spike) = nspike
%                            : (# of pattarn) x (# of spike)% --- Input parameters
% parm.fs        :  sampling rate for observed data [Hz]
% parm.fs_est    : estimation freq [Hz]
% parm.Twin      :  time window length [sec]
% parm.decay     : decay time [sec]
%                  using estimation spike configulation 
% parm.max_spike : Max number of spikes in one window
% parm.max_freq  : Max frequency of spikes 

if nargin<2, print_mode=1; end

plot_mode=0;
if isfield(parm,'plot_mode'),plot_mode = parm.plot_mode;end;
% Max number of spikes in one window
if isfield(parm,'max_spike')
	max_spike = fix(parm.max_spike); 
else
	max_spike = 2;
end

% Max frequency of spikes 
if isfield(parm,'max_freq')
	max_freq = parm.max_freq; 
else
	max_freq = 10;
end

% minimum ISI
MinISI = 1/max_freq;

% parm.fs      =  sampling rate for observed data [Hz]
% parm.fs_est  =  estimation freq [Hz]
% time step
dt      = 1/parm.fs;
dt_est  = 1/parm.fs_est;

%decayRate = 100;
%[tpeak, gpeak, decay] = spike_func_peak(tau, decayRate);
decay = parm.decay;

% parm.Twin      :  time window length [sec]
Twin = parm.Twin; 

%% number of time samples
%N  = fix(Twin*parm.fs);
%% ----- time samples of data in a window
%t  = (1:(N))*dt;

% ----- time window range
tpre = -fix( decay * parm.fs_est); 
%tend =  fix((Twin - dt) * parm.fs_est);
tend =  fix( Twin * parm.fs_est);
% ----- spike timing range for search
Trange = (tpre:tend)*dt_est;

% number of possible spike position in Twin
M = length(Trange);

% ----- start (t1) and end (t2) time for multiple spike packet
[t1,t2] = meshgrid(Trange);
t1 = t1(:);
t2 = t2(:);

% (t2 -t1) should be larger than MinISI = 1/max_freq
ix = find( (t1 + MinISI) < t2 );
% number of possible combination of [t1 t2]
Nlist = length(ix) ;

t1 = t1(ix); % start time
t2 = t2(ix); % end time

% ----- make multiple spike sequence
% - inter spike interval (ISI) are the same in one spike packet

% check variable
tlen_list = zeros(max_spike,Nlist);
ds_list   = zeros(max_spike,Nlist);
tt_list   = zeros(max_spike,Nlist);

% number of states for (spike number) = nspike
Nstate = zeros(max_spike,1); 

% time list of (spike number) = nspike
% Tspike{nspike} : (# of pattarn) x (# of spike)
Tspike = cell(max_spike,1);  

% 1 spike state
Nstate(1) = length(Trange);
Tspike{1} = Trange(:);

% 2 spike state
if max_spike >= 2
	Nstate(2) = Nlist;
	Tspike{2} = [t1(:) t2(:)];
end
% ----- multiple spike state
for nspike=3:max_spike
	for n=1:Nlist
		% period length
		tlen = (t2(n) - t1(n)); 
		% inter spike interval (ISI)
		ds = tlen/(nspike-1); 
		
		tlen_list(nspike,n) = tlen;
		ds_list(nspike,n) = ds;
		
		if ds < dt_est, continue; end;
		
		% spike time list
		tt = 0:(nspike-1);
		tt = tt*ds + t1(n);
		tt = tt( tt <= t2(n) ) ;
		
		% number of spikes
		mspike = length(tt);
		
		tt_list(nspike,n) = mspike;
		
		if mspike > 2
			% Tspike{nspike} : (# of pattarn) x (# of spike)
			Tspike{mspike} = [Tspike{mspike}; tt];
			Nstate(mspike) = Nstate(mspike) + 1;
		end
	end
end

spike_on  = find(Nstate > 0);
max_spike = spike_on(end);
Nstate = Nstate(1:max_spike);
Nall = sum(Nstate);

spike_state.Nstate = Nstate;
spike_state.Tspike = Tspike;
spike_state.Twin = Twin;
spike_state.dt = dt;
spike_state.dt_est = dt_est;
spike_state.tlen_list = tlen_list;
spike_state.ds_list = ds_list;
spike_state.tt_list = tt_list;

if print_mode == 0, return; end

n_ind = 1:max_spike;

	fprintf('\n----- Make multiple spike state\n')
	fprintf('Max number of spike = %d, MinISI = %6.3f\n',max_spike, MinISI)
	fprintf('sampling step = %6.3f, estimate sampling step = %6.3f\n',dt, dt_est)
	fprintf('Max number of spike (result) = %d\n',max_spike)
	fprintf('Number of total configuration = %d\n',Nall)
	fprintf('Nconf = %d (Nspike = %d)\n',[Nstate'; n_ind])


if (plot_mode == 2) & (max_spike > 1)
	figure;
	plot(Nstate)
	title('Number of configuration of spike states')
	ylabel('Number of configuration')
	xlabel('Number of spike')
end
NX = 3;
NY = 2;
nfig = NX*NY+1;

nskip  = 1;
nstate = 1;

if plot_mode >= 1,
	while nstate <= max_spike
		if (nfig > NX*NY)
			figure;
			nfig = 1;
		end
		
		if ~isempty(Tspike{nstate})
			subplot(NY,NX,nfig)
			plot(Tspike{nstate})
			title(sprintf('Nspike=%d',nstate))
			
			nfig = nfig + 1;
		end
		nstate = nstate + nskip;
	end
end

return

