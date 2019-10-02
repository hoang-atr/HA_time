function	post_info = estimate_posterior_state(...
			state,spike_state,overlap,Data)
% --- spike estimation based on posterior for overlap window 
%   post_info = estimated_post_state(state,spike_state,overlap,Data)
% --- Estimated spike info for overlap window
% post_info.spike_num(n)  = estimated spike number of n-th overlap window
% post_info.Pspike(ns,n) = P(ns | yk) of n-th overlap window
% post_info.twin{n} = time index of n-th overlap window
% post_info.tsec = sample time of observed data [sec]
% post_info.fs   = sampling frequency [Hz]
% ----- Input
%	overlap.twin{n} = overlap time index in n-th overlap window
%	overlap.nwin{n} = overlap window index in n-th overlap window
%	overlap.num(n)  = number of overlap window in n-th overlap window
%	overlap.fs   = sampling frequency [Hz]
%	overlap.tsec = sample time of observed data [sec]
% --- Estimated Posterior
% state.logPY = log(P(yk))
%          : (Nwin x 1)
% state.P0 = P(Sk = 0|yk) = P(yk,Sk = 0)/P(yk)
%          : (Nwin x 1)
% state.P  = P(Tk,Sk=ns|yk) = P(yk,Tk,Sk=ns)/P(yk)
%          : (NconfMax x Nspike x Nwin);
% ----- Extractted spike window
% Data.twin(n,:)  : time index of n-th extractted window
% ----- spike_state informative
% spike_state.Tspike{nspike} = spike time onset list of (# of spike) = nspike

% max value for log
logMax = log(1.0e+15);

dt   = 1/overlap.fs;
tsec = overlap.tsec; % time of data [sec]
nwin = length(overlap.twin);

% copy overlap window informative
post_info.twin = overlap.twin;
post_info.tsec = overlap.tsec;
post_info.fs   = overlap.fs;

% post_info.spike_num(n)  = spike number within n-th overlap window
post_info.spike_num  = zeros(1,nwin);

% ----- Spike state prob.
% state.P  =  P(Tk,Sk=ns|yk)  : (NconfMax x Nspike x Nwin);
[NconfMax,Nspike,Nwin] = size(state.P);

% PS = sum_Tk P(Tk,Sk=ns|yk) 
%    = P(Sk=ns|yk)  : (Nspike x Nwin);
% PS = reshape(sum(state.P ,1), [Nspike,Nwin]);

% state.P0 = P(Sk = 0|yk) = P(yk,Sk = 0)/P(yk)

% logPY = log(P(yk))
% PY = P(yk)
LPmax = max(state.logPY);
logPY = state.logPY - LPmax;
PY = zeros(1,Nwin);
ix = find( logPY > - logMax);
PY(ix) = exp( logPY(ix) );

post_info.Pspike = zeros(Nspike,nwin);

% P(Sk , yk) = P(Sk |yk) * P(yk)

for n=1:nwin
	% time of divided period
	t  = overlap.twin{n} ; % time index
	ts = tsec(t); % time [sec]
	
	% window index in this period
	id = overlap.nwin{n};
	M  = length(id);
	
	% spike state prob
	P0 = zeros(M,1);       % prob for ns = 0
	Ps = zeros(M,Nspike);  % prob for ns > 0
	
	% windows in this period
	for m=1:M
		% m-th window index
		iwin = id(m);
		% start time of m-th window
		t0 = tsec(Data.twin(iwin,1));
		
		% state.P(nc,ns,iwin) = P(nc,ns|yk) , nc=1:Nconf, ns=1:Nspike
		% 1 = P(ns=0|yk) + sum_nc sum_ns P(nc,ns|yk)
		%   = state.P0(iwin) + sum_nc sum_ns state.P(nc,ns,iwin)
		
		Pz  = 0;               % prob for ns = 0
		Psm = zeros(1,Nspike); % prob for ns > 0
		
		% Loop for spike number
		for ns=1:Nspike
			% spike time list for ns-spike state
			Tk = spike_state.Tspike{ns} + (t0 - dt); % Nconf x ns
			Nconf = size(Tk,1);
			%P  = zeros(1,Nspike); % spike on  prob
			
			for nc=1:Nconf
				% find spike within this period
				ix = find( (Tk(nc,:) >= (ts(1) - dt)) ...
						&  (Tk(nc,:) < ts(end)));
				if isempty(ix),
					% no spike in this window & configuration
					Pz = Pz + state.P(nc,ns,iwin); 
				else
					%  nspike in this window & configuration
					nspike = length(ix);
					Psm(nspike) = Psm(nspike) + state.P(nc,ns,iwin); 
				end
			end
			
			% P(ns > 0) for m-th window
			%Psm(ns) = sum(P) ;
			%[pmax, imax] = max(P);
			
		end
		
		%[pmax, imax] = max(Psm);
		
		% add spike Prob for outside of this window
		P0m  = state.P0(iwin) + Pz; 
		Psum = sum(Psm) + P0m; % = 1
		
		% normalize prob for each window & multiply P(yk)
		P0(m)   = P0m * PY(iwin) / Psum;
		Ps(m,:) = Psm * PY(iwin) / Psum;
	end
	
	% spike state prob
	% P0 = zeros(M,1);       % prob for ns = 0
	% Ps = zeros(M,Nspike);  % prob for ns > 0
	
	% sum of prob for overlap window
	P0sum  = sum(P0);
	Pspike = sum(Ps,1);
	
	% normalise posterior 
	Psum   = P0sum + sum(Pspike);
	P0sum  = P0sum/Psum;
	Pspike = Pspike/Psum;
	
	post_info.Pspike(:,n) = Pspike(:);
	
	[Pmax,ns_opt] = max(Pspike);
	
	% most probable spike number in this period
%	if sum(Pspike) >= P0sum
	if Pmax >= P0sum
		post_info.spike_num(n)  = ns_opt;
	end
end

