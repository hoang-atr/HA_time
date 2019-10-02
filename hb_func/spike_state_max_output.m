function		out_info = spike_state_max_output(state, spike_state)
% --- estimated spike state for window 
%   out_info = spike_state_max_output(state,spike_state)
% --- Estimated spike info for extracted window
% out_info.spike_num(n)  = spike peak number within n-th window : (1 x Nwin) 
% out_info.spike_time{n} = estimated spike onset time of n-th window [sec] 
%                          relative time when window onset time is t=0
% out_info.Pspike(ns,n)  = P(ns | yk)  : (Nspike x Nwin)
% out_info.PY(n)         = P(yk)  : (1 x Nwin)
% --- possible range of spike state
% spike_state.Tspike{nspike} = time list of (spike number) = nspike
%                            : (# of pattarn) x (nspike)
% --- Estimated Posterior
% state.P  = P(Tk,Sk=ns|yk) = P(yk,Tk,Sk=ns)/P(yk)
%          : (NconfMax x Nspike x Nwin);
% state.P0 = P(Sk = 0|yk) = P(yk,Sk = 0)/P(yk)
%          : (Nwin x 1)

% Model.tau: Estimated tau
%tpeak = spike_func_peak(Model.tau);
%Twin  = spike_state.Twin; % = parm.Twin: window length [sec]

% max value for log
logMax = log(1.0e+15);

[Nconf,Nspike,Nwin] = size(state.P);

% spike_num(n) = spike number of n-th window
spike_num  = zeros(1,Nwin);
spike_flg  = zeros(1,Nwin);
spike_time = cell(1,Nwin);
%entropy    = zeros(1,Nwin);

% ----- Spike state prob.
% state.P  =  P(Tk,Sk=ns|yk) : (NconfMax x Nspike x Nwin)
% PS = sum_Tk P(Tk,Sk=ns|yk) 
%    = P(Sk=ns|yk)           : (Nspike x Nwin)
PS  = reshape(sum(state.P ,1), [Nspike,Nwin]);
% Pspike = P(Sk > 0|yk)      : 1 x Nwin
Pspike = sum(PS,1);

% ----- Max prob. over Nspike states
[PSmax,nspike] = max(PS, [], 1);

% ----- Check spike or non-spike condition
% P(Sk > 0|yk) >= P(Sk = 0|yk) 
ix_on = find( Pspike >= state.P0' );

% number of spike for maximum posterior state
spike_num(ix_on) = nspike(ix_on);

% logPY = log(P(yk))
% PY = P(yk)
LPmax = max(state.logPY);
logPY = state.logPY - LPmax;
PY = zeros(1,Nwin);
ix = find( logPY > - logMax);
PY(ix) = exp( logPY(ix) );

out_info.PY = PY;

% ----- Spike timing list for spike states
for nk=1:Nwin
%	nk = ix_on(n);      % nk-th window
	
	% spike number of nk-th window
	ns = nspike(nk); 
	
	% spike time list for ns-spike state
	Tk = spike_state.Tspike{ns};
	
	% Find max prob. configuration for Tk
	% Ps = P(Tk,Sk = ns|yk) 
	Ps = state.P(:,ns,nk);
	[pmax,imax] = max(Ps);
	
	Ton = Tk(imax,:);
	spike_time{nk} = Ton;
	
%	% find spike within window
%	ix = find( Ton >= 0 & Ton < Twin);
%	if ~isempty(ix)
%		spike_flg(nk) = 1;
%	end
%	
%	Tpeak = Ton + tpeak;
%	
%	if (Ton(end) > 0 ) && ( Ton(1) <= Twin)
%		spike_flg(nk) = 1;
%	end
end

%ix_off = find( spike_flg == 0 );
%spike_num(ix_off) = 0;

out_info.spike_time = spike_time;
out_info.spike_num  = spike_num;
out_info.Pspike     = PS;
%out_info.entropy   = entropy;

return
%hist(state_entropy,0:Nspike+1);
