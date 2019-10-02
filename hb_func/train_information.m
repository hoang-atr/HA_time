function	train_info = train_information(data, parm)
% train_information of spike data
%   train_info = train_information(data, parm)
%
% --- train_info
% train_info.spike_num(n)  = spike number of n-th window
% train_info.spike_time{n} = spike peak time 
%                            of n-th window [sec] 
%                            relative time 
%                            when window onset time is t=0
%                            (# of spike) x 1
% train_info.twin(n,:)    = time index of n-th window
% train_info.twin_sec(n,:) = sample time 
%                            of n-th window [sec]
% train_info.fs : sampling frequency [Hz]
% --- Input
% data     : observed continuous data
% 
% --- sample time for observed continuous data
% data.t          = sample time of data.y [sec]
% --- onset spike timing for observed continuous data
% data.spike_time = time of spike firing [sec] with (fs_raw) resolution 

global DEBUG

% Window size parameters
Twin = ceil(parm.Twin*parm.fs);
Tpre = ceil(parm.Tpre*parm.fs);
Twin = max(Twin,1);

decay = parm.decay;

% sampling step
dt = 1/data.fs;
t  = data.t;
Tmax = length(t);

Nspike = length(data.spike_time);

% spike_time{n} = spike timing of n-th window
spike_win  = zeros(Nspike,Twin);
spike_num  = zeros(1,Nspike);
spike_time = {};

for n=1:Nspike
	% n-th spike time index
	t0 = ceil(data.spike_time(n) * data.fs);
%	t0 = max(t1 - Tpre, 1);
	t2 = t0 + Twin -1;
	
	if t2 > Tmax,
		Nspike = n-1;
		break;
	end
	
	% window time index
	spike_win(n,:) = t0:t2;
	
	% spike onset at t2 give no contribution in this window
	% spike onset at (t0 - decay) 
	%   give contribution in this window
	ix = find( (data.spike_time >= (t(t0) - decay) ) ...
		     & (data.spike_time <   t(t2)) );
	
	% relative time when window onset time (t0 - dt) is t=0
	ts = data.spike_time(ix) - (t(t0) - dt);
	
	% number of spike in n-th window
	spike_num(n)  = length(ix);
	
	% spike onset timing in n-th window
	spike_time{n} = ts(:)';
end

% train_info.spike_num(n) = spike number of n-th window
train_info.spike_num  = spike_num(1:Nspike)  ;
train_info.spike_time = spike_time ;
train_info.twin       = spike_win(1:Nspike,:);
train_info.twin_sec   = t(train_info.twin);
train_info.fs = data.fs;

% rest state
flg = zeros(1,length(t));
flg(train_info.twin(:)) = 1;
trest = find(flg == 0);

train_info.trest = trest;

% 95 percentile of spike number
[spnum, indx] = sort(spike_num);
pos = fix(length(spike_num) * 0.8);
train_info.max_spike = spnum(pos);

%if exist('DEBUG','var') && DEBUG == 1, 
%	plot_train_info(data, train_info, decay)
%end
