function	spike_info = spike_overlap_info(...
			data, overlap, print_mode)
% spike_information for overlap window 
%   spike_info = spike_overlap_info(data, overlap)
%
% --- spike info for overlap window 
% spike_info.spike_num(n)  = spike number of n-th overlap window
% spike_info.spike_time{n} = spike onset time of n-th overlap window [sec] 
%                            absolute time [sec] (# of spike) x 1
% spike_info.twin{n} = overlap time index in n-th overlap window
% spike_info.tsec = data.t; sample time of data.y [sec]
% spike_info.twin_sec{n} = sample time of n-th overlap window [sec]
% spike_info.NSoutside = number of spike outside of window
% spike_info.fs   : sampling frequency [Hz]
% spike_info.hist : spike number histgram
% --- Input
%	overlap.twin{n} = overlap time index in n-th overlap window
% --- sample time for observed continuous data
% data.t          = sample time of data.y [sec]
% --- onset spike timing for observed continuous data
% data.spike_time = time of spike firing [sec] with (fs_raw) resolution 

if nargin <3, print_mode = 1; end;

plot_mode = 0;

Nspike = length(data.spike_time);

tsec = data.t; % time of data [sec]
spike_info.tsec = tsec;

nwin = length(overlap.twin);

spike_info.twin = overlap.twin;
spike_info.twin_sec = cell(1,nwin);

spike_num  = zeros(1,nwin);
spike_time = cell(1,nwin);

% time index of windows
%twin  = data_win.twin; 
% sampling step
dt = 1/data.fs;

% spike count check flag for data.spike_time
spike_flg  = zeros(1,Nspike);
spike_index = [];

for n=1:nwin
	t  = overlap.twin{n} ; % time index
	ts = tsec(t); % time [sec]
	t0 = ts(1);   % start time of n-th window [sec]
	t1 = ts(end); % end time
	
	spike_info.twin_sec{n} = ts;
	
	% signal at t(n) = (n * dt) is representative of t(n-1) < t <= t(n)
	% spike onset at t1 give no contribution in this window
	% spike onset at (t0 - dt) give contribution in this window

	% find spike onset in this window (including pre dt period)
	ix = find( (data.spike_time >  (t0 - dt) ) ...
		     & (data.spike_time <= (t1 ) ) );
	
	if ~isempty(ix)
		spike_index = [spike_index; ix(:)];
		
		% relative time when window onset time (t0 - dt) is t=0
%		ts = data.spike_time(ix) - (t0 - dt);
		ts = data.spike_time(ix) ;
		
		% number of spike in n-th window
		spike_num(n)  = length(ix);
		% spike count check for data.spike_time
		spike_flg(ix) = 1;
		
		% spike onset timing in n-th window
		spike_time{n} = ts;
		
	end
end

spike_index = unique(spike_index);

% spike_info.spike_num(n) = spike number of n-th window
spike_info.spike_num  = spike_num  ;
spike_info.spike_time = spike_time ;
spike_info.spike_index = spike_index;

spike_info.fs = data.fs;

NSmax  = max(spike_num);

spike_info.hist = hist(spike_num, 0:NSmax);

NSall  = length(data.spike_time);% number of true spike
NSsum  = length(spike_index);    % independent spike number in windows
NSwin  = sum(spike_num > 0);     % number of windows with spike
Nfalse = sum(spike_num == 0);    % number of windows without spike
NSrest = NSall - NSsum;          % number of spike outside of window

spike_info.Nspike = NSsum;       % independent spike number in windows
spike_info.Nfalse = Nfalse;      % number of windows without spike
spike_info.Npositive = NSwin;    % number of windows with spike
spike_info.NSoutside = NSrest;   % number of spike outside of window

if print_mode ==1
	fprintf('\n----- Spike information of extracted window\n')
	fprintf('# of windows = %d\n',nwin)
	fprintf('# of windows with spikes = %d\n',NSwin)
	fprintf('# of all spikes = %d\n',NSall)
	fprintf('# of spikes within windows = %d\n',NSsum)
	fprintf('# of false window detection = %3.0f\n',Nfalse)
	fprintf('# of false negative = %d\n',NSrest)
	fprintf('Max number of spike in one window = %d\n',NSmax)

	% spike count check 
	NSerr = sum(spike_flg == 0) - NSrest;
	if NSerr ~= 0
		fprintf('Spike detect error = %d\n',NSerr)
	end;
end

