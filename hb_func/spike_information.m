function	spike_info = spike_information(data, data_win, parm)
% spike_information of extracted window data
%   spike_info = spike_information(data, data_win, parm)
%
% --- spike info for extracted window
% spike_info.spike_num(n)  = spike number of n-th window
% spike_info.spike_time{n} = spike onset time of n-th window [sec] 
%                            absolute time (# of spike) x 1
% spike_info.twin(n,:)     = sample time index of n-th window
% spike_info.twin_sec(n,:) = sample time of n-th window [sec]
% spike_info.fs   : sampling frequency [Hz]
% spike_info.hist : spike number histgram
% --- Input
% data     : observed continuous data
% data_win : extracted window info
% parm.plot_mode = [0], 1  (optional)
% 
% --- sample time for observed continuous data
% data.t          = sample time of data.y [sec]
% --- onset spike timing for observed continuous data
% data.spike_time = time of spike firing [sec] with (fs_raw) resolution 
% data.peaks_time = time of spike peaks  [sec] with (fs_raw) resolution 
% --- extracted window info
% data_win.twin  : time index of windows
%                  (# of windows) x (time sample in one window)

plot_mode=0;
if nargin <3, parm=[];end
if isfield(parm,'plot_mode'),plot_mode = parm.plot_mode;end;

Nspike = length(data.spike_time);

% time index of windows
twin  = data_win.twin; 
spike_info.twin  = twin; 

% sampling step
dt = 1/data.fs;

% (# of window) x (# of time sample)
[Nwin ,Twin]= size(twin);

% spike_num(n)  = spike number of n-th window
% spike_time{n} = spike timing of n-th window
spike_num  = zeros(1,Nwin);
spike_time = cell(1,Nwin);

% spike count check flag for data.spike_time
spike_flg  = zeros(1,Nspike);
spike_index = [];

for n=1:Nwin
	tt = twin(n,:);       % time index of n-th window
	t0 = data.t(tt(1));   % start time of n-th window [sec]
	t1 = data.t(tt(end)); % end time
	
	% find spike onset in this window (including pre dt period)
	% find spike in this window (including pre dt period)
	ix = find( (data.spike_time >  (t0 - dt) ) ...
		     & (data.spike_time <= (t1 ) ) );
	
	if ~isempty(ix)
		spike_index = [spike_index; ix(:)];
		
		ts = data.spike_time(ix);
		
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

spike_info.twin_sec   = data.t(data_win.twin);
spike_info.fs = data.fs;

NSmax  = max(spike_num);

spike_info.hist = hist(spike_num, 0:NSmax);

NSall  = length(data.spike_time);% number of true spike
NSsum  = length(spike_index);    % independent spike number in windows
NSwin  = sum(spike_num > 0);     % number of windows with spike
Nfalse = sum(spike_num == 0);    % number of windows without spike
NSrest = NSall - NSsum;          % number of false negative

spike_info.Nspike = NSsum;       % independent spike number in windows
spike_info.Nfalse = Nfalse;      % number of windows without spike
spike_info.Npositive = NSwin;    % number of windows with spike

if plot_mode >=0
	fprintf('\n----- Spike information of extracted window\n')
	fprintf('# of windows = %d\n',Nwin)
	fprintf('# of windows with spikes = %d\n',NSwin)
	fprintf('# of all spikes = %d\n',NSall)
	fprintf('# of spikes within windows = %d\n',NSsum)
	%fprintf('Spike detect rate = %3.0f\n',100*NSsum/NSall)
	fprintf('# of false window detection = %3.0f\n',Nfalse)
	fprintf('# of false negative = %d\n',NSrest)
	fprintf('Max number of spike in one window = %d\n',NSmax)
end
% spike count check 
NSerr = sum(spike_flg == 0) - NSrest;
if NSerr ~= 0
	fprintf('Spike detect error = %d\n',NSerr)
end;

if plot_mode==0, return; end;

% plot extracted window data according spike number
n_spike  = find( spike_num > 0);
n_single = find( spike_num == 1);
n_multi  = find( spike_num > 1);
n_rest   = find( spike_num == 0);

ymax = max(data.y);

figure;
nfig = 0;

if length(n_single) > 0,
	nfig = nfig+1;
	subplot(2,2,nfig)
	tmin = 0;
	tmax = 0;
	for m=1:length(n_single)
		n = n_single(m);
		
		Ts = spike_info.spike_time{n};
		
		t = data.t(twin(n,:));
		%t = t - t(1);
		t = t - Ts(1);
		
		plot(t, data.y(twin(n,:)));
		hold on
		tmin = min([tmin, t]);
		tmax = max([tmax, t]);
	end
	tmax = max(tmax, Twin*dt);
	xlim([tmin tmax])
	title('Extracted window with single spike')
end	
if length(n_multi) > 0,
	nfig = nfig+1;
	subplot(2,2,nfig)
	tmin = 0;
	tmax = 0;
	for m=1:length(n_multi)
		n = n_multi(m);
		
		Ts = spike_info.spike_time{n};
		
		t = data.t(twin(n,:));
		t = t - t(1);
		%t = t - Ts(1);
		
		plot(t, data.y(twin(n,:)));
		hold on
		tmin = min([tmin, t]);
		tmax = max([tmax, t]);
	end
	tmax = max(tmax, Twin*dt);
	xlim([tmin tmax])
	title('Extracted window with multi-spike')
end
if length(n_rest) > 0,
	nfig = nfig+1;
	subplot(2,2,nfig)
	for m=1:length(n_rest)
		n = n_rest(m);
		t = data.t(twin(n,:));
		t = t - t(1);
		plot(t, data.y(twin(n,:)));
		hold on
	end
	title('Extracted window without spike')
end
if NSmax > 0
	nfig = nfig+1;
	subplot(2,2,nfig)
	hist(spike_num, 0:NSmax)
	xlim([0 NSmax])
	title('Histgram of #(spike) in one window')
end

if plot_mode < 2, return; end;

NX = 3;
NY = 3;
nfig = 0;
Nfig = length(n_spike);

if Nfig==0, return; end;

figure;
for m=1:Nfig
	if nfig >= NX*NY
		mesg = 'wait: ';
		% input char
		Result = input(mesg ,'s');
		
		if strcmp(Result , 'z')==1, 
			fprintf('\n')
			close all
			return; 
		end;
		
		close all
		nfig = 1;
		figure;
	else
		nfig = nfig + 1;
	end

	subplot(NY,NX,nfig)
	
	n = n_spike(m);
	
	t = data.t(twin(n,:));
	t = t - (t(1) -dt);
	
	plot(t, data.y(twin(n,:)));
	hold on
	plot(t, data.y(twin(n,:)), '--');
	
	nspike = spike_num(n);
	Tn = spike_info.spike_time{n};
	plot([Tn; Tn],[zeros(1,nspike); ones(1,nspike)],'r--')
	
	tmin = min( min(Tn) , 0);
	tmax = max( max(Tn) , t(end));
	xlim([tmin,  tmax])
end

return

figure
for n=1:Nwin
	if ~isempty(spike_time{n})
		plot(spike_time{n});
		hold on
	end
end
