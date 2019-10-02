function	Data = extract_spike_window(parm, y, print_mode)
% extract_spike_window
%    Data = extract_spike_window(parm, y)
% ----- parameters for window extraction
% parm.fs        =  sampling rate for observed data [Hz]
% parm.Twin      =  time window length [sec] 
% parm.Tpre      =  time period before peak onset [sec]
% ----- threshold for window extraction
% parm.threshold     :  1st threshold
% parm.threshold_max :  2nd threshold
%  - Condition : y >= mean(y) + std(y)*threshold
%  - if threshold is empty, windows are simply devided to Twin length
%    parm.Tstep    : step size of window [sec]
% ----- Extractted spike window
% Data.twin       : (# of windows) x (time sample in one window)
% Data.twin(n,:)  : time index of 'y' for n-th window data
% Data.y(n,:)     : data im n-th windows
% Data.yrest = data outside of windows (rest period)
% Data.trest = time in rest period
% Data.dt = 1/parm.fs : sampling step [sec]
% Data.ym = mean of data y
% Data.sd = std  of data y

if nargin<3, print_mode=1; end

plot_mode=0;
if isfield(parm,'plot_mode'),plot_mode = parm.plot_mode;end;

if isfield(parm,'threshold'),
	threshold = parm.threshold;
else
	threshold = [];
end;

if ~isempty(threshold)
	threshold_max = parm.threshold_max;
end

% Window size parameters
Twin = ceil(parm.Twin*parm.fs);
Twin = max(Twin,1);
if isfield(parm,'Tpre'),
	Tpre = ceil(parm.Tpre*parm.fs);
else
	Tpre = 0;
end

% # of samples
y = y(:);
T = length(y);

% Mean & STD
if isfield(parm,'sd0')
	ym = parm.y0;
	sd = parm.sd0;
else
	ym = mean(y);
	sd = std(y);
end

if isempty(threshold)
	% Window step parameters
	if isfield(parm,'Tstep') && ~isempty(parm.Tstep)
		Tstep = ceil(parm.Tstep*parm.fs);
	else
		Tstep = Twin;
	end
	
	% simple devision of windows without threshold
	nwin = ceil(T/(Tstep));
	onset = 0:(nwin-1);
	onset = onset*(Tstep) + 1;
	onset = onset( (onset + Twin - 1) <= T);
	nwin  = length(onset);
	th = ym + sd;
	th_max = th;
else
	% Two threshold
	th     = ym + sd*threshold;
	th_max = ym + sd*threshold_max;
	
	% 1. find start and end point over threshold
	ton  = find( ([0; y(1:T-1)] < th) & (y >= th));
	toff = find( ([y(2:T); 0]   < th) & (y >= th));
	Non  = length(ton);
	
	tlen = toff - ton;
	
	if any(tlen < 0), error('extraction error'); end;
	
	ix = find( tlen >= Twin);
	Nadd = sum(ceil( tlen(ix)/(Twin) ));
	Nonset = Non + Nadd;
	
	% window onset time variable
	onset  = zeros(1,Nonset); % onset time of n-th window
	offset = zeros(1,Nonset); % end   time of n-th window
	peak   = zeros(1,Nonset); % peak  time of n-th window
	
	% window index
	tt    = 0:(Twin-1);
	nwin  = 0;
	tlast = 0;
	
	% 2. set window onset time
	for n=1:Non
		% start/end time index of n-th region over threshold
		t1 = ton(n);
		t2 = toff(n);
		
		% check end of last selected window
		if tlast >= t1,	t1 = tlast + 1;	end
		
		while (t1 <= t2) && (t1 + Twin < T)
			% peak value & time in window
			%[ymax, imax] = max(y(tt + t1));
			%ymax >= th_max,
			
			% find onset of threshold_max 
			tst = find( y(tt + t1) > th_max );
			
			if isempty(tst)
				% skip this window 
				t1 = t1 + Twin;
			else
				% select new window if ymax is over th_max
				nwin = nwin + 1;
				
				% onset of threshold_max 
				%peak(nwin)   = t1 + imax - 1;
				peak(nwin)   = tst(1) + t1 -1;
				% onset of selected window
				onset(nwin)  = peak(nwin) - Tpre;
				onset(nwin)  = max(onset(nwin),1);
				
				% end of selected window
				offset(nwin) = onset(nwin) + Twin - 1; 
				tlast = offset(nwin);
				
				t1 = tlast + 1;
			end
		end
	end
	
	peak   = peak(1:nwin);
	onset  = onset(1:nwin);
	offset = offset(1:nwin);
end

% 3. signal for extracted window

% extracted window time index
t    = 0:(Twin-1);
twin = repmat(onset(:), [1 Twin]) + repmat(t, [nwin 1]);

Data.twin = twin; % nwin x Twin
Data.y    = y(twin);
Data.dt   = 1/parm.fs;
Data.ym   = ym;
Data.sd   = sd;

% peak time in window
[ymax, tpeak] = max(Data.y, [], 2);
% time index aligned at peak point
twin_pk = repmat(t(:), [1, nwin]) - repmat((tpeak(:))', [Twin, 1]);

% 5. flag = 1 for samples in extracted window
on_flg = zeros(1,T);
on_flg(twin(:)) = 1;

% Tota number of samples in window
Ntotal  = sum(on_flg > 0);
% number of samples in window with overlap
Nsample = length(twin(:));
% number of overlap samples
Nover = Nsample - Ntotal;

% 6. signal outside of extracted window : rest signal
% rest time index
trest = find(on_flg == 0);
Nrest = length(trest);
% subtract number of overlap sample
%Nrest = max(Nrest - Nover,0);
trest = trest(1:Nrest);

Data.yrest = y(trest);
Data.trest = trest;

if print_mode ==1
	fprintf('\n----- Extract spike window\n')
	fprintf('# of points in one window = %d (%4.0f [ms])\n',...
			Twin,parm.Twin*1000);

	fprintf('Extracted number of windows = %d\n',nwin)
	fprintf('Extracted number of samples in windows = %d\n',Nsample)
	fprintf('Number of rest samples = %d\n',Nrest)
	fprintf('Number of total samples = %d\n',T)
	fprintf('Number of overlap samples = %d\n', Nover)
	
	% --- check signal over threshold outside of extracted windows
	if ~isempty(threshold)
		% over threshold flag
		th_flg = zeros(1,T);
		th_ix = find( y >= th_max);
		th_flg(th_ix) = 1;
		
		N_th  = sum(th_flg > 0);
		N_and = sum((th_flg .* on_flg) > 0);
		Nlost = N_th - N_and;
		fprintf('Non-selected samples over threshold = %d\n',Nlost)
	end
end

if plot_mode == 0, return; end;

NY = 2;
NX = 1;

if plot_mode >= 1
	figure;
	subplot(NY,NX,1)
%	plot(Data.y')
	plot(twin_pk, Data.y')
	hold on
%	t1 = min(twin_pk(:));
%	t2 = min(twin_pk(:));
%	plot([t1; t2], [th; th], 'r-')
%	plot([t1; t2], [ym; ym], 'r--')
	title('Extracted window signal')
	
	if ~isempty(Data.yrest)
		subplot(NY,NX,2)
		plot(Data.yrest)
		hold on
		Trest = length(Data.yrest);
		
		plot([0; Trest], [th_max; th_max], 'r-')
		plot([0; Trest], [ym; ym], 'r--')
		title('Rest signal outside of extracted region')
	end
end

if plot_mode >= 2
	figure;
	plot(y)
	hold on
	
	plot([0; T], [th_max; th_max], 'r-')
	plot([0; T], [ym; ym], 'r--')
	title('Original signal with mean and threshold')
end

return
