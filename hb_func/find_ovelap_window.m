function	overlap = find_ovelap_window(...
			data, Data, spike_state, print_mode)
% find overlap window for post processing
%   overlap = find_ovelap_window(data, Data, spike_state)
% 
%	overlap.twin{n} = overlap time index 
%                     in n-th overlap window
%	overlap.nwin{n} = overlap window index 
%                     in n-th overlap window
%	overlap.num(n)  = number of overlap window 
%                     in n-th overlap window
%	overlap.tlen(n) = length of overlap window 
%                     in n-th overlap window
%	overlap.conn{m} = sets of connected window index
%
%	overlap.fs   = data.fs;
%	overlap.tsec = data.t;
% ----- Input
% data.t          = sample time [sec]
% data.fs         = sampling rate for observed data [Hz]
% ----- Extractted spike window
% Data.twin : (# of windows) x (time sample in one window)
% Data.twin(n,:)  : time index of 'y' for n-th window data
%
% spike_state.Tspike{nspike} = spike time onset list 
%                              of (# of spike) = nspike
%                            : (# of pattarn) x (# of spike)
if nargin<4, print_mode=1; end

if nargin==0
	Tw = 5;
	t0 = [2; 4; 7; 15; 19; 25];
	T  = t0(end) + Tw;
	t  = 1:Tw;
	twin = repmat(t, [size(t0,1) 1]) + repmat(t0, [1 size(t,2)]);
	plot_mode = 1;
else
	T = length(data.t);
	fs = data.fs;
	
	% Window time index
	twin0 = Data.twin;
	[Nw, Nt] = size(twin0);
	% earliest spike time in window (relative time)
	Tspike = min(spike_state.Tspike{1});
	% earliest spike time index
	tpre = fix(Tspike*fs);
	% extended window time index
	twin = zeros(Nw, Nt + abs(tpre) + 1);
	
	for n=1:Nw
		tw = twin0(n,:);
		tp = tpre:0;
		tp = tp + tw(1) - 1;
		twin(n,:) = [tp, tw];
	end
	
	twin = max(twin,1);
end

overlap.fs   = data.fs;
overlap.tsec = data.t;

[Nw, Nt] = size(twin);

% find overlap window
% flag(t) : # of overlap window at t
flag = zeros(1,T);

for n=1:Nw
	flag(twin(n,:)) = flag(twin(n,:)) + 1;
end

%  [0 flag] : [ 0  1  1  2  2  0  0  0  1  2  1  1  0 ]
%  diff     :    [ 1  0  1  0 -2  0  0  1  1 -1  0 -1]
dif_flag = diff([0 flag 0]);
% start index of overlap window 
k1 = find( dif_flag ~= 0 );
% add start index of window 
k1 = [k1(:) ; twin(:,1)];
k1 = unique(k1);
NO = length(k1);

% overlap time & window index 
ovlp_t   = cell(1,NO);
ovlp_win = cell(1,NO);
ovlp_num = zeros(1,NO);

for n=1:NO-1
	% overlap time index 
	t = [k1(n):k1(n+1)-1];
	ovlp_t{n} = t;
	
	tflg = zeros(1,T);
	tflg(t) = 1;
	
	for k=1:Nw
		flg = zeros(1,T);
		flg(twin(k,:)) = 1;
		
		% check k-th window is overlaped with t
		if sum(tflg .* flg) > 0
			% overlap window index 
			ovlp_win{n} = [ovlp_win{n}, k];
			ovlp_num(n) = ovlp_num(n) + 1;
		end
	end
end

% overlap time index 
overlap.twin = [];
% overlap window index 
overlap.nwin = [];
overlap.num  = [];

% delete empty elements
nwin =0;

for n=1:NO
	if ~isempty(ovlp_win{n})
		nwin = nwin + 1;
		overlap.twin{nwin} = ovlp_t{n};
		overlap.nwin{nwin} = ovlp_win{n};
		overlap.num(nwin)  = ovlp_num(n);
	end
end

% connected window

over_conn = cell(1,nwin);
tlast = 0;
nconn = 0;

for n=1:nwin
	t = overlap.twin{n};
	
	% check n-th window is separated with last window
	if t(1) > tlast,
		nconn = nconn + 1;
	end
	
	over_conn{nconn} = [over_conn{nconn}, overlap.nwin{n}];
	tlast = t(end) + 1;
end

overlap.conn = cell(1,nconn);
for n=1:nconn
	overlap.conn{n} = unique(over_conn{n});
end

if print_mode == 0, return; end;

Noverlap = sum( overlap.num > 1 );

fprintf('------- Overlap Window --------\n')
fprintf('Number of overlap window = %d\n',Noverlap)
fprintf('Number of consecutive overlap window = %d\n',length(overlap.conn))

return
% ==================== END =====================
% ==================== END =====================
