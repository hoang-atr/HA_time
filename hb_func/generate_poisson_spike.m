function	spike_time = generate_poisson_spike(parm)
% ----- Poisson spike generation 
%    spike_time = generate_poisson_spike(parm)
%
% spike_time = time of spike firing [sec]  (Nspike x 1)
% ----- Poisson spike generation setting
% parm.Nspike   = number of spikes
% parm.fs_spike = spike freq [Hz]
% ----- multiple spike freq case
% parm.fs_spike(n) = n-th spike freq    [Hz]
% parm.t_period(n) = n-th window period [sec]
% ----- gamma_spike parameter
% parm.Cv = std(ISI)/mean(ISI) : coefficient of variation

%ポアソン過程に従うとき、
%inter-spike interval(ISI)は指数分布に従う
%指数分布に従うサンプルは、一様乱数の対数を取ればよい

fprintf('----- Poisson spike generation\n')

NT = length(parm.fs_spike);

if isfield(parm,'MinISI')
	MinISI = parm.MinISI;
else
	MinISI = [];
end

% gamma_spike parameter
if isfield(parm,'Cv')
	shape_parm = 1/parm.Cv^2;
else
	shape_parm = [];
end

% period length
if isfield(parm,'t_period')
	t_period = parm.t_period;
else
	t_period = [];
end

if NT == 1, 
	fs_spike = parm.fs_spike;
	if isfield(parm,'Nspike')
		Nspike = parm.Nspike; 
	elseif ~isempty(t_period)
		Nspike = ceil(fs_spike .* t_period); 
	else
		Nspike = 20;
	end
else
	if length(t_period)~=NT
		error('# of fs_spike ans t_period do not match')
	end
	% mean firing frequency
	fs_spike = parm.fs_spike;
	% spike number of the period
	Nspike   = ceil(fs_spike .* t_period) * 2; 
end

ISI_all = [];

for n=1:NT
	if fs_spike(n) == 0
		% Add no-spike period
		ISI_all = [ISI_all; t_period(n)];
	else
		% --- Generation of inter-spike interval (ISI)
		% 1. uniform distribution sample generation
		ISI = rand( Nspike(n), 1);
		
		ISI = max(ISI, eps);
		
		if isempty(shape_parm)
			% ISI of Poisson spike
			% 2. exponential distribution sample
			ISI = -(1./fs_spike(n)) .*log(ISI);
			
			if ~isempty(MinISI)
				ISI = max(ISI ,MinISI);
			end
		else
			% ISI of gamma spike
			% 2. gamma distribution sample
			ISI = gammaincinv(ISI,shape_parm)./(fs_spike(n)*shape_parm);
		end
		
		if ~isempty(t_period)
			ix = find(cumsum(ISI) <= t_period(n));
			ISI = ISI(1:ix(end));
		end
		
		ISI_all = [ISI_all; ISI];
	end
end

% 3. firing time = cumulative sum of ISI
% time of spike firing [sec]
spike_time = cumsum(ISI_all);

% parm.fs_spike = spike firing freq [Hz]
if isfield(parm,'t_period')
	for n=1:NT
		fprintf('spike firing freq = %8.3f [Hz], duration = %8.3f [sec]\n', ...
			parm.fs_spike(n) , parm.t_period(n))
	end
else
	fprintf('spike firing freq = %8.3f [Hz]\n', parm.fs_spike)
end
