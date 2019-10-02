function	[Model] = estimate_spike_init(tau, Data)
% ----- Initial Estimation of model parameters 
%        using peak of windows
%	[Model] = estimate_spike_init(tau, Data)
% Model.opt_threshold = opt_th;

plot_mode = -1;

if plot_mode >= 0
	fprintf('--- Estimate_spike initial step\n')
end

Ntau = length(tau);
% tau : spike waveform time constants 
%       (this value is fixed in this step)
Model.tau = tau;
tpeak = spike_func_peak(tau);

% --- threshold list
threshold  = 0.5:0.1:2.5;

% Dim of data : 
% (# of windows) x (# of time sample in one windows)
y = Data.y;
[K,Nt] = size(y);

% time step
dt = Data.dt;
data.dt = dt;
% ----- time samples of data in a window
t  = (1:Nt)*dt;

ym = Data.ym;
sd = Data.sd;

% find peak value & time in each window
[ymax, tmax] = max(y, [], 2);

Nlist = length(threshold);
% model parameter list
switch	Ntau
case	2
	a_list  = zeros(1,Nlist);
case	3
	a_list  = zeros(2,Nlist);
end

b_list  = zeros(Nlist,1);
b0_list  = zeros(Nlist,1);

% model variance list
sx_list = zeros(Nlist,1);
Nspike  = zeros(Nlist,1);

for n=1:Nlist
	% threshold
	th = ym + sd*threshold(n);
	
	ix_spike = find( ymax >= th);
	ix_rest  = find( ymax <  th);
	
	if isempty(ix_spike)
		Nlist = n-1;
		break;
	end
	
	Nspike(n) = length(ix_spike);
	
	% peak time in window
	Tk = t( tmax(ix_spike) );
	% onset time in window
	Dt = cell(Nspike(n),1);
	
	for m=1:Nspike(n)
		Dt{m} = Tk(m) - tpeak; % : 1 x (# of spike)
	end
	
	data.Dt = Dt;
	data.y  = y(ix_spike,:);
	
	% update spike amplitude & bias
	[a ,b] = update_spike_amplitude(tau, data);
	
	% update noise variance
	data.a   = a;
	data.b   = b;
	
	% no spike state data
	data.yrest = y(ix_rest,:);
	
	[sx, b0] = update_model_variance(tau, data);
	
	a_list(:,n) = a;
	b_list(n)   = b;
	b0_list(n)  = b0;
	sx_list(n)  = sx;
end

n_list  = 1:Nlist;
sx_list = sx_list(n_list);

[sxmin,imin] = min(sx_list);
opt_th = threshold(imin);

Model.opt_threshold = opt_th;

Model.sx_list = sx_list;

Model.sx  = sx_list(imin);
Model.a   = a_list(:,imin);
Model.b   = b_list(imin);
Model.b0  = b0_list(imin);
Model.tau = tau;

if ~isempty(Data.yrest)
	% update noise variance
	% threshold
	th = ym + sd*opt_th;
	
	ix_spike = find( ymax >= th);
	ix_rest  = find( ymax <  th);
	% peak time in window
	Tk = t( tmax(ix_spike) );
	% onset time in window
	Dt = cell(Nspike(n),1);
	
	for m=1:Nspike(n)
		Dt{m} = Tk(m) - tpeak; % : 1 x (# of spike)
	end
	
	data.Dt = Dt;
	data.y  = y(ix_spike,:);
	
	data.a   = Model.a;
	data.b   = Model.b;
	
	% no spike state data
	yrest = y(ix_rest,:);
	data.yrest = [yrest(:); Data.yrest(:)];
	
	[Model.sx, Model.b0] = update_model_variance(tau, data);
end

if plot_mode >= 0
	fprintf('# of initial spike window at opt_threshold (%5.2f) =%d\n', ...
		opt_th, Nspike)
end
if plot_mode > 0
	figure;
	subplot(2,1,1)
	plot(n_list,sx_list)
	title('SX')
	xlim([threshold(1) threshold(Nlist)])
	
	subplot(2,1,2)
	plot(n_list,Nspike(n_list))
	title('Nspike')
	xlim([threshold(1) threshold(Nlist)])
end	
return
