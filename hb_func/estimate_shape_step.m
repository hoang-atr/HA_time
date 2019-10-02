function	[Model] = estimate_shape_step(...
			data, parm, spike_state, tau, tau_id, tau_list)
% training shape parameter using test data

% data.y(n,:) : y(n,t)
% data.Dt{n}  : spike onset time of n-th window [sec] 
%              relative time when window onset time is t=0
% data.dt     = sampling step [sec]

print_mode = 0;

Ntau  = length(tau);
Nlist = length(tau_list);

% model parameter list
switch	Ntau
case	2
	a_list  = zeros(1,Nlist);
case	3
	a_list  = zeros(2,Nlist);
end

b_list  = zeros(Nlist,1);
b0_list = zeros(Nlist,1);

% model variance list
sx_list = zeros(Nlist,1);
% threshold list
th_list = zeros(Nlist,1);
py_list = zeros(Nlist,1);
py_hist = cell(Nlist,1);

%fprintf('--- seach start\n')

% search for tau
for n_tau=1:Nlist
%	fprintf('-')
	
	tau(tau_id) = tau_list(n_tau);
	
	parm.threshold = 0.5;
	parm.threshold_max = 0.5;
	% ----- Extract_spike_window;
	Data = extract_spike_window(parm, data.y, print_mode);
	% Exclude non-window rest period for tau search
	Data.yrest = [];
	
	Model_ini = estimate_spike_init(tau, Data);
	Model_ini = set_train_parm(Model_ini);

	th_list(n_tau) = Model_ini.opt_threshold;

	parm.threshold = Model_ini.opt_threshold;
	parm.threshold_max = Model_ini.opt_threshold;
	% ----- Extract_spike_window;
	Data = extract_spike_window(parm, data.y, print_mode);
	% Exclude non-window rest period for tau search
	Data.yrest = [];

	[Model] = estimate_multispike_step(...
				Model_ini,Data,spike_state);

	a_list(:,n_tau) = Model.a;
	b_list(n_tau)  = Model.b;
	b0_list(n_tau) = Model.b0;
	sx_list(n_tau) = Model.sx;
	py_list(n_tau) = Model.LogPostY(end);
	py_hist{n_tau} = Model.LogPostY;
end

%fprintf('\n--- seach end\n')

opt_mode = 2;

switch	opt_mode
case	1
	% find minimum variance
	[sxmin, iopt] = min(sx_list);
case	2
	% find maximum posterior
	[pymax, iopt] = max(py_list);
end

% optimal tau
tau(tau_id) = tau_list(iopt);

Model.tau = tau;

Model.sx_list = sx_list;
Model.LogPostY = py_hist{iopt};

Model.sx  = sx_list(iopt);
Model.a   = a_list(:,iopt);
Model.b   = b_list(iopt);
Model.b0  = b0_list(iopt);
Model.opt_threshold = th_list(iopt);
