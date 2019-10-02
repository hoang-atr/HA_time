function	[Model, train_info] = train_shape(data, parm)
% training shape parameter using training data
%  [Model] = train_shape(data, parm)
%
% data.y   = observed spike data with noise 
% data.t   = sample time [sec]
% data.fs  = sampling rate for observed data [Hz]
% data.spike_time  = time of spike onset [sec] with (fs_raw) resolution 
%
%	parm.fs_est  = 100;   % estimation freq [Hz]
%	parm.Twin    = 1.5;   % (15) sample
%	parm.Tpre    = 0.4;   
%	parm.decay    = 1; % decay time [sec]
%	                   % using estimation spike configulation 
%	
%	parm.threshold = 1;
%	parm.threshold_max = 1.5;
%	
%	parm.max_freq = 10;  % max firing frequency [Hz]
%	parm.max_spike = 2;
%	parm.Ntau = 3;
%
%	Model.tau 
%	Model.a   
%	Model.b   
%	Model.sx  
%	Model.b0
%	Model.sx_list 

% ---- initial search parameter
% Tau.tau_ini  = [0.005; 0.1];
% Tau.max_ini  = max value of tau
% Tau.min_ini  =  grid step
% Tau.step_ini =  grid step
% 
% ---- main search parameter
% Tau.range =  search range of tau
% Tau.step  =  grid step
% Tau.min   =  minimum value
% 

% ----- Get train information
train_info = train_information(data, parm);

% ----- set train data
train_data.dt = 1/data.fs; 
train_data.Dt = train_info.spike_time;
train_data.y  = data.y(train_info.twin);

% ----- set of 1D-grid  seach parameter
if isfield(parm,'Ntau')
	Tau  = set_tau_grid_cs(parm.Ntau);
else
	Tau  = set_tau_grid_cs;
end

Ntau = Tau.Ntau;

sx_list = [];

% ----- pre-train model shape
% initial tau for [??1; ??2]
tau = Tau.tau_ini;

% ------ initial search for [??1; ??2]
%  optimization search for ??2
tau_id = 2;
tau_list = Tau.min_ini(tau_id):Tau.step_ini(tau_id):Tau.max_ini(tau_id);
Model = train_shape_step(train_data, tau, tau_id, tau_list);

tau = Model.tau;
sx_list = [sx_list, Model.sx];

%  optimization search for ??1
tau_id = 1;
tau_list = Tau.min_ini(tau_id):Tau.step_ini(tau_id):Tau.max_ini(tau_id);
tau_list = tau_list( tau_list < tau(2));
Model = train_shape_step(train_data, tau, tau_id, tau_list);

tau = Model.tau;
sx_list = [sx_list, Model.sx];

% ------ optimization search for [??1; ??2]
for itr=1:Tau.Nupdate_ini
	for tau_id = 2:-1:1
	
		tau_list = -Tau.range_ini(tau_id):Tau.step_ini(tau_id): ...
		            Tau.range_ini(tau_id);
		tau_list = tau_list + tau(tau_id);
		
		switch	tau_id
		case	1
			tau_list = tau_list( tau_list < tau(2) ...
			& tau_list <=Tau.max_ini(1) ...
			& tau_list >=Tau.min_ini(1));
		case	2
			tau_list = tau_list( tau_list > tau(1) ...
			& tau_list <=Tau.max_ini(2) ...
			& tau_list >=Tau.min_ini(2));
		end
		
		Model = train_shape_step(train_data, tau, tau_id, tau_list);
		
		tau = Model.tau;
		sx_list = [sx_list, Model.sx];
	end
end

if  Ntau == 3
	%  ----- optimization search for ??3 with [??1; ??2; ??3] 
	tau = [tau; Tau.min_ini(3)];
	
	tau_id = 3;
	tau_list = Tau.min_ini(tau_id):Tau.step_ini(tau_id):Tau.max_ini(tau_id);
	tau_list = tau_list( tau_list > tau(2));
	Model = train_shape_step(train_data, tau, tau_id, tau_list);
	
	tau = Model.tau;
	sx_list = [sx_list, Model.sx];
end

switch	Ntau
case	2
	tau_order = 2:-1:1;
case	3
	tau_order = 1:3;
end

% ------ optimization search for [??1; ??2; ??3] 
for itr=1:Tau.Nupdate
	for tau_id = tau_order
	
		tau_list = -Tau.range(tau_id):Tau.step(tau_id):Tau.range(tau_id);
		tau_list = tau_list + tau(tau_id);
		
		switch	tau_id
		case	1
			tau_list = tau_list( tau_list < tau(2) ...
			& tau_list <=Tau.max(1) ...
			& tau_list >=Tau.min(1));
		case	2
			switch	Ntau
			case	2
				tau_list = tau_list( tau_list > tau(1) ...
				& tau_list <=Tau.max(2) ...
				& tau_list >=Tau.min(2));
			case	3
				tau_list = tau_list( tau_list < tau(3) & tau_list > tau(1) ...
				& tau_list <=Tau.max(2) ...
				& tau_list >=Tau.min(2));
			end
		case	3
			tau_list = tau_list( tau_list > tau(2) ...
			& tau_list <=Tau.max(3) ...
			& tau_list >=Tau.min(3));
		end
		
		Model = train_shape_step(train_data, tau, tau_id, tau_list);
		
		tau = Model.tau;
		sx_list = [sx_list, Model.sx];
	end
end

Model.p = [0 0];        % default nonlinear parameter (i.e, linear)
Model.sx_list = sx_list ;

% set spike amplitude & bias
train_data.a   = Model.a;
train_data.b   = Model.b;

% set no spike state in train_data
train_data.yrest = data.y(train_info.trest);

% update noise variance
[Model.sx, Model.b0] = ...
	update_model_variance(Model.tau, train_data);
