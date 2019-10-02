function	[Model] = train_shape_step(data, tau, tau_id, tau_list)
% training shape parameter using training data
% data.y(n,:) : y(n,t)
% data.Dt{n}  : spike onset time of n-th window [sec] 
%               relative time when window onset time is t=0
% data.dt     = sampling step [sec]

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

% model variance list
sx_list = zeros(Nlist,1);

% search for tau
for n_tau=1:Nlist
	tau(tau_id) = tau_list(n_tau);
	
	% update spike amplitude & bias
	[a ,b] = update_spike_amplitude(tau, data);
	
	% update noise variance
	data.a   = a;
	data.b   = b;
	
	[sx] = update_model_variance(tau, data);
	
	a_list(:,n_tau) = a;
	b_list(n_tau)   = b;
	sx_list(n_tau)  = sx;
end

% find minimum variance
[sxmin, imin] = min(sx_list);

% optimal tau
tau(tau_id) = tau_list(imin);

Model.sx_list = sx_list;
Model.sx  = sx_list(imin);
Model.a   = a_list(:,imin);
Model.b   = b_list(imin);
Model.tau = tau;
