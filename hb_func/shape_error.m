function	err = shape_error(data,Model,parm)

Ntau = length(Model.tau);

decayRate = 100;
[tpeak, gpeak, tdecay] = spike_func_peak(parm.tau, decayRate);
Nt = ceil(tdecay*parm.fs);
dt = 1/parm.fs;

data_eval.dt = dt; 
data_eval.Dt = 0;
data_eval.y  = zeros(1,Nt);

switch	Ntau
case	2
	g1 = spike_func_evaluate2( parm.tau(1:2), data_eval) * data.a(1) + data.b;
	g2 = spike_func_evaluate2( Model.tau, data_eval) * Model.a + Model.b;
case	3
	g0 = spike_func_evaluate2( parm.tau(1:3), data_eval) ;
	g1 = data.a(1) * g0(:,:,1) + data.a(2) * g0(:,:,2) + data.b;
	
	g0 = spike_func_evaluate2(Model.tau(1:3), data_eval) ;
	g2 = Model.a(1) * g0(:,:,1) + Model.a(2) * g0(:,:,2) + Model.b;
end

err = mean((g1 - g2).^2)/mean((g1).^2);
