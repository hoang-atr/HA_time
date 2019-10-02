function	[g, dg1, dg2] = spike_func_evaluate(tau, data)
%  spike_func_evaluate
% --- 
%   [g, dg1, dg2] = spike_func_evaluate(tau, data)
% --- data in spike window
%   data.y(k,:) : y(k,t) : observed data in k-th window (Nwin x Nt)
%   data.Dt{k}  : spike onset in  k-th window [sec] (1 x Nspike)
%   data.dt     : sampling time step [sec]
% --- time constant of spike response
%  tau = [É—1; É—2]   or  [É—1; É—2; É—3] 
%  g   : (Nwin x Nt)  or (Nwin x Nt x 2)
% ---  double exponential function for [É—1; É—2]
% g(t) = ( 1 - exp(-t/É—1) ) * exp(-t/É—2) for t >= 0
%      = ( 1 - f1(t) ) * f2(t) for t >= 0
% ---  double exponential function with two decay constant for [É—1; É—2; É—3] 
% g1(t) = ( 1 - exp(-t/É—1) ) * exp(-t/É—2) for t >= 0
% g2(t) = ( 1 - exp(-t/É—1) ) * exp(-t/É—3) for t >= 0
%
% g1(t) = ( 1 - f1(t) ) * f2(t) 
% g2(t) = ( 1 - f1(t) ) * f3(t) 
% --- gradient
%   f1(t) = exp(-t/É—1) 
%   f2(t) = exp(-t/É—2) 
%
%	dg1 = dg/dÉ—1 =  - f1 *f2*(t/É—1^2)
%	dg2 = dg/dÉ—2 = (1-f1)*f2*(t/É—2^2)

global DEBUG

Ntau = length(tau);

% Nwindow x Ntime
[Nw,Nt] = size(data.y);
% sampling time step
dt  = data.dt;

% sampling time 
t  = (1:(Nt)) * (dt) ;
if Ntau==2
	g  = zeros(Nw,Nt);
elseif Ntau==3
	g  = zeros(Nw,Nt,2);
end

if nargout==3
	dg1  = zeros(Nw,Nt);
	dg2  = zeros(Nw,Nt);
end
% DEBUG variable
if DEBUG == 2,
	gp = cell(Nw,1);
	Dtp= cell(Nw,1);
end

for n=1:Nw
	% data.Dt{n} : spike onset time: 1 x Nspike
	Dt0 = data.Dt{n} ; 
	Ns = length(Dt0);
	
	% Dt : t - Tk : Nt x Ns
	Dt = repmat(t(:), [1 Ns]) - repmat(Dt0, [Nt 1]);
	
	%   It : It(k,n) = 1 for Dt(k,n) >  0
	%                = 0 for Dt(k,n) <= 0
	It = (Dt > 0);
	Dt =  Dt .* It;
	
	% ---  double exponential function
	f1 = exp(- Dt / tau(1)) ;
	f2 = exp(- Dt / tau(2)) ;
	g0 = (1 - f1) .* f2 .* It;
	
	switch	Ntau
	case	2
		% ---  sum over Ns spikes
		g(n,:) = sum(g0,2)';
	case	3
		f3 = exp(- Dt / tau(3)) ;
		g1 = (1 - f1) .* f3 .* It;
		
		% ---  sum over Ns spikes
		g(n,:,1) = sum(g0,2)';
		g(n,:,2) = sum(g1,2)';
	end
	
	if nargout==3
		% dg/dÉ— : sum over multiple spike waveforms
		%	dg/dÉ—1 =  - f1 *f2*(t/É—1^2)
		%	dg/dÉ—2 = (1-f1)*f2*(t/É—2^2)
		dg1(n,:)  = (1/tau(1)^2)*sum(   - f1  .* f2 .* Dt, 2)';
		dg2(n,:)  = (1/tau(2)^2)*sum((1 - f1) .* f2 .* Dt, 2)';
	end
	% DEBUG variable
	if DEBUG == 2,
		gp{n}  = g0;
		Dtp{n} = data.Dt{n}/dt;
	end
end

if DEBUG == 2,
	%check_plot(data.y , g0 , data.Dt);
	check_plot(data.y , gp , Dtp);
end
