function	g = spike_func_evaluate2(tau, data)
%  spike_func_evaluate
% --- 
%   
% --- data in spike window
%   data.y  : y(t)  (1 x Nt)
%   data.Dt : spike onset (# of pattarn) x (# of spike) : (Nconf x Ns)
%   data.dt     : sampling time step [sec]
% --- time constant of spike response
%  tau = [ƒÑ1; ƒÑ2]    or  [ƒÑ1; ƒÑ2; ƒÑ3] 
%  g   : (Nconf x Nt)  or (Nconf x Nt x 2)
% ---  double exponential function for [ƒÑ1; ƒÑ2]
% g(t) = ( 1 - exp(-t/ƒÑ1) ) * exp(-t/ƒÑ2) for t >= 0
%      = ( 1 - f1(t) ) * f2(t) for t >= 0
% ---  double exponential function with two decay constant for [ƒÑ1; ƒÑ2; ƒÑ3] 
% g1(t) = ( 1 - exp(-t/ƒÑ1) ) * exp(-t/ƒÑ2) for t >= 0
% g2(t) = ( 1 - exp(-t/ƒÑ1) ) * exp(-t/ƒÑ3) for t >= 0
%
% g1(t) = ( 1 - f1(t) ) * f2(t) 
% g2(t) = ( 1 - f1(t) ) * f3(t) 

Ntau = length(tau);

% Nwindow x Ntime
[Nw,Nt] = size(data.y);
% sampling time step
dt = data.dt;

% (# of pattarn) x (# of spike)
[Nc,Ns] = size(data.Dt);

% sampling time 
t  = (1:(Nt)) * (dt) ;

% data.Dt : spike onset time: (# of pattarn) x (# of spike)
Dt0 = data.Dt; 
% Dt : t - Tk : Nt x (Nc*Ns)
Dt = repmat(t(:), [1 Nc*Ns]) - repmat(Dt0(:)', [Nt 1]);

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
	gt = reshape(g0, [Nt Nc Ns]);
	% ---  sum over Ns spikes
	g = reshape(sum(gt,3), [Nt Nc])';
case	3
	gt = reshape(g0, [Nt Nc Ns]);
	% ---  sum over Ns spikes
	g(:,:,1) = reshape(sum(gt,3), [Nt Nc])';
	
	f3 = exp(- Dt / tau(3)) ;
	g1 = (1 - f1) .* f3 .* It;
	
	gt = reshape(g1, [Nt Nc Ns]);
	% ---  sum over Ns spikes
	g(:,:,2) = reshape(sum(gt,3), [Nt Nc])';
end
