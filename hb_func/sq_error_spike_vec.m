function	[Evec, Ndata] = sq_error_spike_vec(tau, data)
% square error vector for spike windows
%  [Evec, Ndata] = sq_error_spike_vec(tau, data)
% --- 
%   data.a  : spike amplitude
%   data.b  : bias
%   data.y  : y(k,n)
%   data.Dt : Dt(k,n) = t(k,n) - T(k)
% ---  double exponential function
%   g(t) = ( 1 - exp(-t/T1) ) * exp(-t/T2) for t >= 0
% --- error
% E = sum((y(t) - a*g(t) - b)^2) 

Ntau = length(tau);

% evaluate spike function

if iscell(data.Dt)
% --- time constant of spike response
%  tau = [É—1; É—2]   or  [É—1; É—2; É—3] 
%  g   : (Nwin x Nt)  or (Nwin x Nt x 2)
	g = spike_func_evaluate(tau, data);
	
	% error for data.y
	if Ntau==2
		e = data.y - ( data.a * g + data.b ) ;
	elseif Ntau==3
		e = data.y - ( data.a(1) * g(:,:,1) + data.a(2) * g(:,:,2) + data.b );
	end
else
% --- time constant of spike response
%  tau = [É—1; É—2]    or  [É—1; É—2; É—3] 
%  g   : (Nconf x Nt)  or (Nconf x Nt x 2)
	g = spike_func_evaluate2(tau, data);
	
	% error for data.y
	if Ntau==2
		e  = repadd( data.a * g , data.b - data.y ) ;
	elseif Ntau==3
		e  = repadd(data.a(1) * g(:,:,1) + data.a(2) * g(:,:,2) , ...
				data.b - data.y ) ;
	end
end

Evec  = sum(e.^2, 2);

[Nw,Nt] = size(e);
Ndata = Nw*Nt;
