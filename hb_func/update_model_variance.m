function	[sx, b0] = update_model_variance(tau, data)
%  update variance for spike model
% --- 
%   [sx, b0] = update_model_variance(tau, data);
% --- data in spike window
%   data.a  : spike amplitude
%   data.b  : bias of spike data
%   data.y  : y(k,n) : window data 
%   data.Dt : Dt(k,n) = t(k,n) - T(k)
% --- data in no spike period
%   data.yrest  : window rest data
%   data.Yrest  : outside rest data
% --- time constant of spike response
%  tau = [T1; T2]
% ---  double exponential function
%   g(t) = ( 1 - exp(-t/T1) ) * exp(-t/T2) for t >= 0
% --- error
%  E  = sum((y(t) - a*g(t) - b)^2) + sum((yrest(t) - b0)^2)
%  sx = E/(Ndata + Nrest)

% error of spike state data
[Evec, Ndata] = sq_error_spike_vec(tau, data);
E  = sum(Evec) ;

% error of window rest data
if isfield(data,'yrest') && ~isempty(data.yrest)
	b0 = mean(data.yrest(:));
	Erest = sum((data.yrest(:) - b0).^2);
	Nrest = length(data.yrest(:));
else
	b0 = 0;
	Erest = 0;
	Nrest = 0;
end

sx = (E + Erest)/(Ndata + Nrest);

