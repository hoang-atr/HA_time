function	[a,b] = update_spike_amplitude(tau, data)
%  update_spike_amplitude & bias for spike model
% --- 
%   [a,b]  = update_spike_amplitude(tau, data);
% --- data in spike window
%   data.y(k,:) : y(k,t) : observed data in k-th window
%   data.Dt{k}  : spike onset in  k-th window [sec] (1 x Nspike)
%   data.dt     : sampling time step [sec]
% --- time constant of spike response
%  tau = [T1; T2]
% ---  double exponential function
% g(t) = ( 1 - exp(-t/T1) ) * exp(-t/T2) for t >= 0
% --- spike state amplitude & bias estimation
% E = sum((y(t) - a*g(t) - b)^2) 
% --- spike state parameter
%
% a*<g>    + b      = <y>      : dE/db = 0
% a*<g^2>  + b*<g>  = <y*g>    : dE/da = 0
% --- spike state amplitude & bias estimation
% E = sum((y(t) - a1*g1(t) - a2*g2(t) - b)^2) 
% a1*<g1>    + a2*<g2>    + b      = <y>  
% a1*<g1^2>  + a2*<g1*g2> + b*<g1> = <g1*y>  
% a1*<g1*g2> + a2*<g2^2>  + b*<g2> = <g2*y>  
%
% b = <y> - a*<g>
% a = (<g*y> - <g>*<y>)/( <g^2> - <g>^2 )
%   = <(g-<g>)*(y-<y>)>/<(g-<g>)^2>
% --- non spike state bias estimation
% E0 = sum(y(t) - b0)^2)
% --- spike state amplitude & bias estimation
% E = sum((y(t) - a1*g1(t) - a2*g2(t) - b)^2) 
% a1*<g1>    + a2*<g2>    + b      = <y>  
% a1*<g1^2>  + a2*<g1*g2> + b*<g1> = <g1*y>  
% a1*<g1*g2> + a2*<g2^2>  + b*<g2> = <g2*y>  

Ntau = length(tau);

% evaluate spike function
g = spike_func_evaluate(tau, data);

switch	Ntau
case	2
% < < (y-<y>)^2 >_t >_k 
% = < < ( (y-<y>_t) - (<y> - <y>_t) )^2>_t >_k 
% = < < (y-<y>_t)^2 >_t >_k + <<(<y> - <y>_t) )^2>_t >_k
	
	% --- mean over time & window
	gm  = mean( g(:) );
	ym  = mean( data.y(:) );
	
	% deviation within window
	gd  = g - gm;
	yd  = data.y - ym;
	
	gg = mean( gd(:).^2 );
	yg = mean( yd(:) .* gd(:) );
	
	% --- update spike amplitude & bias
	a = yg/max(gg, eps);
	%b = ym - a * gm;
	b = 0;
    
    if isnan(a), a=1; end

case	3
	% --- spike state amplitude & bias estimation
	% E = sum((y(t) - a1*g1(t) - a2*g2(t) - b)^2) 
	% a1*<g1^2>  + a2*<g1*g2> + b*<g1> = <g1*y>  
	% a1*<g1*g2> + a2*<g2^2>  + b*<g2> = <g2*y>  
	% a1*<g1>    + a2*<g2>    + b      = <y>  
	
	g1 = g(:,:,1);
	g2 = g(:,:,2);
	% --- mean over time & window
	gm1 = mean( g1(:) );
	gm2 = mean( g2(:) );
	
	gg1 = mean( g1(:).^2 );
	gg2 = mean( g2(:).^2 );
	gg3 = mean( g1(:).* g2(:) );
	
	gy1 = mean( g1(:).* data.y(:));
	gy2 = mean( g2(:).* data.y(:));
	ym  = mean( data.y(:) );
	
	gg = [gg1, gg3, gm1; 
	      gg3, gg2, gm2; 
	      gm1, gm2, 1];
	ab = pinv(gg) * [gy1; gy2; ym];
	a  = ab(1:2);
	b  = ab(3);
end
return

%if isfield(data,'yrest') && ~isempty(data.yrest)
%	b0 = mean(data.yrest(:));
%else
%	b0 = 0;
%end

