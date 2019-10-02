function	stat = spike_state_statics(tau, data, p)
%  spike_state_statics for spike window data
%  stat = spike_state_statics(tau, data)
% --- 
%   data.a  : spike amplitude
%   data.b  : bias
%   data.y  : y(k,n)
%   data.Dt : Dt(k,n) = t(k,n) - T(k)
%           : (# of pattarn) x (# of time sample) x (# of spike)
%    Nstate = (# of pattarn)
% ---  double exponential function
%   g(t) = ( 1 - exp(-t/T1) ) * exp(-t/T2) for t >= 0
% --- error
% stat.E  = (Nstate x 1);
%      E  = sum((y(t) - a*g(t) - b)^2) 
% mean over time in nspike state
% stat.g  = (Nstate x 1);
% stat.gg = (Nstate x 1);
% stat.yg = (Nstate x 1);
if nargin<3, p = [0 0]; end
    
Ntau = length(tau);

%  g : (Nstate x Nt)  or (Nstate x Nt x 2)
g0 = spike_func_evaluate2(tau, data);

switch	Ntau
case	2
	% error for data.y
    %e  = repadd( data.a * g0 , data.b - data.y ) ;
    
    aa = data.a * g0;
    e  = repadd( aa+p(1)*(aa.^2-aa)+p(2)*(aa.^3-aa) , data.b - data.y ) ;
	
	
	% sum over time
	stat.E  = sum(e.^2, 2);
	
	% mean over time (within window)
	stat.y  = mean(data.y, 2);
	
	%  g: (# of pattarn) x 1
	% mean over time
	stat.g  = mean(g0, 2);
	
	stat.gg = mean(g0.^2, 2);
	
	yg = repmultiply( g0 , data.y ) ;
	
	stat.yg = mean(yg, 2);
case	3
	g1 = g0(:,:,1);
	g2 = g0(:,:,2);
	
	% error for data.y
	e  = repadd( data.a(1) * g1 + data.a(2) * g2 , data.b - data.y ) ;
	
	% sum over time
	stat.E  = sum(e.^2, 2);
	
	% mean over time (within window)
	stat.y  = mean(data.y, 2);
	
	%  g: (# of pattarn) x 1
	% mean over time
	stat.g  = [mean(g1, 2), mean(g2, 2)];
	
	stat.gg = [mean(g1.^2, 2), mean(g2.^2, 2), mean(g1.*g2, 2)];
	
	yg1 = repmultiply( g1 , data.y ) ;
	yg2 = repmultiply( g2 , data.y ) ;
	
	stat.yg = [mean(yg1, 2), mean(yg2, 2)];
end

return

%%%%%%% Consistency check
% ----- Deviation over time (within window)
%  g0: (# of pattarn) x (# of time sample) 
dy  = repadd(data.y, - stat.y);
dg  = repadd(g0, - stat.g);

% variance over time (within window)
yg = repmultiply( dg , dy ) ;

stat.dyg = mean( yg, 2);
stat.dgg = mean( dg.^2, 2);


