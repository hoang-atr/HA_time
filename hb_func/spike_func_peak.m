function	[tpeak, gpeak, tdecay] = ...
				spike_func_peak(tau, decayRate)
%  peak of double exponential spike function
%  [tpeak, gpeak, tdecay] = spike_func_peak(tau,decayRate)
%  
% --- time constant of spike response
%  tau = [É—1; É—2]    or  [É—1; É—2; É—3]   (T2 < T3)
%  g   : (Nconf x Nt)  or  (Nconf x Nt x 2)
% ---  double exponential function for [É—1; É—2]
% g(t) = ( 1 - exp(-t/É—1) ) * exp(-t/É—2) for t >= 0
%      = ( 1 - f1(t) ) * f2(t) for t >= 0
% ---  double exponential function 
%      with two decay constant for [É—1; É—2; É—3] 
% g1(t) = ( 1 - exp(-t/É—1) ) * exp(-t/É—2) for t >= 0
% g2(t) = ( 1 - exp(-t/É—1) ) * exp(-t/É—3) for t >= 0
%
% g1(t) = ( 1 - f1(t) ) * f2(t) 
% g2(t) = ( 1 - f1(t) ) * f3(t) 
%
%   g(t) represents Ca chemical response
%   or neuron Voltage response
%  1. ( 1 - exp(-t/T1) ) represents Ca molecule 
%                        or Na current inflow
%  2. (a1 * exp(-t/T2) + a2 * exp(-t/T3) )   
%          represents Ca molecule or K current outflow
%
% --- peak amplitude
%  dg/dt = (1/T1)*exp(-t/T1) )*exp(-t/T2) ...
%        - (1/T2)*( 1-exp(-t/T1) )*exp(-t/T2)
%        = [(1/T1)*exp(-t/T1) ) ...
%        - (1/T2)*( 1-exp(-t/T1) ))]*exp(-t/T2)
%        = [- (1/T2) + (1/T1 + 1/T2) *exp(-t/T1)]...
%          *exp(-t/T2)
%        = 0
%  exp(-t/T1) = (1/T2)/(1/T1 + 1/T2) = T1/(T1 + T2)
%  tpeak = T1 * log( (T1 + T2)/T1 )
%  gpeak = g(tpeak)
% --- tdecay : g = gpeak/decayRate  at tdecay
%  decayRate = [100]
%  exp(-tdecay/T2) = gpeak/decayRate
%  tdecay = T2*log(decayRate/gpeak)

tpeak = tau(1,:) .* log( (tau(1,:) + tau(2,:))./tau(1,:) );

if nargout < 2, return; end;

gpeak = ( 1 - exp(-tpeak./tau(1,:)) ) ...
		.* exp(-tpeak./tau(2,:));

if nargout < 3, return; end;

if nargin == 1, decayRate = 100; end;

tdecay  = tau(2,:).*log(decayRate./gpeak);

