function	nhist = spike_hist(fs_spike, Twin, Nmax, Nspike)
% Poisson distribution histogram
%   spike_hist(fs_spike, Twin, Nmax, Nspike)
% ---
% fs_spike : spike firing freq  [Hz]
% Twin     : time window length [sec]
% --- optional
% Nmax     : max number of spikes in Twin
% Nspike   : total number of spikes
% ---
% nhist(n) : histgram for number of spikes in Twin
% --- Prob. n-spike occurs in dt
% P(n) = (fs*dt)^n * exp(-fs*dt)/n!

if nargin < 3, Nmax = 10; end;
if nargin < 4, Nspike = 100; end;

DT = fs_spike * Twin;

n = 0:Nmax;

% Poisson distribution
% P(n) = (fs*dt)^n * exp(-fs*dt) / n!
P = (DT).^n .* exp(-DT) ./ factorial(n);

Psum = sum(n .* P);

nhist = (Nspike) * P / Psum;
