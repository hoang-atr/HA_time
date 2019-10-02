% gamma_spike_generator
%定常過程、Cvとシェイプパラメータκの関係は、
%Cv=1/sqrt(κ) : κ= 1/Cv^2

DEBUG = 1;

fs = 100; % sampling rate [Hz]
Nspike = 100; % # of spikes
firing_rate = 100; % mean firing rate [Hz]
Cv = 0.5;

shape_parm = 1/Cv^2; % shape parameter of gamma process

ISI = rand(Nspike-1,1);
ISI = gammaincinv(ISI,shape_parm)./(firing_rate*shape_parm);

% spike_time
spike_time = [1/firing_rate; cumsum(ISI)]; % output

if DEBUG==1
	fr = 1/mean(ISI) % mean firing rate
	Cv = std(ISI)/mean(ISI) % coefficient of variation
	minISI = min(ISI)
	
	% spike signal
	t = ceil(spike_time*fs);
	T = max(t)+1;
	y = zeros(1,T);
	
	for n=1:length(t)
		y(t(n)) = y(t(n)) + 1;
	end
	
	plot(y)
	ylim([-0.5 (max(y)+0.5)])
end
