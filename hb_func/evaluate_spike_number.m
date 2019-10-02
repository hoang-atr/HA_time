function	[spike_num] =  evaluate_spike_number(spike_info)
% evaluate spike number at evaluation sampling frequency
%   [spike_num] =  evaluate_spike_number(spike_info)
% --- Estimated spike info for overlap window
% spike_info.spike_num(n)  = estimated spike number of n-th overlap window
% spike_info.twin{n} = time index of n-th overlap window
% spike_info.tsec = sample time of observed data [sec]
% spike_info.fs   = sampling frequency [Hz]

fs   = spike_info.fs;
Tmax = length(spike_info.tsec);
Nwin = length(spike_info.twin);

% spike_num(j)  = spike number of j-th sample
spike_num  = zeros(1,Tmax);


for n=1:Nwin
	twin = spike_info.twin{n};
	tlen = length(twin);
	
	% spike number per sec
	spike_num(twin) = spike_info.spike_num(n) * (fs/tlen);
end

% smoothing
spike_num = (spike_num ...
           + [spike_num(2:end)  0 ]*0.5 ...
           + [0 spike_num(1:end-1)]*0.5 ) / 2;
spike_num = (spike_num ...
           + [spike_num(2:end)  0 ]*0.5 ...
           + [0 spike_num(1:end-1)]*0.5 ) / 2;
