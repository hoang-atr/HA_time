function	post_info = estimate_spike_num(post_info)
% --- spike estimation based on posterior for overlap window 
%   post_info = estimated_post_state(state,spike_state,overlap,Data)
% --- Estimated spike info for overlap window
% post_info.Pspike(ns,n) = P(ns | yk) of n-th overlap window

nwin = size(post_info.Pspike,2);

% post_info.spike_num(n)  = spike number within n-th overlap window
post_info.spike_num  = zeros(1,nwin);

% P(Sk , yk) = P(Sk |yk) * P(yk)

for n=1:nwin
	Pspike = post_info.Pspike(:,n);
	P0sum = 1 - sum(Pspike);
	
	[Pmax,ns_opt] = max(Pspike);
	
	% most probable spike number in this period
	if Pmax >= P0sum
		post_info.spike_num(n)  = ns_opt;
	end
end

