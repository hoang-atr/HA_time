function [hit, miss, false_alarm, spike_time_err, spike_time, est_spike_time] = measure_perf(spike_time, est_spike_time, constraint_interval)

%% MAIN
hit = 0;
Nspike = length(spike_time);
Nest_spike = length(est_spike_time);

if Nest_spike>0
    [ste1] = findClosest(est_spike_time', spike_time');
    [ste2] = findClosest(spike_time', est_spike_time');
    spike_time_err = [ste1(:); ste2(:)];
else
    spike_time_err = [];
end

% spike_time_err = zeros(1,1);

while (Nspike>0) && (Nest_spike>0)
    [dt, index]=findClosest(spike_time',est_spike_time');
    [min_dt, min_index] = min(dt);
    
    t = est_spike_time(index(min_index)) - spike_time(min_index);
       
    if min_dt < constraint_interval
        hit = hit + 1;        
%         spike_time_err(hit,1) = abs(t);
        
        spike_time(min_index) = [];        
        Nspike = Nspike - 1;
        
        est_spike_time(index(min_index)) = [];
        Nest_spike = Nest_spike - 1;
    else
        break;
    end
end

miss = Nspike;
false_alarm = Nest_spike;