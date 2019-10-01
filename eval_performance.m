function perf = eval_performance(spike_time, est_spike_time, constraint_interval, cost)
if nargin < 3, constraint_interval = 0.05; end
if nargin < 4, cost = 1; end

perf.spike_time = spike_time;
perf.est_spike_time = est_spike_time;

perf.spk_distance = spkd(spike_time, est_spike_time, cost);
%perf.spk_distance = compute_normalized_dist(spike_time, est_spike_time, cost);

% detection rate
[hit, miss, false_alarm, spk_time_err, spike_time, est_spike_time] = ...
    measure_perf(spike_time, est_spike_time, constraint_interval);
if hit > 0
    perf.sensitivity = hit / (hit + miss);
    perf.precision = hit / (hit + false_alarm);
    perf.f1_score = 2 * (perf.precision * perf.sensitivity) / (perf.precision + perf.sensitivity);
else
    perf.sensitivity = 0;
    perf.precision = 0;
    perf.f1_score = 0;
end  


perf.TP = hit;
perf.FN = miss;
perf.FP = false_alarm;
perf.spk_time_err = spk_time_err;