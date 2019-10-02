function [data] = generate_ca_response(parm)

% generate spike train
spike_train = generate_poisson_spike(parm);
duration = max(spike_train) + 2./parm.fs_spike;

t_raw = 1/parm.fs_raw:1/parm.fs_raw:duration;
t = 1/parm.fs:1/parm.fs:duration;
[~,idxList] = findClosest(t,t_raw);

% elementary Ca transient
model.tau = [parm.tauOn parm.tauOff];
model.a = parm.amp;
ca_transient = dbl_exp_func(t_raw, model);

spkVector = zeros(1,numel(t_raw));
for k = 1:numel(spike_train)
    [~,idx] = min(abs(spike_train(k)-t_raw));
    spkVector(idx) = spkVector(idx)+1;
end

% raw DFF
dffConv = conv(spkVector,ca_transient);
dff = dffConv(1:length(t_raw));

ix = dff>1;     % nonlinearity occurs only when dff>1
dff(ix) = dff(ix).^parm.p;

% add noise to raw DFF
PeakA = parm.amp .* (parm.tauOff./parm.tauOn.*(parm.tauOff./parm.tauOn+1).^-(parm.tauOn./parm.tauOff+1));
sdnoise = PeakA/parm.SNR;
noise = sdnoise.*randn(1,length(t_raw));

noisy_dff = dff + noise;

% DFF at the frame rate (low resolution)
y = noisy_dff(idxList);

% output data
data.t = t;
data.y = y;
data.spike_time = spike_train;
data.fs = parm.fs;
data.Nspike = numel(spike_train);
data.dt = 1/parm.fs;

