function [data] = generate_ca_response(parm, maxd, ifp)
if nargin<2, maxd = 0.01; end;
if nargin<3, ifp = 0; end;
    
spike_train = generate_poisson_spike(parm);

delay = maxd.*rand(parm.Nspike,1);
delayed_spike_train = spike_train + delay;

% calculate noise param
PeakA = parm.amp .* (parm.tauOff./parm.tauOn.*(parm.tauOff./parm.tauOn+1).^-(parm.tauOn./parm.tauOff+1));
sdnoise = PeakA/parm.SNR;

duration = max(delayed_spike_train) + 2./parm.fs_spike;

t_raw = 1/parm.fs_raw:1/parm.fs_raw:duration;
t = 1/parm.fs:1/parm.fs:duration;
[~,idxList] = findClosest(t,t_raw);

% elementary Ca transient
model.tau = [parm.tauOn parm.tauOff];
model.a = parm.amp;
ca_transient = dbl_exp_func(t_raw, model);

spkVector = zeros(1,numel(t_raw));
for k = 1:numel(delayed_spike_train)
    [~,idx] = min(abs(delayed_spike_train(k)-t_raw));
    spkVector(idx) = spkVector(idx)+1;
end

% raw DFF
dffConv = conv(spkVector,ca_transient);
dff = dffConv(1:length(t_raw));

ix = dff>1;     % nonlinearity occurs only when dff>1
dff(ix) = dff(ix).^parm.p;

% add noise to raw DFF
%noise = sdnoise.*randn(1,length(t_raw));
cn = dsp.ColoredNoise('InverseFrequencyPower',ifp,'SamplesPerFrame',numel(t_raw));
noise = sdnoise*cn()';

noisy_dff = dff + noise;

% DFF at the frame rate (low resolution)
y = noisy_dff(idxList);

data.t = t;
data.y = y;
data.spike_time = spike_train;
data.fs = parm.fs;
data.Nspike = numel(spike_train);
data.dt = 1/parm.fs;

