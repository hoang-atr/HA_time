function [X, I, L] = do_sampling(data, score, index, Npre, Npoint) 

Npost = Npoint-Npre-1;
N = length(data.y);

spike_time =  data.spike_time;
Nspike = numel(spike_time);

t = data.t;
dt = data.dt;

spk_flg = zeros(Nspike,1);

index( (index>=N-Npoint) | (index<=Npre) ) = [];
Nsegment = length(index);

X = zeros(Nsegment,2*Npoint);
L = zeros(Nsegment,1);
I = zeros(Nsegment,1);

for i = 1:Nsegment
    str_ix = index(i)-Npre;
    end_ix = index(i)+Npost;
    
    X(i,:) = [data.y(str_ix:end_ix) score(str_ix:end_ix)];
    I(i,1) = str_ix;     % sampling index
        
    if ~isempty(spike_time)
        spk_dist = spike_time - t(index(i)-Npre);
        spk_dist( (spk_dist<0) | (spk_flg>0) ) = inf;            
        [min_dist,min_spk] = min(spk_dist);

        if min_dist <= (Npoint-1)*dt
            L(i) = 1;
            spk_flg(min_spk) = 1;
        else
            L(i) = 0;
        end      
    end                              
end