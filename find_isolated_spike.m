function Segment = find_isolated_spike(data, parm)    
    N = numel(data.y);

    isi_pre = [0; diff(sort(data.spike_time))];
    isi_post = [diff(sort(data.spike_time)); 0];
    ix = isi_pre > parm.MinISI & isi_post > parm.MinISI;
        
    single_spike_time = data.spike_time(ix);
    Nsegment = numel(single_spike_time);
    
    Segment = cell(Nsegment,1);
    
    for i = 1:Nsegment
        [~,samp_n] = min(abs(single_spike_time(i)-data.t));
        
        str_sampling = samp_n - parm.Npre;
        end_sampling = samp_n + parm.Npost;
        if str_sampling < 1, str_sampling = 1; end;
        if end_sampling > N, end_sampling = N; end;              
               
        Segment{i}.t = data.t(str_sampling:end_sampling)-single_spike_time(i);
        Segment{i}.y = data.y(str_sampling:end_sampling);
    end        
end