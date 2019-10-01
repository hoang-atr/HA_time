function data = estimate_spike_time_linear(data, model, parm)

% acquiring data
class_label = data.class_label;     % classified label by SVM

% select spikes detected by SVM
str_index = data.index(class_label==1);

t = data.t;
y = data.y;    

if strcmp(parm.sort_dir, 'ascend') || strcmp(parm.sort_dir, 'descend')     
    data.est_spike_time = estimate_hyperacuity_spike(t, y, str_index, ...
        parm.sort_dir, parm.length0, parm.dt0, parm.Npre0, parm.Nwin0, model);
    data.est_spike_time = sort(data.est_spike_time, 'ascend');
else
    fw = estimate_hyperacuity_spike(t, y, str_index, 'ascend', ...
        parm.length0, parm.dt0, parm.Npre0, parm.Nwin0, model);
    
    bw = estimate_hyperacuity_spike(t, y, str_index, 'descend', ...
        parm.length0, parm.dt0, parm.Npre0, parm.Nwin0, model);
    bw = sort(bw, 'ascend');
    
    data.est_spike_time = (fw + bw)./2;
end

end