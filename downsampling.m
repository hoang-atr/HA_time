function [data] = downsampling(data, Ndown)

if Ndown<=1 
    data.downsampled = 0;
    return; 
end

data.downsampled = 1;
data.t_raw = data.t;
data.y_raw = data.y;

idx = Ndown:Ndown:numel(data.t);

data.t = data.t(idx);
data.y = data.y(idx);
data.dt = mean(diff(data.t(:)));
data.fs = 1/data.dt;