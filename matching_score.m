function score = matching_score(data, Model)
N = numel(data.y);

tmodel = 0:0.001:data.t(end);
ymodel = dbl_exp_func(tmodel, Model);

t = 0:data.dt:tmodel(end);
[~,idx] = findClosest(t, tmodel);
temperate = ymodel(idx);

dy = diff( [data.y data.y(end)] );
dt = diff( temperate );
coin = conv(dy, dt);

score = coin(1:N);