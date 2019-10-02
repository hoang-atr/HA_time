function	post_info = estimate_spike_num_weight(post_info,Model)
% --- 
%   post_info = estimate_spike_num_weight(post_info,Model)
% --- 
% post_info.Pspike(ns,n) = P(ns | yk) of n-th overlap window

X = post_info.Pspike;
X = [X; ones(1,size(X,2))];

Y = Model.W * X;

[ytmp, Yid] = max(Y,[],1);

post_info.spike_num = Yid - 1;

