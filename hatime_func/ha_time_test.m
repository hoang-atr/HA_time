function [ data, feature ] = ha_time_test( data, classifier, model, parm )
% Descriptions for the ha_time_test
% + Input
%     - data: test data
%         - data.t [1 x Nt]: sampling time vector
%         - data.y [1 x Nt]: observed Ca response
%     - classifier: the trained classifier
%     - model: the estimated Ca model
%     - parm: parameter for training HA_time classifier
%         - parm.Ndown: downsampling parameter (default, 1 - no downsample data)
%         - parm.nonlinear: nonlinearity flag ('linear', 'saturation', 'super-linear')
%         - parm.thr: threshold (in term of SD) for selecting spike candidates
%         - parm.Npre: number of preceding "points" before the point that exceeds the threshold
%         - parm.Npost: number of succeeding "points"
%         - parm.Nrefrac: number of refractory "points"
% + Output
%     - data: data with 'est_spike_time' at hyperacuity precision
%     - feature [Nsample x Nfeature]: the feature vector extracted 

% checking the required data
if ~isfield(data,'t') || ~isfield(data,'y')
    print('The required data is lacking!\n')
    return
end

% estimate the model if not provided by users
if isempty(classifier) || isempty(model)  
    print('The trained model is required!\n')
    return
end   

%% Preprocessing the data
% compensate for nonlinearity
if ~strcmp(parm.nonlinear,'linear')
    data = nl2lin_transform(data, model);
end

%% Feature extraction and apply the classifier
[feature, index] = feature_extraction(data, model, parm); 

data.index = index;
data.class_label = predict(classifier, feature);
data = estimate_hyperacuity_spike_time(data, model, parm);

