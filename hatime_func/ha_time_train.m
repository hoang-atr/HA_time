function [ classifier, model, parm, feature, label ] = ha_time_train( data, parm, model )
% Descriptions for the ha_time_train
% + Input
%     - data: training data
%         - data.t [1 x Nt]: sampling time vector
%         - data.y [1 x Nt]: observed Ca response
%         - data.spike_time [1 x Ns]: ground-truth spikes
%     - parm: parameter for training HA_time classifier
%         - parm.Ndown: downsampling parameter (default, 1 - no downsample data)
%         - parm.nonlinear: nonlinearity flag ('linear', 'saturation', 'super-linear')
%         - parm.thr: threshold (in term of SD) for selecting spike candidates
%         - parm.Npre: number of preceding "points" before the point that exceeds the threshold
%         - parm.Npost: number of succeeding "points"
%         - parm.Nrefrac: number of refractory "points"
% + Output
%     - classifer: the trained SVM classifer
%     - model: the Ca spike model estimated by EM algorithm
%     - feature [Nsample x Nfeature]: the feature vector extracted 
%     - label: the label (0: non-spike, 1-spike)

% checking the required data
if ~isfield(data,'t') || ~isfield(data,'y') || ~isfield(data,'spike_time')
    print('The required data is lacking!\n')
    return
end

% setting up the parameter if needed
if nargin<2, parm = set_default_param(); end

% estimate the model if not provided by users
if nargin<3
    % estimate the Ca spike model
    parm.fs = data.fs;
    [model] = train_shape(data, parm);
    
    % compensate for the nonlinearity
    if ~strcmp(parm.nonlinear,'linear')
        % estimate the nonlinearity model
        model = estimate_nonlinear_model(data, model, parm.nonlinear);        
        
        % transform the signal
        data = nl2lin_transform(data, model);
    end
end

% optimize the SVM param if required
if parm.optimize > 0
    parm = optimize_svm_param(data, parm, model);
end
    
% Feature extraction and train the classifier (SVM)
[feature, ~, label] = feature_extraction(data, model, parm); 
classifier = fitcsvm(feature, label, 'Standardize', 1);


