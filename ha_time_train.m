function [ classifier, model, feature, label ] = ha_time_train( data, model, parm )
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

% estimate the model if not provided by users
if isnull(model)
    parm_model.Ntau = 2;
    parm_model.plot_mode = 0;
    parm_model.decay = 1;
    parm_model.Twin    = 1.5;
    parm_model.Tpre    = 0.5; 
    
    % estimate the Ca spike model
    [model] = train_shape(data, parm_model);
    
    % estimate the nonlinearity
    if ~strcmp(parm.nonlinear,'linear')
        model = estimate_nonlinear_model(data, model, parm.nonlinear);
    end
end
    
% setting up the parameter if needed
if nargin<3             % default parameters
    parm.thr = 0;
    parm.Npre = 1;
    parm.Npost = 4;
    parm.Nrefrac = 0;
end
if ~isfield(parm,'thr'), parm.thr = 0; end
if ~isfield(parm,'Npre'), parm.Npre = 1; end
if ~isfield(parm,'Npost'), parm.Npost = 4; end
if ~isfield(parm,'Nrefrac'), parm.Nrefrac = 0; end

%% Preprocessing the data
% compensate for nonlinearity
if ~strcmp(parm.nonlinear,'linear')
    data = nl2lin_transform(data, model);
end

%% Feature extraction and train the classifier (SVM)
[feature, ~, label] = feature_extraction(data, model, parm); 
classifier = fitcsvm(feature, label, 'Standardize', 1);


