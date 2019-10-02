function		[Model, Info, Y] = ...
	train_post_weight(post_info, spike_info,Model)

% Sparse classification by probit Model (Covariance method)
%
%   [Model, Info] = ...
%			probit_sparse_cov(X,Yid,Model_probi,parm)
%
% --- Input
%  X   : Input data  ( M x T )
%  Yid : Class id = {1, ..., N}   ( 1 x T )
%  N  =  # of class
%  M  =  # of input
%  T  =  # of data

%
% --- Set learning parameter
%
% --- Number of training iteration for each method
parm.Ntrain = 500;    % Total iteration number for training parm.Npre_train = parm.Ntrain;
parm.Nskip  = 100;	  
% skip steps for display info

% ARD-Sparse parameters
parm.a_min = 1e-10;% Threshold for pruning input dimension
parm.Fdiff = 1e-10;% Threshold for free energy convergence
parm.Prune = 0;	
% = 1 : Prune irrelevant input dimension in training process
parm.Ta0   = 0.0;	% Degree of sparce prior. Ta0 = 0 : Strong sparse prior
parm.a0    = 1.0;	
% Prior weight variance. If Ta0=0, no effect

X = post_info.Pspike;
X = [X; ones(1,size(X,2))];
Yid = spike_info.spike_num+1;

fprintf('Input dim = %d x %d \n',size(X))
fprintf('Class No = %d\n',max(Yid))

% --- Initialization of Model_probi Parameters
Model_probit = [];

% --- Multiclass Model_probi
[Model_probit, Info] = ...
		probit_regularize(X,Yid,Model_probit,parm);
% [Model_probit, Info] = ...
% 		probit_sparse_cov(X,Yid,Model_probit,parm);

Model.W = Model_probit.W;

Y = predict_probit(X, Model);

Info.err = sum(Y~=Yid)/length(Y);
