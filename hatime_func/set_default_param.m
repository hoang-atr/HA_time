function parm = set_default_param(parm)
if nargin<1, parm = struct; end

%% parameters for estimating the spike model
% number of time constant, 2 or 3
if ~isfield(parm,'Ntau'), parm.Twin = 2; end

% decay rate [s]
if ~isfield(parm,'decay'), parm.decay = 1; end

% sampling window in second
if ~isfield(parm,'Twin'), parm.Twin = 1.5; end
if ~isfield(parm,'Tpre'), parm.Tpre = 0.5; end

%% parameters for sampling spike candidates for SVM
% optimize the parameter by f1-score for training data or not?
if ~isfield(parm,'optimize'), parm.optimize = 0; end

% threshold of coincidence score
if ~isfield(parm,'threshold'), parm.threshold = 0; end

% #point of refractory period
if ~isfield(parm,'Nrefrac'), parm.Nrefrac = 0; end

% #point of sampling window
if ~isfield(parm,'Npre'), parm.Npre = 1; end
if ~isfield(parm,'Nwin'), parm.Nwin = 4; end


%% parameters for calibrating hyperacuity spike time
% acquired hyperacuity interval [s]
if ~isfield(parm,'dt0'), parm.dt0 = 0.002; end

% searching window (from the sampling point) [s]
if ~isfield(parm,'length0'), parm.length0 = 0.1; end  

% sampling point to calculate the residual
if ~isfield(parm,'Npre0'), parm.Npre0 = 2; end        
if ~isfield(parm,'Nwin0'), parm.Nwin0 = 10; end

% searching direction: 'ascend' or 'descend'
if ~isfield(parm,'sort_dir'), parm.sort_dir = 'ascend'; end  

%% others
% nonlinearity: 'linear', 'saturation' or 'superlinear'
if ~isfield(parm,'nonlinear'), parm.nonlinear = 'linear'; end  