function data = nl2lin_transform(data, Model, range, thr)
if nargin<3, range = [1.5 max(data.y)]; end
if nargin<4, thr = max(data.y)*3; end

if isempty(Model.nl_model)
    disp('No nonlinear model found!')
    return
end

ix = data.y>=range(1) & data.y<=range(2);

% transform data
ty = data.y;
ty(ix) = real(feval(Model.nl_model, data.y(ix)));

% just to make sure nothing weird happens
ty(ty>thr) = thr;

% save data
data.y = ty;
