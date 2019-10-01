function [Model] = estimate_nonlinear_model(data, Model, mode, plot_flag)
if nargin<4, plot_flag = 0; end
    
% compute model prediction
t_model = data.t(1):0.001:data.t(end);
y_model = calc_linear_response(data.spike_time, t_model, Model);
[~,idx] = findClosest(data.t,t_model);
y_model = y_model(idx);

% find index of spike segments
%ix = find_spike_ix(data.t, data.spike_time);
ix = data.y>0;

% stack the values
X = y_model(ix);
Y = data.y(ix); 

if strcmp(mode,'saturation')
    fit_options = fitoptions('gauss2', 'Lower',[-Inf,-Inf,-Inf],'Upper',[Inf,Inf,Inf]);
    fit_func = fittype('a/exp(b*X)+c','independent','X','coefficients',{'a','b','c'},'dependent','Y');
elseif strcmp(mode,'super-linear')
    fit_options = fitoptions('gauss2', 'Lower',[0,0],'Upper',[1,1]);
    fit_func = fittype('a*exp(b*X)','independent','X','coefficients',{'a','b'},'dependent','Y');
else
    fit_options = fitoptions('gauss2');
    fit_func = fittype('a*X','independent','X','coefficients',{'a'},'dependent','Y');
end

f = fit(X(:),Y(:),fit_func,fit_options);

if plot_flag>0
    figure
    plot(f,X(:),Y(:))
    hold on        
end

if strcmp(mode,'linear'), return; end

if strcmp(mode,'saturation')
    ff = @(X) f.a/exp(f.b*X)+f.c;
elseif strcmp(mode, 'super-linear')
    ff = @(X) f.a*exp(f.b*X);
end

g = sym(ff);
ffi = finverse(g);
F=matlabFunction(ffi);

Model.nl_model = F;


