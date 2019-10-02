function		Model = update_spike_model(state)
%  update noise variance & spike amplitude/bias for spike model
% --- 
%   	Model = update_spike_model(state)
% --- Model parameters (update parameter)
% Model.a
% Model.b
% Model.sx

% --- amplitude & bias estimation
% E = sum((y(t) - a*g(t) - b)^2) 
%
% y0 = mean(yall)
%
% a*<g>    + b      = y0       : dE/db = 0
% a*<g^2>  + b*<g>  = <y*g>    : dE/da = 0
%
% b = y0 - a*<g>*(NT/Nall)
% a = (<g*y> - <g>*y0)/( <g^2> - <g>^2*(NT/Nall))

% --- Posterior
% state.P = P(Tk,Sk,Ak=1|yk) = P(yk,Tk,Sk,Ak=1)/P(yk)
%         : (NconfMax x Nspike x Nwin);
% --- Statics value of spike state in k-th window
% state.g : (NconfMax x Nspike x Nwin);

[NconfMax,Nspike,Nwin,Ng]	= size(state.g);

switch	Ng
case	1
	% --- Expection of statics over Tk & Sk w.r.t. P(Tk,Sk,Ak=1|yk)
	%  <f(k)> = sum_Sk sum_Tk f(yk,Tk,Sk) P(Tk,Sk,Ak=1|yk)
	g  = sum(sum(state.g  .* state.P ,1), 2);
	gg = sum(sum(state.gg .* state.P ,1), 2);
	yg = sum(sum(state.yg .* state.P ,1), 2);
	
	% --- Expection of statics over window
	g  = mean(g(:) );
	gg = mean(gg(:));
	yg = mean(yg(:));
	
	% --- Mean of y over all windows (rest data are excluded)
	y0 = state.y0;
	%R  = state.NTwin / state.Ndata ;
	
	% --- update spike amplitude & bias based on window data
	Model.a = (yg - y0*g)/max(gg - g^2 , eps);
	Model.b = y0 - Model.a * g ;
case	2
	% --- Expection of statics over Tk & Sk w.r.t. P(Tk,Sk,Ak=1|yk)
	%  <f(k)> = sum_Sk sum_Tk f(yk,Tk,Sk) P(Tk,Sk,Ak=1|yk)
	PT  = sum(sum(state.P ,1), 2);
	g1  = sum(sum(state.g(:,:,:,1)  .* state.P ,1), 2);
	g2  = sum(sum(state.g(:,:,:,2)  .* state.P ,1), 2);
	gg1 = sum(sum(state.gg(:,:,:,1) .* state.P ,1), 2);
	gg2 = sum(sum(state.gg(:,:,:,2) .* state.P ,1), 2);
	gg3 = sum(sum(state.gg(:,:,:,3) .* state.P ,1), 2);
	yg1 = sum(sum(state.yg(:,:,:,1) .* state.P ,1), 2);
	yg2 = sum(sum(state.yg(:,:,:,2) .* state.P ,1), 2);
	
	% --- Expection of statics over window
	PT  = mean(PT(:) );
	g1  = mean(g1(:) );
	g2  = mean(g2(:) );
	gg1 = mean(gg1(:));
	gg2 = mean(gg2(:));
	gg3 = mean(gg3(:));
	yg1 = mean(yg1(:));
	yg2 = mean(yg2(:));
	
	% --- Mean of y over all windows (rest data are excluded)
	y0 = state.y0;
	%R  = state.NTwin / state.Ndata ;
	
	gg = [gg1, gg3, g1; 
	      gg3, gg2, g2; 
	      g1,  g2,  1];
	ab = pinv(gg) * [yg1; yg2; y0];
	% --- update spike amplitude & bias based on window data
	Model.a = ab(1:2);
	Model.b = ab(3);
end

% --- sum of spike state error over K windows
Esum  = sum(sum(sum(state.E  .* state.P ,1), 2));
Esum0 = sum(state.E0 .* state.P0) ;

% sx = (NK * (<E>_s + <(y-b0)^2>_0) + sum(y_rest-b0)^2))/(NK+Nrest)
Model.sx = (Esum + Esum0 + state.Erest)/state.Ndata;

DEBUG = 0;
if DEBUG==0, return; end;

% --- Check prob. sum
EPS = 1.0e-10;

% P = P(Ak=1|yk) = sum_Sk sum_Tk P(Tk,Sk,Ak=1|yk)
P  = sum(sum(state.P ,1), 2);
Perr = max(abs( state.P0 + P(:) - 1 ));

if Perr > EPS
	fprintf('Prob. sum error = %e\n', Perr)
end
