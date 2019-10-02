function	Tau = set_tau_grid_cs(Ntau)
% set_tau_grid_cs
%	parm.tau     = [0.01; 0.5; 1.0];   % [onset time; decay time] [sec]
% ??on was constrained to 2 - 30 ms to avoid steep or shallow onsets. 
% decay constant ??2 was set to a fixed value of 70 ms

if nargin == 0, Ntau = 3; end;

Tau.Ntau = Ntau;

Tau.tau_ini  = [0.01; 0.1];

Tau.Nupdate_ini = 1; % # of iteration of initial search
Tau.Nupdate     = 3; % # of iteration of detailed search

switch	Ntau
case	2
	% initial search parameter
	Tau.max_ini   = [0.1   ; 1.5];   % max value of tau
	Tau.min_ini   = [0.01  ; 0.1];   % minimum value
	Tau.step_ini  = [0.005 ; 0.1];   % grid step
	Tau.range_ini = [0.02  ; 0.5];  % search range of tau
	
	% detailed search parameter
	Tau.max   = [0.1    ; 1.5 ];   % max value of tau
	Tau.min   = [0.01   ; 0.1 ];   % minimum value
	Tau.step  = [0.005  ; 0.05];   % grid step
	Tau.range = [0.01   ; 0.1 ];   % search range of tau
	
case	3
	% initial search parameter
	Tau.max_ini   = [0.1   ; 0.5; 2.0];   % max value of tau
	Tau.min_ini   = [0.01  ; 0.1; 0.5];   % minimum value
	Tau.step_ini  = [0.005 ; 0.1; 0.1];   % grid step
	Tau.range_ini = [0.02  ; 0.5; 0.5 ];  % search range of tau
	
	% detailed search parameter
	Tau.max   = [0.1    ; 0.5 ; 2.0 ];   % max value of tau
	Tau.min   = [0.01   ; 0.1 ; 0.5 ];   % minimum value
	Tau.step  = [0.005  ; 0.05; 0.1 ];   % grid step
	Tau.range = [0.01   ; 0.1 ; 0.5 ];   % search range of tau
end
