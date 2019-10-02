function	spike_time = generate_poisson_spike(parm)
% ----- Poisson spike generation 
%    spike_time = generate_poisson_spike(parm)
%
% spike_time = time of spike firing [sec]  (Nspike x 1)
% ----- Poisson spike generation setting
% parm.Nspike   = number of spikes
% parm.fs_spike = spike freq [Hz]

%�|�A�\���ߒ��ɏ]���Ƃ��A
%inter-spike interval(ISI)�͎w�����z�ɏ]��
%�w�����z�ɏ]���T���v���́A��l�����̑ΐ������΂悢

MinISI = parm.MinISI;
fs_spike = parm.fs_spike;
Nspike = parm.Nspike; 

% --- Generation of inter-spike interval (ISI)
% 1. uniform distribution sample generation
ISI = rand( Nspike, 1);
ISI = max(ISI, eps);

% ISI of Poisson spike
% 2. exponential distribution sample
ISI = -(1./fs_spike) .*log(ISI);
ISI = max(ISI ,MinISI);

% 3. firing time = cumulative sum of ISI
% time of spike firing [sec]
spike_time = cumsum(ISI) + 2./fs_spike;
