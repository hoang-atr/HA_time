function	[nlist]= select_list(tau_list,parm)
% ----- grid search list for tau

Nlist = length(tau_list);
Nloop = fix(Nlist/parm.Ncpu); 
Nrest = Nlist - Nloop*parm.Ncpu;

Nnum = zeros(1,parm.Ncpu);

for n=1:parm.Ncpu
	if n <= Nrest
		Nnum(n) = Nloop + 1;
	else
		Nnum(n) = Nloop;
	end
end

nstart = sum(Nnum(1:(parm.job_num-1))) + 1;
nend   = sum(Nnum(1:(parm.job_num)));
nlist  = nstart:nend;


%for n=1:length(nlist)
%	n_tau = nlist(n);
%	tau = tau_list(:,n_tau);
