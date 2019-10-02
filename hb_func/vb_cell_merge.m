function	[xarray , Ndata] = vb_cell_merge(xcell)
% merge cell array into vector
%  [xarray , Ndata] = vb_cell_merge(xcell)
% --- Input
% xcell : cell array
% ---Output
% xarray : concatenated data vector
% Ndata(n) : number of data in xcell{n}
%
% M. Sato 2006-7-21
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

N = length(xcell);

Ndata = zeros(N,1);

for n=1:N
	Ndata(n) = length(xcell{n});
end

xarray = zeros(sum(Ndata),1);

next = 0;

for n=1:N
	xdata = xcell{n};
	n1 = next + 1;
	n2 = next + Ndata(n);
	xarray(n1:n2) = xdata(:);
	next = n2;
end
