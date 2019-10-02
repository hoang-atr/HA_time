function	idiff = vb_setdiff2(inew,iold)
% Fast calculation of set difference 'setdiff' for interger index
% idiff = vb_setdiff2(inew,iold)
% --- Input
% inew , iold : Integer index array
% --- Output
% idiff : index in 'inew' which is not included in 'iold'
%
% inew �����Ǥ� iold �����ǤǤʤ��ͤ� idiff �˽���
% inew , iold , idiff ������
%
% Ver 1.0 written by M. Sato  2003-3-15
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

inew  = inew(:);
iold  = iold(:);

N	  = max([inew ; iold]);

jflag = zeros(N,1);
imask = ones(N,1);

jflag(inew) = 1; 	% inew �Υ���ǥå����� 1 �Υե饰��Ω�Ƥ�
imask(iold) = 0;	% iold �Υ���ǥå����� 0 �Υޥ���

jflag = jflag.*imask; 

idiff = find(jflag==1);
