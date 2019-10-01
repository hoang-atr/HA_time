function [XX, LL, II] = feature_extraction(data, spk_model, parm)

XX = [];
LL = [];
II = [];

Ndata = numel(data);
for n = 1:Ndata
    [X,L,I] = feature_extraction_linear(data{n}, spk_model, parm.thr_svm, parm.Nrefrac, parm.Npre, parm.Nwin);
    
    XX = [XX; X];
    LL = [LL; L];
    II = [II; I];
end