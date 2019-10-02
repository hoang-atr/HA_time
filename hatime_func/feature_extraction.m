function [X, I, L] = feature_extraction(data, spk_model, parm)

% linear matching score
l_score = matching_score(data, spk_model);
%l_score = [0 diff(data.y)];

% nonlinear matching score
ref_sig = l_score;
svm_threshold = mean(ref_sig)+std(ref_sig)*parm.threshold;

      
flag_refrac = 0;
cand_idx = [];
cand_value = [];
for i = 1:numel(ref_sig)
    if (ref_sig(i) >= svm_threshold) && (flag_refrac <= 0)
        cand_value = [cand_value; ref_sig(i)];
        cand_idx = [cand_idx; i];

        flag_refrac = parm.Nrefrac;
    else
        flag_refrac = flag_refrac - 1;
    end
end

% sort the candidates by amplitude
[~,sort_idx] = sort(cand_value,'descend');
select_idx = cand_idx(sort_idx);

[X, I, L] = do_sampling(data, l_score, select_idx, parm.Npre, parm.Nwin);    

