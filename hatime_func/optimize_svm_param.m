function [param] = optimize_svm_param(data, param, model)

opt_thr = 0;
opt_Nrefrac = 0;
opt_Npre = 1;
opt_Nwin = 4;

max_f1 = 0;
for svm_thr = 0:0.5:2
    for Nrefrac = [0 1]
        for Npre = [1 2 3]
            for Nwin = [4 6 8 10] 
                param.threshold = svm_thr;
                param.Nrefrac = Nrefrac;
                param.Npre = Npre;
                param.Nwin = Nwin;
                
                [X, I, L] = feature_extraction(data, model, param);                 
                if isempty(L), continue, end;
                
                svm_model = fitcsvm(X, L, 'Standardize', 1);
                
                data.index = I;
                data.class_label = predict(svm_model, X);
                
                data = estimate_hyperacuity_spike_time(data, model, param);
                perf = eval_performance(data.spike_time, data.est_spike_time);                                

                if max_f1 < perf.f1_score
                    max_f1 = perf.f1_score;
                    
                    opt_thr = svm_thr;
                    opt_Nrefrac = Nrefrac;
                    opt_Npre = Npre;
                    opt_Nwin = Nwin;
                end
            end
        end
    end
end

% update the optimized svm param
param.svm_thr = opt_thr;
param.Nrefrac = opt_Nrefrac;
param.Npre = opt_Npre;
param.Nwin = opt_Nwin;