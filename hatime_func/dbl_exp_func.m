function x = dbl_exp_func(t, model)
if length(model.tau)==2
    x = model.a*(1-exp(-t/model.tau(1))).*exp(-t/model.tau(2));
else
    x = model.a(1)*(1-exp(-t/model.tau(1))).*(exp(-t/model.tau(2))+model.a(2)*exp(-t/model.tau(3)));
end