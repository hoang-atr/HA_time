function recon_spike = calibrate_spike_time(t, y, str_index, dir, length0, dt0, Npre, Nwin, model)
    Ntp_spike = length(str_index);
    str_index = sort(str_index, dir);
    
    % acquiring spike model
    t0 = 0:dt0:t(end);
    N0 = length(t0);
    
    N = length(y);

    transient = dbl_exp_func(t0, model);

    cur_spkVector = zeros(1,N0);
    recon_spike = zeros(Ntp_spike,1);

    % estimate spike time
    for k = 1:Ntp_spike
        % search length
        t_str = t(str_index(k));    
        t_end = t_str + length0;   

        [~,i_str] = min(abs(t0 - t_str));
        [~,i_end] = min(abs(t0 - t_end));

        % sampling range
        range0 = str_index(k)-Npre;
        range1 = str_index(k)-Npre+Nwin-1;
        if range0 < 1, range0 = 1; end
        if range1 > N, range1 = N; end    
        sampling_range = range0:1:range1;        

        min_err = 1E+100;
        for i = i_str:i_end
            spkVector = cur_spkVector;
            spkVector(i) = spkVector(i)+1;

            ys = conv(spkVector,transient);
            ys = ys(1:N0);

            [~,select_i] = findClosest(t(sampling_range), t0);
            err = sum((y(sampling_range) - ys(select_i)).^2);

            if min_err > err
                min_err = err;
                opt_spkVector = spkVector;
                recon_spike(k) = t0(i);
            end
        end

        cur_spkVector = opt_spkVector;
    end    
end