function sarkka_compare(obj)

    obj = obj.filter;
    obj = obj.smooth;
    Y = obj.y;
    A = obj.par.A;
    Q = obj.par.Q;
    H = obj.par.H;
    R = obj.par.R;
    m0 = obj.par.x0.mu;
    P0 = obj.par.x0.sigma;
    T = obj.d.T;
    
    [kf_m, kf_P] = ds.utils.sarkka_kalman(Y, A, Q, H, R, m0, P0);
    
    [rts_m, rts_P] = ds.utils.sarkka_rts(Y, kf_m, kf_P, A, Q);
    
    delta_filter_m = kf_m - obj.infer.filter.mu;
    delta_smooth_m = rts_m - obj.infer.smooth.mu;
    
    fprintf('Difference in filter mean: norm = %.5e\n', norm(delta_filter_m)./T);
    fprintf('Difference in smooth mean: norm = %.5e\n', norm(delta_smooth_m)./T);
    
    for tt = 1:T
        delta_filter_P = obj.infer.filter.sigma{tt} - squeeze(kf_P(:,:,tt));
        delta_smooth_P = obj.infer.smooth.sigma{tt} - squeeze(rts_P(:,:,tt));
        fprintf('filter P delta: %.7f, smooth P delta %.7f', norm(delta_filter_P), norm(delta_smooth_P));
        evs = sort(eig(obj.infer.filter.sigma{tt}), 'ascend');
        fprintf('   ---   min 2 ev.s: obj.filter.P (%.5f, %.5f)\n',  evs(1), evs(2))
    end
end