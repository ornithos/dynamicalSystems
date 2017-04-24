% Setup globals
    cols          = zeros(3,3);
    cnums         = [2 4 5];
    for kk = 1:3; cols(kk,:) = utils.plot.varyColor2(cnums(kk)); end
    litecols      = utils.plot.colortint(cols, 0.8);
    cols          = utils.plot.colortint(litecols, 1.4);
    pforward      = 150;


%% Main parameter set
    R      = 0.2;
    ss     = 4;        % state space
    os     = 3;        % observation space
    opts   = struct('warnings', true);
    
    % initialise randomly / DA   --- SSID more complex for (a) control, (b)
    % NaN and (c) multiple time series. So not doing it.
    A        = diag(rand(ss, 1));
    H        = rand(os, ss)*2 - 1;
    bias     = nanmean(cell2mat(cellfun(@(x) subsref([x; NaN], struct('type','()','subs',{{1:3}})), cellfun(@(x) nanmean(x(:,5:10), 2), Ys, 'Un', 0), 'Un', 0)'),2);
    
    dsBatch  = ds.dynamicalSystemBatch(ss, os, 'x0', zeros(ss,1), eye(ss)*1e-6, ...
        'evolution', A, 1, 'emission', H, R, bias, ...
        'data', Ys, 'control', Us, true, false, opts);
    %%
    epsilon = 0.01;
    opts = struct('maxiter', 500, 'epsilon', 5e-4, 'sampleStability', 5, ...
           'multistep', 4, 'ssid', false, 'verbose', true, 'stableVerbose', false, ...
           'annealingSchedule', 0.1, 'annealingIter', 100, 'annealingMin', 1e-3, ...
           'diagA', true, 'diagQ', false, 'diagAconstraints', [-1+epsilon, 1-epsilon], 'fixBias2', true, 'strictNegativeCheck', false);
    llhhist = dsBatch.parameterLearningEM(opts);
    dsBatch.filter([],true);
    dsBatch.smooth;
    dsBatch.save('initial');
    
    
    % visualise
    figure
    
    [impy, impPy] = dsLDS.impute_y('variance', true, 'smooth', true);
    futurY        = dsLDS.getPredictFreeRun(dsLDS.d.T, pforward);
    plot(Ys{ii}', ':'); hold on;
%     plot(impy'); plot(dsLDS.d.T+1:dsLDS.d.T+pforward, futurY); hold off;
    plot(impy'); plot(1:dsLDS.d.T+pforward, futurY); hold off;
    stdy          = sqrt(cell2mat(cellfun(@(x) diag(x)', impPy, 'Un', 0)))';
    for jj = 1:os
        nanidx          = find(isnan(dsLDS.y(jj,:)));
        if isempty(nanidx); continue; end
        nanminmax       = nanidx(1):nanidx(end);
        utils.plot.confidenceInterval(nanminmax, impy(jj,nanminmax), stdy(jj,nanminmax), [], 'facecolor', cols(kk,:));
    end
    cmap           = flipud([1 0.86 0.91; 0.86 0.91 1]);
    utils.plot.dataShadeVertical(1:(dsLDS.d.T), pump.target{ii}', cmap, 'edgedetect', 'manual', 'edgemanual', 3);
    
    dsNaNCell{ii} = dsLDS.copy;

