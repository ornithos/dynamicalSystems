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
    % ---> This bias 
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

    %% Use (covariate-based) individual biases in batch learning (NOT LEARNED!!)
    
    % tblDemog should be loaded into workspace from
    % ionlds.dataprep.loadMLHCdata
    covLkp    = tblDemog;
    pxLkp     = ionlds.utils.getPatientLookup;
    covLkp.ID = pxLkp(:,2);
    covLkp    = sortrows(covLkp, 'ID');
    covBiases = zeros(3,40);
    
    % expected BP for Iranian adults according to Hosseini et al. 2015
    for nn = 1:40
        pxCovs = covLkp(nn,:);
        bp     = ionlds.utils.getBPGivenAgeWeightHeight(pxCovs.age, pxCovs.weight, pxCovs.height, pxCovs.gender{1});
        covBiases(1:2,nn) = [bp.sys(1); bp.dia(1)];
    end
    
    % appears to be significant difference in diastolic BP vs. empirical
    % data --> dia is much lower, particularly given the elevated physiological
    % state assumed in some patients.
    covBiases(2,:) = covBiases(2,:) - (mean(covBiases(2,:)) - bias(2));
    covBiases(3,:) = 95;
    
    dsBatchBiasCov       = dsBatch.copy;
    dsBatchBiasCov.par.c = mat2cell(covBiases, 3, ones(1,40));
    llhhist = dsBatchBiasCov.parameterLearningEM(opts);
    
%% Compare predictions
batchTable = table;
batchTable.orig_1  = ds.utils.getKstepSMSE(dsBatch, 1, 'average', true)';
batchTable.orig_10 = ds.utils.getKstepSMSE(dsBatch, 10, 'average', true)';
batchTable.orig_20 = ds.utils.getKstepSMSE(dsBatch, 20, 'average', true)';
batchTable.orig_fr = ds.utils.getKstepSMSE(dsBatch, Inf, 'average', true)';
batchTable.bias_1  = ds.utils.getKstepSMSE(dsBatchBiasCov, 1, 'average', true)';
batchTable.bias_10 = ds.utils.getKstepSMSE(dsBatchBiasCov, 10, 'average', true)';
batchTable.bias_20 = ds.utils.getKstepSMSE(dsBatchBiasCov, 20, 'average', true)';
batchTable.bias_fr = ds.utils.getKstepSMSE(dsBatchBiasCov, Inf, 'average', true)';
%%
% These are better than in abProjectWork1-6 (e.g.) table, since we
% normalise by 'allVar'. However, even using remVar, things are
% considerably better for 10-step, and a little better for 20-step.

normaliseFactor = 'allVar';
% normaliseFactor = 'remVar';
invdlTable = table;
invdlTable.pred_1 = (1:40)'; invdlTable.pred_10 = (1:40)'; invdlTable.pred_20 = (1:40)'; invdlTable.pred_fr = (1:40)';
for nn = 1:40
    invdlTable.pred_1(nn) = ds.utils.getKstepSMSE(dsNaNCell{nn}, 1, 'standardize', normaliseFactor);
    invdlTable.pred_10(nn) = ds.utils.getKstepSMSE(dsNaNCell{nn}, 10, 'standardize', normaliseFactor);
    invdlTable.pred_20(nn) = ds.utils.getKstepSMSE(dsNaNCell{nn}, 20, 'standardize', normaliseFactor);
    invdlTable.pred_fr(nn) = ds.utils.getKstepSMSE(dsNaNCell{nn}, Inf, 'standardize', normaliseFactor);
end

    %% Visualise difference between batch parameters, and those with individual biases.
    % The individual bias SMSE is better, but also qualitatively very
    % encouraging for use with nonlinearity. Nevertheless, one imagines
    % emission parameters may be required to be very different.
    figh1 = figure;
    figh2 = figure;
    
for nn = 1:40
    ds.plot.predictions(dsBatchBiasCov, [1,10,20,Inf], pump.target, 'nn', nn, 'titlePrefix', sprintf('Patient %d',nn), 'verbose', false, 'figure', figh1); 
    ds.plot.predictions(dsBatch, [1,10,20,Inf], pump.target, 'nn', nn, 'titlePrefix', sprintf('Patient %d',nn), 'verbose', false, 'figure', figh2);
    pause
end
