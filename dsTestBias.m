%% test data
if false
    load 20170327NaNTestY   % gives us yy (test data)
    % test data generated from PK model for patient 15
    ss     = 3;
    os     = 2;
    A      = diag(rand(ss, 1));
    H      = rand(os, ss)*2 - 1;

    dsLDS  = ds.dynamicalSystem(ss, os, 'x0', [1;0;0;0], eye(ss)*100, ...
        'evolution', A, 1e-2, 'emission', H, 0.01, ...
        'data', Ys{12}, 'control', Us{12}, true, false);

    opts = struct('maxiter', 1600, 'epsilon', 5e-4, 'sampleStability', 5, ...
               'multistep', 4, 'ssid', false, 'verbose', true, 'stableVerbose', false, ...
               'annealingSchedule', 0.1, 'annealingIter', 100, 'annealingMin', 1e-3, ...
               'diagA', true, 'diagQ', false);

    dsLDS.parameterLearningEM(opts);
    dsLDS.filter([],true);
    dsLDS.smooth;


    % Visualise

    [impy, impPy] = dsLDS.impute_y('variance', true, 'smooth', true);
    plot(impy');
    stdy          = sqrt(cell2mat(cellfun(@(x) diag(x)', impPy, 'Un', 0)))';
    cols          = {'blue', 'red'};
    for jj = 1:os
        nanidx          = find(isnan(dsLDS.y(jj,:)));
        nanminmax       = nanidx(1):nanidx(end);
        utils.plot.confidenceInterval(nanminmax, impy(jj,nanminmax), stdy(jj,nanminmax), [], 'facecolor', cols{jj});  %#ok
    end
end

%% Now for annotations of MLHC Data
% Note: must have run loadMLHCData first

if ~exist('dsNaNCell', 'var')
    dsNaNCell = cell(40,1);
end
cols          = zeros(3,3);
cnums         = [2 4 5];
for kk = 1:3; cols(kk,:) = utils.plot.varyColor2(cnums(kk)); end
litecols      = utils.plot.colortint(cols, 0.8);
cols          = utils.plot.colortint(litecols, 1.4);
pforward      = 150;

for ii = 1:13
    fprintf('*********** Patient %d ******************* \n', ii);
    R      = 0.2;
    ss     = 4;
    os     = size(Ys{ii},1);
    opts   = struct('warnings', false);
    
    uu     = Us{ii};
    yy     = Ys{ii};
    rmObs  = logical(annotations{ii});
    rmObs(1:2,:)  = repmat(rmObs(1,:) | rmObs(2,:), 2, 1);
    yy(rmObs) = NaN;
    yy        = yy(:,1:end);
    uu        = uu(:,1:end);
    
    % initialise randomly / DA
    A      = diag(rand(ss, 1));
    H      = rand(os, ss)*2 - 1;
    dsLDS  = ds.dynamicalSystem(ss, os, 'x0', zeros(ss,1), eye(ss)*1e-6, ...
        'evolution', A, 1, 'emission', H, R, nanmean(yy(:,1:10),2), ...
        'data', yy, 'control', uu, true, false, opts);
    
    epsilon = 0.01;
    opts = struct('maxiter', 1000, 'epsilon', 5e-4, 'sampleStability', 5, ...
           'multistep', 4, 'ssid', false, 'verbose', true, 'stableVerbose', false, ...
           'annealingSchedule', 0.1, 'annealingIter', 100, 'annealingMin', 1e-3, ...
           'diagA', true, 'diagQ', false, 'diagAconstraints', [-1+epsilon, 1-epsilon], 'fixBias2', true, 'strictNegativeCheck', false);
    llhhist = dsLDS.parameterLearningEM(opts);
    dsLDS.filter([],true);
    dsLDS.smooth;
    dsLDS.save('initial');
    
    
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
end


%% visualise predictions
figure;
plthdl = plot(1,1);
testpredahead = [1, 10, 20, Inf];

cols          = zeros(3,3);
cnums         = [2 4 5];
for kk = 1:3; cols(kk,:) = utils.plot.varyColor2(cnums(kk)); end
litecols      = utils.plot.colortint(cols, 0.8);
cols          = utils.plot.colortint(litecols, 1.4);

for ii = 1:4
    os = dsNaNCell{ii}.d.y;
    for jj = 1:4
        subplot(2,2,jj);

        [impy, impPy] = dsNaNCell{ii}.impute_y('variance', true, 'smooth', true);
        for kk = 1:os
            plot(Ys{ii}(kk,:)', ':', 'Color', litecols(kk,:)); hold on;
            plot(impy(kk,:)', 'Color', litecols(kk,:));
        end
        hold off;
        stdy          = sqrt(cell2mat(cellfun(@(x) diag(x)', impPy, 'Un', 0)))';
        
        for kk = 1:os
            nanidx          = find(isnan(dsNaNCell{ii}.y(kk,:)));
            if isempty(nanidx); continue; end
            nanminmax       = nanidx(1):nanidx(end);
            utils.plot.confidenceInterval(nanminmax, impy(kk,nanminmax), stdy(kk,nanminmax), [], 'facecolor', cols(kk,:));
        end
        if jj == 4
            predvals = dsNaNCell{ii}.getPredictFreeRun(1);
        else
            prdahead = testpredahead(jj);
            predvals = dsNaNCell{ii}.getPredictedValues(prdahead);
        end
        T         = dsNaNCell{ii}.d.T;
        ix        = T-size(predvals,2)+1:T;
        smse      = mean(nansum((dsNaNCell{ii}.y(:,ix) - predvals).^2,2)./(nanvar(dsNaNCell{ii}.y(:,ix), [], 2).*sum(~isnan(dsNaNCell{ii}.y(:,ix)),2)));
        hold on; 
        for kk = 1:os
            plot(ix, predvals(kk,:), '-', 'Color', cols(kk,:)); 
        end
        hold off;
        cmap           = flipud([1 0.86 0.91; 0.86 0.91 1]);
        utils.plot.dataShadeVertical(1:(dsNaNCell{ii}.d.T), pump.target{ii}', cmap, 'edgedetect', 'manual', 'edgemanual', 3);
        title(sprintf('Patient %d: (%d-step). Average SMSE %.2f%%', ii, testpredahead(jj), smse*100));
    end
%     if os == 2
%         legend({'BPsys', 'BPdia'});
%     else
%         legend({'BPsys', 'BPdia', 'BIS'});
%     end
    pause
end

    %% visualise return to steady state
    figure
    cols          = utils.plot.colortint(litecols, 1.2);
    for ii = 1:40
        dsLDS = dsNaNCell{ii}.copy;
        [impy, impPy] = dsLDS.impute_y('variance', true, 'smooth', true);
        futurY        = dsLDS.getPredictFreeRun(dsLDS.d.T, pforward);
        plot(Ys{ii}(1,:)', ':', 'Color', cols(1,:)); 
        hold on;
        for kk = 2:dsLDS.d.y
            plot(Ys{ii}(kk,:)', ':', 'Color', cols(kk,:)); 
        end
        for kk = 1:dsLDS.d.y
            plot(impy(kk,:)', 'Color', cols(kk,:));
            plot(dsLDS.d.T+1:dsLDS.d.T+pforward, futurY(kk,:), 'Color', litecols(kk,:));
        end
        
        hold off;
        stdy          = sqrt(cell2mat(cellfun(@(x) diag(x)', impPy, 'Un', 0)))';
        for jj = 1:dsLDS.d.y
            nanidx          = find(isnan(dsLDS.y(jj,:)));
            if isempty(nanidx); continue; end
            nanminmax       = nanidx(1):nanidx(end);
            utils.plot.confidenceInterval(nanminmax, impy(jj,nanminmax), stdy(jj,nanminmax), [], 'facecolor', cols(kk,:));
        end
        cmap           = flipud([1 0.86 0.91; 0.86 0.91 1]);
        utils.plot.dataShadeVertical(1:(dsLDS.d.T), pump.target{ii}', cmap, 'edgedetect', 'manual', 'edgemanual', 3);
        pause
    end