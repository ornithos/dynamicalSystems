%% test data
if false
    load 20170327NaNTestY   % gives us yy (test data)
    % test data generated from PK model for patient 15
    ss     = 3;
    os     = 2;
    A      = diag(rand(ss, 1));
    H      = rand(os, ss)*2 - 1;

    dsLDS  = ds.dynamicalSystem(ss, os, 'x0', [1;0;0], eye(ss)*100, ...
        'evolution', A, 1e-2, 'emission', H, 0.01, ...
        'data', yy, 'control', Us{15}, true, false);

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


for ii = 14:40
    fprintf('*********** Patient %d of %02d ******************* \n', ii, maxPx);
    R      = 0.2;
    ss     = 4;
    os     = size(Ys{ii},1);
    opts   = struct('warnings', false);
    
    yy     = Ys{ii};
    rmObs  = logical(annotations{ii});
    rmObs(1:2,:)  = repmat(rmObs(1,:) | rmObs(2,:), 2, 1);
    yy(rmObs) = NaN;
    
    % initialise randomly / DA
    A      = diag(rand(ss, 1));
    H      = rand(os, ss)*2 - 1;
    dsINL  = ds.dynamicalSystem(ss, os, 'x0', pars{ii}.x0, pars{ii}.V0, ...
        'evolution', A, 1, 'emission', H, R, 'data', yy, 'control', Us{ii}, true, false, opts);
    
    opts = struct('maxiter', 3000, 'epsilon', 5e-4, 'sampleStability', 5, ...
           'multistep', 4, 'ssid', false, 'verbose', true, 'stableVerbose', false, ...
           'annealingSchedule', 0.1, 'annealingIter', 100, 'annealingMin', 1e-3, ...
           'diagA', true, 'diagQ', false);
    llhhist = dsINL.parameterLearningEM(opts);
    dsINL.filter([],true);
    dsINL.smooth;
    dsINL.save('initial');
    
    % visualise
    figure
    
    [impy, impPy] = dsINL.impute_y('variance', true, 'smooth', true);
    plot(Ys{ii}', ':'); hold on;
    plot(impy'); hold off;
    stdy          = sqrt(cell2mat(cellfun(@(x) diag(x)', impPy, 'Un', 0)))';
    cols          = {'blue', 'red', 'green'};
    for jj = 1:os
        nanidx          = find(isnan(dsINL.y(jj,:)));
        if isempty(nanidx); continue; end
        nanminmax       = nanidx(1):nanidx(end);
        utils.plot.confidenceInterval(nanminmax, impy(jj,nanminmax), stdy(jj,nanminmax), [], 'facecolor', cols{jj});
    end
    cmap           = flipud([1 0.86 0.91; 0.86 0.91 1]);
    utils.plot.dataShadeVertical(1:(dsINL.d.T), pump.target{ii}', cmap, 'edgedetect', 'manual', 'edgemanual', 3);
    dsNaNCell{ii} = dsINL.copy;
end


%% visualise predictions
figure;
h = plot(1,1);
testpredahead = [1, 10, 20, Inf];

cols          = zeros(3,3);
cnums         = [2 4 5];
for kk = 1:3; cols(kk,:) = utils.plot.varyColor2(cnums(kk)); end
litecols      = utils.plot.colortint(litecols, 0.8);
cols          = utils.plot.colortint(litecols, 1.4);

for ii = 1:13
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
