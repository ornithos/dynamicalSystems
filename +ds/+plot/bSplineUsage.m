function bSplineUsage(obj, predT, dim, varargin)
    assert(isa(obj, 'ds.ionlds'), 'object is not an IONLDS model');
    assert(isnumeric(predT) && all(utils.is.int(predT)), 'predT must be a vector of integers');
    assert(isnumeric(dim) && all(utils.is.int(dim, 1)) && all(dim <= obj.d.y), 'dim must be a vector of integers in 1, ..., %d', obj.d.y);
    
    optsDefault.showY  = false;
    optsDefault.useX   = false;
    optsDefault.showLines = false;
    optsDefault.labelPts = false;
    opts               = utils.base.processVarargin(varargin, optsDefault);
    
    nPreds = numel(predT);
    nDims  = numel(dim);
    
    isBatch = isa(obj, 'ds.dynamicalSystemBatch');
    if isBatch
        N = obj.d.n;
    else
        N = 1;
    end
    
    predvals = cell(nPreds, 1);
    X = cell(nPreds,1);
    y = obj.y;
    if ~isBatch; y = {y}; end
    
    for jj = 1:nPreds
        [predvals{jj}, X{jj}] = obj.getPredictedValues(predT(jj), struct('alpha', 1, 'beta', 0, 'kappa', 0));
        if ~isBatch; predvals{jj} = {predvals{jj}}; X{jj} = {X{jj}}; end
        X{jj} = cellfun(@(x) bsxfun(@plus, obj.par.emiNLParams.C * x, obj.par.emiNLParams.bias), X{jj}, 'Un', 0);
    end
        
    spdims = utils.plot.subplotdims(nPreds);
    
    figure;
    
    if nPreds > 1 || N > 1
        % outer loop over dimensions
        for ii = 1:nDims
            cDim = dim(ii);
            for nn = 1:N

                % find min-max range over all prediction windows
                xrng = cellfun(@(x) [min(x{nn}(cDim,:),[],2), max(x{nn}(cDim,:),[],2)], X, 'Un', 0);
                xrng = cell2mat(xrng);
                xrng = [min(xrng(:,1)), max(xrng(:,2))];
                
                for jj = 1:nPreds
                    subplot(spdims(1), spdims(2), jj);

                    fvals = xrng(1):.1:xrng(2);
                    plot(fvals, obj.par.emiNLParams.bspl.functionEval(fvals, obj.par.emiNLParams.eta(cDim,:)), 'k-');
                    hold on;
                    plot(X{jj}{nn}(cDim,:), predvals{jj}{nn}(cDim,:), '*'); 
                    if opts.labelPts
                        text(X{jj}{nn}(cDim,:), predvals{jj}{nn}(cDim,:), cellstr(num2str((1:size(X{jj}{nn},2))')), 'VerticalAlignment', 'top'); 
                    end
                    if opts.showY
                        yy = y{nn}(cDim,:);
                        yy = yy(numel(yy)-size(X{jj}{nn},2)+1:end);
                        plot(X{jj}{nn}(cDim,:), yy, 'ro');
                        if opts.showLines
                            line([X{jj}{nn}(cDim,:); X{jj}{nn}(cDim,:)], [yy; predvals{jj}{nn}(cDim,:)]);
                        end
                    end
                    hold off
                    xlim([xrng(1), xrng(2)]);
                    titletxt = sprintf('Dimension %d, Prediction window: %d', cDim, predT(jj));
                    if N > 1; titletxt = sprintf('SERIES %d --- %s', nn, titletxt); end
                    title(titletxt);
                end
                if nn < N || ii < nDims
                    pause
                end
            end
        end
    else
        % dimensions on same subplot
        dimsOnSamePlot(obj, dim, X, predvals, predT, opts)
    end 
end

function dimsOnSamePlot(obj, dim, x, predvals, predT, opts)
    nDims  = numel(dim);
    spdims = utils.plot.subplotdims(nDims);
    
    x      = x{1}{1};
    if opts.useX; x = obj.par.emiNLParams.C * obj.x(:,2:end) + obj.par.emiNLParams.bias; end
    predvals = predvals{1}{1};
    
    for ii = 1:nDims
        subplot(spdims(1), spdims(2), ii);
        cDim = dim(ii);
        xrng = [min(x(cDim,:)), max(x(cDim,:))];
        fvals = xrng(1):.1:xrng(2);
        plot(fvals, obj.par.emiNLParams.bspl.functionEval(fvals, obj.par.emiNLParams.eta(cDim,:)), 'k-');
        hold on;
        plot(x(cDim,:), predvals(cDim,:), '*');
        if opts.labelPts; text(x(cDim,:), predvals(cDim,:), cellstr(num2str((1:size(x,2))')), 'VerticalAlignment', 'top'); end
        if opts.showY
            yy = obj.y(cDim,:);
            yy = yy(numel(yy)-size(x,2)+1:end);
            plot(x(cDim,:), yy, 'ro');
            if opts.showLines
                line([x(cDim,:); x(cDim,:)], [yy; predvals(cDim,:)]);
            end
        end
        hold off
        xlim([xrng(1), xrng(2)]);
        titletxt = sprintf('Dimension %d, Prediction window: %d', cDim, predT(1));
        title(titletxt);
    end
end