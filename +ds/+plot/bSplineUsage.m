function bSplineUsage(obj, predT, dim, varargin)
    assert(isa(obj, 'ds.ionlds'), 'object is not an IONLDS model');
    assert(isnumeric(predT) && all(utils.is.int(predT)), 'predT must be a vector of integers');
    assert(isnumeric(dim) && all(utils.is.int(dim, 1)) && all(dim <= obj.d.y), 'dim must be a vector of integers in 1, ..., %d', obj.d.y);
    
    optsDefault.showY  = false;
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
%     figure
    for ii = 1:nDims
        cDim = dim(ii);
        for nn = 1:N
            
            xrng = cellfun(@(x) [min(x{nn}(cDim,:),[],2), max(x{nn}(cDim,:),[],2)], X, 'Un', 0);
            xrng = cell2mat(xrng);
            xrng = [min(xrng(:,1)), max(xrng(:,2))];
            for jj = 1:nPreds
%                 subplot(spdims(1), spdims(2), jj);
                
                fvals = xrng(1):.1:xrng(2);
                plot(fvals, obj.par.emiNLParams.bspl.functionEval(fvals, obj.par.emiNLParams.eta(cDim,:)), 'k-');
                hold on;
                plot(X{jj}{nn}(cDim,:), predvals{jj}{nn}(cDim,:), '*'); 
                if opts.showY
                    yy = y{nn}(cDim,:);
                    plot(X{jj}{nn}(cDim,:), yy(numel(yy)-size(X{jj}{nn},2)+1:end), 'ro');
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
end